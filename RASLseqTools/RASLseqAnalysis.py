
import pandas as pd
import os,sys
import numpy as np
import argparse

from RASLseqAligner import *


class RASLseqAnalysis_obj(object):
    '''
    This class creates a pandas DataFrame for RASLseq fastq sequences.

    Attributes of this class annotate RASLseq fastq sequences.

    The get_attributes method allows the user to specify what additional
    information to add to the pandas DataFrame.


    Attributes
    ----------

    fastq_file: str
        path to Fastq file

    probes_path: str
        path to Probe file

    blastdb: str
        path to write a blast database

    blastn_dir: str
        path to directory holding the blastn executable

    well_annot: str
        path to Barcode Annotations file

    write: str
        path to write directory

    print_on: boolean, default=False
        Whether to print information during data processing

    read_df: pandas DataFrame
        index: (PlateBarcode, sequence)
        columns: ['PlateBarcode', 'seq', 'seq_count']
            PlateBarcode - Index read from fastq file
            seq - fastq read sequence
            seq_count - number of occurrences of seq in fastq


    APRIL 14 2014: ToDos:
        2) Refactor external methods ( raslblast)
        3) Package everything
        4) Consider GUI

    '''

    def __init__(self, fastq_path, sequencer_id, probes_path, blastdb_path, blastn_dir, well_annot, write_path, print_on=False ):
        #create a probe with name name

        self.RASLseqReads_obj = RASLseqReads.RASLseqReads(fastq_path, sequencer_id, print_on)

        self.RASLseqProbes_obj = RASLseqProbes.RASLseqProbes(probes_path, blastdb_path, blastn_dir)

        self.RASLseqBCannot_obj = RASLseqBCannot.RASLseqBCannot(well_annot)

        self.RASLseqBCannot_obj.well_bc = self.RASLseqBCannot_obj.well_bc

        self.sequencer_id = sequencer_id

        self.print_on = print_on

        self.write_path = write_path

        #self.write_obj = open(write_path,'w')

        self.read_df = self.RASLseqReads_obj.read_df

        self.fastq_read_count = self.read_df.seq_count.sum()  #number of total reads found in fastq input


        #filter thresholds
        self.bc_edit_dist_filter = 2

        self.blast_results_filter = {'length':30, 'qstart':6, 'obs_wellbc_len_max':10, 'obs_wellbc_len_min':6}



    def _get_probe_well_read_counts(self, collapsed_read_counts):
        '''
        This function aggregates probe-specific read counts for each plate and well

        Parameters
        ----------
        collapsed_read_counts: pandas dataframe
            must possess cols: ['plate_barcode','mapped_bc','probe']

        Returns
        -------
        Pandas Dataframe
            index: ['plate_barcode','WellBarcode']
            columns: Probe-specific read count sum (sum of counts across reads
            mapping by BLAST to probe)
        '''

        #grouping reads by plate barcode, well barcode, and probe
        collapsed_read_counts_group = collapsed_read_counts.groupby(['plate_barcode','mapped_bc','probe'])

        #aggregating probe counts for each well
        counts = collapsed_read_counts_group.seq_count.aggregate(np.sum)
        counts_df = counts.unstack(level=2)  #creating matrix of aggregated probe counts indexed on 'plate_barcode','mapped_bc'


        #id and removal of off-target ligation hits found by ! character in blast sseq name
        on_target_col = [i for i in counts_df.columns if "!" not in i]
        counts_df = counts_df[on_target_col]  #removing off-target ligation hits from df
        counts_df.index.names= ['PlateBarcode', 'WellBarcode']
        return counts_df



    def _merge_plate_well_annot(self, probe_counts_df, well_annot_df):
        '''
        This function merges gene_counts_df with well annotations

        Parameters
        ----------
        probe_counts_df: Pandas DataFrame
            Requires pandas index: ('plate_barcode','WellBarcode')

        well_annot_path: Pandas DataFrame
            Requires pandas index: ('plate_barcode','WellBarcode')

        Returns
        -------
        Pandas DataFrame
            well_annot_df right joined to gene_counts_df
            index: ('plate_barcode','WellBarcode')

        '''
        return well_annot_df.join(probe_counts_df,how='right')




    def get_target_counts_df(self):
        '''
        This method uses a functional programming design to
        produce a pandas DataFrame describing the RASLseq probe-specific
        read counts prepended with well specific annotations.


        '''

        #Add observed wellbc
        self.read_df = RASLseqWellbc.get_rasl_probe_and_wellbc_exact_df(self.read_df, print_on=self.print_on)

        #add fuzzy matched wellbc (mapped_bc)
        self.read_df = RASLseqWellbc.get_rasl_probe_and_wellbc_fuzzy_df(self.read_df, self.RASLseqBCannot_obj.well_bc, \
                            print_on=self.print_on)

        #filter fuzzy wellbc mappings
        self.read_df = self.read_df[(self.read_df.bc_edit_dist.astype(float) < self.bc_edit_dist_filter)]

        #add observed rasl_probe_seq
        self.read_df = RASLseqSeq.get_rasl_probe_exact_df(self.read_df)

        #blast alignment
        self.read_df = RASLseqAlign.get_rasl_blast_df(self.read_df, self.RASLseqProbes_obj.blastn_dir, \
                            self.RASLseqProbes_obj.blastdb_file, print_on=self.print_on)

        #filtering blast results
        self.read_df = self.read_df[(self.read_df.length > self.blast_results_filter['length']) & \
                                    (self.read_df.qstart < self.blast_results_filter['qstart']) & \
                                    (self.read_df.observed_wellbc.map(len) < \
                                            self.blast_results_filter['obs_wellbc_len_max']) & \
                                    (self.read_df.observed_wellbc.map(len) > \
                                            self.blast_results_filter['obs_wellbc_len_min'])]

        #SUM PROBE-SPECIFIC READ COUNTS
        self.read_df_counts = self._get_probe_well_read_counts(self.read_df)

        self.read_count_mapped = self.read_df_counts.sum().sum()


        #MERGING WELL ANNOTATIONS AND PROBE READ COUNTS
        #return self.read_df_counts, self.RASLseqBCannot_obj.well_annot_df
        self.RASLseqAnalysis_df = self._merge_plate_well_annot(self.read_df_counts, self.RASLseqBCannot_obj.well_annot_df)

        #WRITING DATAFRAMES TO FILE
        self.RASLseqAnalysis_df.to_csv(self.write_path, sep="\t", index=False)  #final analysis df
        self.read_df_counts.to_csv(self.write_path+'.read_df_counts', sep="\t")  #probe reads counts by well df
        self.read_df.to_csv(self.write_path+'.read_df', sep="\t")  #all unique sequences with blast alignments df



parser = argparse.ArgumentParser()
parser.add_argument('-f','--fastq', type=str,help='Specifies the input fastq file, /path/to/CG_data/RASLseq.fq')
parser.add_argument('-s','--sequencer_id', type=str,help='Specifies the sequencer identifier in the fastq index lines, e.g. @HISEQ')
parser.add_argument('-p','--probes', type=str,help='Specifies the input probes file containing the following columns: AcceptorProbeSequence DonorProbeSequence AcceptorAdaptorSequence DonorAdaptorSequence ProbeName, /path/to/probes.txt')
parser.add_argument('-w','--well_annot', type=str,help='Specifies the input well annotations file containing the following columns: PlateBarcode and WellBarcode, /path/to/well/annotations.txt')
parser.add_argument('-d','--blastdb_path', type=str,help='Specifies the directory for writing blast database, /path/to/write/blastdb/')
parser.add_argument('-b','--blastn_bin', type=str,help='Specifies the path to directory holding blastn executable, /path/to/blastn_dir/')
parser.add_argument('-P','--print_on', action='store_true', default=False,help='Specifies whether to print summary stats during alignment')
parser.add_argument('-o','--output_path', type=str,help='Specifies the output file path, e.g. /path/to/output/counts_df.txt')


#opts = parser.parse_known_args()
#fastq_path, sequencer_id, probes_path, blastdb_path, blastn_dir = opts[0].fastq, opts[0].sequencer_id, opts[0].probes, opts[0].blastdb_path, opts[0].blastn_bin
#well_annot, write_path, print_on_bool = opts[0].well_annot, opts[0].output_path, opts[0].print_on




#if __name__ == '__main__':
#
#
#    rasl_analysis = RASLseqAnalysis_obj(fastq_path, sequencer_id, probes_path, blastdb_path, blastn_dir, well_annot, write_path, print_on=print_on_bool)
#
#    rasl_analysis.get_target_counts_df()
#
#    print 'Complete'



"""
TEST
python ~/git/raslseq/dev/RASLseqAligner/RASLseqAnalysis.py -f /gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/fastq/P1_CGGTTCT_L001_R1_001_100k.fastq.gz -p /gpfs/home/erscott/Datasets_raw/temp/probes.txt -w /gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/Bcell_exp_2013_10_bc_key.txt -d /gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/blastdb/ -b /gpfs/home/erscott/Tools/blast/ncbi-blast-2.2.26+/bin/ -P -o /gpfs/home/erscott/temp/2014_04_21_output.txt
"""

