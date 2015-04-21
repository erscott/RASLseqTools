import pandas as pd
import os,sys
import numpy as np
import argparse

import inspect
source_dir =  os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
source_dir = '/'.join(source_dir.split('/')[:-1])
sys.path.append(source_dir)
from RASLseqTools import *



class RASLseqAnalysis_BLAST(object):
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

    def __init__(self, fastq_path, sequencer_id, probes_path, blastdb_path, blastn_dir, well_annot, write_path, print_on=False,\
                 offset_5p=24, offset_3p=22, wellbc_start=0, wellbc_end=8):
        #create a probe with name name
        
        self.RASLseqReads_obj = RASLseqReads.RASLseqReads(fastq_path, sequencer_id, print_on)
        
        self.RASLseqProbes_obj = RASLseqProbes.RASLseqProbes(probes_path, blastdb_path, blastn_dir, aligner='blast')   
        
        self.RASLseqBCannot_obj = RASLseqBCannot.RASLseqBCannot(well_annot)
        
        self.RASLseqBCannot_obj.well_bc = self.RASLseqBCannot_obj.well_bc
        
        self.sequencer_id = sequencer_id
        
        self.print_on = print_on
        
        self.write_path = write_path
        
        self.offset_5p = int(offset_5p)
        
        self.offset_3p = int(offset_3p)
        
        self.wellbc_start = int(wellbc_start)
        
        self.wellbc_end = int(wellbc_end)
        
        self.read_df = self.RASLseqReads_obj.get_blast_read_df()
        
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
        #self.read_df = RASLseqWellbc.get_rasl_probe_and_wellbc_exact_df(self.read_df, print_on=self.print_on)
        self.read_df['observed_wellbc'] = self.read_df['seq'].str[self.wellbc_start:self.wellbc_end]
         
        #add fuzzy matched wellbc (mapped_bc)
        self.read_df = RASLseqWellbc.get_rasl_probe_and_wellbc_fuzzy_df(self.read_df, self.RASLseqBCannot_obj.well_bc, \
                            print_on=self.print_on)         
        
        #filter fuzzy wellbc mappings
        self.read_df = self.read_df[(self.read_df.bc_edit_dist.astype(float) < self.bc_edit_dist_filter)]
        
        #add observed rasl_probe_seq
        #self.read_df = RASLseqSeq.get_rasl_probe_exact_df(self.read_df)
        self.read_df['rasl_probe'] = self.read_df['seq'].str[self.offset_5p : -self.offset_3p]
        
        #blast alignment
        self.read_df = RASLseqAlign.get_rasl_aligned_df(self.read_df, self.RASLseqProbes_obj.aligner_dir, \
                            self.RASLseqProbes_obj.probedb_file, print_on=self.print_on, aligner='blast')
        
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
        
        
        self.RASLseqAnalysis_df.to_csv(self.write_path, sep="\t")
        self.read_df_counts.to_csv(self.write_path+'.read_df_counts', sep="\t")
        self.read_df.to_csv(self.write_path+'.read_df', sep="\t")
        
        os.system('gzip ' + self.write_path+'.read_df')
        
        return

parser = argparse.ArgumentParser()
parser.add_argument('-f','--fastq', type=str,help='Specifies the input fastq file, /path/to/CG_data/RASLseq.fq')

parser.add_argument('-s','--sequencer_id', type=str,help='Specifies the sequencer identifier in the fastq index lines, e.g. @HISEQ')

parser.add_argument('-p','--probes', type=str,help='Specifies the input probes file containing the following columns: AcceptorProbeSequence DonorProbeSequence AcceptorAdaptorSequence DonorAdaptorSequence ProbeName, /path/to/probes.txt') 

parser.add_argument('-w','--well_annot', type=str,help='Specifies the input well annotations file containing the following columns: PlateBarcode and WellBarcode, /path/to/well/annotations.txt')

parser.add_argument('-d','--blastdb_path', type=str,help='Specifies the directory for writing blast database, /path/to/write/blastdb/')

parser.add_argument('-b','--blastn_bin', type=str,help='Specifies the path to directory holding blastn executable, /path/to/blastn_dir/')

parser.add_argument('-P','--print_on', action='store_true', default=False,help='Specifies whether to print summary stats during alignment')

parser.add_argument('-o','--output_path', type=str,help='Specifies the output file path, e.g. /path/to/output/counts_df.txt')

parser.add_argument('-o5','--offset_5p', type=int, default=24, help='Specifies the number of bases to clip from 5-prime end of read to isolate probe sequence, default 24')

parser.add_argument('-o3','--offset_3p', type=int, default=22, help='Specifies the number of bases to clip from 3-prime end of read to isolate probe sequence, default 22')

parser.add_argument('-ws','--wellbc_start', type=int, default=0, help='Specifies the index position of the wellbc start base, default 0')

parser.add_argument('-we','--wellbc_end', type=int, default=8, help='Specifies the index position of the wellbc start base, default 8')


opts = parser.parse_known_args()
fastq_path, sequencer_id, probes_path, blastdb_path, blastn_dir = opts[0].fastq, opts[0].sequencer_id, opts[0].probes, opts[0].blastdb_path, opts[0].blastn_bin 
well_annot, write_path, print_on_bool = opts[0].well_annot, opts[0].output_path, opts[0].print_on
offset_5p, offset_3p, wellbc_start, wellbc_end = opts[0].offset_5p, opts[0].offset_3p, opts[0].wellbc_start, opts[0].wellbc_end
            
        




if __name__ == '__main__':
    
    rasl_analysis = RASLseqAnalysis_BLAST(fastq_path, sequencer_id, probes_path, blastdb_path, blastn_dir, well_annot, write_path, print_on=print_on_bool, \
                                          offset_5p=offset_5p, offset_3p=offset_3p, wellbc_start=wellbc_start, wellbc_end=wellbc_end)
    
    rasl_analysis.get_target_counts_df()
    
    print 'Complete'
    



"""
TEST
python /path/to/RASLseqAnalysis.py -s '@HISEQ' -f /path/to/your.fastq.gz -p /paht/to/RASL.probes -w /path/to/annotations.bc -d /path/to/blastdb/write_dir/ -b /path/to/blast/ncbi-blast-2.2.26+/bin/ -P -o /path/to/output.txt
"""