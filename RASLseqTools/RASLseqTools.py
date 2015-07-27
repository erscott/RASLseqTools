from dataframes.RASLseqProbes import RASLseqProbes

import numpy as np
import argparse
import raven
import uuid
import sys
import os


PROJECT_PATH = os.path.realpath(os.path.join(os.path.dirname(__file__), os.path.pardir))
client = raven.Client(dsn='http://4939f4e75211487ab3821e1af9d9f3d2:c8569ad1ff844e6f8b7284fca7c0bffe@sentry.sulab.org/7', release='1.0.0')


class RASLseqTools(object):

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

    def __init__(self,
            fastq_path, probes_path, well_annot_path,
            aligner_dir,

            offset_5p=24, offset_3p=22,
            wellbc_start=0, wellbc_end=8,

            blastdb_path='',
            verbose=False,

            n_job=1, write_alignments=False):

        # (TODO) Validate these paths now
        self.fastq_path = fastq_path
        self.probes_path = probes_path
        self.well_annot_path = well_annot_path

        # (TODO) Validate this path now
        self.aligner_dir = aligner_dir

        # (TODO) Validate these ints. Ex. A > B
        self.offset_5p = int(offset_5p)
        self.offset_3p = int(offset_3p)
        self.wellbc_start = int(wellbc_start)
        self.wellbc_end = int(wellbc_end)

        self.blastdb_path = blastdb_path
        self.verbose = bool(verbose)

        # Setup path for saving results
        self.results_dir = PROJECT_PATH + '/results/'
        try:
            os.mkdir(self.results_dir)
        except OSError as e:
            if not e.errno == 17:
                pass
            else:
                client.captureException()
        self.run_results_dir = self.results_dir + str(uuid.uuid1())
        os.mkdir(self.run_results_dir)


        if self.aligner == 'blast':
            # (TODO) Does this happen with STAR (? why was it excluded ?)
            self.RASLseqReads_obj = RASLseqReads.RASLseqReads(fastq_path, sequencer_id, print_on)
            self.read_df = self.RASLseqReads_obj.get_blast_read_df()
            self.fastq_read_count = self.read_df.seq_count.sum()  # Number of total reads found in fastq input

            self.sequencer_id = sequencer_id

            # >> Filter thresholds
            self.bc_edit_dist_filter = 2
            self.blast_results_filter = {'length': 30, 'qstart': 6, 'obs_wellbc_len_max': 10, 'obs_wellbc_len_min': 6}

        if self.aligner == 'star':
            self.n_jobs = int(n_jobs)
            self.write_alignments = write_alignments
            self.bc_edit_dist_filter = 2

        # Dataframe assignment
        self.RASLseqProbes_obj = RASLseqProbes(self.probes_path, self.aligner_dir, self.run_results_dir, aligner=self.aligner)

        '''
        self.RASLseqBCannot_obj = RASLseqBCannot.RASLseqBCannot(well_annot)
        self.RASLseqBCannot_obj.well_bc = self.RASLseqBCannot_obj.well_bc
        '''

    @property
    def aligner(self):
        return 'star'

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

        # Grouping reads by plate barcode, well barcode, and probe
        collapsed_read_counts_group = collapsed_read_counts.groupby(['plate_barcode', 'mapped_bc', 'probe'])

        # Aggregating probe counts for each well
        counts = collapsed_read_counts_group.seq_count.aggregate(np.sum)
        counts_df = counts.unstack(level=2)  # Creating matrix of aggregated probe counts indexed on 'plate_barcode','mapped_bc'

        # ID and removal of off-target ligation hits found by ! character in blast sseq name
        on_target_col = [i for i in counts_df.columns if "!" not in i]
        counts_df = counts_df[on_target_col]  # Removing off-target ligation hits from df
        counts_df.index.names = ['PlateBarcode', 'WellBarcode']
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
        return well_annot_df.join(probe_counts_df, how='right')

    def get_target_counts_df(self):
        '''
            This method uses a functional programming design to
            produce a pandas DataFrame describing the RASLseq probe-specific
            read counts prepended with well specific annotations.
        '''

        # Add observed wellbc
        self.read_df = RASLseqWellbc.get_rasl_probe_and_wellbc_exact_df(self.read_df)

        # Add fuzzy matched wellbc (mapped_bc)
        self.read_df = RASLseqWellbc.get_rasl_probe_and_wellbc_fuzzy_df(self.read_df,
            self.RASLseqBCannot_obj.well_bc)

        # Filter fuzzy wellbc mappings
        self.read_df = self.read_df[(self.read_df.bc_edit_dist.astype(float) < self.bc_edit_dist_filter)]

        # Add observed rasl_probe_seq
        self.read_df = RASLseqSeq.get_rasl_probe_exact_df(self.read_df)

        # Blast alignment
        self.read_df = RASLseqAlign.get_rasl_blast_df(self.read_df,
            self.RASLseqProbes_obj.blastn_dir,
            self.RASLseqProbes_obj.blastdb_file)

        # Filtering blast results
        self.read_df = self.read_df[
            (self.read_df.length > self.blast_results_filter['length']) &
            (self.read_df.qstart < self.blast_results_filter['qstart']) &
            (self.read_df.observed_wellbc.map(len) < self.blast_results_filter['obs_wellbc_len_max']) &
            (self.read_df.observed_wellbc.map(len) > self.blast_results_filter['obs_wellbc_len_min'])]

        # SUM PROBE-SPECIFIC READ COUNTS
        self.read_df_counts = self._get_probe_well_read_counts(self.read_df)
        self.read_count_mapped = self.read_df_counts.sum().sum()

        # MERGING WELL ANNOTATIONS AND PROBE READ COUNTS
        # return self.read_df_counts, self.RASLseqBCannot_obj.well_annot_df
        self.RASLseqAnalysis_df = self._merge_plate_well_annot(self.read_df_counts, self.RASLseqBCannot_obj.well_annot_df)

        # WRITING DATAFRAMES TO FILE
        self.RASLseqAnalysis_df.to_csv(self.write_path, sep="\t", index=False)  # Final analysis df
        self.read_df_counts.to_csv(self.write_path + '.read_df_counts', sep="\t")  # Probe reads counts by well df
        self.read_df.to_csv(self.write_path + '.read_df', sep="\t")  # All unique sequences with blast alignments df


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fastq', type=str, help='Specifies the input fastq file, e.g. /path/to/CG_data/RASLseq.fq')
# parser.add_argument('-f', '--fastq', type=str, help='Specifies the input fastq file can be series delimitted list of fq, e.g. /path/to/CG_data/RASLseq.fq')

parser.add_argument('-p', '--probes', type=str, help='Specifies the input probes file containing the following columns: AcceptorProbeSequence DonorProbeSequence AcceptorAdaptorSequence DonorAdaptorSequence ProbeName, e.g. /path/to/probes.txt')
parser.add_argument('-w', '--well_annot', type=str, help='Specifies the input well annotations file containing the following columns: PlateBarcode and WellBarcode, e.g. /path/to/well/annotations.txt')

# parser.add_argument('-b', '--blastn_bin', type=str, help='Specifies the path to directory holding blastn executable, /path/to/blastn_dir/')
# parser.add_argument('-b', '--blastn_bin', type=str, help='Specifies the path to directory holding blastn executable, /path/to/blastn_dir/')
parser.add_argument('-a', '--aligner_bin', type=str, help='Specifies the path to directory holding Aligner (STAR or blastn) executable, /path/to/aligner_dir/')

parser.add_argument('-o5', '--offset_5p', type=int, default=24, help='Specifies the number of bases to clip from 5-prime end of read to isolate probe sequence, default 24')
parser.add_argument('-o3', '--offset_3p', type=int, default=22, help='Specifies the number of bases to clip from 3-prime end of read to isolate probe sequence, default 22')
parser.add_argument('-ws', '--wellbc_start', type=int, default=0, help='Specifies the index position of the wellbc start base, default 0')
parser.add_argument('-we', '--wellbc_end', type=int, default=8, help='Specifies the index position of the wellbc start base, default 8')

parser.add_argument('-s', '--sequencer_id', type=str, help='Specifies the sequencer identifier in the fastq index lines, e.g. @HISEQ')

parser.add_argument('-d', '--blastdb_path', type=str, help='Specifies the directory for writing blast database, /path/to/write/blastdb/')
# parser.add_argument('-P', '--print_on', action='store_true', default=False, help='Specifies whether to print summary stats during alignment')
# parser.add_argument('-o', '--output_path', type=str, help='Specifies the output file path, e.g. /path/to/output/counts_df.txt')
# parser.add_argument('-d', '--output_dir', type=str, help='Specifies the output directory path, e.g. /path/to/output/')
parser.add_argument('-o', '--output_file', type=str, help='Specifies the output file path, e.g. /path/to/output/STAR_reads.txt')

parser.add_argument('-n', '--n_jobs', type=int, default=1, help='Specifies the number of processors to use, default 1')
parser.add_argument('-A', '--write_alignments', action='store_true', default=False, help='Specifies whether to write STAR alignments to disk, default=False')

opts = parser.parse_known_args()[0]

fastq_path, probes_path, well_annot_path, aligner_dir = opts.fastq, opts.probes, opts.well_annot, opts.aligner_bin
offset_5p, offset_3p, wellbc_start, wellbc_end = opts.offset_5p, opts.offset_3p, opts.wellbc_start, opts.wellbc_end
n_jobs, write_alignments = opts.n_jobs, opts.write_alignments

if __name__ == '__main__':

    rasl_analysis = RASLseqTools(
            fastq_path, probes_path, well_annot_path,

            aligner_dir,

            offset_5p=offset_5p, offset_3p=offset_3p,
            wellbc_start=wellbc_start, wellbc_end=wellbc_end,

            # sequce_id
            # blastdb_path

            n_job=n_jobs, write_alignments=write_alignments
            )

    '''
    rasl_analysis.get_target_counts_df()
    '''
    print 'Complete'



'''
BLAST Notes
    $ python /path/to/RASLseqAnalysis.py -s '@HISEQ' -f /path/to/your.fastq.gz -p /paht/to/RASL.probes -w /path/to/annotations.bc -d /path/to/blastdb/erite_dir/ -b /path/to/blast/ncbi-blast-2.2.26+/bin/ -P -o /path/to/output.txt

STAR Notes

if __name__ == '__main__':

    # Handles multiple fastq files in serial
    if ',' in fastq_path:
        fq_files = fastq_path.split(',')

        rasl_analysis = RASLseqAnalysis_STAR(
            fq_files[0], probes_path, aligner_dir,
            well_annot, write_path, write_file,
            print_on=False, n_jobs=n_jobs, offset_5p=offset_5p,
            offset_3p=offset_3p, wellbc_start=wellbc_start,
            wellbc_end=wellbc_end, write_alignments=write_alignments)

        rasl_analysis.get_target_counts_df()

        annot_cols = list(rasl_analysis.RASLseqBCannot_obj.well_annot_df.columns)

        master_df = count_df(rasl_analysis.RASLseqAnalysis_df.copy(), annot_cols)

        # Iterating through FASTQ files
        for fastq in fq_files[1:]:
            rasl_analysis = RASLseqAnalysis_STAR(
                fastq, probes_path, aligner_dir,
                well_annot, write_path, write_file,
                print_on=False, n_jobs=n_jobs, offset_5p=offset_5p,
                offset_3p=offset_3p, wellbc_start=wellbc_start,
                wellbc_end=wellbc_end, write_alignments=write_alignments)

            rasl_analysis.get_target_counts_df()

            master_df = master_df.add(
                count_df(rasl_analysis.RASLseqAnalysis_df, annot_cols),
                fill_value=0)

            print
            print 'Demultiplexing, Alignment, & Counting Complete:', fastq_path

        master_df = rasl_analysis.RASLseqBCannot_obj.well_annot_df.join(master_df)
        master_df.to_csv(write_path + 'Aggregated_counts_STAR_alignment.txt', sep='\t')

        print 'All Files Complete:'
        print fastq_path

        # GZip STAR alignment file
        os.system('gzip ' + rasl_analysis.alignment_write_file)

    else:
        # Handles single fastq file
        rasl_analysis = RASLseqAnalysis_STAR(
            fastq_path, probes_path, aligner_dir,
            well_annot, write_path, write_file,
            print_on=False, n_jobs=n_jobs, offset_5p=offset_5p,
            offset_3p=offset_3p, wellbc_start=wellbc_start,
            wellbc_end=wellbc_end, write_alignments=write_alignments)

        rasl_analysis.get_target_counts_df()

        rasl_analysis.RASLseqAnalysis_df.to_csv(rasl_analysis.write_file, sep='\t')

        # GZip STAR alignment file
        os.system('gzip ' + rasl_analysis.write_file)
        print
        print 'Demultiplexing, Alignment, & Counting Complete:', fastq_path

        (TEST)
        python RASLseqAnalysis_STAR.py -f ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/truth_sets/SeqRun1/lane1_Undetermined_L001_R1_001_truth_set_reads.fastq.gz -p ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/on_target_probes_Bcell_2014.probes -a ~/Dropbox/RASLseq/Bcell_exp2/STAR_bin/ -w ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/20131203_Rasl-Seq_bioactive_cmp-Table1.tsv -d ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/temp/ -o ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/temp/STAR_testing.txt -n 1 -o5 24 -o3 22 -ws 0 -we 8


'''

