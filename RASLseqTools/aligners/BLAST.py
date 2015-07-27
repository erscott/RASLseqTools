from RASLseqTools import *
from ..RASLseqAnalysis import RASLseqAnalysis

import numpy as np
import argparse
import inspect
import sys
import os


class BLAST(RASLseqAnalysis):
    name = 'blast'
    version = 0.1

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

    def __init__(self, fastq_path, sequencer_id,
            probes_path, blastdb_path, blastn_dir,
            well_annot, write_path, print_on=False,
            offset_5p=24, offset_3p=22, wellbc_start=0,
            wellbc_end=8):
        # Create a probe with name name
        pass

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


    def get_target_counts_df(self):
        '''
            This method uses a functional programming design to
            produce a pandas DataFrame describing the RASLseq probe-specific
            read counts prepended with well specific annotations.
        '''

        # Add observed wellbc
        # self.read_df = RASLseqWellbc.get_rasl_probe_and_wellbc_exact_df(self.read_df, print_on=self.print_on)
        self.read_df['observed_wellbc'] = self.read_df['seq'].str[self.wellbc_start:self.wellbc_end]

        # Add fuzzy matched wellbc (mapped_bc)
        self.read_df = RASLseqWellbc.get_rasl_probe_and_wellbc_fuzzy_df(self.read_df,
            self.RASLseqBCannot_obj.well_bc,
            print_on=self.print_on)

        # Filter fuzzy wellbc mappings
        self.read_df = self.read_df[(self.read_df.bc_edit_dist.astype(float) < self.bc_edit_dist_filter)]

        # Add observed rasl_probe_seq
        # self.read_df = RASLseqSeq.get_rasl_probe_exact_df(self.read_df)
        self.read_df['rasl_probe'] = self.read_df['seq'].str[self.offset_5p:-self.offset_3p]

        # Blast alignment
        self.read_df = RASLseqAlign.get_rasl_aligned_df(self.read_df,
            self.RASLseqProbes_obj.aligner_dir,
            self.RASLseqProbes_obj.probedb_file,
            print_on=self.print_on,
            aligner='blast')

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
        # Rreturn self.read_df_counts, self.RASLseqBCannot_obj.well_annot_df
        self.RASLseqAnalysis_df = self._merge_plate_well_annot(self.read_df_counts, self.RASLseqBCannot_obj.well_annot_df)

        self.RASLseqAnalysis_df.to_csv(self.write_path, sep="\t")
        self.read_df_counts.to_csv(self.write_path + '.read_df_counts', sep="\t")
        self.read_df.to_csv(self.write_path + '.read_df', sep="\t")

        os.system('gzip ' + self.write_path + '.read_df')

        return

