from ..RASLseqAnalysis import RASLseqAnalysis

import multiprocessing as mp
import pandas as pd
import numpy as np
import argparse
import inspect
import raven
import sys
import os


class STAR(RASLseqAnalysis):
    name = 'star'
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

        aligner_dir: str
        path to directory holding the blastn or STAR executable

        well_annot: str
        path to Barcode Annotations file

        write_path: str
        path to write directory

        offset_5p: int, optional, default=24
        Index position for RASLprobe start,
        Number of bases to clip from 5' end of read
        to isolate probe sequence

        offset_3p: int, optional, default=22
        Number of bases from 3' end of read to clip
        in order to isolate probe sequence in read

        wellbc_start: int, default=0
        Expected wellbarcode start position

        wellbc_end: int, default=8
        Expected wellbarcode end position


        Returns
        -----------
        read_df: pandas DataFrame
        index: (PlateBarcode, sequence)
        columns: ['PlateBarcode', 'seq', 'seq_count']
        PlateBarcode - Index read from fastq file
        seq - fastq read sequence
        seq_count - number of occurrences of seq in fastq
    '''

    def __init__(self, fastq_path, probes_path,
            aligner_dir, well_annot, write_path,
            write_file, print_on=False, n_jobs=1,
            offset_5p=24, offset_3p=22, wellbc_start=0,
            wellbc_end=8, write_alignments=False):
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
        collapsed_read_counts['count'] = 1
        collapsed_probe_counts = collapsed_read_counts.groupby(['PlateBarcode', 'WellBarcode', 'ProbeName'])['count'].aggregate(np.sum)

        # Aggregating probe counts for each well
        counts_df = collapsed_probe_counts.unstack('ProbeName')  # Creating matrix of aggregated probe counts indexed on 'plate_barcode','mapped_bc'

        # ID and removal of off-target ligation hits found by ! character in blast sseq name
        on_target_col = [i for i in counts_df.columns if "!" not in i]
        counts_df = counts_df[on_target_col]  # Removing off-target ligation hits from df
        counts_df.index.names = ['PlateBarcode', 'WellBarcode']
        return counts_df


    def get_demultiplexed_id(self):
        '''
            This function maps wellbarcodes and platebarcodes
        '''

        # Parsing barcodes from fastq data
        self.read_df['wellbarcode'] = self.read_df['seq'].str[self.wellbc_start:self.wellbc_end]  # Parsing wellbc from read seq
        self.read_df['platebarcode'] = self.read_df['id_line'].str.split(' ').str[-1].str.split(':').str[-1].str.rstrip('N')  # Parsing platebc from id_line

        self.read_df = self.read_df[['id_line', 'wellbarcode', 'platebarcode']]  # Removing unnecessary data

        # Demultiplexing platebarcode
        imperfect_platebcs = self.read_df[~self.read_df['platebarcode'].isin(self.RASLseqBCannot_obj.plate_bc)][['platebarcode']].drop_duplicates()
        imperfect_platebcs['PlateBarcode'] = imperfect_platebcs['platebarcode'].apply(RASLseqWellbc.map_bc, args=[self.RASLseqBCannot_obj.plate_bc])
        perfect_platebcs = pd.DataFrame([[i, i] for i in self.RASLseqBCannot_obj.plate_bc], columns=['PlateBarcode', 'platebarcode'])
        imperfect_platebcs = imperfect_platebcs.append(perfect_platebcs)

        # Demultiplexing wellbarcode
        imperfect_wellbcs = self.read_df[~self.read_df['wellbarcode'].isin(self.RASLseqBCannot_obj.well_bc)][['wellbarcode']].drop_duplicates()
        imperfect_wellbcs['WellBarcode'] = imperfect_wellbcs['wellbarcode'].apply(RASLseqWellbc.map_bc, args=[self.RASLseqBCannot_obj.well_bc])
        perfect_wellbcs = pd.DataFrame([[i, i] for i in self.RASLseqBCannot_obj.well_bc], columns=['WellBarcode', 'wellbarcode'])
        imperfect_wellbcs = imperfect_wellbcs.append(perfect_wellbcs)

        # Combining demultiplexed platebarcodes with read data
        self.read_df = self.read_df.merge(imperfect_platebcs, on='platebarcode', how='left')
        del self.read_df['platebarcode']

        # Combining demultiplexed wellbarcodes with read data
        self.read_df = self.read_df.merge(imperfect_wellbcs, on='wellbarcode', how='left')
        del self.read_df['wellbarcode']

        return True

    def get_target_counts_df(self):
        '''
            This function demultiplexes and aligns RASL-seq fastq reads

            Expects the PlateBarcode in the read id-line
            PlateBarcode as the last element in the read id-line, delimitted by " "

            Expects the following read structure:
            WellBarcode - Adaptor - RASL probe RASL probe - Adaptor

            Processes 250000 reads per chunk
        '''

        self.fastq_df_chunk = pd.read_table(self.fastq_path, sep='\t', chunksize=1000000, header=None, compression='gzip' if '.gz' in self.fastq_path else None)

        self.master_df = []
        for n, read_df in enumerate(self.fastq_df_chunk):

            # READ FASTQ CHUNK INTO DF
            self.read_df = FastqCollapse.get_fastq_unstack(read_df)

            # WRITE TO TEMP FILE
            random_str = self.RASLseqProbes_obj.random_str
            self.temp_fastq = self.write_path + 'temp_' + random_str + '.fastq'
            FastqCollapse.write_fastq(self.read_df, self.temp_fastq)

            print
            print 'Starting Asynchronous Alignment Process'
            pool = mp.Pool(1)
            task = (self.temp_fastq, self.aligner_path.rstrip('/'), self.RASLseqProbes_obj.probedb_path,
                    self.print_on, 'star', 1,
                    self.offset_5p, self.offset_3p)
            self.aligned_df = pool.apply_async(RASLseqAlign.get_rasl_aligned_df, task)

            # BARCODE DEMULTIPLEXING
            print 'Demultiplexing Reads'
            self.get_demultiplexed_id()
            print 'Demultiplexing Complete'

            # FORMATTING ID_LINE FOR MERGING
            self.read_df['id_line'] = self.read_df['id_line'].str.split(' ').str[0].str.lstrip('@')

            pool.close()
            pool.join()
            self.aligned_df = self.aligned_df.get()
            pool.terminate()
            print 'Alignment Complete'

            if self.write_alignments:
                print 'Writing STAR Alignments To Disk'
                self.alignment_write_file = self.write_file + '.Aligned.out'
                self.alignment_write_file_header = self.aligned_df.columns
                self.aligned_df.to_csv(self.alignment_write_file, sep='\t', mode='a', header=False, index=False)
                print 'Writing STAR Alignments To Disk Complete'
                print

            # MERGING BARCODES AND ALIGNED READS
            self.read_df = self.read_df.merge(self.aligned_df, on='id_line', how='inner')
            del self.aligned_df

            # REMOVING MISMATCH WELLBARCODE AND PLATEBARCODE READS
            self.read_df = self.read_df[(self.read_df.WellBarcode != 'mismatch') & (self.read_df.PlateBarcode != 'mismatch')]

            # SUM PROBE-SPECIFIC READ COUNTS
            self.read_df = self._get_probe_well_read_counts(self.read_df)

            self.master_df.append(self.read_df)
            del self.read_df

            print (n + 1) * 250000, 'Reads Processed'

            os.system('rm ' + self.temp_fastq)

        # AGGREGATING THE PROBE COUNTS FOR EACH SAMPLE
        self.RASLseqAnalysis_df = self.master_df[0]
        for df in self.master_df[1:]:
            self.RASLseqAnalysis_df = self.RASLseqAnalysis_df.add(df, fill_value=0)
        del self.master_df

        self.probe_columns = self.RASLseqAnalysis_df.columns
        self.read_count_mapped = self.RASLseqAnalysis_df.sum().sum()

        # MERGING WELL ANNOTATIONS AND PROBE READ COUNTS
        self.annotation_columns = self.RASLseqBCannot_obj.well_annot_df.columns
        self.RASLseqAnalysis_df = self._merge_plate_well_annot(self.RASLseqAnalysis_df, self.RASLseqBCannot_obj.well_annot_df)

        # DELETING TEMP FILES
        os.system('rm -r ' + self.RASLseqProbes_obj.probedb_path)

        return

    def load_target_counts_df(self, file_path):
        '''
            This function loads a RASLseqAnalysis_df from the file_path

            Parameters
            ----------------
            file_path: str
                Specifies the path to the RASLseqAnalysis_df

                Expects tab-separated sample by count dataframe with
                PlateBarcode, WellBarcode columns

            Returns
            ----------------
            self.RASLseqAnalysis_df, pandas DataFrame

        '''

        compression = ''
        if 'gz' in file_path:
            compression = 'gzip'

        df = pd.read_table(file_path, sep='\t', index_col=['PlateBarcode', 'WellBarcode'], compression=compression)
        self.RASLseqAnalysis_df = df

        self.probe_columns = list(set(self.RASLseqProbes_obj.probe_columns) & set(self.RASLseqAnalysis_df.columns))

        self.annot_columns = list(set(self.RASLseqBCannot_obj.annot_columns) & set(self.RASLseqAnalysis_df.columns))

        return


def count_df(df, annot_cols):
    '''
        This function aggregates dataframes

        Parameters
        --------------
        df: pandas dataframe

        count_cols: list of probe count columns


        Returns
        --------------
        pandas df with aggregated counts

    '''

    def get_count_cols(df, annot_cols):
        return list(set(df.columns) - set(annot_cols))

    aggregated_counts = df[get_count_cols(df, annot_cols)]
    aggregated_counts.fillna(value=0, inplace=True)

    return aggregated_counts


