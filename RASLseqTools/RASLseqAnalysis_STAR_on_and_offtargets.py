import pandas as pd
import os,sys
import numpy as np
import argparse
import subprocess
import multiprocessing as mp
import sys,os
import time

import inspect
source_dir =  os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
source_dir = '/'.join(source_dir.split('/')[:-1])
sys.path.append(source_dir)
from RASLseqTools import *



class RASLseqAnalysis_STAR(object):
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
        
        print_on: boolean, default=False
        Whether to print information during data processing
        
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
    
    def __init__(self, fastq_path, probes_path, aligner_dir, well_annot, \
                 write_path, write_file, print_on=False, n_jobs=1, offset_5p=24, \
                 offset_3p=22, wellbc_start=0, wellbc_end=8, write_alignments=False ):
        
        
        
        
        
        self.RASLseqProbes_obj = RASLseqProbes.RASLseqProbes(probes_path, write_path, aligner_dir, aligner='star')
        
        self.RASLseqBCannot_obj = RASLseqBCannot.RASLseqBCannot(well_annot)
        
        self.RASLseqBCannot_obj.well_bc = self.RASLseqBCannot_obj.well_bc
        
        self.aligner = 'star'
        
        self.aligner_path = aligner_dir
        
        self.fastq_path = fastq_path
        
        self.n_jobs = n_jobs
        
        self.offset_5p = int(offset_5p)
        
        self.offset_3p = int(offset_3p)
        
        self.wellbc_start = int(wellbc_start)
        
        self.wellbc_end = int(wellbc_end)
        
        self.print_on = print_on
        
        self.write_alignments = write_alignments
        
        self.write_path = write_path
        
        self.write_file = write_file
        
        self.bc_edit_dist_filter = 2
    
    
    
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
        collapsed_read_counts['count'] = 1
        collapsed_probe_counts = collapsed_read_counts.groupby(['PlateBarcode','WellBarcode','ProbeName'])['count'].aggregate(np.sum)
        
        #aggregating probe counts for each well
        counts_df = collapsed_probe_counts.unstack('ProbeName')  #creating matrix of aggregated probe counts indexed on 'plate_barcode','mapped_bc'
        
        #id and removal of off-target ligation hits found by ! character in blast sseq name
        on_target_col = [i for i in counts_df.columns if "!" not in i]
        counts_df = counts_df[on_target_col]  #removing off-target ligation hits from df
        counts_df.index.names= ['PlateBarcode', 'WellBarcode']

        offtarget_cols = [i for i in counts_df.columns if "!" in i]
        off_taget_counts_df = counts_df[offtarget_cols] 
        off_taget_counts_df.index.names = [‘PlateBarcode’, ‘WellBarcode’]
        

        return counts_df, off_taget_counts_df
    
    
    def _merge_plate_well_annot(self, probe_counts_df, well_annot_df):
        '''
            This function merges gene_counts_df with well annotations
            
            Parameters
            ----------
            probe_counts_df: Pandas DataFrame
            Requires pandas index: ('PlateBarcode','WellBarcode')
            
            well_annot_path: Pandas DataFrame
            Requires pandas index: ('PlateBarcode','WellBarcode')
            
            Returns
            -------
            Pandas DataFrame
            well_annot_df right joined to gene_counts_df
            index: ('plate_barcode','WellBarcode')
            
            '''
        return well_annot_df.join(probe_counts_df,how='right')
    
    
    
    def get_demultiplexed_id(self):
        '''
        This function maps wellbarcodes and platebarcodes
        '''
        
        #parsing barcodes from fastq data
        self.read_df['wellbarcode'] = self.read_df['seq'].str[self.wellbc_start:self.wellbc_end]  #parsing wellbc from read seq
        self.read_df['platebarcode'] = self.read_df['id_line'].str.split(' ').str[-1].str.split(':').str[-1].str.rstrip('N')  #parsing platebc from id_line
        
        self.read_df = self.read_df[['id_line', 'wellbarcode', 'platebarcode']]  #removing unnecessary data
        
        
        #demultiplexing platebarcode
        imperfect_platebcs = self.read_df[~self.read_df['platebarcode'].isin(self.RASLseqBCannot_obj.plate_bc)][['platebarcode']].drop_duplicates()
        imperfect_platebcs['PlateBarcode'] = imperfect_platebcs['platebarcode'].apply(RASLseqWellbc.map_bc, args=[self.RASLseqBCannot_obj.plate_bc])
        perfect_platebcs = pd.DataFrame([[i,i] for i in self.RASLseqBCannot_obj.plate_bc], columns=['PlateBarcode', 'platebarcode'])
        imperfect_platebcs = imperfect_platebcs.append(perfect_platebcs)
        
        #demultiplexing wellbarcode
        imperfect_wellbcs = self.read_df[~self.read_df['wellbarcode'].isin(self.RASLseqBCannot_obj.well_bc)][['wellbarcode']].drop_duplicates()
        imperfect_wellbcs['WellBarcode'] = imperfect_wellbcs['wellbarcode'].apply(RASLseqWellbc.map_bc, args=[self.RASLseqBCannot_obj.well_bc])
        perfect_wellbcs = pd.DataFrame([[i,i] for i in self.RASLseqBCannot_obj.well_bc], columns=['WellBarcode', 'wellbarcode'])
        imperfect_wellbcs = imperfect_wellbcs.append(perfect_wellbcs)
        
        #combining demultiplexed platebarcodes with read data
        self.read_df = self.read_df.merge(imperfect_platebcs, on='platebarcode', how='left')
        del self.read_df['platebarcode']
        
        #combining demultiplexed wellbarcodes with read data
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
        
        
        if '.gz' in self.fastq_path:
            self.fastq_df_chunk = pd.read_table(self.fastq_path, sep='\t', chunksize=1000000, header=None, compression='gzip')
        
        else:
            self.fastq_df_chunk = pd.read_table(self.fastq_path, sep='\t', chunksize=1000000, header=None)
        
        
        self.master_df = []
        self.master_offtarget_df = []
        for n,read_df in enumerate(self.fastq_df_chunk):
            
            #READ FASTQ CHUNK INTO DF
            self.read_df = FastqCollapse.get_fastq_unstack(read_df)
            
            
            #WRITE TO TEMP FILE
            random_str = self.RASLseqProbes_obj.random_str
            self.temp_fastq = self.write_path + 'temp_' + random_str + '.fastq'
            FastqCollapse.write_fastq(self.read_df, self.temp_fastq)
            
            
            #LAUNCH STAR ALIGNMENT PROCESS, STAR ALIGN STDOUT
            #             self.aligned_df = RASLseqAlign.get_rasl_aligned_df(self.temp_fastq, aligner_path.rstrip('/'), \
            #                              self.RASLseqProbes_obj.probedb_path, print_on=self.print_on, aligner='star', \
            #                              n_jobs=1, offset_5p=self.offset_5p, offset_3p=self.offset_3p)
            
            #ASYNCHRONOUS STAR ALIGNMENT
            #print time.gmtime(), 1
            print
            print 'Starting Asynchronous Alignment Process'
            pool = mp.Pool(1)
            task = (self.temp_fastq, self.aligner_path.rstrip('/'), self.RASLseqProbes_obj.probedb_path, \
                    self.print_on, 'star', 1, self.offset_5p, self.offset_3p)
            self.aligned_df = pool.apply_async(RASLseqAlign.get_rasl_aligned_df, task)
            
            
            #BARCODE DEMULTIPLEXING
            print 'Demultiplexing Reads'
            #print time.gmtime(), 2
            self.get_demultiplexed_id()
            #print time.gmtime(), 3
            print 'Demultiplexing Complete'
                    
                    
            #FORMATTING ID_LINE FOR MERGING
            self.read_df['id_line'] = self.read_df['id_line'].str.split(' ').str[0].str.lstrip('@')
                    
                    
            #print time.gmtime(), 4
            pool.close()
            pool.join()
            self.aligned_df = self.aligned_df.get()
            pool.terminate()
            print 'Alignment Complete'
            #print time.gmtime(), 5
            
            if self.write_alignments:
                print 'Writing STAR Alignments To Disk'
                self.alignment_write_file = self.write_file + '.Aligned.out'
                self.alignment_write_file_header = self.aligned_df.columns
                self.aligned_df.to_csv(self.alignment_write_file, sep='\t', mode='a', header=False, index=False)
                print 'Writing STAR Alignments To Disk Complete'  
                print   
            
            #MERGING BARCODES AND ALIGNED READS
            
            
            self.read_df = self.read_df.merge(self.aligned_df, on='id_line', how='inner')
            del self.aligned_df
                    
                    
            #REMOVING MISMATCH WELLBARCODE AND PLATEBARCODE READS
            self.read_df = self.read_df[(self.read_df.WellBarcode != 'mismatch') & (self.read_df.PlateBarcode != 'mismatch')]
                    
                    
            #SUM PROBE-SPECIFIC READ COUNTS
            self.read_df, self.offtarget_df = self._get_probe_well_read_counts( self.read_df )
                    
                    
            self.master_df.append(self.read_df)
            self.master_offtarget_df.append(self.offtarget_df)
            del self.read_df
            del self.offtarget_df
            
            print (n+1)*250000, 'Reads Processed'
            
            os.system('rm ' + self.temp_fastq)
            
        
        
        #AGGREGATING THE PROBE COUNTS FOR EACH SAMPLE
        self.RASLseqAnalysis_df = self.master_df[0]
        self.RASLseqAnalysis_offtarget_df = self.master_offtarget_df[0]

        for df in self.master_df[1:]:
            self.RASLseqAnalysis_df = self.RASLseqAnalysis_df.add(df, fill_value=0)
        del self.master_df

        for df in self.master_offtarget_df[1:]:
            self.RASLseqAnalysis_offtarget_df = self.RASLseqAnalysis_offtarget_df.add(df, fill_value=0)
        del self.master_offtarget_df
    
        self.probe_columns = self.RASLseqAnalysis_df.columns
        self.read_count_mapped = self.RASLseqAnalysis_df.sum().sum()

        self.offtarget_probe_columns = self.RASLseqAnalysis_offtarget_df.columns
        self.read_count_offtarget_mapped = self.RASLseqAnalysis_offtarget_df.sum().sum()
        
        
        #MERGING WELL ANNOTATIONS AND PROBE READ COUNTS
        self.annotation_columns = self.RASLseqBCannot_obj.well_annot_df.columns
        self.RASLseqAnalysis_df = self._merge_plate_well_annot(self.RASLseqAnalysis_df, self.RASLseqBCannot_obj.well_annot_df)

        
        #DELETING TEMP FILES
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
            compression='gzip'
        
        df = pd.read_table(file_path, sep='\t', index_col=['PlateBarcode', 'WellBarcode'], compression=compression)
        self.RASLseqAnalysis_df = df
        
        self.probe_columns = list( set(self.RASLseqProbes_obj.probe_columns) & set(self.RASLseqAnalysis_df.columns) )
        
        self.annot_columns = list( set(self.RASLseqBCannot_obj.annot_columns) & set(self.RASLseqAnalysis_df.columns) )
        
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
        return list( set(df.columns) - set(annot_cols) )
    
    aggregated_counts = df[get_count_cols(df, annot_cols)]
    aggregated_counts.fillna(value=0, inplace=True)
    
    return aggregated_counts








parser = argparse.ArgumentParser()

parser.add_argument('-f','--fastq', type=str,help='Specifies the input fastq file can be series delimitted list of fq, e.g. /path/to/CG_data/RASLseq.fq')

parser.add_argument('-p','--probes', type=str,help='Specifies the input probes file containing the following columns: AcceptorProbeSequence DonorProbeSequence AcceptorAdaptorSequence DonorAdaptorSequence ProbeName, e.g. /path/to/probes.txt')

parser.add_argument('-a','--aligner_bin', type=str,help='Specifies the path to directory holding STAR executable, /path/to/aligner_dir/')

parser.add_argument('-w','--well_annot', type=str,help='Specifies the input well annotations file containing the following columns: PlateBarcode and WellBarcode, e.g. /path/to/well/annotations.txt')

parser.add_argument('-d','--output_dir', type=str,help='Specifies the output directory path, e.g. /path/to/output/')

parser.add_argument('-o','--output_file', type=str,help='Specifies the output file path, e.g. /path/to/output/STAR_reads.txt')

parser.add_argument('-P','--print_on', action='store_true', default=False,help='Specifies whether to print summary stats during alignment, default=False')

parser.add_argument('-A','--write_alignments', action='store_true', default=False, help='Specifies whether to write STAR alignments to disk, default=False')

parser.add_argument('-n','--n_jobs', type=int, default=1, help='Specifies the number of processors to use, default 1')

parser.add_argument('-o5','--offset_5p', type=int, default=24, help='Specifies the number of bases to clip from 5-prime end of read to isolate probe sequence, default 24')

parser.add_argument('-o3','--offset_3p', type=int, default=22, help='Specifies the number of bases to clip from 3-prime end of read to isolate probe sequence, default 22')

parser.add_argument('-ws','--wellbc_start', type=int, default=0, help='Specifies the index position of the wellbc start base, default 0')

parser.add_argument('-we','--wellbc_end', type=int, default=8, help='Specifies the index position of the wellbc start base, default 8')

opts = parser.parse_known_args()

fastq_path, probes_path, aligner_dir, well_annot, write_path  = opts[0].fastq, opts[0].probes, opts[0].aligner_bin, opts[0].well_annot, opts[0].output_dir

n_jobs, offset_5p, offset_3p, wellbc_start, wellbc_end, write_file = opts[0].n_jobs, opts[0].offset_5p, opts[0].offset_3p, opts[0].wellbc_start, opts[0].wellbc_end, opts[0].output_file

print_on, write_alignments = opts[0].print_on, opts[0].write_alignments





if __name__ == '__main__':
    
    if ',' in fastq_path:  #Handles multiple fastq files in serial
        fq_files = fastq_path.split(',')
        rasl_analysis = RASLseqAnalysis_STAR(fq_files[0], probes_path, aligner_dir, well_annot, \
                                                         write_path, write_file, print_on=False,n_jobs=n_jobs, \
                                                         offset_5p=offset_5p, offset_3p=offset_3p, \
                                                         wellbc_start=wellbc_start, wellbc_end=wellbc_end, \
                                                         write_alignments=write_alignments)
        rasl_analysis.get_target_counts_df()
        
        annot_cols = list(rasl_analysis.RASLseqBCannot_obj.well_annot_df.columns)
        
        master_df = count_df(rasl_analysis.RASLseqAnalysis_df.copy(), annot_cols)
        master_offtarget_df = rasl_analysis.RASLseqAnalysis_offtarget_df.copy()
        for fastq in fq_files[1:]:  #IRERATING THROUGH FASTQ FILES
            rasl_analysis = RASLseqAnalysis_STAR(fastq, probes_path, aligner_dir, well_annot, write_path, \
                                                 write_file, print_on=False,n_jobs=n_jobs, offset_5p=offset_5p, \
                                                 offset_3p=offset_3p, wellbc_start=wellbc_start, \
                                                 wellbc_end=wellbc_end, write_alignments=write_alignments)
         
            rasl_analysis.get_target_counts_df()
            
            master_df = master_df.add( count_df(rasl_analysis.RASLseqAnalysis_df, annot_cols), fill_value=0 )
            offtarget_master_df = offtarget_master_df.add(rasl_analysis.RASLseqAnalysis_offtarget_df, fill_value=0)
            
            print 
            print 'Demultiplexing, Alignment, & Counting Complete:', fastq_path 
             
        master_df = rasl_analysis.RASLseqBCannot_obj.well_annot_df.join(master_df)
        master_df.to_csv(write_path + 'Aggregated_counts_STAR_alignment.txt', sep='\t')
        offtarget_master_df.to_csv(write_path + 'Aggregated_off_target_counts_STAR_alignment.txt', sep='\t')
            
        print 'All Files Complete:'
        print fastq_path
        
        os.system('gzip ' + rasl_analysis.alignment_write_file)  #gzip STAR alignment file
         
    else:  #Handles single fastq file
        rasl_analysis = RASLseqAnalysis_STAR(fastq_path, probes_path, aligner_dir, well_annot, write_path, \
                                             write_file, print_on=False,n_jobs=n_jobs, offset_5p=offset_5p, \
                                             offset_3p=offset_3p, wellbc_start=wellbc_start, \
                                             wellbc_end=wellbc_end, write_alignments=write_alignments)
     
        rasl_analysis.get_target_counts_df()
        
        rasl_analysis.RASLseqAnalysis_df.to_csv(rasl_analysis.write_file, sep='\t')

        rasl_analysis.RASLseqAnalysis_offtarget_df.to_csv(rasl_analysis.write_file + '.offtarget', sep='\t')
        
        os.system('gzip ' + rasl_analysis.write_file)  #gzip STAR alignment file
        print 
        print 'Demultiplexing, Alignment, & Counting Complete:', fastq_path 
     



"""
TEST
python RASLseqAnalysis_STAR.py -f ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/truth_sets/SeqRun1/lane1_Undetermined_L001_R1_001_truth_set_reads.fastq.gz -p ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/on_target_probes_Bcell_2014.probes -a ~/Dropbox/RASLseq/Bcell_exp2/STAR_bin/ -w ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/20131203_Rasl-Seq_bioactive_cmp-Table1.tsv -d ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/temp/ -o ~/Dropbox/RASLseq/Bcell_exp2/ipynb/data/temp/STAR_testing.txt -n 1 -o5 24 -o3 22 -ws 0 -we 8
"""

