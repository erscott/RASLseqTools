

# Libraries


import pandas as pd
import editdist
import mmap
import os
import sys
import time
import itertools
import gzip
import numpy as np
from collections import Counter
import random



# ALIGNER BLAST


def rasl_probe_blast(read_file_path, blastn_path, db_path):
    """
    This function returns blast results from queries optimized for rasl-seq conditions and outputs
    a custom format
    
    Parameters
    ----------
    read_file_path: str, path to temp blast input file
        Specifies where to write temporary blast input file
        
    blastn_path: str, path to blastn executable
        Specifies full path to blastn executable
    
    db_path: str, path to target BLAST database
        Specifies path to on_off_target BLAST database
    
    Returns
    -------
    Pandas Blast Results dataframe
        index = qseqid, rasl_probe sequence
        cols = ['qseqid','mismatch','evalue','sseqid','qlen', 'length','qseq','qstart','sseq','sstart','send']
    """
    
    #SETTING WRITE PATH FOR BLAST OUTPUT
    write_path = read_file_path.rstrip(".txt") + "_blast_output.txt"
    
    
    db = " -db " + db_path #/gpfs/home/erscott/Datasets_raw/UCSC_tables/commonSnp135_masked_hg19/commonSnp135_masked_hg19_transcript_coding_seq" # -db /gpfs/home/erscott/Datasets_raw/UCSC_tables/masked_hg19_snp135/masked_hg19_transcript_coding_seq"
    
    #SETTING BLAST WORD SIZE
    wordsize = "-word_size 8 "
    
    #PRINTING BLAST COMMAND
    print "/gpfs/home/erscott/Tools/blast/ncbi-blast-2.2.26+/bin/blastn -task blastn-short -query "+ read_file_path +" -evalue 20 "+ wordsize + db +" -outfmt 5"

    #blast_run = subprocess.Popen("/gpfs/home/erscott/Tools/blast/ncbi-blast-2.2.26+/bin/blastn -task megablast -query "+ read_file_path +" -evalue 20 "+ wordsize + db +" > new_seq_blast_results.txt",shell=True,stdout=subprocess.PIPE).stdout
    blast_run = os.system(blastn_path + " -task blastn-short -query " \
                          + read_file_path +" -evalue 1e-6 "+ wordsize + db +" -max_target_seqs 1" +" -strand 'plus' -xdrop_gap 7" \
                          + " -outfmt '6 -outfmt qseqid -outfmt mismatch -outfmt evalue -outfmt sseqid -outfmt qlen -outfmt length -outfmt qseq -outfmt qstart -outfmt sseq -outfmt sstart -outfmt send' >"+ write_path)
    # Evlaue of 20 was used to account for a cutoff of 13mer matching
    
    #CONVERTING BLAST OUTPUT INTO DATAFRAME
    bl_results = pd.read_table(write_path, header=None)
    bl_results.columns=['qseqid','mismatch','evalue','sseqid','qlen', 'length','qseq','qstart','sseq','sstart','send']
    bl_results['probe'] = bl_results.sseqid.apply(lambda x: "_".join(x.split("_")))
    bl_results.ix[0]
    bl_results.set_index('qseqid',drop=False,inplace=True)
    
    #REMOVING INPUT/OUTPUT BLAST FILES
    os.system('rm ' + read_file_path)
    os.system('rm '+ write_path)
    
    return bl_results

# <headingcell level=2>

# COLLAPSING DICTIONARY

# <codecell>


def fastq_collapsing_Counter(fastq_file_obj):
    """
    This function creates a dictionary of FASTQ reads sequences (plus 3' barcode) by collapsing identical reads
    
    Parameters
    ----------
    fastq_file_obj: unread python file object
    
    
    Returns
    -------
    Counter object (python Collections)
        key: (FASTQ read sequence, bc_3')  (bc_3' is found in the header sequence)
        value: int, count of FASTQ read sequence
    """
    seq_library_counter = Counter()
    header_bc_2 = ""
    next_line_seq = 0
    for fq_read in fastq_file_obj:
        if next_line_seq == 1:  #occurs on second line of sequence
            seq_library_counter.update([(header_bc_2, fq_read.rstrip())])
            next_line_seq = 0
        if fq_read.startswith("@HISEQ"):  #save barcode_3' read
            header_bc_2 = fq_read.rstrip().split(":")[-1]  #grab 3' barcode
            next_line_seq = 1
        
    return seq_library_counter



# ADAPTOR AND BARCODE MATCHING

# <codecell>

def rasl_probe_seq_extraction(line, AD1='GGAAGCCTTGGCTTTTG', AD2='AGATCGGAAGAGCACAC'):
    ''' 
    This function returns the ligation sequence between the two adaptors
    if there is a perfect string match with adaptor
        
    Parameters
    ----------
    line: Pandas dataframe line (Series)
        requires column: 'seq' possessing fastq read sequence
        
    AD1: str, optional, default='GGAAGCCTTGGCTTTTG'
        Specifies the first adaptor sequence expected in the fastq read
    
    AD2: str, optional, default='AGATCGGAAGAGCACAC'
        Specifies the second adaptor sequence expected in the fastq read
    
    Returns
    -------
    Fastq nucleotide sequence between AD1 and AD2 
      or
    '0' if AD1 or AD2 is not in seq sequence
    '''
    
    seq = line['seq']
    if AD1 in seq and AD2 in seq:
        return seq[seq.index(AD1)+17: seq.index(AD2)]
    else: return '0'


    
def observed_wellbc_extraction(line, AD1='GGAAGCCTTGGCTTTTG'):
    ''' 
    This function returns the sequence (observed wellbc) before the first adaptor sequence 
    ToDo, fix alignments to well bc
    
    Parameters
    ----------
    line: Pandas dataframe line (Series)
        requires column: 'seq' possessing fastq read sequence
        
    AD1: str, optional, default='GGAAGCCTTGGCTTTTG'
        Specifies the first adaptor sequence expected in the fastq read
        
    Returns
    -------
    Fastq sequence before first occurent of AD1
      or
    '0' if AD1 not in seq sequence
    '''
    
    seq = line['seq']
    if AD1 in seq:
        return seq[ : seq.index(AD1)]
        
    else: return '0'



# <codecell>

def get_file_obj(file_path):
    '''
    This function returns a file object
    given a gzip or uncompressed file path
    '''
    if file_path.endswith('.gz'):
        return gzip.open(file_path,'r')
    else:
        return open(file_path, 'r')
    

# <headingcell level=2>

# ALIGNER PIPELINE START

# <codecell>

def collapsed_read_df(fastq_path, print_on=False):
    '''
    This function collapses identical fastq read sequences using a Counter
    object and returns a pandas dataframe:
        index: (plate_barcode, sequence)
        columns: ['plate_barcode', 'seq', 'seq_count']
    
    Parameters
    ----------
    fastq_path: str
        Specifies the full path to the fastq file, handles gzip and non-gzip fastq files
    
    print_on: boolean, default=False
        If print_on == True prints "number of unique bc/reads to map total"
    
    Returns
    -------
    Pandas dataframe
        multi-index: (plate_barcode, sequence)
        columns: ['plate_barcode', 'seq', 'seq_count']
    
    plate_barcode - Index read from fastq file
    seq - fastq read sequence
    seq_count - number of occurrences of seq in fastq
    '''
    
    #GET FILE OBJECT
    fastq_obj = get_file_obj(fastq_path)
    
    #CREATE COLLAPSED DICTIONARY OF READS
    seq_lib = fastq_collapsing_Counter(fastq_obj)  #~4-fold compression
    
    
    #CREATING seq_count DATAFRAME THAT CAN BE ACCESSED BY MULTI-INDEX
    df = pd.DataFrame.from_dict(seq_lib,orient='index')  #index = (barcode, sequence) ; value = seq_count from collapsing dictionary
    df.columns = ['seq_count']
    
    
    #CREATING plate_barcode AND seq DF
    df_plate_seq = pd.DataFrame.from_records(seq_lib.keys())
    df_plate_seq.columns = ['plate_barcode','seq']
    df_plate_seq.set_index(['plate_barcode','seq'], drop=False, inplace=True)
    df_plate_seq = df_plate_seq.join(df)  #merged dataframe; index=(plate_barcode, sequence); columns = ['plate_barcode', 'seq', 'seq_count'] 

    
    if print_on:
        print len(seq_lib), "number of unique bc/reads to map total"
        
    
    return df_plate_seq
    

# <codecell>

def get_rasl_probe_and_wellbc_exact(collapsed_read_df, AD1='GGAAGCCTTGGCTTTTG', AD2='AGATCGGAAGAGCACAC', perfect_matches=True, print_on=False):   
    '''
    This function identifies the rasl_probe (~40mer) sequence and
    well barcode sequence within the fastq read using an exact
    match to the adaptor sequences (AD1 & RCAD2)
    
    Parameters
    ----------
    collapsed_read_df: pandas dataframe
        Requires 'seq' column (fastq read sequence)
    
    AD1: str, optional, default='GGAAGCCTTGGCTTTTG'
        Specifies the first adaptor sequence expected in the fastq read
    
    AD2: str, optional, default='AGATCGGAAGAGCACAC'
        Specifies the second adaptor sequence expected in the fastq read
    
    perfect_matches: boolean, optional, default=False
        If True, returns a pandas dataframe containing fastq reads with
        an exact match to AD1 and AD2
    
    
    Returns
    -------
    Pandas Dataframe
        multi-index: (plate_barcode, sequence)
        columns: ['plate_barcode', 'seq', 'seq_count', 'rasl_probe', 'observed_wellbc']    
    
    '''
    
    #CREATING PERFECT MATCH SEQUENCE COLUMN ('rasl_probe')
    collapsed_read_df['rasl_probe'] = collapsed_read_df.apply(rasl_probe_seq_extraction, args=[AD1, AD2], axis=1)
    
    #CREATING OBSERVED WELL BARCODE ('observed_wellbc') COLUMN
    collapsed_read_df['observed_wellbc'] = collapsed_read_df.apply(observed_wellbc_extraction, args=[AD1], axis=1)
    
    
    if print_on:
        print len(collapsed_read_df[collapsed_read_df.rasl_probe != '0']), 'number of perfect matches to adaptor seq'
        print len(collapsed_read_df[collapsed_read_df.observed_wellbc != '0']), 'number of perfect matches to ad1'
        print len(collapsed_read_df.rasl_probe.unique()), 'number of unique wellbarcodes that need to be mapped'
    
    if perfect_matches:
        #ONLY CONSIDER PERFECT MATCHES TO ADAPTORS
        return collapsed_read_df[(collapsed_read_df.rasl_probe !='0') & (collapsed_read_df.observed_wellbc !='0')]     
    
    return collapsed_read_df
    

# <codecell>

def get_blast_alignments(collapsed_read_df, blastn_path, db_path, blast_write_path):
    '''
    This function returns the rasl_probe sequence blast alignments
    
    Parameters
    ----------
    collapsed_read_df: pandas dataframe, must contain columns ['rasl_probe']
        collapsed_read_df['rasl_probe'] - fastq sequence observed between rasl adaptor sequences
        
    blastn_path: str, path to blastn executable
        Specifies full path to blastn executable
    
    db_path: str, path to target BLAST database
        Specifies path to on_off_target BLAST database
    
    blast_write_path: str, path to scratch space
        Specifies where to write temporary blast input/output files
        

    Returns
    -------
    Pandas Blast Results dataframe
        index = qseqid  #which is the ~40 nt rasl_probe sequence
        cols = ['qseqid','mismatch','evalue','sseqid','qlen', 'length','qseq','qstart','sseq','sstart','send']

    '''
    
    random_file_handle = str(random.randrange(0,1000000))  #setting random seed for temp file writing
    blast_write_path = blast_write_path + 'temp_blast_' + random_file_handle +".txt" 
    
    #WRITING TEMP BLAST INPUT FILE
    blast_input = open(blast_write_path,"w")
    for items in collapsed_read_df.rasl_probe.unique():  #only blasting unique perfect matches
        blast_input.write(">"+items + "\n" + items + "\n")
        blast_input.flush()
    blast_input.close()
    
    
    #ARGUMENTS FOR BLAST FUNCTION
    bl_results = rasl_probe_blast(blast_write_path, blastn_path, db_path)
    
    
    #collapsed_read_df.set_index('rasl_probe',inplace=True,drop=False)
    #collapsed_read_df = collapsed_read_df.join(bl_results,how='inner')
    

    return bl_results


# <codecell>

#TESTING
# a = get_collapsed_read_df('/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/fastq/P1_CGGTTCT_L001_R1_001_100k.fastq.gz',print_on=True)     
# b = find_get_rasl_probe_and_wellbc_exact(a,print_on=True)

# bl_write_path = "/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/blastdb/"
# blastn_path1 = '/gpfs/home/erscott/Tools/blast/ncbi-blast-2.2.26+/bin/blastn'
# db_path = "/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/blastdb/on_off_targetseq"
# bl_results = get_blast_alignments(b, blastn_path1, db_path, bl_write_path)


# b.set_index('rasl_probe',inplace=True,drop=False)
# collapsed_read_df = b.join(bl_results,how='inner')
# collapsed_read_df = collapsed_read_df[(collapsed_read_df.length >30) & (collapsed_read_df.qstart < 6) & (collapsed_read_df.observed_wellbc.map(len) < 10) & (collapsed_read_df.observed_wellbc.map(len) > 6)]     


# <headingcell level=2>

# 
# 
# Well and Plate Barcodes

# <codecell>

#plate_barcode = index read
#WellBarcode = first ~8nt in read
path = '/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/'
barcodes = pd.read_table(path+'SampleID_Barcodes.txt',sep="\t")
wellbarcode = barcodes.WellBarcode

# <codecell>

def fuzzy_wellbc_match(obs_wellbc, well_barcodes, start_pos, end_pos):
    '''
    This function takes a read and searches for supplied barcode sequences.
    
    Parameters
    ----------
    obs_wellbc: str, fastq sequence before the first instance of AD1
        e.g. ATGCATG
    
    well_barcodes: list, expected well barcodes
        e.g. ATGCATG
    
    start_pos: int, limits the string search space of the obs_wellbc
    
    end_position: int, limits the string search space of the obs_wellbc
    
    Returns
    -------
    The expected barcode found in the obs_wellbc OR 'mismatch' (if no barcode is found) 
    
    
    '''
    
    assert type(start_pos) == int, 'start_pos should be an int'
    assert type(end_pos) == int, 'end_pos should be an int'
    
    #initializing obs_wellbc variable and set of expected well barcodes
    FASTQ, bc_set = obs_wellbc.upper(), set(well_barcodes) 
    
    #limit search for exact matches to [: end_pos]
    matches = set(FASTQ[n:n+8] for n in range(0, end_pos) if FASTQ[n:n+8] in bc_set)
    
    #DEPRECATED
    #[matches.add(FASTQ[n:n+8]) for n in range(0, end_pos) if FASTQ[n:n+8] in bc_set]
    
    
    #RETURNS EXPECTED BARCODE SEQUENCE IF UNIQUE MATCH FOUND
    if len(matches) == 1: return (0, list(matches)[0])  #return the single best match

    #BRUTE FORCE SEARCH OF obs_wellbc subsequence
    else: 
        matches = set()
        for bars in well_barcodes:  #Differentiate between more than 1 exact match, or find fuzzy matches
            BARS = bars.upper()  #added a bit of ADAPTOR1 to make the mappings more stringent
            edist = editdist.distance(FASTQ, BARS)
            delta_8 = abs(len(FASTQ) - 8)  #correcting for difference in string seq lengths
            if edist < 2:
                return (edist, BARS)
            else:
                if edist + delta_8 < 2:  #looser thresholds performed poorly and only added maybe a couple hundred reads out of a million
                    matches.add((edist, bars))
        if len(matches) >0: return (len(matches), ";".join([i[1] for i in matches]) )
        return (8, "mismatch")

# <codecell>

#Fuzzy matching observed wellbc with expected wellbc    
wellbc_mappings = dict((obs_bc, fuzzy_wellbc_match(obs_bc, wellbarcode, 0, 9)) for obs_bc in collapsed_read_df['observed_wellbc'].unique())

#converting wellbc mappings into a dataframe
wellbc_mappings_df = pd.DataFrame.from_dict(wellbc_mappings,orient='index')  #index is observed_wellbc
wellbc_mappings_df = wellbc_mappings_df[[0,1]]
wellbc_mappings_df.columns = ['bc_edit_dist','mapped_bc']  #index is wellbarcode


# <codecell>

#MERGING COLLAPSED_READ_DF WITH FUZZY WELLBC MAPPINGS
collapsed_read_df.index = collapsed_read_df.observed_wellbc  #setting index to observed_wellbc
collapsed_read_df1 = collapsed_read_df.join(wellbc_mappings_df)  #Merging wellbarcode mappings with blast aligned fastq reads
##df_plate_seq_blast_reliable.head()
print 'done'

# <codecell>

#collapsed_read_df1[['observed_wellbc','mapped_bc','bc_edit_dist','seq_count','seq']].to_csv('/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/analysis/P1_well_barcode_mappings.txt',sep="\t", index=False)        



# <headingcell level=3>

# Well Annotations

# <codecell>

def well_annot_df(well_annot_path, delim):
    '''
    This function creates a well_annot_df from a specified file
    
    Parameters
    ----------
    well_annot_path: str
        Specifies path to well annotations file possessing
        columns headers 'plate_barcode' and 'WellBarcode'
        
    delim: str
        Specifies delimmiter to use when parsing well annotation file
        
    Returns
    -------
    Pandas DataFrame
        index: ('plate_barcode','WellBarcode')
    
    '''
    bc_condition_key = pd.read_table(well_annot_path, sep=delim)
    bc_condition_key.set_index(['plate_barcode','WellBarcode'],inplace=True)
    return bc_condition_key

# <codecell>

well_annot_df = well_annot_df('/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/Bcell_exp_2013_10_bc_key.txt', '\t')     


# <headingcell level=3>

# Summing Probe-specific Read Counts

# <codecell>

def get_probe_well_read_counts(collapsed_read_counts):
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
    
    collapsed_read_counts_group = collapsed_read_counts.groupby(['plate_barcode','mapped_bc','probe'])
    
    counts = collapsed_read_counts_group.seq_count.aggregate(np.sum)
    counts_df = counts.unstack(level=2)  #creating matrix of aggregated probe counts indexed on 'plate_barcode','mapped_bc'
    
    #ID AND REMOVAL OF OFF-TARGET LIGATION HITS FOUND BY ! CHARACTER IN BLAST SSEQ NAME
    on_target_col = [i for i in counts_df.columns if "!" not in i]  
    counts_df = counts_df[on_target_col]  #removing off-target ligation hits from df
    
    return counts_df

# <codecell>

#Filtering using bc_edit_dist
collapsed_read_df1 = collapsed_read_df1[(collapsed_read_df1.bc_edit_dist.astype(float)<2)]


#Final Well Read Count dataframe
collapsed_read_df_counts = get_probe_well_read_counts(collapsed_read_df1)



# <headingcell level=3>

# Merging Probe-specific Read Counts and Well Annotations

# <codecell>

def merge_plate_well_annot(probe_counts_df, well_annot_df):
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
    

# <codecell>

#MERGING WELL ANNOTATIONS AND PROBE READ COUNTS
counts_df = merge_plate_well_annot(collapsed_read_df_counts, well_annot_df)

# <headingcell level=2>

# ENTIRE PROCESS

# <codecell>

class RASLseqAnalysis(object):
    def __init__(self, fastq_path, blastdb_path, blastn_path, write_path, print_on=False):
        '''
        Initialize RASLseqAnalysis object
        
        Parameters
        ----------
        fastq_path: str
            Specifies the full path to the fastq file, handles gzip and non-gzip fastq files
        
        '''
        
        self.fastq_file = fastq_path
        
        self.blastdb = blastdb_path
        
        self.blastn = blastn_path
        
        self.write = write_path
        
        self.read_df = collapsed_read_df(self.fastq_file, print_on = self.print_on)
        
    
    def get_collapsed_reads_df(self):
        '''
        This function collapses identical fastq read sequences and extracts
        the 3' barcode read to return a pandas dataframe:
            index: (plate_barcode, sequence)
            columns: ['plate_barcode', 'seq', 'seq_count']
        
        Parameters
        ----------
        print_on: boolean, default=False
            If print_on == True prints "number of unique bc/reads to map total"
        
        Returns
        -------
        Pandas dataframe
            multi-index: (plate_barcode, sequence)
            columns: ['plate_barcode', 'seq', 'seq_count']
        
                plate_barcode - Index read from fastq file
                seq - fastq read sequence
                seq_count - number of occurrences of seq in fastq
        '''
        
        self.collapsed_read_df = FastqCollapse.get_collapsed_fastq_df(self.fastq_file, print_on=self.print_on)
        
        return True
    
    
    
    def get_rasl_probe_and_wellbc_exact(self):
        
        return get_rasl_probe_and_wellbc_exact(self.get_collapsed_reads_df(), print_on=self.print_on)
    
    
    
        

# <codecell>

bl_write_path = "/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/blastdb/"
blastn_path = '/gpfs/home/erscott/Tools/blast/ncbi-blast-2.2.26+/bin/blastn'
db_path = "/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/blastdb/on_off_targetseq"
fq_path = '/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/fastq/P1_CGGTTCT_L001_R1_001_100k.fastq.gz'

test = RASLseqAnalysis(fq_path, db_path, blastn_path, bl_write_path)

# <codecell>

c = test.get_collapsed_reads_df()

# <codecell>

if __name__ == '__main__':
    ##RUN FROM FASTQ TO COUNTS_DF##
    
    
    #CREATE FASTQ COLLAPSED DATAFRAME
    collapsed_read_df = get_collapsed_read_df('/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/fastq/P1_CGGTTCT_L001_R1_001_100k.fastq.gz',print_on=True)     
    
    
    #ID RASL_PROBE AND OBSERVED WELL BARCODE (BC) SEQUENCES USING ADAPTOR SEQUENCES AND ADD TO COLLAPSED_READ_DF
    collapsed_read_df = get_rasl_probe_and_wellbc_exact(collapsed_read_df, print_on=True)
    
    
    #SET BLAST PARAMETERS
    bl_write_path = "/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/blastdb/"
    blastn_path = '/gpfs/home/erscott/Tools/blast/ncbi-blast-2.2.26+/bin/blastn'
    db_path = "/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/blastdb/on_off_targetseq"
    
    #BLAST RASL_PROBE SEQ AGAINST ALL COMBINATIONS OF ACCEPTOR AND DONOR PROBES
    bl_results = get_blast_alignments(collapsed_read_df, blastn_path, db_path, bl_write_path)
    
    #JOINING BLAST RESULTS WITH COLLAPSED_READ_DF
    collapsed_read_df.set_index('rasl_probe',inplace=True,drop=False)
    collapsed_read_df = collapsed_read_df.join(bl_results,how='inner')
    
    #FILTERING BLAST RESULTS
    collapsed_read_df = collapsed_read_df[(collapsed_read_df.length >30) & (collapsed_read_df.qstart < 6) & (collapsed_read_df.observed_wellbc.map(len) < 10) & (collapsed_read_df.observed_wellbc.map(len) > 6)]     
    
    
    
    #SETTING wellbarcode SEQUENCES
    path = '/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/'
    barcodes = pd.read_table(path+'SampleID_Barcodes.txt',sep="\t")
    wellbarcode = barcodes.WellBarcode
    
    #FUZZY MATCHING OBSERVED WELL BC WITH EXPECTED WELL BC    
    wellbc_mappings = dict((obs_bc, fuzzy_wellbc_match(obs_bc, wellbarcode, 0, 9)) for obs_bc in collapsed_read_df['observed_wellbc'].unique())
    
    #CONVERTING WELL BC MAPPINGS INTO A DF
    wellbc_mappings_df = pd.DataFrame.from_dict(wellbc_mappings,orient='index')  #index is observed_wellbc
    wellbc_mappings_df = wellbc_mappings_df[[0,1]]
    wellbc_mappings_df.columns = ['bc_edit_dist','mapped_bc']  #index is wellbarcode
    
    #MERGING COLLAPSED_READ_DF WITH FUZZY WELLBC MAPPINGS
    collapsed_read_df.index = collapsed_read_df.observed_wellbc  #setting index to observed_wellbc
    collapsed_read_df = collapsed_read_df.join(wellbc_mappings_df)  #Merging wellbarcode mappings with blast aligned fastq reads
    
    #FILTERING COLLAPSED_READ_DF USING bc_edit_dist
    collapsed_read_df = collapsed_read_df[(collapsed_read_df.bc_edit_dist.astype(float)<2)]
    
    
    #WELL ANNOTATION DATAFRAME
    well_annot_df = get_well_annot_df('/gpfs/home/erscott/Datasets_raw/RASLseq/Bcell_exp_2013_10/Bcell_exp_2013_10_bc_key.txt', '\t')     
    
    #SUM PROBE-SPECIFIC READ COUNTS
    collapsed_read_df_counts = get_probe_well_read_counts(collapsed_read_df)
    
    #MERGING WELL ANNOTATIONS AND PROBE READ COUNTS
    counts_df = merge_plate_well_annot(collapsed_read_df_counts, well_annot_df)
    
    print
    print 'FASTQ to COUNTS_DF COMPLETE'
    







