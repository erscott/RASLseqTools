
import pandas as pd
from collections import Counter
import gzip


def get_file_obj(file_path):
    '''
    This function returns a file object
    given a gzip or uncompressed file path
    '''
    if file_path.endswith('.gz'):
        return gzip.open(file_path,'r')
    else:
        return open(file_path, 'r')



def fastq_collapsing_Counter(fastq_file_obj, sequencer_id):
    """
    Implementation of a Counter (Collections library) object filled with FASTQ reads sequences 
    (plus 3' barcode) by iterating through a fastq file and collapsing identical reads
    
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
            seq_library_counter[(header_bc_2, fq_read.rstrip())] +=1  #faster than Counter.update
            next_line_seq = 0
        if fq_read.startswith(sequencer_id):  #save barcode_3' read
            header_bc_2 = fq_read.rstrip().split(":")[-1][:7]  #grab 3' barcode
            next_line_seq = 1
        
    return seq_library_counter



def get_collapsed_fastq_df(fastq_path, sequencer_id, print_on=False):
    '''
    This function collapses identical fastq read sequences and extracts
    the 3' barcode read to return a pandas dataframe:
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
    
    #CREATE COUNTER OBJ OF READS
    seq_lib = fastq_collapsing_Counter(fastq_obj, sequencer_id)  #~6.5-fold compression
    
    
    #CREATING seq_count DATAFRAME THAT CAN BE ACCESSED BY MULTI-INDEX
    df = pd.DataFrame.from_dict(seq_lib,orient='index')  #index = (barcode, sequence) ; value = seq_count from collapsing dictionary
    df.columns = ['seq_count']
    df.index = pd.MultiIndex.from_tuples(df.index)
    df.index.names = ['plate_barcode', 'seq']
    
    
    #CREATING plate_barcode AND seq DF
    df_plate_seq = pd.DataFrame.from_records(seq_lib.keys())
    df_plate_seq.columns = ['plate_barcode','seq']
    df_plate_seq.set_index(['plate_barcode','seq'], drop=False, inplace=True)


    df_plate_seq = df_plate_seq.join(df)  #merged dataframe; index=(plate_barcode, sequence); columns = ['plate_barcode', 'seq', 'seq_count'] 

    
    if print_on:
        print len(seq_lib), "number of unique bc/reads to map total"
        
    
    return df_plate_seq














