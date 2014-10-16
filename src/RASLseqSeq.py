
'''
These functions extract RASLseq Probe Sequences from a fastq read
by finding the sequence between 2 adaptor strings
'''



def rasl_probe_seq_extraction(df, AD1='GGAAGCCTTGGCTTTTG', AD2='AGATCGGAAGAGCACAC'):
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
    initial_df_len = len(df)
    
    df['AD1_in'] = df.seq.str.count(AD1)
    df['AD2_in'] = df.seq.str.count(AD2)
    
    df_perfect_match = df[(df['AD1_in']==1) & (df['AD2_in']==1)]
    df_perfect_match['rasl_probe'] = df_perfect_match.seq.str.split(AD1).str[1].str.split(AD2).str[0]
    
    df = df[(df['AD1_in']!=1) | (df['AD2_in']!=1)]
    
    df['rasl_probe'] = '0'
    df = df.append(df_perfect_match)
    assert initial_df_len == len(df)
    return df
    
#     seq = line['seq']
#     if AD1 in seq and AD2 in seq:
#         return seq[seq.index(AD1)+17: seq.index(AD2)]
#     else: return '0'



def get_rasl_probe_exact_df(collapsed_read_df, AD1='GGAAGCCTTGGCTTTTG', AD2='AGATCGGAAGAGCACAC', perfect_matches=True, print_on=False):   
    '''
    This function identifies the rasl_probe (~40mer) sequence  
    within the fastq read using an exact match to the adaptor 
    sequences (AD1 & RCAD2)
    
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
    #collapsed_read_df['rasl_probe'] = collapsed_read_df.apply(rasl_probe_seq_extraction, args=[AD1, AD2], axis=1)
    collapsed_read_df = rasl_probe_seq_extraction(collapsed_read_df, AD1, AD2)
    
    
    if print_on:
        print len(collapsed_read_df[collapsed_read_df.rasl_probe != '0']), 'number of perfect matches to adaptor seq'
        
    
    if perfect_matches:
        #ONLY CONSIDER PERFECT MATCHES TO ADAPTORS
        return collapsed_read_df[(collapsed_read_df.rasl_probe !='0')]     
    
    return collapsed_read_df



    
    
    
    