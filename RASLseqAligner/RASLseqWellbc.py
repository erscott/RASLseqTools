


import editdist
import pandas as pd



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
    
    
    



def get_rasl_probe_and_wellbc_exact_df(collapsed_read_df, AD1='GGAAGCCTTGGCTTTTG', AD2='AGATCGGAAGAGCACAC', perfect_adaptor_matches=True, print_on=False):   
    '''
    This function identifies the observed well barcode
    sequence within the fastq read using an exact
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
    
    
    #CREATING OBSERVED WELL BARCODE ('observed_wellbc') COLUMN
    collapsed_read_df['observed_wellbc'] = collapsed_read_df.apply(observed_wellbc_extraction, args=[AD1], axis=1)
    
    
    if print_on:
        print len(collapsed_read_df[collapsed_read_df.observed_wellbc != '0']), ':number of observed well barcodes using perfect matches to AD1'
        
    
    if perfect_adaptor_matches:
        #ONLY CONSIDER PERFECT MATCHES TO ADAPTORS
        return collapsed_read_df[(collapsed_read_df.observed_wellbc !='0')]     
    
    return collapsed_read_df






def get_rasl_probe_and_wellbc_fuzzy_df(collapsed_read_df, wellbarcodes, AD1='GGAAGCCTTGGCTTTTG', AD2='AGATCGGAAGAGCACAC', print_on=False):   
    '''
    This function identifies the observed well barcode
    sequence within the fastq read using an exact
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
    
        
    
    #Fuzzy matching observed wellbc with expected wellbc    
    wellbc_mappings = dict((obs_bc, fuzzy_wellbc_match(obs_bc, wellbarcodes, 0, 9)) for obs_bc in collapsed_read_df['observed_wellbc'].unique())
    
    #converting wellbc mappings into a dataframe
    wellbc_mappings_df = pd.DataFrame.from_dict(wellbc_mappings,orient='index')  #index is observed_wellbc
    wellbc_mappings_df = wellbc_mappings_df[[0,1]]
    wellbc_mappings_df.columns = ['bc_edit_dist','mapped_bc']  #index is wellbarcode
    

    
    #MERGING COLLAPSED_READ_DF WITH FUZZY WELLBC MAPPINGS
    collapsed_read_df.index = collapsed_read_df.observed_wellbc  #setting index to observed_wellbc
    collapsed_read_df = collapsed_read_df.join(wellbc_mappings_df)  #Merging wellbarcode mappings with blast aligned fastq reads
    
    if print_on:
        print len(collapsed_read_df[collapsed_read_df.mapped_bc != 'mismatch' ]), ':number of well barcodes with edit distance < 2 using perfect matches to AD1'
    
    
    return collapsed_read_df













    



