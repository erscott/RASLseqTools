
import random
import pandas as pd
import os
import sys

def rasl_probe_blast(read_file_path, blastn_path, db_path, print_on=False):
    """
    This function returns blast results from queries optimized for rasl-seq conditions and outputs
    a custom format
    
    Parameters
    ----------
    read_file_path: str, path to temp blast input file
        Specifies where temporary blast input file has been written
        
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
    
    
    db = " -db " + db_path 
    
    
    #SETTING BLAST WORD SIZE
    wordsize = "-word_size 8 "
    
    
    
    blast_run = os.system(blastn_path+"/blastn" + " -task blastn-short -query " \
                          + read_file_path +" -evalue 1e-6 "+ wordsize + db +" -max_target_seqs 1" +" -strand 'plus' -xdrop_gap 7 -num_threads 12" \
                          + " -outfmt '6 -outfmt qseqid -outfmt mismatch -outfmt evalue -outfmt sseqid -outfmt qlen -outfmt length -outfmt qseq -outfmt qstart -outfmt sseq -outfmt sstart -outfmt send' >"+ write_path)
    # Evlaue of 20 was used to account for a cutoff of 13mer matching
    
    if print_on:
        print "BLAST COMMAND\n" + blastn_path + " -task blastn-short -query " \
                              + read_file_path +" -evalue 1e-6 "+ wordsize + db +" -max_target_seqs 1" +" -strand 'plus' -xdrop_gap 7 -num_threads 12" \
                              + " -outfmt '6 -outfmt qseqid -outfmt mismatch -outfmt evalue -outfmt sseqid -outfmt qlen -outfmt length -outfmt qseq -outfmt qstart -outfmt sseq -outfmt sstart -outfmt send' >"+ write_path
    #CONVERTING BLAST OUTPUT INTO DATAFRAME
    bl_results = pd.read_table(write_path, header=None)
    bl_results.columns=['qseqid','mismatch','evalue','sseqid','qlen', 'length','qseq','qstart','sseq','sstart','send']
    bl_results['probe'] = bl_results.sseqid.apply(lambda x: "_".join(x.split("_")))
    bl_results.ix[0]
    bl_results.set_index('qseqid',drop=False,inplace=True)
    
    #REMOVING INPUT/OUTPUT BLAST FILES
    os.system('rm ' + read_file_path)
    os.system('rm '+ write_path)
    os.system('rm ' + db_path + ".*")
    
    return bl_results




def get_blast_alignments(collapsed_read_df, blastn_path, db_path, print_on=False):
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
    

    Returns
    -------
    Pandas Blast Results dataframe
        index = qseqid  #which is the ~40 nt rasl_probe sequence
        cols = ['qseqid','mismatch','evalue','sseqid','qlen', 'length','qseq','qstart','sseq','sstart','send']

    '''
    #CREATING BLAST INPUT FILE
    random_file_handle = str(random.randrange(0,1000000))  #setting random seed for temp file writing
    blast_write_path = db_path + 'temp_blast_' + random_file_handle +".txt" 
    
    #WRITING TEMP BLAST INPUT FILE
    blast_input = open(blast_write_path,"w")
    for items in collapsed_read_df.rasl_probe.unique():  #only blasting unique perfect matches
        blast_input.write(">"+items + "\n" + items + "\n")
        blast_input.flush()
    blast_input.close()
    
    
    #ARGUMENTS FOR BLAST FUNCTION
    bl_results = rasl_probe_blast(blast_write_path, blastn_path, db_path, print_on=False)
    
    
    return bl_results





##PASS AN OBJECT (RASLseqReads) WITH ALL OF THESE ATTRIBUTES OR PASS THE INDIVIDUAL OBJECTS?
def get_rasl_blast_df(collapsed_read_df, blastn_path, db_path, print_on=False):
    
    #PASS OBJECT
    #BLAST RASL_PROBE SEQ AGAINST ALL COMBINATIONS OF ACCEPTOR AND DONOR PROBES
    bl_results = get_blast_alignments(collapsed_read_df, blastn_path, db_path, print_on=False)
    
    #JOINING BLAST RESULTS WITH COLLAPSED_READ_DF
    collapsed_read_df.set_index('rasl_probe',inplace=True,drop=False)
    collapsed_read_df = collapsed_read_df.join(bl_results,how='inner')
    
    
    return collapsed_read_df


    
    










    
