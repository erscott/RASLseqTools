
import random
import pandas as pd
import os
import sys

def rasl_probe_blast(read_file_path, aligner_path, db_path, print_on=False):
    """
    This function returns blast results from queries optimized for rasl-seq conditions and outputs
    a custom format
    
    Parameters
    ----------
    read_file_path: str, path to temp blast input file
        Specifies where temporary blast input file has been written
        
    aligner_path: str, path to blastn or STAR executable
        Specifies full path to blastn or STAR executable
    
    db_path: str, path to target BLAST or STAR database
        Specifies path to on_off_target BLAST or STAR database
    
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
    
    
    #blast_run = subprocess.Popen("/gpfs/home/erscott/Tools/blast/ncbi-blast-2.2.26+/bin/blastn -task megablast -query "+ read_file_path +" -evalue 20 "+ wordsize + db +" > new_seq_blast_results.txt",shell=True,stdout=subprocess.PIPE).stdout
    blast_run = os.system(aligner_path+"/blastn" + " -task blastn-short -query " \
                          + read_file_path +" -evalue 1e-6 "+ wordsize + db +" -max_target_seqs 1" +" -strand 'plus' -xdrop_gap 7 -num_threads 12" \
                          + " -outfmt '6 -outfmt qseqid -outfmt mismatch -outfmt evalue -outfmt sseqid -outfmt qlen -outfmt length -outfmt qseq -outfmt qstart -outfmt sseq -outfmt sstart -outfmt send' >"+ write_path)
    # Evlaue of 20 was used to account for a cutoff of 13mer matching
    
    if print_on:
        print "BLAST COMMAND\n" + aligner_path + " -task blastn-short -query " \
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




def get_blast_alignments(collapsed_read_df, aligner_path, db_path, print_on=False):
    '''
    This function returns the rasl_probe sequence blast alignments
    
    Parameters
    ----------
    collapsed_read_df: pandas dataframe, must contain columns ['rasl_probe']
        collapsed_read_df['rasl_probe'] - fastq sequence observed between rasl adaptor sequences
        
    aligner_path: str, path to blastn or STAR executable
        Specifies full path to STAR executable
    
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
    bl_results = rasl_probe_blast(blast_write_path, aligner_path, db_path, print_on=False)

    return bl_results
    
    
   
   
        
    

def get_star_alignments(collapsed_read_df, aligner_path, db_path, print_on=False, n_jobs=1, offset_5p=16, offset_3p=16):  
    
        #/Users/ers_vader/Tools/STAR_2.3.0e.OSX_x86_64/STAR   --genomeDir ../../STAR_RRASL/ --readFilesIn lane1_Undetermined_L001_R1_001_truth_set_reads.fastq --runThreadN 4 
        #--seedSearchStartLmax 1 --seedPerReadNmax 15000 --scoreDelOpen -1 --scoreDelBase -2 --scoreInsOpen -1 --scoreInsBase -2 --outFilterMultimapNmax 100 
        #--outFilterMatchNminOverLread 0.6 --scoreGapNoncan -1 --outFilterScoreMinOverLread 0.6 --clip5pNbases 16 --clip3pNbases 16

    if 'gz' in collapsed_read_df:
        os.system('gunzip ' + collapsed_read_df)
        collapsed_read_df.replace('.gz', '')
        
    random_file_handle = str(random.randrange(0,1000000))
    star_write_path = db_path + 'temp_star_' + random_file_handle +"_"
    cmd = [aligner_path + '/STAR', '--genomeDir', db_path, '--readFilesIn', collapsed_read_df, '--runThreadN', str(n_jobs), \
           '--seedSearchStartLmax 1', '--seedPerReadNmax', str(15000), '--scoreDelOpen -1', '--scoreDelBase -2', '--scoreInsOpen -1', \
           '--scoreInsBase -2', '--outFilterMultimapNmax 100', '--outFilterMatchNminOverLread 0.6', '--scoreGapNoncan -1', \
           '--outFilterScoreMinOverLread 0.6', '--clip5pNbases', str(offset_5p), '--clip3pNbases', str(offset_3p), \
           '--outStd SAM']
    
    #reading STAR ALGINMENTS FROM STDOUT into pandas df
    star_results = pd.DataFrame.from_records([i.split('\t') for i in os.popen(' '.join(cmd)) if i[0] !='@'])[[0,1,2,4,5,9,10]]
    star_results.rename(columns={2:'ProbeName', 0:'id_line'},inplace=True)
    star_results = star_results[star_results[1] !=256]  #removing multimapper
    
    if print_on:
        print 'STAR COMMAND'
        print ' '.join(cmd)
        
        
    return star_results
    
    



##PASS AN OBJECT (RASLseqReads) WITH ALL OF THESE ATTRIBUTES OR PASS THE INDIVIDUAL OBJECTS?
def get_rasl_aligned_df(collapsed_read_df, aligner_path, db_path, print_on=False, aligner='star', n_jobs=1, offset_5p=16, offset_3p=16):
    
    #PASS OBJECT
    if aligner == 'blast':
        #BLAST RASL_PROBE SEQ AGAINST ALL COMBINATIONS OF ACCEPTOR AND DONOR PROBES
        bl_results = get_blast_alignments(collapsed_read_df, aligner_path, db_path, print_on=False)
        
        #JOINING BLAST RESULTS WITH COLLAPSED_READ_DF
        collapsed_read_df.set_index('rasl_probe',inplace=True,drop=False)
        collapsed_read_df = collapsed_read_df.join(bl_results,how='inner')
        return collapsed_read_df
    
    
    if aligner == 'star':
        return get_star_alignments(collapsed_read_df, aligner_path, db_path, print_on=False, n_jobs=n_jobs, offset_5p=offset_5p, offset_3p=offset_3p )




def mp_get_rasl_aligned_df(args):
    collapsed_read_df, aligner_path, db_path, print_on, aligner, n_jobs, offset_5p, offset_3p = args
    return get_star_alignments(collapsed_read_df, aligner_path, db_path, print_on=False, n_jobs=n_jobs, offset_5p=offset_5p, offset_3p=offset_3p )
        












