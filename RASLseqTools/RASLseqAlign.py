
import random
import pandas as pd
import os
import sys
import Levenshtein as editdist
import swalign
import numpy as np
import itertools
from scipy.spatial import distance


#BLAST ALIGNMENT

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
    
    
    db = " -db " + db_path 
    
    
    #SETTING BLAST WORD SIZE
    wordsize = "-word_size 8 "
    
    
    #blast_run = subprocess.Popen("/path/to/bin/blastn -task megablast -query "+ read_file_path +" -evalue 20 "+ wordsize + db +" > new_seq_blast_results.txt",shell=True,stdout=subprocess.PIPE).stdout
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




#STAR ALGINMENT




#STAR ALIGNMENT
def get_star_alignments(read_fq_path, aligner_path, db_path, print_on=False, n_jobs=1, offset_5p=16, offset_3p=16):  
    '''
    This function aligns fastq file using the STAR aligner
    
    Parameters
    ----------
    read_fq_path: str, path to fastq file
        
        
    aligner_path: str, path to blastn or STAR executable
        Specifies full path to STAR executable
    
    db_path: str, path to target STAR database
        Specifies path to on_off_target STAR database
        
    print_on: bool, default=False
        Print STAR Command 
        
    n_jobs: int, default=1
        Number of threads to use for alignment
        
    offset_5p: int, default=16
        Index position for RASLprobe start,
        Number of bases to clip from 5' end of read
        to isolate probe sequence
        
    offset_3p: int, defaults=16
        Number of bases from 3' end of read to clip
        in order to isolate probe sequence in read
        
    
    Returns
    -------
    star_results, pandas dataframe
        columns: SAM columns [0,1,2,4,5,9,10]
        Drops multi-mapping reads

    
    EXAMPLE COMMAND
        #/path/to/STAR_2.3.0e.OSX_x86_64/STAR   --genomeDir ../../STAR_RRASL/ --readFilesIn lane1_Undetermined_L001_R1_001_truth_set_reads.fastq --runThreadN 4 
        #--seedSearchStartLmax 1 --seedPerReadNmax 15000 --scoreDelOpen -1 --scoreDelBase -2 --scoreInsOpen -1 --scoreInsBase -2 --outFilterMultimapNmax 100 
        #--outFilterMatchNminOverLread 0.6 --scoreGapNoncan -1 --outFilterScoreMinOverLread 0.6 --clip5pNbases 16 --clip3pNbases 16
    '''
    if 'gz' in read_fq_path:
        os.system('gunzip ' + read_fq_path)
        read_fq_path.replace('.gz', '')
        
    random_file_handle = str(random.randrange(0,1000000))
    star_write_path = db_path + 'temp_star_' + random_file_handle +"_"
    cmd = [aligner_path + '/STAR', '--genomeDir', db_path, '--readFilesIn', read_fq_path, '--runThreadN', str(n_jobs), \
           '--seedSearchStartLmax 1', '--seedPerReadNmax', str(15000), '--scoreDelOpen -1', '--scoreDelBase -2', '--scoreInsOpen -1', \
           '--scoreInsBase -2', '--outFilterMultimapNmax 100', '--outFilterMatchNminOverLread 0.7', '--scoreGapNoncan -1', \
           '--outFilterScoreMinOverLread 0.7', '--clip5pNbases', str(offset_5p), '--clip3pNbases', str(offset_3p), \
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
    '''
    This function coordinates aligment of collapsed_read_df
    '''
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
        


#LEVENSHTEIN EDIT DISTANCE
def edist_mp(seq1_seq2):
    '''
    Pairwise Levenshtein Edit Distances
    '''
    return editdist.distance(seq1_seq2[0], seq1_seq2[1])


#PAIRWISE LEVENSHTEIN EDIT DISTANCE DataFrame
def pairwise_distances(sequences):
    '''
    This function returns a dataframe of pair-wise distances between sequences
    
    Parameters
    --------------
    sequences: list or set of barcode sequences
    
    
    Returns
    --------------
    sequences distance matrix, pandas dataframe
        squareform
        
    '''
    barcode_count = len(sequences)
    sequences = sorted(list(sequences))
    levenshtein_dist = np.array( map(edist_mp, itertools.combinations(sequences, 2)) )
    levenshtein_dist = distance.squareform(levenshtein_dist).reshape(barcode_count,barcode_count)
    levenshtein_dist = pd.DataFrame(levenshtein_dist, columns=sequences, index=sequences)

    return levenshtein_dist



def cigar_to_align(sw_cigar):
    '''
    Parse swalign CIGAR
    '''
    
    cigar_trans = {'M':'|','D':'D','I':'I'}
    l = [i[0] * cigar_trans[i[1]]  for i in sw_cigar]
    return ''.join(l)
    
    

def swalign_df(ref, query):
    '''
    This function returns swalign info:
    ref, query, r_pos, r_end, q_pos, q_end, score, matches, mismatches, identity, cigar
    '''
    match = 2
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    
    sw = swalign.LocalAlignment(scoring, gap_extension_penalty=-2, prefer_gap_runs=False)  # you can also choose gap penalties, etc...
    aligned = sw.align(ref, query)  #ref, query
    #return alignment
    
    align_series = pd.Series([ref, query, aligned.r_pos, aligned.r_end, aligned.q_pos, aligned.q_end,\
             aligned.score, aligned.matches, aligned.mismatches, \
             aligned.identity, cigar_to_align(aligned.cigar)])
    align_series.index = 'ref, query, r_pos, r_end, q_pos, q_end, score, matches, mismatches, identity, cigar'.split(', ')
    return  align_series



def mp_swalign(align_probes):
    '''
    This function uses smith waterman alignment
    to map a ref sequence against a set of probes
    
    Parameters
    ------------
    align_probes: tuple(ref, query_probes), required
        ref is the reference string to align against
        query_probes is a set of probes to align against the reference
        
    Returns
    -------------
    max alignment score
    
    max alignment frequency
    
    alignment objects 
    

    Usage:
    #results = map(mp_swalign, alignment_tasks)
    #ontarget_align_df = [n for i in results for n in i[-1]]  #flatten results
    #ontarget_align_df = pd.DataFrame(ontarget_align_df)  #aggregated dataframe
    '''
    ref, q_probes = align_probes[0], align_probes[1]  #str and list
    
    min_dist = 0  #min distance
    min_dist_count = 0  #count of non-self probes with min distance
    best_alignments = []
    for q_probe in q_probes:  #iterate through reference probes
        if q_probe == ref: continue  # only consider non-self alignments
        temp_align = swalign_df(ref, q_probe)  #sw alignment
        temp_align_score = temp_align.score 
        
        if temp_align_score == min_dist:  #if new alignment has highest score
            min_dist_count += 1
            best_alignments.append(temp_align)
        
        else:
            if temp_align_score > min_dist:  #new best sw score
                min_dist = temp_align_score  
                min_dist_count = 1  #start at 1 alignment
                best_alignments = [temp_align]
            else:
                continue
    return min_dist, min_dist_count, best_alignments










