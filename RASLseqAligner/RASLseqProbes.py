
import itertools
import random
import pandas as pd
import os

class RASLseqProbes(object):
    '''
    This class creates a pandas DataFrame for RASLseq Probes and
    offers functions to create a fasta file and blast database
    using the cartesian product of the Donor and Acceptor probes.
    
    
    Parameters
    ----------
    probe_file: str
        Specifies the path to the tab seperated RASLseq probe file
            format: AcceptorProbeSequence DonorProbeSequence AcceptorAdaptorSequence DonorAdaptorSequence ProbeName
    
    blastdb_path: str
        Specifices the path for writing the blastdb
            e.g. /path/to/blastdb.fa
        The resulting blastdb will can be found at /path/to/blastdb

    blastn_dir: str
        Specifies path to the directory holding the makeblastdb executable
        
    
    Attributes
    ----------
    self.probes: pandas DataFrame
        columns = [AcceptorProbeSequence, DonorProbeSequence, AcceptorAdaptorSequence,
            DonorAdaptorSequence, acceptor_no_adaptor, donor_no_adaptor, on_target_seq,
            on_target_name, ProbeName]
    
    self.blastdb_path: str
        Specifices the path for writing the blastdb
            e.g. /path/to/blastdb.fa
        The resulting blastdb will can be found at /path/to/blastdb
    
    self.blastn_dir: str
        Specifies path to the makeblastdb executable
    
    '''
    
    import itertools
    
    
    def remove_adaptor(self, line, probe_col, adaptor_col):
        '''
        Convenience function to remove adaptor sequences
        '''
        probe = line[probe_col].replace('r','').replace('/5Phos/','')
        return probe.replace(line[adaptor_col], '')
    
    
    
    def _on_off_target_probes_df(self):
        '''
        This function returns a dataframe with all combinations of 
        Acceptor and Donor probes (excluding on-target probe pairs)
        '''
        
        acc_adap = self.probes['AcceptorAdaptorSequence'].unique()[0]
        don_adap = self.probes['DonorAdaptorSequence'].unique()[0]
        
        #constructing off_target list
        probe_cartesian_product = itertools.product(self.probes['acceptor_no_adaptor'], self.probes['donor_no_adaptor'])
        probe_name_cartesian_product = [i for i in itertools.product(self.probes['ProbeName'], self.probes['ProbeName'])]   
        
        
        #off_target probes only
        on_off_target_df_data = []
        for i,probes in enumerate(probe_cartesian_product):
            
            #if probe_name_cartesian_product[i][0] == probe_name_cartesian_product[i][1]: #removing on-targets
                #continue
            
            off_target_seq = "".join(probes)
            
            on_off_target_df_data.append([probe_name_cartesian_product[i][0], probe_name_cartesian_product[i][1], off_target_seq, acc_adap, don_adap])
        
        on_off_target_df = pd.DataFrame.from_records(on_off_target_df_data)
        on_off_target_df.columns = ['acceptor_probe_id', 'donor_probe_id', 'off_target_seq', 'AcceptorAdaptorSequence', 'DonorAdaptorSequence']
        on_off_target_df.index = on_off_target_df['acceptor_probe_id'] + "///" + on_off_target_df['donor_probe_id']
        return on_off_target_df
    
    
    
    def make_fasta(self):
        '''
        This function writes a fasta file using the cartesian product of the 
        Acceptor and Donor Probe Name/Sequences (adaptor sequences removed)
            e.g. 
            >DonorName///AcceptorName
            DonorSeqAcceptorSeq
        
        The fasta file can be found at the user-specified blastdb_path with a .fa suffix
        
        '''
        
        #OPENING BLAST DATABASE
        random_file_handle = self.blastdb_path + str(random.randrange(0,1000000)) + ".fa"
        self.fasta = random_file_handle
        blast_fasta_input = open(self.fasta,"w")
        
        #writing on_target probes
        seq_tracker = set()  #tracking on-target probe sequences
        for i in self.probes.index:  #only blasting unique perfect matches
            acceptor = self.probes.ix[i]['acceptor_no_adaptor']
            donor = self.probes.ix[i]['donor_no_adaptor']
            blast_fasta_input.write(">"+self.probes.ix[i]['on_target_name'] + "\n" + acceptor+donor + "\n")
            blast_fasta_input.flush()
            seq_tracker.add(acceptor+donor)
        
        #constructing off_target list
        probe_cartesian_product = itertools.product(self.probes['acceptor_no_adaptor'], self.probes['donor_no_adaptor'])
        probe_name_cartesian_product = [i for i in itertools.product(self.probes['ProbeName'], self.probes['ProbeName'])]   
        
        #writing off_target probes
        for i,probes in enumerate(probe_cartesian_product):
            if "".join(probes) in seq_tracker: continue
            if probe_name_cartesian_product[i][0] == probe_name_cartesian_product[i][1]: 
                continue  #on_target_probe
            
            blast_fasta_input.write(">"+"///".join(probe_name_cartesian_product[i])+"!\n" + "".join(probes) + "\n")    
            blast_fasta_input.flush()
            seq_tracker.add("".join(probes))
        
        blast_fasta_input.close()
        

        return True



    def make_blastdb(self):
        '''
        This function creates a blastdb for RASLseq fastq read matching
        
        The blastdb path is set by the user and can be accessed at the
        Attribute 'blastdb_path'
        
        '''
        self.blastdb_file = self.fasta.rstrip(".fa")
        
        cmd = [self.blastn_dir+'/makeblastdb', '-in', self.fasta, '-dbtype', 'nucl', '-max_file_sz', '1GB', '-out', self.blastdb_file, '-title', self.blastdb_path]    
        
        if self.print_on:
            print "maskblastdb Command:"
            print " ".join(cmd)
            print "blastdb stats"
        
        blastdb = os.system(" ".join(cmd))   
        
        
        return blastdb

        
    
    
    def __init__(self, probe_file, blastdb_path, blastn_dir, print_on=False):
        
        
        
        self.probes = pd.read_table(probe_file, sep="\t")
        
        
        if len(set(['AcceptorProbeSequence', 'DonorProbeSequence', 'AcceptorAdaptorSequence', \
                    'DonorAdaptorSequence', 'ProbeName']) & set(self.probes.columns)) != 5:
            print "!!!!"
            print "!!!!"
            print '!!!!  Probes file Missing Column(s) or Column Header(s) AcceptorProbeSequence, DonorProbeSequence, AcceptorAdaptorSequence, DonorAdaptorSequence, ProbeName  !!!!!"'
            print "!!!!"
            print "!!!!"
            assert False
        
        self.probes['acceptor_no_adaptor'] = self.probes.apply(self.remove_adaptor, args=['AcceptorProbeSequence','AcceptorAdaptorSequence'], axis=1)
        
        self.probes['donor_no_adaptor'] = self.probes.apply(self.remove_adaptor, args=['DonorProbeSequence','DonorAdaptorSequence'], axis=1)
        
        self.probes['on_target_seq'] = self.probes['acceptor_no_adaptor'] + self.probes['donor_no_adaptor']
        
        self.probes['on_target_seq_adaptor'] = self.probes['AcceptorAdaptorSequence'] + self.probes['acceptor_no_adaptor'] + self.probes['donor_no_adaptor'] + self.probes['DonorAdaptorSequence']
        
        self.probes['on_target_name'] = self.probes['ProbeName'] +"///"+self.probes['ProbeName']
        
        self.on_off_target_probes_df = self._on_off_target_probes_df()
        
        self.print_on = print_on
        
        self.blastdb_path = blastdb_path
        
        self.blastn_dir = blastn_dir
        
        self.make_fasta()
        
        self.make_blastdb()
 

   


