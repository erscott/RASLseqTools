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
    
    probedb_path: str
    Specifices the path for writing the blastdb
    e.g. /path/to/blastdb.fa
    The resulting blastdb will can be found at /path/to/blastdb
    
    aligner_dir: str
    Specifies path to the directory holding the makeblastdb executable
    
    
    Attributes
    ----------
    self.probes: pandas DataFrame
    columns = [AcceptorProbeSequence, DonorProbeSequence, AcceptorAdaptorSequence,
    DonorAdaptorSequence, acceptor_no_adaptor, donor_no_adaptor, on_target_seq,
    on_target_name, ProbeName]
    
    self.probedb_path: str
    Specifices the path for writing the blastdb
    e.g. /path/to/blastdb.fa
    The resulting blastdb will can be found at /path/to/blastdb
    
    self.aligner_dir: str
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
        on_off_target_df.columns = ['acceptor_probe_id', 'donor_probe_id', 'probe_seq', 'AcceptorAdaptorSequence', 'DonorAdaptorSequence']
        on_off_target_df.index = on_off_target_df['acceptor_probe_id'] + "///" + on_off_target_df['donor_probe_id']
        return on_off_target_df
    
    
    
    def make_fasta(self):
        '''
        This function writes a fasta file using the cartesian product of the
        Acceptor and Donor Probe Name/Sequences (adaptor sequences removed)
        e.g.
        >DonorName///AcceptorName
        DonorSeqAcceptorSeq
        
        The fasta file can be found at the user-specified probedb_path with a .fa suffix
        
        '''
        
        #OPENING BLAST DATABASE
        
        random_file_handle = self.probedb_path + self.random_str + ".fa"
        
        self.fasta = random_file_handle
        fasta_input = open(self.fasta,"w")
        
        #writing on_target probes
        counter = 0
        seq_tracker = set()  #tracking on-target probe sequences
        for i in self.probes.index:  #only blasting unique perfect matches
            acceptor = self.probes.ix[i]['acceptor_no_adaptor']
            donor = self.probes.ix[i]['donor_no_adaptor']
            if acceptor+donor in seq_tracker: continue #no duplicate probes even if in probes file
            fasta_input.write(">"+self.probes.ix[i]['on_target_name'] + "\n" + acceptor+donor + "\n")
            fasta_input.flush()
            seq_tracker.add(acceptor+donor)
            counter += 1
        
        #constructing off_target list
        probe_cartesian_product = itertools.product(self.probes['acceptor_no_adaptor'], self.probes['donor_no_adaptor'])
        probe_name_cartesian_product = [i for i in itertools.product(self.probes['ProbeName'], self.probes['ProbeName'])]
        
        #writing off_target probes
        
        for i,probes in enumerate(probe_cartesian_product):
            if "".join(probes) in seq_tracker: continue
            if probe_name_cartesian_product[i][0] == probe_name_cartesian_product[i][1]:
                continue  #on_target_probe
            
            fasta_input.write(">"+"///".join(probe_name_cartesian_product[i])+"!\n" + "".join(probes) + "\n")
            fasta_input.flush()
            seq_tracker.add("".join(probes))
            counter += 1
        
        fasta_input.close()
        self.probe_count = counter

        return True



    def make_blastdb(self):
        '''
        This function creates a blastdb for RASLseq fastq read matching
        
        The blastdb path is set by the user and can be accessed at the
        Attribute 'probedb_path'
        
        '''
        self.probedb_file = self.fasta.rstrip(".fa")

                
        cmd = [self.aligner_dir.rstrip('/')+'/makeblastdb', '-in', self.fasta, '-dbtype', 'nucl', '-max_file_sz', '1GB', '-out', self.probedb_file, '-title', self.probedb_path]
                    
        if self.print_on:
            print "maskblastdb Command:"
            print " ".join(cmd)
            print "blastdb stats"
                        
        blastdb = os.system(" ".join(cmd))
                            
        self.blastdb_cmd = ' '.join(cmd)           
        return blastdb
    
    

    def make_stardb(self):
        '''
        This function creates a STAR probe database
        '''
        #make genome file:/Users/ers_vader/Dropbox/RASLseq/Bcell_exp2/ipynb/data/STAR_genome/
        #~/Tools/STAR_2.3.0e.OSX_x86_64/STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ./on_off_probes.fa --genomeChrBinNbits 6 --genomeSAindexNbases 4
        
        os.system('mkdir ' + self.probedb_path.rstrip('/') + '/' + self.random_str)
        self.probedb_file = self.probedb_path.rstrip('/') + '/' + self.random_str
        cmd = [self.aligner_dir+'/STAR', '--runMode genomeGenerate', '--genomeDir', self.probedb_path, '--genomeFastaFiles', self.fasta, '--genomeChrBinNbits 6', '--genomeSAindexNbases 4' ]
        stardb = os.system(' '.join(cmd))
                
        if self.print_on:
            print 'make stardb command:'
            print ' '.join(cmd)
                    
        return stardb



    def __init__(self, probe_file, probedb_path, aligner_dir, aligner='star', print_on=False):
        
        
        
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
        
        self.random_str = str(random.randrange(0,1000000))
        
        self.aligner_dir = aligner_dir
        
        self.probe_columns = list( self.on_off_target_probes_df.index )
        
        
        if aligner == 'blast':
            self.probedb_path = probedb_path
            self.blastn_dir = aligner_dir
            self.make_fasta()
            self.make_blastdb()
            
    
        if aligner == 'star':
            
            os.system('mkdir ' + probedb_path.rstrip('/') + '/' + self.random_str + '/')
            
            self.probedb_path = probedb_path.rstrip('/') + '/' + self.random_str + '/'
            self.make_fasta()
            self.make_stardb()


