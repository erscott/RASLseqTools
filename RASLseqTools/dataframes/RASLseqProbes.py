import pandas as pd
import itertools
import random
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
    def __init__(self, probe_file, aligner_dir, run_results_dir, aligner='star', verbose=False):

        self.probes = pd.read_table(probe_file, sep='\t')
        self.aligner_dir = aligner_dir
        self.run_results_dir = run_results_dir
        self.verbose = verbose

        if len(set(['AcceptorProbeSequence', 'DonorProbeSequence', 'AcceptorAdaptorSequence', 'DonorAdaptorSequence', 'ProbeName']) & set(self.probes.columns)) != 5:
                print '!!!!\n'*2
                print '!!!!  Probes file Missing Column(s) or Column Header(s) AcceptorProbeSequence, DonorProbeSequence, AcceptorAdaptorSequence, DonorAdaptorSequence, ProbeName  !!!!!"'
                print '!!!!\n'*2
                assert False

        self.probes['acceptor_no_adaptor'] = self.probes.apply(self.remove_adaptor, args=['AcceptorProbeSequence', 'AcceptorAdaptorSequence'], axis=1)
        self.probes['donor_no_adaptor'] = self.probes.apply(self.remove_adaptor, args=['DonorProbeSequence', 'DonorAdaptorSequence'], axis=1)
        self.probes['on_target_seq'] = self.probes['acceptor_no_adaptor'] + self.probes['donor_no_adaptor']
        self.probes['on_target_seq_adaptor'] = self.probes['AcceptorAdaptorSequence'] + self.probes['acceptor_no_adaptor'] + self.probes['donor_no_adaptor'] + self.probes['DonorAdaptorSequence']
        self.probes['on_target_name'] = self.probes['ProbeName'] + "///" + self.probes['ProbeName']
        self.on_off_target_probes_df = self._on_off_target_probes_df()

        self.probe_columns = list(self.on_off_target_probes_df.index)

        self.make_fasta()

        if aligner == 'star':
            '''
                This function creates a STAR probe database
            '''
            cmd = '{aligner_dir}/STAR --runMode genomeGenerate --genomeDir {probedb_path} --genomeFastaFiles {fasta} --genomeChrBinNbits 6 --genomeSAindexNbases 4 --outFileNamePrefix {star_log}'.format(
                aligner_dir=self.aligner_dir,
                probedb_path=self.run_results_dir,
                fasta=self.fasta,
                star_log=self.run_results_dir+'/star_'
            )


        elif aligner == 'blast':
            '''
                This function creates a blastdb for RASLseq fastq read matching

                The blastdb path is set by the user and can be accessed at the
                Attribute 'probedb_path'
            '''
            cmd = '{aligner_dir}/makeblastdb -in {fasta} -dbtype nucl -max_file_sz 1GB -out {probedb_file} -title {probedb_path}'.format(
                aligner_dir=self.aligner_dir,
                probedb_path=self.run_results_dir,
                fasta=self.fasta,
                probedb_file=self.fasta.rstrip('.fa')
            )

        else:
            # (TODO) Throw error
            pass

        if self.verbose:
            print "Make DB Command:", cmd

        os.system(cmd)

    def remove_adaptor(self, line, probe_col, adaptor_col):
        '''
            Convenience function to remove adaptor sequences
        '''
        probe = line[probe_col].replace('r', '').replace('/5Phos/', '')
        return probe.replace(line[adaptor_col], '')

    def _on_off_target_probes_df(self):
        '''
            This function returns a dataframe with all combinations of
            Acceptor and Donor probes (excluding on-target probe pairs)
        '''

        acc_adap = self.probes['AcceptorAdaptorSequence'].unique()[0]
        don_adap = self.probes['DonorAdaptorSequence'].unique()[0]

        # Constructing off_target list
        probe_cartesian_product = itertools.product(self.probes['acceptor_no_adaptor'], self.probes['donor_no_adaptor'])
        probe_name_cartesian_product = [i for i in itertools.product(self.probes['ProbeName'], self.probes['ProbeName'])]

        # Off_target probes only
        on_off_target_df_data = []
        for i, probes in enumerate(probe_cartesian_product):
            off_target_seq = ''.join(probes)
            on_off_target_df_data.append([probe_name_cartesian_product[i][0], probe_name_cartesian_product[i][1], off_target_seq, acc_adap, don_adap])

        on_off_target_df = pd.DataFrame.from_records(on_off_target_df_data)
        on_off_target_df.columns = ['acceptor_probe_id', 'donor_probe_id', 'probe_seq', 'AcceptorAdaptorSequence', 'DonorAdaptorSequence']
        on_off_target_df.index = on_off_target_df['acceptor_probe_id'] + '///' + on_off_target_df['donor_probe_id']
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

        # OPENING BLAST DATABASE
        self.fasta = self.run_results_dir + '/BLAST_DATABASE.fa'
        fasta_input = open(self.fasta, 'w')

        # Writing on_target probes
        counter = 0
        seq_tracker = set()  # Tracking on-target probe sequences
        for i in self.probes.index:  # Only blasting unique perfect matches
            acceptor = self.probes.ix[i]['acceptor_no_adaptor']
            donor = self.probes.ix[i]['donor_no_adaptor']
            if acceptor + donor in seq_tracker:
                continue  # No duplicate probes even if in probes file
            fasta_input.write(">" + self.probes.ix[i]['on_target_name'] + "\n" + acceptor + donor + "\n")
            fasta_input.flush()
            seq_tracker.add(acceptor + donor)
            counter += 1

        # Constructing off_target list
        probe_cartesian_product = itertools.product(self.probes['acceptor_no_adaptor'], self.probes['donor_no_adaptor'])
        probe_name_cartesian_product = [i for i in itertools.product(self.probes['ProbeName'], self.probes['ProbeName'])]

        # Writing off_target probes
        for i, probes in enumerate(probe_cartesian_product):
            if ''.join(probes) in seq_tracker:
                continue
            if probe_name_cartesian_product[i][0] == probe_name_cartesian_product[i][1]:
                continue  # on_target_probe

            fasta_input.write('>' + '///'.join(probe_name_cartesian_product[i]) + '!\n' + ''.join(probes) + '\n')
            fasta_input.flush()
            seq_tracker.add("".join(probes))
            counter += 1

        fasta_input.close()
        self.probe_count = counter

        return True

