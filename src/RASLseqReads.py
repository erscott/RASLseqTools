

import FastqCollapse



class RASLseqReads(object):
    '''
    This class creates a pandas DataFrame for RASLseq fastq sequences.
    
    Attributes of this class annotate RASLseq fastq sequences.
    
    The get_attributes method allows the user to specify what additional
    information to add to the pandas DataFrame.
    
    
    Attributes
    ----------
    
    fastq_file: str 
        path to fastq_file

    
    print_on: boolean, default=False
        Whether to print information during data processing
    
    read_df: pandas DataFrame
        index: (PlateBarcode, sequence)
        columns: ['PlateBarcode', 'seq', 'seq_count']
            PlateBarcode - Index read from fastq file
            seq - fastq read sequence
            seq_count - number of occurrences of seq in fastq
    
    
    APRIL 14 2014: ToDos:  
        2) Refactor external methods (raslblast vs STAR) 
        3) Package everything
        4) Consider GUI
    
    '''

    def __init__(self, fastq_path, fastq_header, print_on=False):
        #create a probe with name name
        self.fastq_file = fastq_path
        
        self.print_on = print_on
        
        self.read_df = FastqCollapse.get_collapsed_fastq_df(self.fastq_file, fastq_header, print_on= self.print_on)
        
        #self.read_df = collapsed_read_df(self.fastq_file, print_on = self.print_on)
        
        


