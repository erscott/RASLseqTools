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

    '''

    def __init__(self, fastq_path, fastq_header, print_on=False):
        # Create a probe with name name
        self.fastq_file = fastq_path

        self.fastq_header = fastq_header

        self.print_on = print_on

    def get_blast_read_df(self):

        return FastqCollapse.get_collapsed_fastq_df(self.fastq_file, self.fastq_header, print_on=self.print_on)

