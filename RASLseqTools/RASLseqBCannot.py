import pandas as pd


class RASLseqBCannot(object):
    '''
    This class creates a pandas DataFrame containing
    RASLseq annotations for unique combinations of
    RASLseq plate and well barcodes.

    Parameters
    ----------
    plate_well_annot_file: str
        path to tab separated annotation file for plate, well barcodes
        file format: plate_bc well_bc annotation1 annotationN...

    Attributes
    ----------
    well_annot_df: pandas DataFrame
        index: (plate_bc, well_bc)
        columns: [annotations]

    '''

    def __init__(self, plate_well_annot_file):

        self.well_annot_df = pd.read_table(plate_well_annot_file, sep="\t")

        try:

            self.well_bc = set(self.well_annot_df.WellBarcode)

            self.plate_bc = set(self.well_annot_df.PlateBarcode)

            self.well_annot_df.set_index(['PlateBarcode', 'WellBarcode'], inplace=True)

            self.annot_columns = self.well_annot_df.columns

        except AttributeError:
            print "!!!!"
            print "!!!!"
            print '!!!!  Well Barcode Annotation File Missing Column(s) or Column Header(s) WellBarcode or PlateBarcode  !!!!!'
            print "!!!!"
            print "!!!!"
            assert False
