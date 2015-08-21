
import os, sys
import pandas as pd
import Levenshtein as editdist
import itertools
import numpy as np 
from scipy.spatial import distance

from RASLseqTools import RASLseqAlign
import matplotlib
import pylab as plt


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
        
            self.well_annot_df.set_index(['PlateBarcode','WellBarcode'],inplace=True)
            
            self.annot_columns = self.well_annot_df.columns

        except AttributeError:
            print "!!!!"
            print "!!!!"
            print '!!!!  Well Barcode Annotation File Missing Column(s) or Column Header(s) WellBarcode or PlateBarcode  !!!!!'
            print "!!!!"
            print "!!!!"
            assert False
    
    
    
    def get_pairwise_barcode_distances(self, barcodes='well'):
        
        '''
        This function creates a pandas df containing Levenshtein
        edit distances between all pairs of well or plate barcodes
        
        Parameters
        ------------
        barcodes: str, 'well' or 'plate'
            default: 'well' 
            
        Returns
        ------------
        pandas DataFrame
        
        if barcodes=='well':
            self.wellbc_levenshtein_dist
        
        if barcodes=='plate'
            self.platebc_levenshtein_dist
        '''
        
        if barcodes == 'well':
            #pairwise well barcode edit distances
            self.wellbc_levenshtein_dist = RASLseqAlign.pairwise_distances(self.well_bc)
        
        else:
            #pairwise plate barcode edit distances
            self.platebc_levenshtein_dist = RASLseqAlign.pairwise_distances(self.plate_bc)
            
        return
        
    
    
    def plot_distances(self, barcodes='well'):
        '''
        This function plotes the pairwise edit distances
        between barcodes, masking the diagonal
        
        Parameters
        ------------
        barcodes: str, 'well' or 'plate'
            default: 'well' 
            
        Returns
        ------------
        fig, ax
            Matplotlib Subplot Object 
                heatmap of pairwise edit distances
        
        '''
        
        
        if barcodes == 'well':
            barcode_plt_name = 'WellBarcode'
            try:
                plot_df = self.wellbc_levenshtein_dist
            
            except AttributeError:
                self.get_pairwise_barcode_distances('well')
                plot_df = self.wellbc_levenshtein_dist
                
        else:
            barcode_plt_name = 'PlateBarcode'
            try:
                plot_df = self.platebc_levenshtein_dist
            
            except AttributeError:
                self.get_pairwise_barcode_distances('plate')
                plot_df = self.platebc_levenshtein_dist      
        
        

        #PLOTTING
        matplotlib.rcParams['xtick.labelsize']='xx-large'
        matplotlib.rcParams['ytick.labelsize']='xx-large'
        
        fig1, ax1 = plt.subplots(1, 1, dpi=300, figsize=(8,5))
        
        mask = np.diag(np.ones(len(plot_df)))  #masking diagonal
        masked_data = np.ma.masked_array(plot_df, mask)  #masked diagonal array
        matplotlib.cm.spectral.set_bad(color='white', alpha=None)  #setting color map
        heatmap = ax1.pcolor(masked_data, cmap='spectral')  #plotting heatmap
        ax1.invert_yaxis()
        ax1.set_aspect('equal')
        ax1.axis('tight')
        ax1.set_title(barcode_plt_name + ' Levenshtein Edit Distances', fontsize=20)
        ax1.set_xlabel(barcode_plt_name, fontsize=20)
        ax1.set_ylabel(barcode_plt_name, fontsize=20)
        cbar = plt.colorbar(heatmap, ax=ax1, pad=0.005)
        plt.show()
        plt.close()
                
        
        return fig1, ax1

    
    
    def get_min_edist_counts(self, barcodes='well'):

        '''
        Calculates the number of barcodes possessing the 
        minimum Levenshtein edit distance.
        
        Parameters
        ------------
        barcodes: str, 'well' or 'plate'
            default: 'well' 
            specifies which barcode distance dataframe to use
            
        '''
        
        if barcodes == 'well':
            barcode_name = 'WellBarcode'
            try:
                barcode_dist_df = self.wellbc_levenshtein_dist
            
            except AttributeError:
                self.get_pairwise_barcode_distances('well')
                barcode_dist_df = self.wellbc_levenshtein_dist
                
        else:
            barcode_name = 'PlateBarcode'
            try:
                barcode_dist_df = self.platebc_levenshtein_dist
            
            except AttributeError:
                self.get_pairwise_barcode_distances('plate')
                barcode_dist_df = self.platebc_levenshtein_dist    
        
        
        
        np.fill_diagonal(barcode_dist_df.values, 100)  #setting self-alignments to 100
        
        barcode_min_edist_count = []
        for i in barcode_dist_df.index:
            minimum = barcode_dist_df.ix[i].min()
            barcode_min_edist_count.append([ i, minimum, barcode_dist_df.ix[i].value_counts().ix[minimum] ] )
        barcode_min_edist_count = pd.DataFrame(barcode_min_edist_count, columns=[barcode_name, 'min_edist', 'min_edist_count'])
        barcode_min_edist_count.set_index(barcode_name, inplace=True)
        
        self.barcode_min_edist_count = barcode_min_edist_count
        
        return 
    
    
    