
import pandas as pd
import os, sys, time, random
import numpy as np
from scipy import stats

#NORMALIZATION FUNCTIONS
def LibSize_norm(df, normalization_probes=[]):
    '''
    Expects sample (rows) by gene (cols) pandas DataFrame
    '''
    
    if normalization_probes==[]:
        normalization_probes = df.columns
        
    return df.divide(df[normalization_probes].apply(np.sum, axis=1), axis=0) * 10000


def hkgene_norm(df, hk_gene):
    
    '''
    Housekeeping Gene Normalization
    
    
    '''

    return df.divide(df[hk_gene].apply(np.mean, axis=1), axis=0)



def MGMR_norm(df, normalization_probes=[], by_plate=False):
    '''
    Median Geometric Mean Ratio (MGMR)
    
    Expects sample by probe count dataframe, no annotations
    
    Anders Hubel 2010: http://www.genomebiology.com/content/pdf/gb-2010-11-10-r106.pdf
    
    scale factor for sample j = median (for each gene i: ( kij / geometric mean gene i ) ); where kij is sample j gene i       
    
    R code:
    data = gene by sample matrix 
    lib.size <- colSums(x)
    RLE = .calcFactorRLE(x)/lib.size
    
    .calcFactorRLE <- function (data)
    {
    gm <- exp(rowMeans(log(data)))  #geometric mean of log gene counts
    
    apply(data, 2, function(u) median((u/gm)[gm > 0]))  #divide each gene count by geometric mean for gene and take median value
    }

    '''
    def gmean(array_vals):
        return stats.gmean([i for i in array_vals if i > 0])
    
    def get_median(array_vals):
        return np.median([i for i in array_vals if i > 0])
    
    
    if len(normalization_probes) == 0:
        normalization_probes = df.columns
        
        
    if by_plate:
        return pd.concat([MGMR_norm(n, normalization_probes) for i,n in df.groupby(level=0)])
        
    
    #Removing wells with zero counts
    df_nonzero = df[df.apply(np.sum, axis=1)>100] #minimum of 100 reads
    #df_sample_sums = df_sample_sums[df_sample_sums>100].index 
    #df_nonzero = df[df.index.isin(df_sample_sums)]
    
    
    #calculating the pseudo-reference probe count using the geometric mean
    gmean_genei = df_nonzero[normalization_probes].apply(gmean, axis=0)  #geometric mean for each probe
    
    
    #calculating ratio of probe count to pseudo-reference probe count (geo mean)
    df_ratio = df_nonzero[normalization_probes].divide(gmean_genei, axis=1)  #counts divide by geomean by probe

    
    #calculating the scale factor by taking the median across probe ratios for each sample
    size_factor = df_ratio[normalization_probes].apply(get_median, axis=1)  #pick the median gene as the size factor for each sample
    
    
    #returning the scaled counts
    return df_nonzero.divide(size_factor, axis=0)   #dividing counts by size_factor
     

def quantile_norm(df):
    df = df.copy()
    df_rank = df.T.rank(method='min', axis=0)
    
    
    data = []
    for i in df.index:
        df_sorted = df.ix[i].values
        df_sorted.sort()
        data.append(df_sorted)
        
    quantile_values = pd.DataFrame(data).T.apply(np.mean, axis=1).to_dict() 
    
    quantile_df = df_rank.replace(quantile_values)
    
    return quantile_df


#STANDARDIZATION FUNCTIONS

def standardize_df(df, by_plate=False):
    
    if by_plate==False:
        df_min = df.apply(np.min, axis=0)
        df_max = df.apply(np.max, axis=0)
        return df.subtract(df_min, axis=1).divide(df_max, axis=1)
    
    return pd.concat([standardize_df(n) for i,n in df.groupby(level=0)])



def mean_center(values):
    '''
    input vector of values
    '''
    mean_val = np.mean(values)
    return [i-mean_val for i in values]



def mean_center_by_plate(df):
    '''
    input df
    '''
    data = []
    for i,n in  df.groupby(level=0):
        data.append(n.apply(mean_center, axis=0))
        
    return pd.concat(data)    



def std_scale(df, by_plate=False, level=0, with_mean=True, with_std=True):
    '''
    Convenience function to standard scale by pandas index level
    
    '''
    from sklearn.preprocessing import StandardScaler
    
    if by_plate:
        return pd.concat([pd.DataFrame(StandardScaler(with_mean=with_mean, with_std=with_std).fit_transform(n), columns=n.columns, index=n.index)  for i,n in df.groupby(level=level)])    
    
    return pd.DataFrame(StandardScaler(with_mean=with_mean, with_std=with_std).fit_transform(df), columns=df.columns, index=df.index)


