# -*- coding: utf-8 -*-
"""
Spyder Editor


@david angeles
dangeles@caltech.edu
"""

import tissue_enrichment_analysis as tea #the library to be used
import pandas as pd
import os
import importlib as imp
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import matplotlib as mpl

sns.set_context('paper')

pd.set_option('display.float_format', lambda x:'%f'%x) 


#this script generates a few directories. 
dirOutput= '../output/'
dirSummaries= '../output/SummaryInformation/'
dirHGT25_any= '../output/HGT25_any_Results/'
dirHGT33_avg= '../output/HGT33_avg_Results/'
dirHGT33_any= '../output/HGT33_any_Results/'
dirHGT50_any= '../output/HGT50_any_Results/'
dirHGT100_any= '../output/HGT100_any_Results/'

DIRS= [dirOutput, dirSummaries, dirHGT25_any, dirHGT33_avg,
       dirHGT50_any, dirHGT100_any]
#open the relevant file
path_sets= '../input/genesets_golden/'
path_dicts= '../input/WS252AnatomyDictionary/'


#Make all the necessary dirs if they don't already exist
for d in DIRS:
    if not os.path.exists(d):
        os.makedirs(d)

#Make the file that will hold the summaries and make the columns. 
with open(dirSummaries+'ExecutiveSummary.csv', 'w') as fSum:
    fSum.write('#Summary of results from all benchmarks\n')
    fSum.write('NoAnnotations,Threshold,Method,EnrichmentSetUsed,TissuesTested,GenesSubmitted,TissuesReturned,GenesUsed,AvgFold,AvgQ,GenesInDict\n')

#==============================================================================
#==============================================================================
# # Perform the bulk of the analysis, run every single dictionary on every set
#==============================================================================
#==============================================================================
i= 0
#look in the dictionaries
for folder in os.walk(path_dicts):
    #open each one
    for f_dict in folder[2]:
        if f_dict == '.DS_Store':
            continue
        tissue_df= pd.read_csv(path_dicts+f_dict)
        
        #tobedropped when tissue dictionary is corrected
        annot, thresh= re.findall(r"[-+]?\d*\.\d+|\d+", f_dict)    
        annot= int(annot); thresh= float(thresh) #typecasting
        method= f_dict[-7:-4]
        
        ntiss= len(tissue_df.columns)
        ngenes= tissue_df.shape[0]
        
        #open each enrichment set
        for fodder in os.walk(path_sets):
            for f_set in fodder[2]:
                df= pd.read_csv(path_sets + f_set)
                test= df.gene.values
                ntest= len(test)
                short_name= f_set[16:len(f_set)-16]
                
                df_analysis, unused= \
                tea.enrichment_analysis(test, tissue_df, alpha= 0.05, show= False)
                
                nana= len(df_analysis) #len of results
                nun= len(unused) #number of genes dropped
                avf= df_analysis['Fold Change'].mean()
                avq= df_analysis['Q value'].mean()
                s= '{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10}'.format(
                annot, thresh, method, f_set, ntiss, ntest, nana, ntest-nun, avf,avq,ngenes)
                with open(dirSummaries+'ExecutiveSummary.csv', 'a+') as fSum:
                    fSum.write(s)
                    fSum.write('\n')
    
#Print summary to csv
df_summary= pd.read_csv(dirSummaries+'ExecutiveSummary.csv', comment= '#')

#some entries contain nulls. before I remove them, I can inspect them
df_summary.isnull().any()
indexFold = df_summary['AvgFold'].index[df_summary['AvgFold'].apply(np.isnan)]
indexQ = df_summary['AvgQ'].index[df_summary['AvgQ'].apply(np.isnan)]
df_summary.ix[indexFold[0]]
df_summary.ix[indexQ[5]]


#kill all nulls!
df_summary.dropna(inplace= True)
#calculate fraction of tissues that tested significant in each run
df_summary['fracTissues']= df_summary['TissuesReturned']/df_summary['TissuesTested']

df_summary.sort_values(['NoAnnotations', 'Threshold', 'Method'], inplace= True)
#==============================================================================
#==============================================================================
# # Plot summary graphs
#==============================================================================
#==============================================================================
sel= lambda x, y, z: (df_summary.NoAnnotations == x) & (df_summary.Threshold == y) & (df_summary.Method == z)


#KDE of the fraction of all tissues that tested significant
cols= ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e'] #used with varying colors
ls= ['-', '--', ':'] # used with varying thresh
thresh= df_summary.Threshold.unique()
NoAnnotations= df_summary.NoAnnotations.unique()

def resplot(column, method= 'any'):
    """
    A method to quickly plot all combinations of cutoffs, thresholds.
    All cutoffs are same color
    All Thresholds are same line style
    """
    for j, annots in enumerate(NoAnnotations):
        for i, threshold in enumerate(thresh):
            if threshold == 1:
                continue
            df_summary[sel(annots, threshold, method)][column].plot('kde', 
                            color= cols[j], ls= ls[i], lw= 4, 
                            label= 'Annotation Cut-off: {0}, Threshold: {1}'.format(annots, threshold))



resplot('fracTissues')
plt.xlabel('Fraction of all tissues that tested significant')
plt.xlim(0, 1)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend()
plt.savefig(dirSummaries+'fractissuesKDE_method=any.png')
plt.close()

resplot('AvgQ', method= 'avg')
plt.xlabel('Fraction of all tissues that tested significant')
plt.xlim(0, 0.05)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend()
plt.savefig(dirSummaries+'avgQKDE_method=avg.png')
plt.close()


resplot('AvgQ')
plt.xlabel('AvgQ value')
plt.xlim(0,.05)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend()
plt.savefig(dirSummaries+'avgQKDE_method=any.png')
plt.close()

#KDE of the fraction of avgFold
resplot('AvgFold')
plt.xlabel('Avg Fold Change value')
plt.xlim(0, 15)
plt.title('KDE Curves for all dictionaries, benchmarked on all gold standards')
plt.legend()
plt.savefig(dirSummaries+'avgFoldChangeKDE.png')
plt.close()



def line_prepender(filename, line):
    """
    Given a filename, opens it and prepends the line 'line' 
    at the beginning o fthe file
    """
    with open(filename, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)
        
        
#==============================================================================
#==============================================================================
# # Detailed analysis of 25 and 50 genes per node dictionaries
#==============================================================================
#==============================================================================
def walker(tissue_df, directory, save= True):
    
    with open(directory+'empty.txt', 'w') as f:
        f.write('Genesets with no enrichment:\n')
    for fodder in os.walk(path_sets):
        for f_set in fodder[2]:
            df= pd.read_csv(path_sets + f_set)
            
            short_name= f_set
            test= df.gene.values
            df_analysis, unused= tea.enrichment_analysis(test, tissue_df, show= False)
            
            if df_analysis.empty == False & save:
                df_analysis.to_csv(directory+short_name+'.csv')
                line= '#' + short_name+'\n'
                line_prepender(directory+short_name+'.csv', line)
                tea.plot_enrichment_results(df_analysis, title= short_name, dirGraphs= directory)
                plt.close()
            if df_analysis.empty:
                with open(directory+'empty.txt', 'a+') as f:
                    f.write(short_name+'\n')
               
               
tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/cutoff25_threshold0.95_methodany.csv')
walker(tissue_df, dirHGT25_any)
               
tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/cutoff50_threshold0.95_methodany.csv')
walker(tissue_df, dirHGT50_any)


tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/cutoff100_threshold0.95_methodany.csv')
walker(tissue_df, dirHGT100_any)


tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/cutoff33_threshold0.95_methodavg.csv')
walker(tissue_df, dirHGT33_any)
tissue_df= pd.read_csv('../input/WS252AnatomyDictionary/cutoff33_threshold0.95_methodavg.csv')
walker(tissue_df, dirHGT33_avg)

grouped= df_summary.groupby(['NoAnnotations','Threshold', 'Method'])
with open('../output/SummaryInformation/TissueNumbers.csv', 'w') as f:
    f.write('No. Of Annotations,Threshold,Method,No. Of Tissues in Dictionary\n')
    for key, group in grouped:
        f.write('{0},{1},{2},{3}\n'.format(key[0], key[1], key[2], group.TissuesTested.unique()[0]))





tissue_data= pd.read_csv('../output/SummaryInformation/TissueNumbers.csv')
sel= lambda y, z: (tissue_data.iloc[:, 1] == y) & (tissue_data.iloc[:, 2] == z)

#KDE of the fraction of all tissues that tested significant
cols= ['#1b9e77', '#d95f02', '#7570b3'] #used with varying colors
thresh= df_summary.Threshold.unique()
NoAnnotations= df_summary.NoAnnotations.unique()


def resplot(column, cutoff= 25, method= 'any'):
    """
    A method to quickly plot all combinations of cutoffs, thresholds.
    All cutoffs are same color
    All Thresholds are same line style
    """
    for i, threshold in enumerate(thresh):
        ax= plt.gca()
        ax.grid(False)
        if threshold == 1:
            continue
        tissue_data[sel(threshold, method)].plot(x= 'No. Of Annotations',
                    y= 'No. Of Tissues in Dictionary', kind= 'scatter', 
                    color= cols[i], ax= ax, s= 50, alpha= .7)
    ax.set_xlim(20, 110)
    ax.set_xscale('log')
    ax.set_xticks([25, 33, 50, 100])
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    ax.set_ylim(25, 1000)
    ax.set_yscale('log')
    ax.set_yticks([50, 100, 250, 500])
    ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())


resplot('No. Of Tissues in Dictionary')