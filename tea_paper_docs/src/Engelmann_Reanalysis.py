# -*- coding: utf-8 -*-
"""
Created on Mon Feb 15 14:07:04 2016
A script to analyze some of the genes identified by Engelmann et al 2011,
A Comprehensive Analysis of Gene Expression Changes Provoked by Bacterial
and Fungal Infection in C elegans.

@author: dangeles
"""
import tissue_enrichment_analysis as tea  # the library to be used
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# tissue dictionary to use
tissue_df = pd.read_csv('../input/final_cutoff33_threshold0.95_methodany.csv')

# rna-seq datasets to reanalyze
#Bacterial
dfLum = pd.read_csv('../input/Engelmann/luminescens_Engelmann_2011.csv')
dfMarc = pd.read_csv('../input/Engelmann/marcescens_Engelmann_2011.csv')
dfFaec = pd.read_csv('../input/Engelmann/faecalis_Engelmann_2011.csv')

# Fungal
dfHarp = pd.read_csv('../input/Engelmann/harposporium_Engelmann_2011.csv')
dfCon = pd.read_csv('../input/Engelmann/coniospora_Engelmann_2011.csv')
# dfSin= pd.read_csv('../input/Engelmann/sinha_2012.csv')

# dictionary of gene names
names = pd.read_csv('../input/Engelmann/c_elegans.PRJNA13758.WS241.livegeneIDs.unmaprm.txt',
                    sep='\t', comment='#')

# make the directories to place the analyses
dirEngelmann = '../output/Engelmann/'
dirAnalysis = '../output/Engelmann/Analysis'
dirGraphs = '../output/Engelmann/Graphs'
DIRS = [dirEngelmann, dirAnalysis, dirGraphs]

# Make all the necessary dirs if they don't already exist
for d in DIRS:
    if not os.path.exists(d):
        os.makedirs(d)


# Selector functions to draw WBIDs from
f = lambda x, y: (names.HumanReadable.isin(x[y]))
g = lambda x, y: (names.GeneName.isin(x[y]))

# ==============================================================================
# ==============================================================================
# # Engelmann Analysis
# ==============================================================================
# ==============================================================================


class receptacle(object):

    def __init__(self, name):
        self.name = name
        self.result_dict = {}
        self.n = 0

    def add_result(self, name, df):
        self.result_dict[name] = df
        self.n += 1


Lnames = ['Otorhabdus luminescens', 'Serratia marcescens',
          'Enterococcus faecalis', 'Harposporium sp', 'Drechmeria coniospora']
Ldf = [dfLum, dfMarc, dfFaec, dfHarp, dfCon]
Ldirection = ['Infection_upregulated', 'Infection_downregulated']

# direction specific
Lreceptacles = {}
n_genes = []
for i, df in enumerate(Ldf):
    fname = Lnames[i]
    obj = receptacle(fname)

    for direction in Ldirection:
        ind = g(df[df[direction] == 1.0], 'SequenceNameGene')
        x = names[ind].WBID

        print('---------')
        print(fname + ' ' + direction)
        print('Number of genes submitted for analysis ', len(x))
        y = tissue_df[tissue_df.wbid.isin(x)].wbid.unique().shape[0]
        print('Number of genes used for analysis ', y)
        print('\n')

        df_res, unused = tea.enrichment_analysis(x, tissue_df,
                                                 show=False, alpha=0.1)

        print(df_res.empty)
        if len(df_res) > 0:
            obj.add_result(direction, df_res)

        n_genes.append(['{0}, {1}'.format(fname, direction), y, len(unused)])

    Lreceptacles[fname] = obj


# make
for species in Lreceptacles:
    current = Lreceptacles[species]

    n = current.n  # number of dfs
    keys = current.result_dict.keys()

    if n == 0:
        next

    i = 0
    if n > 1:
        fig, ax = plt.subplots(nrows=n, figsize=(8, 8))
        fig.subplots_adjust(top=2)
        fig.suptitle(species, fontsize=15, y=1.02)

        tea.plot_enrichment_results(current.result_dict['Infection_downregulated'],
                                    title='Name', save=False, fig=fig, ax=ax[0])

        # suppress xlabel
        ax[0].set_xlabel('')
        ax[0].set_ylabel('Down-Regulated Tissues')
        ax[0].yaxis.set_label_position('right')
        tea.plot_enrichment_results(current.result_dict['Infection_upregulated'],
                                    title='Name', save=False, fig=fig, ax=ax[1])
        ax[1].set_ylabel('Up-Regulated Tissues')
        xlabel = ax[1].set_xlabel('Enrichment Fold Change - {0}'.format(species))
        ax[1].yaxis.set_label_position('right')
        fig.tight_layout()
        plt.savefig('../output/Engelmann/Graphs/'+species+'Enrichment.pdf',
                    rect=[0, 0.03, 1, 0.95], bbox_extra_artists=[xlabel],
                    bbox_inches='tight')
#        plt.close()

    else:
        fig, ax = plt.subplots(nrows=n, figsize=(8, 4))
        fig.subplots_adjust(top=0.85)
        fig.suptitle(species, fontsize=15, y=1.05)

        if 'Infection_downregulated' in current.result_dict:
            ax = tea.plot_enrichment_results(current.result_dict['Infection_downregulated'],
                                             title='Name', save=False, fig=fig,
                                             ax=ax)
            ax.set_ylabel('Down-Regulated Tissues')
        else:
            ax = tea.plot_enrichment_results(current.result_dict['Infection_upregulated'],
                                             title='Name', save=False, fig=fig,
                                             ax=ax)
            ax.set_ylabel('Up-Regulated Tissues')

        xlabel = ax.set_xlabel('Enrichment Fold Change - {0}'.format(species))
        ax.yaxis.set_label_position('right')
        plt.tight_layout()
        plt.savefig('../output/Engelmann/Graphs/'+species+'Enrichment.pdf',
                    rect=[0, 0.03, 1, 0.95], bbox_extra_artists=[xlabel],
                    bbox_inches='tight')
    i += 1
