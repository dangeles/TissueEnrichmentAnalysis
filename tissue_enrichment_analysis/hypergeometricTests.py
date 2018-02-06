"""
A script to implement a hypergeometric test procedure.

Author: David Angeles
Date: May 26, 2015
Requires Python > 3.5
Needs:
A tissue dictionary
A control list of gene names
An experimental list of gene names
"""
# -*- coding: utf-8 -*-
import pandas as pd
from scipy import stats
import numpy as np
import os
import sys
from urllib.request import urlopen
import contextlib


def pass_list(user_provided, tissue_dictionary):
    """
    A function to check which genes in provided list are in the dictionary.
    """
    ind = tissue_dictionary.wbid.isin(user_provided)
    present = tissue_dictionary[ind].wbid
    return present


def hgf(gene_list, tissue_df):
    """
    Given a list, returns the p-value for each tissue tested.

    Given a list of tissues and a gene-tissue dictionary,
    returns a p-dictionary for the enrichment of every tissue
    (a p-dictionary is a vector of length equal to the number
    of tissues in the tissue_dictionary, sorted by value).

    The entries of the p-vector are p-values not corrected
    for multiple hypothesis testing.
    gene_list should be a list or list-like
    tissue_dictionary should be a pandas df
    """
    # figure out what genes are in the user provided list
    # wanted = pass_list(gene_list, tissue_df)

    # re-index the dictionary s.t. the wbid is the index
    tissue_df = tissue_df.set_index('wbid')

    # number of balls per tissue in the dictionary
    sums_of_tissues = tissue_df.sum()

    # total labels in the dictionary
    total_balls = tissue_df.sum().sum()

    # slice out the rows from tissue_dictionary that came from the list
    wanted_dictionary = tissue_df.loc[gene_list]

    # get the total number of labels from each tissue
    wanted_sum = wanted_dictionary.sum()

    # get the total number of balls provided by the user that are in dictionary
    picked = wanted_sum.sum()

    # make a hash with the p-values for enrichment of each tissue.
    p_hash = {}
    exp_hash = {}
    for i, name in enumerate(tissue_df.columns.values):
        # if the total number of genes is zero, return p= 1 for all tissues
        if picked == 0:
            p_hash[name] = 1
            continue
        # if a certain tissue has never been called, don't test it
        if wanted_sum[name] == 0:
            p_hash[name] = 1
            continue
        # no. of balls of color name picked
        # total number of balls in urn
        # total number of balls of color name in urn
        # total number of balls picked out
        n_obs = wanted_sum[name]
        s_tissue = sums_of_tissues[name]
        p_hash[name] = stats.hypergeom.sf(n_obs, total_balls, s_tissue,
                                          picked)
        exp_hash[name] = stats.hypergeom.mean(total_balls, s_tissue,
                                              picked)

    # return the p-values, the genes associated with each tissue and the user
    # provided genes associate with each tissue.
    return p_hash, exp_hash, wanted_dictionary


def benjamini_hochberg_stepup(p_vals):
    """
    Given a list of p-values, apply FDR correction and return the q values.
    """
    # sort the p_values, but keep the index listed
    index = [i[0] for i in sorted(enumerate(p_vals), key=lambda x:x[1])]

    # keep the p_values sorted
    p_vals = sorted(p_vals)
    q_vals = [None]*len(p_vals)  # initialize an empty list
    prev_q = 0

    # BH Step Up begins here.
    for i, p in enumerate(p_vals):
        q = len(p_vals)/(i+1)*p  # calculate the q_value for the current point
        q = min(q, 1)  # if q >1, make it == 1
        q = max(q, prev_q)  # preserve monotonicity
        q_vals[i] = q  # store the q_value
        prev_q = q  # update the previous q_value

    # prevent the lowest q value from going to zero
    if np.sum(q_vals == 0) > 0:
        # set the min q-value to 10x less than the smallest non-zero value
        q_vals[np.where(q_vals == 0)] = np.min(q_vals[np.where(q_vals != 0)])/10

    # return q_vals and the index so we can match up each q-value to its index
    return q_vals, index


def return_enriched_tissues(p_hash, alpha):
    """
    A function index p-values and call the FDR function.
    """
    # initialize a list, a hash and a counter
    p_values = list(p_hash.values())
    keys = list(p_hash.keys())

    # FDR
    q_values, index = benjamini_hochberg_stepup(p_values)
    q_hash = {keys[index[pair[0]]]: pair[1] for pair in enumerate(q_values)}
    return q_hash


def enrichment_analysis(gene_list, tissue_df, alpha=0.05, aname='', save=False,
                        show=False):
    """
    Execute complete enrichment analysis (hypergeometric test, BH correction).

    ------
    Params:
    gene_list: a list of non-redundant WBIDs
    tissue_df: as provided by WormBase (use fetch_dictionary)
    alpha: significance threshold, defaults to 0.05
    aname= filename to use to save results
    show= Whether to print results or not.

    -------
    output:
    df_final - the final dataframe containing enriched tissues
    """
    if show:
        print('Executing script\n')

    # always make iterable
    if type(gene_list) in [str]:
        gene_list = [gene_list]

    if len(gene_list) == 0:
        raise ValueError('gene_list is empty!')

    # calculate the enrichment
    p_hash, exp_hash, wanted_dic = hgf(gene_list, tissue_df)

    # FDR correct
    q_hash = return_enriched_tissues(p_hash, alpha)

    # TODO: is there a better way to do this?
    def get_observed(x):
        """A function to find the number of observations of a tissue x."""
        return wanted_dic[x].sum()

    # slight modification
    # make a dataframe, index will be tissues column
    df_final = pd.DataFrame.from_dict(exp_hash, orient='index')
    # make the tissues their own column:
    df_final.reset_index(level=0, inplace=True)
    df_final.columns = ['Term', 'Expected']
    df_final['Observed'] = df_final.Term.apply(get_observed)  # v. slow
    df_final['Enrichment Fold Change'] = df_final.Observed/df_final.Expected
    df_final['P value'] = df_final.Term.map(p_hash)
    df_final['Q value'] = df_final.Term.map(q_hash)

    df_final.dropna(inplace=True)
    df_final.Observed = df_final.Observed.astype(int)

    df_final.sort_values('Q value', inplace=True)

    df_final = df_final[df_final['Q value'] < alpha]

    if show:
        if len(df_final) == 0:
            print('Analysis returned no enriched tissues.')
        else:
            print(df_final)  # print statement for raymond

    if save:
        df_final.to_csv(aname)

    return df_final


def plot_enrichment_results(df, y='logq', title='', analysis='tissue',
                            n_bars=15, save=False, **kwargs):
    """
    A plot function for TEA.

    df: dataframe as output by implement_hypergmt_enrichment_tool
    y: One of 'Fold Change', 'Q value' or a user generated column
    title - Title for the graph, also file name
    analysis - one of `tissue`, `phenotype` or `go`
    n_bars: number of bars to be shown, defaults to 15
    dirGraps: directory to save figures to. if not existent,
    generates a new folder

    ------
    output:
    ax - an axis object holding the graph that was generated
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    # sns.choose_colorbrewer_palette('sequential', as_cmap=False)
    if df.empty:
        print('dataframe is empty!')
        return

    if analysis.lower() not in ['tissue', 'phenotype', 'go']:
        raise ValueError('analysis variable must be one of' +
                         '`tissue`, `phenotype` or `go`')

    analysis = analysis.lower()

    ax = kwargs.pop('ax', None)
    ftype = kwargs.pop('ftype', 'svg')
#
    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 8))

    # sort by q value change
    df.sort_values(['Q value', 'Enrichment Fold Change'],
                   ascending=[True, False], inplace=True)

    # make a logq bar:
    logq = -df['Q value'].apply(np.log10)

    # added August 26 2016:
    tissue_ID = 11
    pheno_ID = 19
    go_ID = 10

    if analysis == 'phenotype':
        yvals = df.Term.str[:-pheno_ID-1]
    elif analysis == 'tissue':
        yvals = df.Term.str[:-tissue_ID-1]
    elif analysis == 'go':
        yvals = df.Term.str[:-go_ID-1]

    # plot first n_bars
    with sns.axes_style('whitegrid'):
        ax = sns.barplot(x=logq[:n_bars], y=yvals[:n_bars], ax=ax)

    # fix the plot to prettify it
    ax.set_ylabel('Terms', fontsize=15)
    if y.lower() != 'logq':
        ax.set_xlabel(y, fontsize=15)
    else:
        ax.set_xlabel('$\log_{10}{q}$', fontsize=15)
    ax.tick_params(axis='x', labelsize=13)
    ax.tick_params(axis='y', labelsize=13)
    plt.tight_layout()

    # save
    if save:
        plt.savefig('{0}.{1}'.format(title, ftype), dpi=1200)

    return ax


def fetch_dictionary(analysis='tissue'):
    """
    Fetch the dictionary we want.

    If analysis isn't specified, fetches the tissue dictionary.

    Params:
    ------
    analysis - one of `tissue`, `phenotype` or `go`

    Output:
    data - a dataframe containing the dictionary of interest
    """
    analysis = analysis.lower()
    if analysis not in ['tissue', 'phenotype', 'go']:
        raise ValueError('analysis must be one of `tissue`, `phenotype`' +
                         ' or `go`')

    url_tissue = 'http://caltech.wormbase.org/TissueEnrichmentAnalysis/'

    if analysis == 'tissue':
        url_tissue += 'anatomy_dict.csv'
    elif analysis == 'phenotype':
        url_tissue += 'phenotype_dict.csv'
    elif analysis == 'go':
        url_tissue += 'go_dict.csv'

    try:
        with contextlib.closing(urlopen(url_tissue)) as conn:
            data = pd.read_csv(conn)
            return data
    except:
        print('Cannot fetch dictionary. Please check internet connection.')


# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================

if __name__ == '__main__':

    import re
    import argparse
    import matplotlib
    matplotlib.use('Agg')
    import seaborn as sns

    sns.set_context('paper')
    sns.set_style('whitegrid')

    path = './'
    os.chdir(path)

    defQ = 0.1

    parser = argparse.ArgumentParser(description='Run EA.')
    parser = argparse.ArgumentParser()
    parser.add_argument("gene_list",
                        help='The full path to the gene list (WBIDs) you would\
                         like to analyse in .csv format')
    parser.add_argument('title', help='Title for your analysis (shouldn\'t\
                        include file extension)',
                        type=str)
    parser.add_argument('kind', help='What kind of analysis will be ' +
                        'performed. One of `tissue`, `phenotype` or `go`',
                        type=str)
    parser.add_argument("-d", '--dictionary', nargs='?', help='Provide a\
                        dictionary to test. If none given, WormBase URL \
                        will be used to download the corresponding file')
    parser.add_argument("-q", help='Qvalue threshold for significance. \
                        Default is {0} if not provided'.format(defQ),
                        type=float)
    parser.add_argument('-p', '--print', help='Indicate whether you would like \
                        to print results', action='store_true')
    parser.add_argument('-s', "--save", help='Indicate whether to save your \
                        plot.', action='store_true')
    args = parser.parse_args()

    gl_name = args.gene_list
    title = args.title

    # optional args
    if args.dictionary:
        dict_name = args.tissue_dictionary
        dictionary = pd.read_csv(dict_name)
    else:
        dictionary = fetch_dictionary(analysis=args.kind)

    if args.q:
        q = args.q
    else:
        q = defQ

    if args.print:
        prnt = True
    else:
        prnt = False

    if args.save:
        save = True

    else:
        save = False

    gene_list = []
    with open(gl_name, 'r') as f:
        gene_list = [x.strip() for x in f.readlines()]

    df_results = enrichment_analysis(gene_list, dictionary, alpha=q,
                                     show=False)

    dfname = title + '.tsv'
    df_results.to_csv(dfname, index=False, sep='\t')

    if prnt:
        with open(dfname, 'r') as f:
            printer = f.readlines()

        for value in printer:
            value = value.split(',')
            for val in value:
                if re.findall("\d+\.\d+", val):
                    ind = value.index(val)
                    x = float(val)
                    # if x < 10**-64:
                    #     value[ind] = '<10^-65'
                    # else:
                    value[ind] = '{0:.2g}'.format(x)

            value = '\t'.join(value)
            print(value)

    if save:
        plot_enrichment_results(df_results, title=title, save=save,
                                analysis=args.kind)

    sys.exit()
