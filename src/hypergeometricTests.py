 # -*- coding: utf-8 -*-
"""
Author: David Angeles
Date: May 26, 2015
A script to implement a hypergeometric test
Needs:
A tissue dictionary
A control list of gene names
An experimental list of gene names
"""

import os
import sys, getopt
os.environ['MPLCONFIGDIR']='/tmp'
# os.environ['HOME']='/tmp/'

import matplotlib
matplotlib.use('Agg') 

# from __future__ import division, print_function, absolute_import
import pandas as pd
from scipy import stats

# for arg in sys.argv:
#     print arg

q_threshold=0.1
# q_threshold= 0.018
# in the case of anatomical expression patterns, \
# there are ~5000 genes total in the dictionary whereas a tissue has \
# minimally 100 genes. Thus, if the input contains a singleton gene, \
# the p can be as low as 0.02 for tissues that does not express that \
# gene at all. It makes no sense to say that the tissue is enriched.


path= "./"
os.chdir(path)
if len(sys.argv) == 2:
    gene_file = sys.argv[1]
#    q_threashold = sys.argv[2]
else:
    exit("Must specify a gene list. Bye")
genes1= pd.read_csv(gene_file) #this file should be only wormbase ID's!
print "using file ", gene_file

#extract the first column of the file as a list of strings for analysis
gene_list1= genes1[genes1.columns[0]].values

#read in the dictionary
# tissue_df= pd.read_csv("../input/smalldictionary.txt")
tissue_df= pd.read_csv("/home/raymond/local/src/git/tissue_enrichment_tool_hypergeometric_test/input/dictionary_anatomy.csv")


#==============================================================================
# a function that checks which genes given by the user are in the dictionary. 
#==============================================================================
def pass_list(user_provided, tissue_dictionary):
    """
    A function to check which genes are in the 
    provided list of user provided names and to
    return a dataframe of presence or fail
    
    Takes two vectors, user_provided and tissue_dictionary in that order
    tissue_dictionary is a pandas dataframe that you should know about
    user_provided is a list of gene names
    """
    length= tissue_dictionary.shape[0] #how many genes in the dict
    
    #make an empty dataframe
    present= pd.DataFrame(index= range(length), columns= ['wbid','provided'])
    
    #fill wbid column with all the gene names
    present.wbid= tissue_dictionary.wbid
    
    #go through and pass attendance -- 1 if present, 0 otherwise
    for item in user_provided:
        present.provided[present.wbid== item]= 1
        
    #remove NA's and make them 0
    present.provided= present.provided.fillna(0)
    
    #return df
    return present
    
#==============================================================================
#hgf is short for hypergeometric function
#this function implements a test for enrichment (not depletion!!!)
#==============================================================================
def hgf(gene_list, tissue_dictionary):
    """
    Given a list of tissues and a gene-tissue dictionary,
    returns a p-dictionary for the enrichment of every tissue
    (a p-dictionary is a vector of length equal to the number
    of tissues in the tissue_dictionary, sorted by value)
    
    The entries of the p-vector are p-values not corrected
    for multiple hypothesis testing. 
    
    gene_list should be a list or list-like
    tissue_dictionary should be a pandas df
    """    
    
    #figure out what genes are in the user provided list
    present= pass_list(gene_list, tissue_dictionary)  
    #slice out only the genes that were present from the user-provided list    
    observed= present.wbid[present.provided==1]

    #re-index the dictionary s.t. the wbid is the index
    tissue_dictionary= tissue_dictionary.set_index('wbid')
    
    #figure out the number of labels for each tissue in tissue_df
    sums_of_tissues= tissue_dictionary.sum()[0:] #this object can be identified by 
    #the column names of tissues and excludes gene IDs

    total_genes= tissue_dictionary.shape[0] #total genes in the dictionary

    #slice out the rows from tissue_dictionary that came from the user-provided list
    observed_dictionary= tissue_dictionary.loc[observed]
    
    #get the total number of labels from each tissue
    observed_sum= observed_dictionary.sum()
    #get the total number of genes provided by the user that are in the dictionary
    total_observed= observed.shape[0]    
    
    #make a hash with the p-values for enrichment of each tissue. 
    p_hash= {}
    for i, name in enumerate(tissue_dictionary.columns.values): 
		if total_observed == 0:   # no testable gene in list at all
			p_hash[name]= 1
		elif observed_sum[name] >= total_observed:
			p_hash[name]= 0
		elif observed_sum[name] == 0:     # tissue has no queried gene
			p_hash[name]= 1
		else:
#			print(i, "--", name)
#			print(sums_of_tissues[name])
    
        # sf: survival function = (1 - CDF)
        # CDF(cumulative density function): the probability that a \
        # real-valued random variable X with a given probability \
        # distribution will be found to have a value less than or equal to x.
        # function params = x, M, n, N
        # x: Genes of focus tissue in user list (ball of type in sample)
        # M: Total genes in dictionary (pool)
        # n: Genes positive for focus tissue in dict (ball of type in pool)
        # N: Total (valid) genes in user list (sample)
        
        ##function params= k, K, n, N
        ##k is the number of observed number of genes associated with that tissue
        ##K is number of genes associated with the current tissue
        ##n is the length of the list provided by the user
        ##N is the total number of associations in the dictionary
        
			p_hash[name]= stats.hypergeom.sf(observed_sum[name],total_genes,\
							sums_of_tissues[name],total_observed)
#		print(name, p_hash[name], observed_sum[name],total_genes,\
#							sums_of_tissues[name],total_observed)
        
    #return the p-values. 
    return p_hash
    
    
#==============================================================================
#FDR correction
#==============================================================================
def benjamin_hochberg_stepup(p_vals):
    """
    Given a list of p-values, apply a BH
    FDR correction and return the pertaining q values
    
    alpha is a scalar FDR value
    p_vals should be iterable
    """
    #sort the p_values, but keep the index listed
    index= [i[0] for i in sorted(enumerate(p_vals), key=lambda x:x[1])]
    
    #keep the p_values sorted
    p_vals= sorted(p_vals)
    q_vals= [None]*len(p_vals) #initialize an empty list
    prev_q= 0
    
    #BH Step Up begins here. 
    for i, p in enumerate(p_vals):
        q= len(p_vals)/(i+1)*p #calculate the q_value for the current point
        q= min(q, 1) #if q >1, make it == 1
        q= max(q, prev_q) #preserve monotonicity, since at endpoints the procedure sometimes makes it so that prev_q > new_q
        q_vals[i]= q #store the q_value
        prev_q= q #update the previous q_value
        
    #return q_vals and the index so we can match up each 
    #q_value with its appropriate tissue. 
    return q_vals, index
    
    
#==============================================================================
# format the results from benjamin_hochberg_stepup and hgf
# and return any results that are have q_values below alpha
#==============================================================================
def return_enriched_tissues(p_hash, alpha):
    """
    Given a hash of p-values
    (tissue -> p-values)
    apply an FDR correction and
    return the q_values that pass
    significance along with the 
    tissue type. 
    """
    
    #initialize a list, a hash and a counter
    p_values= []
    index_to_tissue= {}
    k= 0
    #go through the p_hash
    for key, value in p_hash.iteritems():
        #place each value in the list
        p_values.append(value)
        #and the index of the value in the hash with value equal to the tissue
        index_to_tissue[k]= key
        #add 1 to the index
        k+=1
    
    #apply FDR, using BH stepup procedure
    q_values, index= benjamin_hochberg_stepup(p_values)
    
    #place everything in a hash
    q_hash= {}
    for i, value in enumerate(q_values):
        j= index_to_tissue[index[i]]
        q_hash[j]= value
    
    #print the results. This will likely be modified to return a 
    #file or some such.
    print("q-values less than alpha = {0:.2} are considered statistically significant".format(alpha))
    print("\n\n")
    print("------------------------")
    print 'tissue,adjusted p value (<{})'.format(q_threshold)
    for key, value in q_hash.iteritems():
        if value < alpha:
			# there may be value=0, integer
           if value == 0:
              print("{0},{1}".format(key, value))
           else:
              print("{0},{1:.2e}".format(key, value))
    print("------------------------\n\n")

        
    return q_hash
#==============================================================================
# call everything
#==============================================================================    

def implement_hypergmt_enrichment_tool(gene_list, tissue_df, alpha= 0.01):
    """
    Calls all the above functions
    
    gene_list: a list of non-redundant gene names
    tissue_df
    alpha: significance threshold, defaults to 0.01
    """

    print('Executing script\n')
    #calculate p_vals
    p_hash= hgf(gene_list, tissue_df)
    #calculate q_vals
    q_hash= return_enriched_tissues(p_hash, alpha)
    
    return q_hash
#==============================================================================
# end of functions    
#==============================================================================

#Run the whole thing:
q1= implement_hypergmt_enrichment_tool(gene_list1, tissue_df, alpha= q_threshold)
