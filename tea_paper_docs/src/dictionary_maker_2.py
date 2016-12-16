# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 14:39:03 2016

@author: dangeles
"""

from urllib.request import urlopen
import simplejson
import json
import numpy as np
import pandas as pd
import contextlib

class solr_query(object):
    """
    A solr_query class that stores URLs
    
    Attbts:
    solr_url - the main solr_url
    """
    def __init__(self, solr_url, query):
        self.solr_url= solr_url
        self.query= query
        
    def set_solr_url(self, url):
        self.solr_url= url
    def add_query_url(self, url):
        self.query= url
        
    def open_query(self):
        """
        Given a query, append it to the main url. Open URL and
        use simplejson to load the results
        """
        try:
            with contextlib.closing(urlopen(self.solr_url+ self.query)) as conn:
                return simplejson.load(conn)
        except:
            raise Warning('URL is invalid or may have timed out')
        

class node(object):
    """
    A node is intended to be a single ontology term
    
    Attributes:
    name - wbbt id
    parents
    genes
    similarity - no. of genes assoc. with this node divided by the set of genes of its sister set
    drop -- whether to drop or not
    good_name -- human readable plus wbbt
    
    
    QUERIES FOR RELATIONS AND GENES ARE LAMBDA FUNCTIONS
    query_relation(x) -- gets families of tissue x
    query_genes(x) -- gets genes assoc with x
    query_readable

    """
    def __init__(self, name):
        self.name= name
        self.daughters= []
        self.parents= []
        self.genes= []
        self.similarity= 0
        self.drop= False
        self.good_name= ''
    
    def get_name(self, human_readable):        
        self.good_name= human_readable + ' '+ self.name
        
    def add_daughter(self, daughter):
        self.daughters.append(daughter)
    
    def add_parent(self, parent):
        self.parents.append(parent)
    
    def add_annotation(self, gene):
        self.gemes.append(gene)
    
    def throw_away(self):
        self.drop= True
    
    def calc_similarity(self, sim):
        self.similarity= sim

    def find_family(self, solr_url, query_relation):
        """
        query_relation(x) --lambda function
        """
        #get the json object
        rsp_rlshp= solr_query(solr_url, query_relation(self.name))  
            
            
        #extract the array with all the right information
        array_of_rlshps= rsp_rlshp.open_query()['response']['docs'][0]
    
        #go through the array, turning each line into a dictionary
        #these mini-dictionaries contain the edges between nodes
        for j in json.loads(array_of_rlshps['topology_graph_json'])['edges']:        
            #if the object isn't the same as the wbbt, object is parent to wbbt        
            #if object is same as wbbt, wbbt is parent to subject
            if self.name != j['obj']:
                self.add_parent(j['obj'])
            else:
                self.add_daughter(j['sub'])
        
    def find_genes(self, solr_url, query_genes):
        """
        For a given wbbt, find the genes associated with it
        query_genes(x) -- lambda function!
        """
        rsp_genes= solr_query(solr_url, query_genes(self.name))
        #extract the array with all the right information
        array_of_genes= rsp_genes.open_query()['response']['docs']
        #go through the array, turning each line into a dictionary
        for entry in array_of_genes:
            self.genes.append(entry['id'][3:]) #remove WB: from the string
            
class sisters(object):
    """
    A sister object that is meant to contain a set of terms that are related
    Sisters are defined as a set of nodes that share a single parent
    If a node is multiparent, it can have as many different sister sets as parents
    
    Attributes:
    parent -- the parent for this set
    sisters -- set of \'node\' objects that are related by the same parent
    geneset -- total set of genes associated with these sisters
    threshold -- similarity threshold that specifies above which similarity sisters must be killed
    dropsisters -- boolean
    dropped -- an array that keeps track of all sisters ever dropped
    
    """
    def __init__(self, parent, threshold):
        self.parent= parent
        self.sisters= []
        self.geneset= []
        self.threshold= threshold
        self.dropsisters= 0
        self.dropped= []
        
    def add_sister(self, sister):
        if self.sisters:
            self.sisters.append(sister)
        else:
            self.sisters= [sister]
        
        self.geneset= list(set(self.geneset+(sister.genes)))
    
    def add_sisters(self, sisters):
        
        self.sisters= list(set(self.sisters+sisters))
        
        for sister in sisters:
            self.geneset= self.geneset+sister.genes
        self.geneset= list(set(self.geneset))
        
    def calc_similarity(self, method):
        """
        A method to calculate the similarity of a set of sisters to each other
        by finding the cardinality of the total gene set and the cardinality of
        the gene set for each node
        
        Depending on the method, the sisters.dropsisters value is modified if 
        the sisters are too similar to each other
        """
        if len(self.sisters) == 0:
            return 0
        if self.geneset == 0:
            return 1
        
        if method not in ['avg', 'any']:
            raise ValueError('method must be one of \'avg\' or \'any\'')
        
        avg= 0
        for sister in self.sisters:
            sim= len(sister.genes)/len(self.geneset)
            sister.calc_similarity(sim)
            
            if method == 'any':
                if sim > threshold:
                    self.dropsisters= 1
            avg+= sim
        
        avg= avg/len(self.sisters)
        
        if method == 'avg':
            if avg > threshold:
                self.dropsisters= 1
    
    def kill(self):
        if self.dropsisters == 1:
            self.dropped= self.sisters            
            self.sisters= []
            
    
    def trim(self, val):
        
        if len(self.sisters) == 0:
            return
            
        for sister in self.sisters:
            if len(sister.genes) < val:
                self.dropped.append(sister)
                self.sisters.pop(self.sisters.index(sister))
        













#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
# # # # # # # #         
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
#==============================================================================
if __name__ == '__main__':
    #Raymond:
    #I have split up the URLs into 2 different variables to make life easier
    #solr_url contains the first part
    #query_xxx contains the second. However, query_xxx can be a lambda function.
    #basically, at a point in the string, I have written something like...
    #'select?qt=standard&indent={0}'.format(x) -- the {0} is replaced by x
    #this allows me to modify the query in predictable ways.
    #hope this is clear.
    
    
    path= '/Users/dangeles/WormFiles/hgf_benchmarking/input/'
    #main solr url
    solr_url = 'http://wobr.caltech.edu:8082/solr/anatomy/';

    #queries must be lambda functions
    #query for terms. Finds terms that have x or more annotating genes
    query_terms= lambda x: 'select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.mincount={0}&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22'.format(x)
    
    #query for relationships. given a wbbt ID, find the nodes connected to it. 
    query_relation= lambda x: "select?qt=standard&fl=topology_graph_json&version=2.2&wt=json&indent=on&rows=1&q=id:%22{0}%22&fq=document_category:%22ontology_class%22".format(x)
    
    #query for number of genes. given a wbbt ID, find genes assoc. with it. 
    query_genes= lambda x: "select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=10000&q=document_category:bioentity&fq=source:%22WB%22&fq=-qualifier:%22not%22&fq=regulates_closure:%22{0}%22".format(x)
    
    #query for readable names
    query_readable= "select?qt=standard&fl=id,annotation_class_label&version=2.2&wt=json&indent=on&rows=100000&q=id:*&fq=document_category:ontology_class&fq=-is_obsolete:true"
    
    threshold= .3
    cutoff= 100
    rd= solr_query(solr_url, query_readable)
    readable= rd.open_query()
    method= 'any'
    min_annot= 5
    
    fname= path + 'tissue_dictionary_cutoff{0}_threshold{1}_method{2}.csv'.format(cutoff, threshold, method)

    
    if threshold > 1-1/cutoff:
        threshold= 1-1.5/cutoff

    #get all terms
    sq= solr_query(solr_url, query_terms(min_annot))
    rsp_terms= sq.open_query()

#==============================================================================
#   cutoff terms
#==============================================================================
    #find all terms with more than cutoff annotations
    wbbts= {}
    sister_dict= {}
    i= 0
    for k in enumerate(rsp_terms['facet_counts']['facet_fields']['regulates_closure']):
        if i%2 == 0:
            n= node(k[1])
            n.find_family(solr_url, query_relation)
            n.find_genes(solr_url, query_genes)
            wbbts[n.name]= n 
            sister_dict[n.name]= sisters(n.name, threshold)
        i+=1

    print('No. of tissues in wbbts {0}'.format(len(wbbts)))
                    
    k= 0
    for parent in sister_dict:
        #count the sisters, find their similarity and kill them if they don't pass
        if parent in wbbts:
            for daughter in wbbts[parent].daughters:
                if daughter in wbbts:
                    sister_dict[parent].add_sister(wbbts[daughter])
            
            sister_dict[parent].calc_similarity(method)
            sister_dict[parent].kill()
            sister_dict[parent].trim(cutoff)
    
        #remove all sisters that were killed or trimmed from wbbts
        for dropped in sister_dict[parent].dropped:
            if dropped.name in wbbts:
                del wbbts[dropped.name]
    print('No. of tissues in wbbts after kill and trim: {0}'.format(len(wbbts)))
        
#==============================================================================
#   reference terms
#==============================================================================
    #find all terms with more than min_annot genes and use it as a reference
    ref= {}
    ref_sisters={}
    i= 0
    for k in enumerate(rsp_terms['facet_counts']['facet_fields']['regulates_closure']):
        if i%2 == 0:
            n= node(k[1])
            n.find_family(solr_url, query_relation)
            ref[n.name]= n 
            ref_sisters[n.name]= sisters(n.name, 1)
        i+=1
    print('No. of tissues in reference {0}'.format(len(ref)))   
    
    #find sisters in ref
    k= 0
    for s in ref_sisters:
        #count the sisters, find their similarity and kill them if they don't pass
        if s in ref:
            for daughter in ref[s].daughters:
                if daughter in ref:
                    ref_sisters[s].add_sister(ref[daughter])
                    
#==============================================================================
#   check completeness of sisters and pop parents (ceiling)           
#==============================================================================
    to_pop= []
    for parent in wbbts:    
        if len(sister_dict[parent].sisters) == len(ref_sisters[parent].sisters):
            to_pop.append(parent)

    to_pop= list(set(to_pop))
    for p in to_pop:
        wbbts.pop(p)
    print('No. of tissues in wbbts after ceiling: {0}'.format(len(wbbts)))
#==============================================================================
#   get human readable names
#==============================================================================
    wbbts_listform= []
    genes= []
    for thingy in readable['response']['docs']:
        #annotate human readable
        if thingy['id'] in wbbts.keys():
            wbbts[thingy['id']].get_name(thingy['annotation_class_label'])
            
            genes= genes+wbbts[thingy['id']].genes
            wbbts_listform.append(wbbts[thingy['id']].good_name)

    genes= list(set(genes))

#==============================================================================
#   build dictionary
#==============================================================================
    def build_dictionary(wbbts, tissue_array, genes):
        #given a list of tissues, find the genes associated with each tissue and 
        #place them in a vector..... 
    
        mat= np.zeros(shape= (len(genes), len(wbbts)))
        for i, gene in enumerate(genes):
            for j, tissue in enumerate(wbbts):
                if gene in wbbts[tissue].genes:
                    mat[i, j]= 1

    
        cols= tissue_array
        df= pd.DataFrame(mat, columns= cols)
        df.insert(0, 'wbid', genes)
    
        #drop the root term, for some reason it causes problems with hgt
        if 'C. elegans Cell and Anatomy WBbt:0000100' in df.columns:
            df.drop('C. elegans Cell and Anatomy WBbt:0000100', axis= 1, inplace= True)    
    
        return df
        
    
    #tocsv
    df= build_dictionary(wbbts, wbbts_listform, genes)
    df.to_csv(fname, index= False)


        
        
        
       
       
       
       
















