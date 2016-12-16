# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 14:24:04 2016

A script to generate 'clean' dictionaries that exclude
golden standards. 

This script 

@author: dangeles
"""
import dictionary_generator as dg
from dictionary_generator import ontology
import sys
x_array= [100, 50, 33, 25]
thresh_array= [1, .95, .9]
methods= ['avg', 'any']
path= '../input/WS252AnatomyDictionary/'


 #main solr url
solr_url = 'http://caprica.caltech.edu:8080/solr/anatomy/'
#queries must be lambda functions
#query for terms. Finds terms that have x or more annotating genes
query_terms= lambda x: 'select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.mincount={0}&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22'.format(x)
#query for relationships. given a wbbt ID, find the nodes connected to it. 
query_relation= lambda x: "select?qt=standard&fl=topology_graph_json&version=2.2&wt=json&indent=on&rows=1&q=id:%22{0}%22&fq=document_category:%22ontology_class%22".format(x)
#query for number of genes. given a wbbt ID, find genes assoc. with it. 
query_genes= lambda x: "select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=10000&q=document_category:bioentity&fq=source:%22WB%22&fq=-qualifier:%22not%22&fq=regulates_closure:%22{0}%22".format(x)
#query for readable names
query_readable= "select?qt=standard&fl=id,annotation_class_label&version=2.2&wt=json&indent=on&rows=100000&q=id:*&fq=document_category:ontology_class&fq=-is_obsolete:true"
queries= [query_terms, query_relation, query_genes, query_readable]
min_annot= 2

for x in x_array:
    for threshold in thresh_array:
        for method in methods:    
            print(x, threshold, method)
            sys.stdout.flush()
            print(solr_url)
            trial1= ontology('tissue_ontology', x, threshold, method, solr_url)
            trial1.set_min_cutoff(5)
            trial1.add_nodes(query_terms, query_readable)
            trial1.find_node_annotations(query_genes)
            trial1.find_node_family(query_relation)
            trial1.find_families()
            trial1.calculate_similarities()
            print(len(trial1.dropped))
            trial1.kill()
            print(len(trial1.dropped))
            trial1.ceiling()
            trial1.find_good()
            print(len(trial1.dropped))
            print('final {0}'.format(len(trial1.good)))
            #extract keys
            tissues= []
            genes= []
            for n in trial1.good:
                tissues.append(n)
                genes= genes+trial1.good[n].genes    
                genes= list(set(genes))
                
            df= dg.build_dictionary(trial1.good, tissues, genes)
            name= 'cutoff{0}_threshold{1}_method{2}.csv'.format(x, threshold, method)
            df.to_csv(path+name, index= False)
            


#generate a 'good' final dictionary not cleaned
#main solr url
solr_url = 'http://wobr.caltech.edu:8082/solr/anatomy/'
#queries must be lambda functions
#query for terms. Finds terms that have x or more annotating genes
query_terms= lambda x: 'select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=0&q=document_category:bioentity&facet=true&facet.field=regulates_closure&facet.limit=-1&facet.mincount={0}&facet.sort=count&fq=source:%22WB%22&fq=-qualifier:%22not%22'.format(x)
#query for relationships. given a wbbt ID, find the nodes connected to it. 
query_relation= lambda x: "select?qt=standard&fl=topology_graph_json&version=2.2&wt=json&indent=on&rows=1&q=id:%22{0}%22&fq=document_category:%22ontology_class%22".format(x)
#query for number of genes. given a wbbt ID, find genes assoc. with it. 
query_genes= lambda x: "select?qt=standard&indent=on&wt=json&version=2.2&fl=id&start=0&rows=10000&q=document_category:bioentity&fq=source:%22WB%22&fq=-qualifier:%22not%22&fq=regulates_closure:%22{0}%22".format(x)
#query for readable names
query_readable= "select?qt=standard&fl=id,annotation_class_label&version=2.2&wt=json&indent=on&rows=100000&q=id:*&fq=document_category:ontology_class&fq=-is_obsolete:true"
queries= [query_terms, query_relation, query_genes, query_readable]
min_annot= 2
x_array= [33]
thresh_array= [.95]
methods= ['any']
path= '../input/'

for x in x_array:
    for threshold in thresh_array:
        for method in methods:    
            print(x, threshold, method)
            sys.stdout.flush()
            print(solr_url)
            trial1= ontology('tissue_ontology', x, threshold, method, solr_url)
            trial1.set_min_cutoff(min_annot)
            trial1.add_nodes(query_terms, query_readable)
            trial1.find_node_annotations(query_genes)
            trial1.find_node_family(query_relation)
            trial1.find_families()
            trial1.calculate_similarities()
            print(len(trial1.dropped))
            trial1.kill()
            print(len(trial1.dropped))
            trial1.ceiling()
            trial1.find_good()
            print(len(trial1.dropped))
            print('final {0}'.format(len(trial1.good)))
            #extract keys
            tissues= []
            genes= []
            for n in trial1.good:
                tissues.append(n)
                genes= genes+trial1.good[n].genes    
                genes= list(set(genes))
                
            df= dg.build_dictionary(trial1.good, tissues, genes)
            name= 'final_cutoff{0}_threshold{1}_method{2}.csv'.format(x, threshold, method)
            df.to_csv(path+name, index= False)