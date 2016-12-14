"""Contains classes solr_query, node and sisters."""
# -*- coding: utf-8 -*-

from urllib.request import urlopen
import simplejson
import json
import numpy as np
import pandas as pd
import contextlib
# import copy


class solr_query():
    """
    A solr_query class that stores URLs.

    Attbts:
    solr_url - the main solr_url
    """

    def __init__(self, solr_url, query):
        """Initialize the solr_query object."""
        self.solr_url = solr_url
        self.query = query

    def set_solr_url(self, url):
        """Assign a url to the solr_query object."""
        self.solr_url = url

    def add_query_url(self, url):
        """Add a query url to the solr_query object."""
        self.query = url

    def open_query(self):
        """
        Given a query, append it to the main url.

        Open URL and use simplejson to load the results
        """
        try:
            with contextlib.closing(urlopen(self.solr_url +
                                    self.query)) as conn:
                return simplejson.load(conn)
        except:
            raise Warning('URL is invalid or may have timed out')


class node():
    """
    A node is intended to be a single ontology term.

    Attributes:
    name - wbbt id
    parents
    genes
    similarity - no. of genes assoc. with this node divided
                 by the set of genes of its sister set
    drop -- whether to drop or not
    good_name -- human readable plus wbbt


    QUERIES FOR RELATIONS AND GENES ARE LAMBDA FUNCTIONS
    query_relation(x) -- gets families of tissue x
    query_genes(x) -- gets genes assoc with x
    query_readable
    """

    def __init__(self, name):
        """Initialization function."""
        self.name = name
        self.daughters = []
        self.parents = []
        self.genes = []
        self.similarity = 0
        self.drop = False
        self.good_name = ''

    def get_name(self, human_readable):
        """Generate a good name (human readable + WBid)."""
        if human_readable == '':
            print('warning, empty human readable name')
        self.good_name = human_readable + ' ' + self.name

    def add_daughter(self, daughter):
        """Add a daughter to this node."""
        self.daughters.append(daughter)
        self.daughters = list(set(self.daughters))  # prevent redundancy

    def add_parent(self, parent):
        """Add a parent to this node."""
        self.parents.append(parent)
        self.parents = list(set(self.parents))

    def add_annotation(self, gene):
        """Add annotation to this node."""
        self.genes.append(gene)
        self.genes = list(set(self.genes))

    def throw_away(self):
        """Set the `drop` variable to True."""
        self.drop = True

    def calc_similarity(self, sim):
        """Calculate similarity."""
        self.similarity = sim

    def find_family(self, solr_url, query_relation):
        """
        Find the family for this node by using solr_url and query_relation.

        query_relation(x) --lambda function
        """
        # get the json object
        rsp_rlshp = solr_query(solr_url, query_relation(self.name))

        # extract the array with all the right information
        array_of_rlshps = rsp_rlshp.open_query()['response']['docs'][0]

        # go through the array, turning each line into a dictionary
        # these mini-dictionaries contain the edges between nodes
        for j in json.loads(array_of_rlshps['topology_graph_json'])['edges']:
            # if the object isnt the same as the wbbt, object is parent to wbbt
            # if object is same as wbbt, wbbt is parent to subject
            if self.name != j['obj']:
                self.add_parent(j['obj'])
            else:
                self.add_daughter(j['sub'])

    def find_genes(self, solr_url, query_genes):
        """
        For a given wbbt, find the genes associated with it.

        query_genes(x) -- lambda function!
        """
        rsp_genes = solr_query(solr_url, query_genes(self.name))
        # extract the array with all the right information
        array_of_genes = rsp_genes.open_query()['response']['docs']
        # go through the array, turning each line into a dictionary
        for entry in array_of_genes:
            self.genes.append(entry['id'][3:])  # remove WB: from the string

        self.genes = list(set(self.genes))


class sisters(object):
    """
    A sister object  that contains related terms.

    A sister object that is meant to contain a set of terms that are related
    Sisters are defined as a set of nodes that share a single parent
    If a node is multiparent, it can have as many different sister sets as
    parents.

    Attributes:
    parent -- the parent for this set
    sisters -- set of `node` objects that are related by the same parent
    geneset -- total set of genes associated with these sisters
    threshold -- similarity threshold that specifies above which similarity
                 sisters must be killed
    dropsisters -- boolean
    dropped -- an array that keeps track of all sisters ever dropped
    """

    def __init__(self, parent, threshold):
        """Initialize function."""
        self.parent = parent
        self.sisters = []
        self.geneset = []
        self.threshold = threshold
        self.dropsisters = 0
        self.dropped = []

    def add_sister(self, sister):
        """Add a sister."""
        if self.sisters:
            self.sisters.append(sister)
        else:
            self.sisters = [sister]

        self.geneset = list(set(self.geneset+(sister.genes)))

    def add_sisters(self, sisters):
        """Add multiple sisters."""
        self.sisters = list(set(self.sisters+sisters))

        for sister in sisters:
            self.geneset = self.geneset+sister.genes
        self.geneset = list(set(self.geneset))

    def add_dropped(self, sister):
        """Add a sister to the `dropped` list."""
        if sister not in list:
            self.dropped.append(sister)

        else:
            self.dropped = self.dropped+sister

    def calc_similarity(self, method):
        """
        Calculate the family wise similarity for this object.

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

        avg = 0
        for sister in self.sisters:
            sim = len(sister.genes)/len(self.geneset)
            sister.calc_similarity(sim)

            if method == 'any':
                if sim > self.threshold:
                    self.dropsisters = 1
            avg += sim

        avg = avg/len(self.sisters)

        if method == 'avg':
            if avg > self.threshold:
                self.dropsisters = 1

    def kill(self):
        """If dropsister variable is 1, set `dropped` = 'sisters'."""
        if self.dropsisters == 1:
            self.dropped = self.sisters

    def trim(self, val):
        """If sister doesn't have `val` genes assoc. with it, drop it."""
        if len(self.sisters) == 0:
            return

        for sister in self.sisters:
            if len(sister.genes) < val:
                self.dropped.append(sister)


class ontology():
    """An ontological object."""

    def __init__(self, name, cutoff, threshold, method, solr_url):
        """Initialization function."""
        self.name = name
        self.threshold = threshold
        self.method = method
        self.nodes = {}
        self.family = {}
        self.solr_url = solr_url
        self.query_min_cutoff = 5
        self.cutoff = cutoff
        self.dropped = {}
        self.good = {}

    def set_min_cutoff(self, x):
        """Set minimum gene cutoff below which nodes are not fetched."""
        self.query_min_cutoff = x

    def add_nodes(self, query_terms, query_readable):
        """Add nodes from solr database."""
        sq = solr_query(self.solr_url, query_terms(self.query_min_cutoff))
        rsp_terms = sq.open_query()
        sd = solr_query(self.solr_url, query_readable)
        rsp_read = sd.open_query()

        i = 0
        for k in enumerate(rsp_terms['facet_counts']
                           ['facet_fields']['regulates_closure']):
            if i % 2 == 0:
                n = node(k[1])
                if n.name not in self.nodes:
                    self.nodes[n.name] = n
                    self.nodes[n.name].get_name(query_readable)
                if n.name not in self.family:
                    self.family[n.name] = sisters(n.name, self.threshold)
            i += 1

        for k, val in enumerate(rsp_read['response']['docs']):
            if val['id'] not in self.nodes:
                continue
            self.nodes[val['id']].get_name(val['annotation_class_label'])

    def find_node_family(self, lambda_query_rlshp):
        """Find the nodes that are related to this one."""
        for n in iter(self.nodes):
            self.nodes[n].find_family(self.solr_url, lambda_query_rlshp)

    def find_node_annotations(self, lambda_query_genes):
        """Fetch the annotations for this node."""
        for n in iter(self.nodes):
            self.nodes[n].find_genes(self.solr_url, lambda_query_genes)
            if len(self.nodes[n].genes) < self.cutoff:
                self.dropped[self.nodes[n].name] = self.nodes[n]

    def annotate_nodes(self, lambda_query_rlshp, lambda_query_genes):
        """Annotate this node with a family and with annotations."""
        self.find_node_family(lambda_query_rlshp)
        self.find_node_annotations(lambda_query_genes)

    def find_families(self):
        """Figure out the family structure for each node."""
        for node in self.nodes:
            n = self.nodes[node]
            for daughter in n.daughters:

                if daughter not in self.nodes:
                    continue
                # if 'WBbt:0002367' == daughter:
                #     print('hi')
                if len(self.nodes[daughter].genes) < self.threshold:

                    # add sister
                    self.family[n.name].add_sister(self.nodes[daughter])
                    # place it in sister.dropped
                    self.family[n.name].add_dropped(self.nodes[daughter])
                    # but also in self.dropped
                    self.dropped[n.name] = n

                else:
                    self.family[n.name].add_sister(self.nodes[daughter])

    def calculate_similarities(self):
        """Calculate the family-wise similarity."""
        for parent in self.family:
            self.family[parent].calc_similarity(self.method)

    def kill(self):
        """Remove whatever nodes fulfill the sisters.kill criterion."""
        for parent in self.family:
            self.family[parent].kill()

            for killed in self.family[parent].dropped:
                if killed.name in self.nodes:
                    self.dropped[killed.name] = killed

    def ceiling(self):
        """If a node has all its complement of daughters, kill it."""
        for parent in self.family:
            if parent not in self.nodes:
                continue

            if len(self.family[parent].sisters) == 0:
                continue

            if len(self.family[parent].dropped) == 0:
                self.dropped[self.nodes[parent].name] = self.nodes[parent]

    def find_good(self):
        """Fetch the surviving nodes."""
        for node in self.nodes:
            if node not in self.dropped:
                self.good[self.nodes[node].good_name] = self.nodes[node]


def build_dictionary(wbbts, tissue_array, genes):
    """Build the dictionary from a list of terms and wbbts."""
    # given a list of tissues, find the genes associated with each tissue and
    # place them in a vector.....

    mat = np.zeros(shape=(len(genes), len(wbbts)))
    d = {}
    for i, gene in enumerate(genes):
        d[gene] = i
        # for j, tissue in enumerate(wbbts):
        #     if gene in wbbts[tissue].genes:
        #         mat[i, j] = 1

    for j, tissue in enumerate(wbbts):
        for gene in wbbts[tissue].genes:
            mat[d[gene], j] = 1

    cols = tissue_array
    df = pd.DataFrame(mat, columns=cols)
    df.insert(0, 'wbid', genes)

    # drop the root term, for some reason it causes problems with hgt
    if 'C. elegans Cell and Anatomy WBbt:0000100' in df.columns:
        df.drop('C. elegans Cell and Anatomy WBbt:0000100', axis=1,
                inplace=True)
    return df

# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
#  # # # # # # #
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
# ==============================================================================
if __name__ == '__main__':
    # Raymond:
    # I have split up the URLs into 2 different variables to make life easier
    # solr_url contains the first part
    # query_xxx contains the second. However, query_xx can be a lambda function
    # basically, at a point in the string, I have written something like...
    # 'select?qt=standard&indent={0}'.format(x) -- the {0} is replaced by x
    # this allows me to modify the query in predictable ways.
    # hope this is clear.

    import argparse
    import sys

    parser = argparse.ArgumentParser(description='Run Dictionary Maker')
    parser.add_argument("threshold", help='The redundancy threshold',
                        type=float)
    parser.add_argument('cutoff', help='The annotation cutoff for each term',
                        type=int)
    parser.add_argument("fname",
                        help='Filename (complete with path) to save to',
                        type=str)
    parser.add_argument("-m", '--method',
                        help='method - defaults to \'any\' if not specified',
                        type=str)
    parser.add_argument("-mc", '--mincutoff',
                        help='The minimum cutoff to fetch. Defaults to 2.',
                        type=int)
    parser.add_argument("-su", '--solrurl',
                        help='The main body of the solr url.', type=str)
    parser.add_argument("-o", "--ontology",
                        help='One of `phenotype`, `tissue` or `gene`. Only\
                        works if --solrurl has not been specified',
                        type=str, default='anatomy',
                        choices=['anatomy', 'phenotype', 'go'])

    args = parser.parse_args()

    # main solr url
    if args.solrurl:
        solr_url = args.solrurl
    else:
        # solr_url = 'http://wobr.caltech.edu:8082/solr/anatomy/'
        s = 'http://wobr.caltech.edu:8082/solr/{0}/'
        solr_url = s.format(args.ontology)

    # queries must be lambda functions
    # query for terms. Finds terms that have x or more annotating genes
    def query_terms(x, ontology=args.ontology):
        """Search solr for terms (nodes) in the ontology."""
        if ontology != 'go':
            s = 'select?qt=standard&indent=on&wt=json&version=2.2&fl=' +\
                'id&start=0&rows=0&q=document_category:bioentity' +\
                '&facet=true&facet.field=regulates_closure&' +\
                'facet.limit=-1&facet.mincount={0}&facet.sort' +\
                '=count&fq=source:%22WB%22&fq=-qualifier:%22not%22'
            return s.format(x)
        else:
            s = 'select?qt=standard&indent=on&wt=json&version=2.2&fl=' +\
                'id&start=0&rows=1&q=document_category:bioentity&facet=' +\
                'true&facet.field=regulates_closure&facet.limit=-1&' +\
                'facet.mincount={0}&facet.sort=count&fq=source:%22WB' +\
                '%22&fq=taxon:%22NCBITaxon:6239%22&fq=-qualifier:%22not%22'
            return s.format(x)

    def query_relation(x, ontology=args.ontology):
        """
        query for relationships between nodes.

        given a wbbt ID `x`, find the nodes connected to it.
        Links are slightly different for [anatomy, phenotype] and GO, because
        in WormBase, the GO solr database includes all other worm species as
        well.
        """
        if ontology != 'go':
            s = "select?qt=standard&fl=topology_graph_json&" +\
                "version=2.2&wt=json&indent=on&rows=1&q=id:" +\
                "%22{0}%22&fq=document_category:%22ontology_class%22"
            return s.format(x)
        else:
            s = 'select?qt=standard&indent=on&wt=json&version=2.2&fl' +\
                '=id&start=0&rows=0&q=document_category:bioentity' +\
                '&facet=true&facet.field=regulates_closure&facet' +\
                '.limit=-1&facet.mincount={0}&facet.sort=count&' +\
                'fq=source:%22WB%22&fq=-qualifier:%22not%22'
            return s.format(x)

    def query_genes(x):
        """
        find the genes associated with every node.

        given a wbbt ID `x`, open URL that contains genes assoc. with it.
        """
        s = "select?qt=standard&indent=on&wt=json&version=2.2&" +\
            "fl=id&start=0&rows=10000&q=document_category:bioentity" +\
            "&fq=source:%22WB%22&fq=-qualifier:%22not%22&" +\
            "fq=regulates_closure:%22{0}%22"
        return s.format(x)

    # query for readable names
    query_readable = "select?qt=standard&fl=id,annotation_class_label" +\
                     "&version=2.2&wt=json&indent=on&rows=100000&q=id:" +\
                     "*&fq=document_category:ontology_class&" +\
                     "fq=-is_obsolete:true"

    queries = [query_terms, query_relation, query_genes, query_readable]
    threshold = args.threshold
    cutoff = args.cutoff

    if args.method:
        method = args.method
    else:
        method = 'any'

    if args.mincutoff:
        min_annot = args.mincutoff
    else:
        min_annot = 2

    trial1 = ontology('tissue_ontology', cutoff, threshold, method, solr_url)
    print('Object made')

    print('Min cutoff set at: {0}....'.format(min_annot))
    sys.stdout.flush()
    trial1.set_min_cutoff(min_annot)

    print('Fetching nodes.....')
    sys.stdout.flush()
    trial1.add_nodes(query_terms, query_readable)

    print('Annotating nodes')
    sys.stdout.flush()
    trial1.find_node_annotations(query_genes)

    print('Finding node families...')
    sys.stdout.flush()
    trial1.find_node_family(query_relation)

    print('Generating node family representation...')
    sys.stdout.flush()
    trial1.find_families()

    message = 'Calculating similarities and \
    removing nodes with more than {0:.2} similarity...'
    print(message.format(threshold))
    sys.stdout.flush()
    trial1.calculate_similarities()

    message = 'killing nodes that have less than {0} annotations...'
    print(message.format(cutoff))
    sys.stdout.flush()
    trial1.kill()

    print('Applying ceiling...')
    sys.stdout.flush()
    trial1.ceiling()

    print('Generating final list of terms...')
    trial1.find_good()
    print('No. of terms in dictionary: {0}'.format(len(trial1.good)))
    # extract keys

    print('Generating file at {0}'.format(args.fname))
    tissues = []
    genes = []

    for n in trial1.good:
        tissues.append(n)
        print(n)
        genes = genes+trial1.good[n].genes
    genes = list(set(genes))

    df = build_dictionary(trial1.good, tissues, genes)
    df.to_csv(args.fname, index=False)
