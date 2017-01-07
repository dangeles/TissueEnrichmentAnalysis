# -*- coding: utf-8 -*-
"""Benchmarking PEA."""
import tissue_enrichment_analysis as tea
import pandas as pd
import matplotlib.pyplot as plt
import os

# focus diseases for the paper. Only do certain analyses on these:
diseases = ['Systemic lupus erythematosus', 'Obesity-related traits',
            'Rheumatoid arthritis']

# import dictionaries
phenotype_df = pd.read_csv('../input/phenotype_ontology.csv')
go_df = pd.read_csv('../input/go_dictionary.csv')
tissue_df = tea.fetch_dictionary()

# figure out the shapes of each thing and print
print("Tissue Dictionary:")
print("Genes: {0}\nTerms: {1}".format(tissue_df.shape[0], tissue_df.shape[1]))
print("\nPhenotype Dictionary:")
print("Genes: {0}\nTerms: {1}".format(phenotype_df.shape[0],
                                      phenotype_df.shape[1]))
print("\nGO Dictionary:")
print("Genes: {0}\nTerms: {1}".format(go_df.shape[0], go_df.shape[1]))

# TODO: fix this path:
# from the gwas catalog, identify diseases with more than 100 associated genes
gwas_catalog = '~/Downloads/gwas_catalog_v1.0-associations_e87_r2016-12-12.tsv'

# SNP/GWAS data
# snp = pd.read_csv('../input/known_snps.csv', sep=',', comment='#')
# load the homolog data for the SNP/GWAS data
homologs = pd.read_csv('../input/gwas_homologs_worm.txt', sep='\t')

# all wormbase IDs
wbids = pd.read_csv('../input/wormbase_ids.txt', sep='\t', header=None,
                    names=['wbid', 'human_read', 'target_id'])

# Ciliary Neuron Transcriptome (from Juan Wang and Maureen Barr)
ciliated_df = pd.read_excel('~/Downloads/mmc3.xlsx')

# the minimum number of observations to keep for the GWAS enrich. analysis
# this ensures that only high-scoring results are kept.
n_min_obs = 4

# melt the dictionaries and save the melted dataframes. These melted
# dataframes will be useful later
pheno_traits = pd.melt(phenotype_df, id_vars='wbid', var_name='term')
pheno_traits = pheno_traits[pheno_traits.value == 1]
go_traits = pd.melt(go_df, id_vars='wbid', var_name='term')
go_traits = go_traits[go_traits.value == 1]
tissue_traits = pd.melt(tissue_df, id_vars='wbid', var_name='term')
tissue_traits = tissue_traits[tissue_traits.value == 1]


# define a function to run e.a.s quickly
def run(path, f, output_path, dictionary, column='gene', **kwargs):
    """Run and save e.a.."""
    genes = pd.read_csv(path+f, **kwargs)
    df = tea.enrichment_analysis(genes[column], dictionary,
                                 show=False)
    df = df[df.Observed > 2]
    df.to_csv(output_path + f + '.csv', index=False)

###############################################################################
###############################################################################
# run the benchmark analyses
print('\nRunning benchmarks...')
for f in os.listdir('../input/genesets_golden'):
    if f != '':
        run('../input/genesets_golden/', f,
            '../output/pea_goldenset_results/', phenotype_df)
        run('../input/genesets_golden/', f,
            '../output/go_goldenset_results/', go_df)
        run('../input/genesets_golden/', f,
            '../output/tea_goldenset_results/', tissue_df)

# run the fog-2 analysis
print('\nAnalyzing fog-2 data...')
for f in os.listdir('../input/fog2'):
    if f != '':
        run('../input/fog2/', f, '../output/fog2/pea_', phenotype_df,
            column='ens_gene')
        run('../input/fog2/', f, '../output/fog2/go_', go_df,
            column='ens_gene')
        run('../input/fog2/', f, '../output/fog2/tea_', tissue_df,
            column='ens_gene')

# analyze the hypoxia datasets
print('\nAnalyzing hypoxia datasets...')
for f in os.listdir('../input/hypoxia'):
    if f != '':
        run('../input/hypoxia/', f, '../output/hypoxia/pea_',
            phenotype_df, column='ens_gene')
        run('../input/hypoxia/', f, '../output/hypoxia/tea_',
            tissue_df, column='ens_gene')
        run('../input/hypoxia/', f, '../output/hypoxia/go_',
            go_df, column='ens_gene')


###############################################################################
###############################################################################
all_snps = pd.read_csv(gwas_catalog, sep='\t', encoding='ISO-8859-1')
array_of_genes = []

# find all the genes in the gwas catalog, and place the set of genes into
# a single list file so you can use DIOPT to find homologs easily
for g in all_snps['REPORTED GENE(S)'].values:
    if type(g) is not str:
        continue
    list_of_genes = [x.strip() for x in g.split(',')]
    for gi in list_of_genes:
        array_of_genes += [gi]

# save the array of genes to a file for translation into DIOPT
array_of_genes = list(set(array_of_genes))
with open('gene_homologs.csv', 'w+') as F:
    for g in array_of_genes:
        F.write(g+'\n')

# find all the diseases with more than 200 associated genes and save the
# list of traits that meet criterion
print('\nCounting gene-disease associations:')
candidates = []
for trait in all_snps.TRAIT.unique():
    if all_snps[all_snps.TRAIT == trait].shape[0] > 200:
        print(trait, all_snps[all_snps.TRAIT == trait].shape[0])
        candidates += [trait]

# what are the pubmedids associated with the list of traits?
# print('studies associated with candidates:')
# for trait in candidates:
#     print(trait)
#     print(all_snps[all_snps.TRAIT == trait].PUBMEDID.unique())

# for each GWAS disease, find the worm homologs with moderate or high scores
# and use those to run e.a. to identify phenologues of diseases in the worm
sel2 = (homologs.Rank == 'moderate') | (homologs.Rank == 'high')


def get_n(df, gene_list):
    """
    Find the number of genes from gene_list that appear at least once in df.

    Given a dataframe and a gene list, find how many genes in the list
    are represented at least once in the dataframe
    """
    return len(df[df.wbid.isin(gene_list)].wbid.unique())

print('\nBeginning Phenolog Search...')

for trait in candidates:
    # get genes associated with trait:
    g = all_snps[all_snps.TRAIT == trait]['REPORTED GENE(S)'].values
    current = []
    for i, gi in enumerate(g):
        if type(gi) is not str:
            continue
        list_of_genes = [x.strip() for x in gi.split(',')]
        for x in list_of_genes:
            current += [x]

    current = list(set(current))
    sel1 = (homologs.HID.isin(current))
    worm_genes = homologs[sel1 & sel2].WormBaseID

    # only run E.A. if there are enough homologs:
    if worm_genes.shape[0] > 100:
        print(trait, 'has', worm_genes.shape[0],
              ' homolog candidates in the worm')

        # certain analyses only perform on focus disease:
        if trait in diseases:
            print('Extra analyses for focus disease: ', trait)
            print('-----------------------------------')
            print('Genes with annotated phenotype terms: ',
                  get_n(pheno_traits, worm_genes))
            print('Genes with annotated tissue terms: ',
                  get_n(tissue_traits, worm_genes))
            print('Genes with annotated go terms: ',
                  get_n(go_traits, worm_genes))
            print('-----------------------------------\n')

        # one of the traits foolishly has a '/'
        if '/' in trait:
            # rename the trait with this character or it breaks your code
            trait = 'post bronchodilator fev1 fevc ratio'

        # enrichment analyses:
        df = tea.enrichment_analysis(worm_genes, phenotype_df, show=False)
        df = df[df.Observed > n_min_obs]
        df.to_csv('../output/phenologues_2/pea_' + trait + '.csv', index=False)

        df = tea.enrichment_analysis(worm_genes, tissue_df, show=False)
        df = df[df.Observed > n_min_obs]
        df.to_csv('../output/disease_tissues_2/tea_' + trait + '.csv',
                  index=False)

        df = tea.enrichment_analysis(worm_genes, go_df, show=False)
        df = df[df.Observed > n_min_obs]
        df.to_csv('../output/disease_go_2/gea_' + trait + '.csv',
                  index=False)


for trait in diseases:
    sel1 = homologs.HID.isin(current)
    worm_genes = homologs[sel1 & sel2].WormBaseID

    if 'Lupus' in trait:
        print(pheno_traits[pheno_traits.term.str.contains(aneuploidy)].wbid)

###############################################################################
###############################################################################
# Begin analysis of ciliary neuron transcriptome
# only analyse stat. sig. genes:

print('Beginning analysis of ciliary neuron transcriptome')
ind = ciliated_df['pval'] < 0.05
sel = (wbids.human_read.isin(ciliated_df[ind].gene_name))
sel2 = (wbids.target_id.isin(ciliated_df[ind].gene_name))
cilia_genes = wbids[sel | sel2].wbid

path = '../output/ciliated_neurons/'
df = tea.enrichment_analysis(cilia_genes, tissue_df, show=False)
df.to_csv(path + 'tea.csv', index=False)

df = tea.enrichment_analysis(cilia_genes, phenotype_df, show=False)
df.to_csv(path + 'pea.csv', index=False)

df = tea.enrichment_analysis(cilia_genes, go_df, show=False)
df.to_csv(path + 'gea.csv', index=False)

###############################################################################
###############################################################################
print('Beginning phenotype-tissue mapping analysis')
egl_genes = pheno_traits[pheno_traits.term.str.contains('egg laying')].wbid

df = tea.enrichment_analysis(egl_genes, tissue_df, show=False)

df



egl_genes = pheno_traits[pheno_traits.term.str.contains('body posture')]
egl_genes.shape
egl_genes.head()

df = tea.enrichment_analysis(egl_genes.wbid, tissue_df, show=False)
df


df = tea.enrichment_analysis(egl_genes.wbid, go_df, show=False)
df


df = tea.enrichment_analysis(egl_genes.wbid, phenotype_df, show=False)
df
