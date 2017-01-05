# -*- coding: utf-8 -*-
"""Benchmarking PEA."""
import tissue_enrichment_analysis as tea
import pandas as pd
# import matplotlib.pyplot as plt
import os

# import data
phenotype_df = pd.read_csv('../input/phenotype_ontology.csv')
go_df = pd.read_csv('../input/gene_ontology.csv')
tissue_df = tea.fetch_dictionary()

tissue_df.shape
phenotype_df.shape

snp = pd.read_csv('../input/known_snps.csv', sep=',', comment='#')
linker_cell = pd.read_csv('../input/linker_cell.txt', sep='\t')
lc_screen = pd.read_csv('../input/lc_screen.csv')


wbids = pd.read_csv('../input/wormbase_ids.txt', sep='\t', header=None,
                    names=['wbid', 'human_read', 'target_id'])

ciliated_df = pd.read_excel('~/Downloads/mmc3.xlsx')


# define a function to run e.a.s quickly
def run(path, f, output_path, dictionary, column='gene', **kwargs):
    """Run and save e.a.."""
    genes = pd.read_csv(path+f, **kwargs)
    df = tea.enrichment_analysis(genes[column], dictionary,
                                 show=False)
    df = df[df.Observed > 2]
    df.to_csv(output_path + f + '.csv', index=False)

# run the benchmark analyses
for f in os.listdir('../input/genesets_golden'):
    if f != '':
        run('../input/genesets_golden/', f,
            '../output/pea_goldenset_results/', phenotype_df)
        run('../input/genesets_golden/', f,
            '../output/go_goldenset_results/', go_df)
        run('../input/genesets_golden/', f,
            '../output/tea_goldenset_results/', tissue_df)

# run the fog-2 analysis
for f in os.listdir('../input/fog2'):
    if f != '':
        run('../input/fog2/', f, '../output/fog2/pea_', phenotype_df,
            column='ens_gene')
        run('../input/fog2/', f, '../output/fog2/go_', go_df,
            column='ens_gene')
        run('../input/fog2/', f, '../output/fog2/tea_', tissue_df,
            column='ens_gene')

# analyze the hypoxia datasets
for f in os.listdir('../input/hypoxia'):
    if f != '':
        run('../input/hypoxia/', f, '../output/hypoxia/pea_',
            phenotype_df, column='ens_gene')
        run('../input/hypoxia/', f, '../output/hypoxia/tea_',
            tissue_df, column='ens_gene')
        run('../input/hypoxia/', f, '../output/hypoxia/go_',
            go_df, column='ens_gene')


# for disease in snp.Disease_Trait.unique():
#     print(disease, snp[snp.Disease_Trait == disease].shape[0])


# depression = ['Major depressive disorder', 'Major depressive disorder (broad)']
# celiac = ['Celiac disease', 'Celiac disease and Rheumatoid arthritis']
# macular_degeneration = ['Age-related macular degeneration',
#                         'Age-related macular degeneration (CNV vs. GA)',
#                         'Age-related macular degeneration (CNV)',
#                         'Age-related macular degeneration (GA)',
#                         'Age-related macular degeneration (extreme sampling)',
#                         'Age-related macular degeneration (wet)']
# alcoholism = ['Alcohol and nictotine co-dependence', 'Alcohol consumption',
#               'Alcohol dependence', 'Alcoholism' +
#               '(12-month weekly alcohol consumption)',
#               'Alcoholism (alcohol dependence factor score)',
#               'Alcoholism (alcohol use disorder factor score)',
#               'Alcoholism (heaviness of drinking)']
# alzheimers = ['Alzheimer\'s disease', 'Alzheimer\'s disease (age of onset)',
#               'Alzheimer\'s disease (cognitive decline)',
#               'Alzheimer\'s disease (late onset)',
#               'Alzheimer\'s disease (neuritic plaque pathology)',
#               'Alzheimer\'s disease biomarkers']
# adhd = ['Attention deficit hyperactivity disorder',
#         'Attention deficit hyperactivity disorder (combined symptoms)',
#         'Attention deficit hyperactivity disorder'
#         '(hyperactivity-impulsivity symptoms)',
#         'Attention deficit hyperactivity disorder (inattention symptoms)',
#         'Attention deficit hyperactivity disorder (time to onset)',
#         'Attention deficit hyperactivity disorder and conduct disorder',
#         'Attention deficit hyperactivity disorder motor coordination',
#         'Attention deficit hyperactivity disorder symptoms (interaction)']
# bipolar = ['Bipolar disorder']
# heart_disease = ['Coronary heart disease']
# crohns = ['Crohn\'s disease']
# mult_scler = ['Multiple sclerosis',
#               'Multiple sclerosis (OCB status)',
#               'Multiple sclerosis (age of onset)',
#               'Multiple sclerosis (severity)',
#               'Multiple sclerosis--Brain Glutamate Levels']
# hiv = ['AIDS', 'AIDS progression', 'HIV (mother-to-child transmission)',
#        'HIV-1 control', 'HIV-1 progression', 'HIV-1 susceptibility',
#        'HIV-1 viral setpoint', 'HIV-associated dementia']
# obesity_related = ['Obesity', 'Obesity (early onset extreme)',
#                    'Obesity (extreme)', 'Obesity and blood pressure',
#                    'Obesity and osteoporisis', 'Obesity-related traits']
# ovarian_cancer = ['Ovarian cancer',
#                   'Ovariance cancer in BRCA1 mutation carriers']
# diabetes = ['Type 2 diabetes']
# lupus = ['Systemic lupus erythematosus',
#          'Systemic lupus erythematosus and Systemic sclerosis']
# prostate_cancer = ['Prostate cancer',
#                    'Prostate cancer (gene x gene interaction)',
#                    'Prostate cancer mortality']
# hash_of_list = {'prostate_cancer': prostate_cancer,
#                 'lupus': lupus,
#                 'diabetes': diabetes,
#                 'ovarian_cancer': ovarian_cancer,
#                 'obesity_related': obesity_related,
#                 'hiv': hiv, 'mult_scler': mult_scler,
#                 'crohns': crohns, 'heart_disease': heart_disease,
#                 'bipolar': bipolar, 'adhd': adhd,
#                 'alzheimers': alzheimers, 'alcoholism': alcoholism,
#                 'macular_degeneration': macular_degeneration,
#                 'celiac': celiac, 'depression': depression}
#
# path1 = '../input/disease_gwas/'
# path2 = '../input/disease_gwas_links/'
#
# for key, value in hash_of_list.items():
#     ind = (snp.Disease_Trait.isin(value))
#     ind2 = (snp.Mapped_gene != 'Intergenic')
#     genes = snp[ind & ind2].dropna().Mapped_gene.unique()
#     links = snp[ind & ind2][['Link', 'Disease_Trait']].dropna().copy()
#     links.drop_duplicates(inplace=True)
#     if len(genes) < 30:
#         continue
#     if os.path.isfile('../input/disease_gwas/links.csv') is False:
#         links.to_csv('../input/disease_gwas/links.csv', mode='w', index=False)
#     else:
#         links.to_csv('../input/disease_gwas/links.csv', mode='a', index=False,
#                      header=False)
#
#     with open(path1 + key + '.csv', 'w') as f:
#         for gene in genes:
#             f.write(gene + '\n')

# # analyze the SNP data from T-cell paper:
# path = '../input/worm_disease_gwas/'
# output_path = '../output/phenologues/'
# output_path2 = '../output/disease_tissues/'
# for f in os.listdir(path):
#     # run(path, f, output_path, phenotype_df, column='WormBaseID', sep='\t')
#     run(path, f, output_path2, tissue_df, column='WormBaseID', sep='\t')

# half of crohns genes have no orthologs
# half of mult scler have no orthologs

# analyze the linker cell transcriptome
# lc = linker_cell[linker_cell.ratio > 10].Gene.str[:14]
# print(len(lc))
# df, _ = tea.enrichment_analysis(lc, phenotype_df, show=False)
# df.to_csv('../output/linker_cell.csv', index=False)
#
# df, _ = tea.enrichment_analysis(lc, tissue_df, show=False)
# df.to_csv('../output/tea_linker_cell.csv', index=False)
#
# df, _ = tea.enrichment_analysis(lc, go_df, show=False)
# df.to_csv('../output/gea_linker_cell.csv', index=False)
#
# lc_deplete = linker_cell[linker_cell.ratio < 10**-1].Gene.str[:14]
# print(len(lc_deplete))
# df, _ = tea.enrichment_analysis(lc_deplete, phenotype_df, show=False)
# df.to_csv('../output/linker_cell_deplete.csv', index=False)
#
# df, _ = tea.enrichment_analysis(lc_deplete, tissue_df, show=False)
# df.to_csv('../output/tea_linker_cell_deplete.csv', index=False)
#
# df, _ = tea.enrichment_analysis(lc_deplete, go_df, show=False)
# df.to_csv('../output/gea_linker_cell_deplete.csv', index=False)
#
# # analyze the linker cell screen
# df, _ = tea.enrichment_analysis(lc_screen.gene, phenotype_df, show=False)
# df.to_csv('../output/pea_lc_screen.csv', index=False)
#
# sel = lc_screen.phenotype == 'nonwt_RNAi'
# df, _ = tea.enrichment_analysis(lc_screen[sel].gene, phenotype_df, show=False)
# df.to_csv('../output/pea_lc_screen_hits.csv', index=False)


# from the gwas catalog, identify diseases with more than 100 associated genes
f = '~/Downloads/gwas_catalog_v1.0-associations_e87_r2016-12-12.tsv'
all_snps = pd.read_csv(f, sep='\t', encoding='ISO-8859-1')
all_snps.head()
all_snps.columns
array_of_genes = []

# find all the genes in the gwas catalog, and place the set of genes into
# a single list file so you can use DIOPT to find homologs easily
for g in all_snps['REPORTED GENE(S)'].values:
    if type(g) is not str:
        continue
    list_of_genes = [x.strip() for x in g.split(',')]
    for gi in list_of_genes:
        array_of_genes += [gi]

array_of_genes = list(set(array_of_genes))
with open('gene_homologs.csv', 'w+') as F:
    for g in array_of_genes:
        F.write(g+'\n')

# find all the diseases with more than 200 associated genes
candidates = []
for trait in all_snps.TRAIT.unique():
    if all_snps[all_snps.TRAIT == trait].shape[0] > 200:
        print(trait, all_snps[all_snps.TRAIT == trait].shape[0])
        candidates += [trait]

# what are the pubmedids associated with the list of traits?
for trait in candidates:
    print(all_snps[all_snps.TRAIT == trait].PUBMEDID.unique())

# load the homolog data
homologs = pd.read_csv('../input/gwas_homologs_worm.txt', sep='\t')

# for each GWAS disease, find the worm homologs with moderate or high scores
# and use those to run e.a. to identify phenologues of diseases in the worm
sel2 = (homologs.Rank == 'moderate') | (homologs.Rank == 'high')
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
        # one of the traits foolishly has a '/'
        if '/' in trait:
            trait = 'post bronchodilator fev1 fevc ratio'
        df = tea.enrichment_analysis(worm_genes, phenotype_df, show=False)
        # print(_.head(), _.shape)
        # break
        df = df[df.Observed > 10]
        df.to_csv('../output/phenologues_2/pea_' + trait + '.csv', index=False)

        df = tea.enrichment_analysis(worm_genes, tissue_df, show=False)
        df = df[df.Observed > 10]
        df.to_csv('../output/disease_tissues_2/tea_' + trait + '.csv',
                  index=False)

        df = tea.enrichment_analysis(worm_genes, go_df, show=False)
        df = df[df.Observed > 10]
        df.to_csv('../output/disease_go_2/gea_' + trait + '.csv',
                  index=False)


ind = ciliated_df['pval'] < 0.05
sel = (wbids.human_read.isin(ciliated_df[ind].gene_name))
sel2 = (wbids.target_id.isin(ciliated_df[ind].gene_name))

wbids.shape
cilia_genes = wbids[sel | sel2].wbid
df = tea.enrichment_analysis(cilia_genes, tissue_df, show=False)

df = tea.enrichment_analysis(cilia_genes, phenotype_df, show=False)

df = tea.enrichment_analysis(cilia_genes, go_df, show=False)
