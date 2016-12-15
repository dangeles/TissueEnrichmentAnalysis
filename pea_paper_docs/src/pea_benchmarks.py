import tissue_enrichment_analysis as tea
import pandas as pd
import matplotlib.pyplot as plt
import os

phenotype_df = pd.read_csv('../input/phenotype_ontology.csv')
go_df = pd.read_csv('../input/gene_ontology.csv')

snp = pd.read_csv('../input/known_snps.csv', sep=',', comment='#')


def run(path, f, output_path, dictionary, **kwargs):
    """Run and save e.a.."""
    genes = pd.read_csv(path+f, **kwargs)
    df, _ = tea.enrichment_analysis(genes.gene, dictionary,
                                    show=False)
    df.to_csv(output_path + f + '.csv', index=False)

for f in os.listdir('../input/genesets_golden'):
    if f != '':
        run('../input/genesets_golden/', f,
            '../output/pea_goldenset_results/', phenotype_df)
        run('../input/genesets_golden/', f,
            '../output/go_goldenset_results/', go_df)

for f in os.listdir('../input/fog2'):
    if f != '':
        run('../input/fog2/', f, '../output/fog2/pea_', phenotype_df)
        run('../input/fog2/', f, '../output/fog2/go_', go_df)

for f in os.listdir('../input/hypoxia'):
    if f != '':
        run('../input/hypoxia/', f, '../output/hypoxia/pea_', phenotype_df)
        run('../input/hypoxia/', f, '../output/hypoxia/go_', go_df)

for disease in snp.Disease_Trait.unique():
    print(disease, snp[snp.Disease_Trait == disease].shape[0])

depression = ['Major depressive disorder', 'Major depressive disorder (broad)']
celiac = ['Celiac disease', 'Celiac disease and Rheumatoid arthritis']
macular_degeneration = ['Age-related macular degeneration',
                        'Age-related macular degeneration (CNV vs. GA)',
                        'Age-related macular degeneration (CNV)',
                        'Age-related macular degeneration (GA)',
                        'Age-related macular degeneration (extreme sampling)',
                        'Age-related macular degeneration (wet)']
alcoholism = ['Alcohol and nictotine co-dependence', 'Alcohol consumption',
              'Alcohol dependence', 'Alcoholism' +
              '(12-month weekly alcohol consumption)',
              'Alcoholism (alcohol dependence factor score)',
              'Alcoholism (alcohol use disorder factor score)',
              'Alcoholism (heaviness of drinking)']
alzheimers = ['Alzheimer\'s disease', 'Alzheimer\'s disease (age of onset)',
              'Alzheimer\'s disease (cognitive decline)',
              'Alzheimer\'s disease (late onset)',
              'Alzheimer\'s disease (neuritic plaque pathology)',
              'Alzheimer\'s disease biomarkers']
adhd = ['Attention deficit hyperactivity disorder',
        'Attention deficit hyperactivity disorder (combined symptoms)',
        'Attention deficit hyperactivity disorder'
        '(hyperactivity-impulsivity symptoms)',
        'Attention deficit hyperactivity disorder (inattention symptoms)',
        'Attention deficit hyperactivity disorder (time to onset)',
        'Attention deficit hyperactivity disorder and conduct disorder',
        'Attention deficit hyperactivity disorder motor coordination',
        'Attention deficit hyperactivity disorder symptoms (interaction)']
bipolar = ['Bipolar disorder']
heart_disease = ['Coronary heart disease']
crohns = ['Crohn\'s disease']
mult_scler = ['Multiple sclerosis',
              'Multiple sclerosis (OCB status)',
              'Multiple sclerosis (age of onset)',
              'Multiple sclerosis (severity)',
              'Multiple sclerosis--Brain Glutamate Levels']
hiv = ['AIDS', 'AIDS progression', 'HIV (mother-to-child transmission)',
       'HIV-1 control', 'HIV-1 progression', 'HIV-1 susceptibility',
       'HIV-1 viral setpoint', 'HIV-associated dementia']
obesity_related = ['Obesity', 'Obesity (early onset extreme)',
                   'Obesity (extreme)', 'Obesity and blood pressure',
                   'Obesity and osteoporisis', 'Obesity-related traits']
ovarian_cancer = ['Ovarian cancer',
                  'Ovariance cancer in BRCA1 mutation carriers']
diabetes = ['Type 2 diabetes']
lupus = ['Systemic lupus erythematosus',
         'Systemic lupus erythematosus and Systemic sclerosis']
prostate_cancer = ['Prostate cancer',
                   'Prostate cancer (gene x gene interaction)',
                   'Prostate cancer mortality']
hash_of_list = {'prostate_cancer': prostate_cancer,
                'lupus': lupus,
                'diabetes': diabetes,
                'ovarian_cancer': ovarian_cancer,
                'obesity_related': obesity_related,
                'hiv': hiv, 'mult_scler': mult_scler,
                'crohns': crohns, 'heart_disease': heart_disease,
                'bipolar': bipolar, 'adhd': adhd,
                'alzheimers': alzheimers, 'alcoholism': alcoholism,
                'macular_degeneration': macular_degeneration,
                'celiac': celiac, 'depression': depression}

path1 = '../input/disease_gwas/'
path2 = '../input/disease_gwas_links/'
for key, value in hash_of_list.items():
    ind = (snp.Disease_Trait.isin(value))
    ind2 = (snp.Mapped_gene != 'Intergenic')
    genes = snp[ind & ind2].dropna().Mapped_gene.unique()
    links = snp[ind & ind2][['Link', 'Disease_Trait']].dropna().copy()
    links.drop_duplicates(inplace=True)
    if len(genes) < 30:
        continue
    if os.path.isfile('../input/disease_gwas/links.csv') is False:
        links.to_csv('../input/disease_gwas/links.csv', mode='w', index=False)
    else:
        links.to_csv('../input/disease_gwas/links.csv', mode='a', index=False,
                     header=False)

    with open(path1 + key + '.csv', 'w') as f:
        for gene in genes:
            f.write(gene + '\n')

path = '../input/worm_disease_gwas/'
output_path = '../output/phenologues/'
for f in os.listdir(path):
    run(path, f, output_path, phenotype_df, sep='\t')

# half of crohns genes have no orthologs
# half of mult scler have no orthologs
