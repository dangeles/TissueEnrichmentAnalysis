import tissue_enrichment_analysis as tea
import pandas as pd


phenotype_df = pd.read_csv('../input/phenotype_ontology.csv')
go_df = pd.read_csv('../input/gene_ontology.csv')
tissue_df = tea.fetch_dictionary()

sfari = pd.read_excel('../input/sfari.xlsx')
name_df = pd.read_excel('../input/sfari_name_converter.xlsx')

sfari.head()
df, _ = tea.enrichment_analysis(sfari.Gene, tissue_df, show=False)

df.to_csv('../output/tea_sfari.csv', index=False)


df, _ = tea.enrichment_analysis(sfari.Gene, phenotype_df, show=False)

df.to_csv('../output/pea_sfari.csv', index=False)

df, _ = tea.enrichment_analysis(sfari.Gene, go_df, show=False)

df.to_csv('../output/goa_sfari.csv', index=False)


melt_pheno = pd.melt(phenotype_df, id_vars='wbid', var_name='phenotype')
melt_pheno = melt_pheno[melt_pheno.value == 1]

def convert(x):
    return name_df[name_df.wbid == x].gene_name.values[0]


melt_pheno = melt_pheno[melt_pheno.wbid.isin(sfari.Gene.unique())]
melt_pheno['gene_name'] = melt_pheno.wbid.apply(convert)

melt_pheno.wbid.unique().shape
fname = '../output/annotated_genes_with_phenotype.csv'
melt_pheno[['gene_name', 'wbid', 'phenotype']].to_csv(fname, index=False)
