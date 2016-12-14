import tissue_enrichment_analysis as tea
import pandas as pd
import matplotlib.pyplot as plt
import os

phenotype_df = pd.read_csv('../input/phenotype_ontology.csv')


phenotype_df.columns

for f in os.listdir('../input/genesets_golden'):
    if f != '':
        print(f)
        genes = pd.read_csv('../input/genesets_golden/'+f)
        df, _ = tea.enrichment_analysis(genes.gene, phenotype_df,
                                        show=False)
        df.to_csv('../output/pea_goldenset_results/' + f + '.csv', index=False)
