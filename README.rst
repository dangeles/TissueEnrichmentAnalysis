Tissue Enrichment Analysis (TEA)
================================

This repository holds the scripts to calculate enrichment of a set of
labels using hypergeometric statistics. Results are output to a
dataframe. A standard plotting function is provided, but it is really
just a thin wrapper around a seaborn plot.

Requirements
================================

This library has been developed in Python >= 3.5, using the Anaconda
distribution. Requirements include ``pandas``, ``matplotlib``, ``numpy``,
``scipy`` and ``seaborn``.

Installation
================================
Use ``pip install tissue_enrichment_analysis``

Basic Usage
================================

Web usage
----------------------

Go to `www.wormbase.org/tea <http://www.wormbase.org/tea>`_, input your gene list
and enjoy the results!


Within a Python Script
----------------------

There are really just two main functions that are provided in TEA:
``enrichment_analysis`` and ``plot_enrichment_results``.

A standard call to this library would be as follows:

``import tissue_enrichment_analysis as tea``

``gene_list= some_gene_list``

``tissue_df= tea.fetch_dictionary()``

``df_results= tea.enrichment_analysis(tissue_df, gene_list, aname= 'FileName')``

``tea.plot_enrichment_results(df_results, title= 'FileName')``



Calling from Terminal
---------------------

Gene enrichment analysis can be generated easily by calling the program via terminal using:
``tea tissue_dictionary your_gene_list -[OPTIONS]``

Type
``tea -h`` for help and full documentation.



Future Work
================================
We may try to add support for other model organisms!



Contact
================================

If you find any bugs, have suggestions or just want to say hi, feel free
to contact me at dangeles@caltech.edu

Good luck!

David Angeles-Albores

Author:
=======
David Angeles-Albores

Contributors:
================================

Raymond Y. Lee, Juancarlos Chan, Paul W. Sternberg

Acknowledgements
================

With special thanks to the entire worm community!
