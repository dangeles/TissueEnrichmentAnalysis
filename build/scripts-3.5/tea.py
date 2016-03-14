#!/Users/dangeles/anaconda3/envs/testenv/bin/python

import hypergeometricTessts as hgt


if __name__ == '__main__':
    
    import re
    
    path= './'
    os.chdir(path)
    
    import argparse
    import matplotlib
    matplotlib.use('Agg')
#    from matplotlib.pylab import * 
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set_context('paper')
    sns.set_style('whitegrid')
#    sns.despine(left= True)
#    sns.despine(trim= True)

    defQ= 0.1    
    
    parser = argparse.ArgumentParser(description='Run TEA.')
    parser = argparse.ArgumentParser()
    parser.add_argument("tissue_dictionary", help= 'The full path to the tissue dictionary provided by WormBase')
    parser.add_argument("gene_list", help= 'The full path to the gene list (WBIDs) you would like to analyse in .csv format')
    parser.add_argument('-t', "--title", help= 'Title for your analysis (shouldn\'t include file extension', type= str)
    parser.add_argument("-q", help= 'Qvalue threshold for significance. Default is {0} if not provided'.format(defQ), type= float)
    parser.add_argument('-p','--print', help= 'Indicate whether you would like to print results', action= 'store_true')    
    parser.add_argument('-pl', "--plot", help= 'Indicate whether you would like to plot results.', action= 'store_true')
    parser.add_argument('-s', "--save", help= 'Boolean variable indicating whether to save your plot or not.', action= 'store_true')
    args = parser.parse_args()
    
    tdf_name= args.tissue_dictionary
    gl_name= args.gene_list
    title= args.title
    
    if args.q:   
        q= args.q
    else:
        q= defQ
    
    if args.print:
        show= True
    else:
        show= False
        
    if args.plot:
        plot= True
    else:
        plot= False
    
    if args.save:
        save= True
                    
    else:
        save= False
            
    tissue_df= pd.read_csv(tdf_name)
    
    gene_list= []
    with open(gl_name, 'r') as f:
        gene_list= [x.strip() for x in f.readlines()]
    
    df_results, unused= enrichment_analysis(gene_list, tissue_df, alpha= q, show= False)
    
    dfname= title+'.csv'
    df_results.to_csv(dfname, index= False)
    
    if show:        
        with open(dfname, 'r') as f:
            printer= f.readlines()
            
        for value in printer:
            value= value.split(',')
            for val in value:
                if re.findall("\d+\.\d+", val):
                    ind= value.index(val)
                    x= float(val)
                    if x < 10**-6:
                        value[ind]= '<10^-6'
                    else:
                        value[ind]= '{0:.2f}'.format(x)

            value= '\t'.join(value)
            print(value)
    
    if plot:
        plot_enrichment_results(df_results, title= title, save= save)
        plt.show()

    if save:
        plot_enrichment_results(df_results, title= title, save= save)


    sys.exit()
