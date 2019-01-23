"""
MicArTools
A toolset for navigating and analyzing microarray datasets
alias: micartools

Copyright (C) 2019  Jordan A. Berg
jordan <dot> berg <at> biochem <dot> utah <dot> edu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <https://www.gnu.org/licenses/>.
"""

"""
DEPENDENCIES
"""
from .utils import custom_list
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
sns.set(font='arial')
jakes_cmap = sns.diverging_palette(212, 61, s=99, l=77, sep=1, n=16, center='dark') #Custom aesthetics
from sklearn.decomposition import PCA
import numpy as np
import scipy.stats as stats
from scipy.stats import linregress

"""
DESCRIPTION: Plots heatmap for prep_data formatted data

METHODS: We use centroid method as default as it looks at the gene set for each sample holistically to determine a 'centroid' then used for clustering. Other valid options are available (see VARIABLES)

VARIABLES:
data_scaled= Scaled dataframe as created with the MICARtools prep_data function
info= MICARtools formatted info matrix
palette= Dictionary of labels and colors for plotting, or valid seaborns.clustermap col_colors option
gene_list= List of genes (either as list or as .csv file path and name with list of genes) to plot (IMPORTANT: Gene names are case-sensitive)
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_to_anchor= Format saved figure (often useful for making sure no text is cut off)
font_scale= Scaling factor for font
cmap= A valid seaborns.clustermap cmap option (default: Rutter Lab colorblind-friendly scale)
Please see seaborns.heatmap documentation for descriptions of the following options:
center, metric, method, xticklabels, linewidths, linecolor, col_cluster, row_cluster, figsize

USAGE:
import micartools as mat
col_colors = {'adenocarcinoma': 'grey',
        'adenoma': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        'normal': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}
mat.heatmap(data_scaled, data_labeled, color_dict=color_dict, gene_list='/path/to/gene_list.csv', figsize=(40,8))

ASSUMPTIONS:
Data has been scaled and labeled with the MICARtools prep_data function
"""
def heatmap(data_scaled, info, palette=None, gene_list=None, save_fig=None, dpi=600, bbox_to_anchor='tight', font_scale=.8,
    cmap=jakes_cmap, center=0, metric='euclidean', method='centroid', xticklabels=True, linewidths=.03, linecolor='#DCDCDC', col_cluster=True,
    row_cluster=False, figsize=(16,6.5)):

    scaled = data_scaled.copy()

    #Set colors for plotting
    if palette != None:
        if type(palette) is dict:
            info = info.T
            info.columns = info.iloc[0]
            info = info.reindex(info.index.drop(0))
            info = info.rename({1: 'samples'})
            labels = info.iloc[0]
            color_map = labels.map(palette)

        else:
            print('Error: a dictionary was not provided as palette')
            return

    #Custom panel heatmap
    if gene_list != None:

        #Check file formats
        if type(gene_list) is list:
            custom_data = data_scaled.reindex(labels=gene_list, axis=0)
        elif type(gene_list) is str:
            genes = custom_list(str(gene_list))
            custom_data = data_scaled.reindex(labels=genes, axis=0)
        else:
            print('Incorrect gene_list type provided')
            return

        #Prep data for plotting
        plot_data = custom_data.dropna(axis=0)

    else:
        plot_data = data_scaled.dropna(axis=0)


    #Generate clustermap
    sns.set(font_scale=float(font_scale))
    if palette is None:
        sns.clustermap(plot_data,
                        cmap=cmap,
                        center=float(center),
                        metric=str(metric),
                        method=str(method),
                        xticklabels=xticklabels,
                        linewidths=float(linewidths),
                        linecolor=str(linecolor),
                        col_cluster=col_cluster,
                        row_cluster=row_cluster,
                        figsize=figsize
                       )
    else:
        sns.clustermap(plot_data,
                        cmap=cmap,
                        center=float(center),
                        metric=str(metric),
                        method=str(method),
                        xticklabels=xticklabels,
                        linewidths=float(linewidths),
                        linecolor=str(linecolor),
                        col_cluster=col_cluster,
                        row_cluster=row_cluster,
                        col_colors=color_map,
                        figsize=figsize
                       )

    #Save figure
    if save_fig is not None:
        plt.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=str(bbox_to_anchor))

"""
DESCRIPTION: Create violin plots of subset of gene expressions or all gene expression by sample

VARIABLES:
data_labeled= Unscaled dataframe with sample labels as created with the MICARtools prep_data function
gene_list= List of genes (either as list or as .csv file path and name with list of genes) to plot (IMPORTANT: Gene names are case-sensitive)
order= List of samples in order to plot
palette= Dictionary of matplotlib compatible colors for samples
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_to_anchor= Format saved figure (often useful for making sure no text is cut off)
title= Provide title for figure and saved file if save_fig option used

ASSUMPTIONS:
Data has been scaled and labeled with the MICARtools prep_data function
"""
def sample_overview(data_scaled, info , gene_list=None, order=None, palette=None, save_fig=None, dpi=600, bbox_to_anchor='tight', title=None, grid=False, whitegrid=False):

    if whitegrid == True:
        sns.set_style("whitegrid")

    #For some reason, function is returning data_scaled with labels, even though not specified to return anything
    #This solves the issue, but still throws the error
    scaled = data_scaled.copy()
    scaled = scaled.dropna(axis=0)

    #Prep data_scaled by adding labels from info
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    scaled.loc['label'] = scaled.columns.map(labels.get)

    #Output collapsed dataframe
    newIndex = ['label'] + [ind for ind in scaled.index if ind != 'label']
    scaled = scaled.reindex(index=newIndex)

    scaled = scaled.T
    scaled = scaled.set_index('label', drop=True)
    scaled = scaled.T
    scaled = scaled.dropna(axis=0)

    data_unstacked = scaled.unstack().reset_index()
    data_unstacked = data_unstacked.rename(columns={data_unstacked.columns[0]: 'type', data_unstacked.columns[1]: 'gene', data_unstacked.columns[2]: 'expr'})

    #Custom gene list to grab distributions from
    if gene_list != None:

        #Check file formats
        if type(gene_list) is list:
            genes = gene_list

        elif type(gene_list) is str:
            genes = custom_list(str(gene_list))

        else:
            print('Incorrect gene_list type provided')
            return

        swarm = data_unstacked.loc[(data_unstacked.iloc[:,1].isin(genes))]

    else:
        swarm = data_unstacked

    #Final formatting check
    swarm[["expr"]] = swarm[["expr"]].apply(pd.to_numeric)

    #Plot
    if order is None or palette is None:
        ax = sns.violinplot(x="type", y="expr", data=swarm)
    elif type(order) is list or type(palette) is list:
        ax = sns.violinplot(x="type", y="expr", data=swarm, order=order, palette=palette)
    else:
        return

    ax.set_xlabel('')
    ax.set_ylabel('Expression')

    if grid == False:
        ax.grid(False)

    #Save fig
    if save_fig is not None:
        if title is not None:
            plt.savefig('./MICARtools_violin_plot.pdf', dpi=dpi, bbox_to_anchor=bbox_to_anchor)
        else:
            ax.set_title(str(title))
            plt.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=bbox_to_anchor)

    #Reset aesthetics
    sns.set_style("darkgrid")

"""
DESCRIPTION: Plot a 2-D PCA with confidence intervals or a 3-D PCA with no confidence intervals

METHODS: 3-D PCA -- option to output as interactive plot by providing plotly credentials, or as a static plot without these credentials (default)

VARIABLES:
data_labeled= Unscaled dataframe with sample labels as created with the MICARtools prep_data function -- can be a dataframe prepared using the prep_data() function or log_scale() function
info= MICARtools formatted sample info dataframe
palette= Dictionary of matplotlib compatible colors for samples; for plotly 3-D PCA,
grouping= Perform PCA sample-wise (default) or gene-wise (grouping='genes')
gene_list= List of genes (either as list or as .csv file path and name with list of genes) to plot (IMPORTANT: Gene names are case-sensitive) (only functional when grouping='samples')
gene_labels= Option for grouping='genes', not currently implemented
ci= For 2-D PCA, confidence interval to display for each sample type (i.e. 1 == CI1 == 68%, 2 == CI2 == 95%, 3 == CI3 == 99%)
principle_components= Principle components to plot (default: [1,2] for 2-D PCA, or [1,2,3] for 3-D PCA)
n_components= Number of components to calculate in dataframe (more applicable if you want to perform a deeper survey of the principle components -- values returned if return_pca_dataframe=True)
_3d_pca= Plot 3 principle components (default: False)
plotly_login= ['userid','api key'], usage of this option creates a plotly interactive plot that can be viewed online using your login credentials
title= Provide title for figure and saved file if save_fig option used
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_to_anchor= Format saved figure (often useful for making sure no text is cut off)
order_legend= List of integers to reorder samples in legend (i.e. if samples are displayed 1:Sample_A, 2:Sample_C, 3:Sample_B, provide the list [1,3,2]) (Not currently compatible with 3-D PCA options)
grid= For non-plotly options, remove gridlines from plot
fig_size= Option not used in function currently
size= Marker size

palette examples:
colors = {'adenocarcinoma': (0.5725490196078431, 0.5843137254901961, 0.5686274509803921),
        'adenoma': (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        'normal': (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)}
colors = {'adenocarcinoma': 'gray',
        'adenoma': 'orange',
        'normal': 'green'}
colors = {'adenocarcinoma': '#59656d',
        'adenoma': '#f0833a',
        'normal': '#0a5f38'}

USAGE:

ASSUMPTIONS:
Plotly and api key generation, addition
Data has been scaled and labeled with the MICARtools prep_data function

FEATURES TO ADD:
Allow for compatibility with adding labels for gene classes for plotting
Option to order legend
Move legend outside plot
Add options to vary marker size and opacity
"""
def pca(data_scaled, info, palette, grouping='samples', gene_list=None, gene_labels=False,
    ci=2, principle_components=[1,2], n_components=10, _3d_pca=False, plotly_login=None,
    scree_only=False, save_scree=None, return_pca_dataframe=False, size=10, whitegrid=False,
    title=None, save_fig=None, dpi=600, bbox_to_anchor='tight', order_legend=None, grid=False, fig_size=(10,10)):

    if whitegrid == True:
        sns.set_style("whitegrid")

    #Initial variable checks
    if len(principle_components) != 2 and _3d_pca == False:
        print('Incompatible options provided for principle_components and _3d_pca')
        return

    if plotly_login != None and _3d_pca == False:
        print('Only provide plotly login for _3d_pca')
        return

    if _3d_pca == True:
        if principle_components == [1,2]:
            principle_components = [1,2,3]
        elif len(principle_components) != 3:
            print('Incompatible options provided for principle_components and _3d_pca')
            return
        else:
            print('Not sure if this will ever be triggered, but if it is need to check case')
            return

        if plotly_login != None and type(plotly_login) is not list:
            print('Provide proper plotly login information in form of list')
            print("['userid','api key']")
            return

    scaled = data_scaled.copy()
    scaled = scaled.dropna(axis=0)

    if gene_list != None:

        #Check file formats and get dataframe with genes of interest
        if type(gene_list) is list:
            scaled = scaled.reindex(labels=gene_list, axis=0)


        elif type(gene_list) is str:
            genes = custom_list(str(gene_list))
            scaled = scaled.reindex(labels=genes, axis=0)

        else:
            print('Incorrect gene_list type provided')
            return

    #Clean data
    scaled = scaled.dropna(axis=0)

    #Format data for sample-wise or gene-wise PCA analysis
    if str(grouping) == 'samples':
        scaled = scaled.T

    elif str(grouping) == 'genes':
        scaled = scaled

    else:
        print('Incorrect grouping variable provided')
        return

    #Prep PCA
    pca = PCA(n_components=n_components)
    pca_result = pca.fit_transform(scaled.values)

    #Record PCs
    x = 1
    while x <= n_components:
        pc = 'PC' + str(x)
        component = x - 1
        scaled[pc] = pca_result[:,component]
        x += 1

    #Scree
    vari = 'Explained variation per principal component: {}'.format(np.round(pca.explained_variance_ratio_, decimals=4)*100)
    scree = np.round(pca.explained_variance_ratio_, decimals=4)*100

    if save_scree != None:

        #Plot scree
        sing_vals = np.arange(n_components) + 1

        ax = sns.lineplot(x=sing_vals, y=scree, color="red")
        ax.set(xlabel='Principal Component', ylabel='Proportion of Variance Explained', title='Scree Plot')
        plt.savefig(str(save_fig[:-4]) + '_scree.pdf', dpi=dpi, bbox_to_anchor=bbox_to_anchor)

        if grid == False:
            ax.grid(False)

        #Remove scree from memory to prevent plot bleeding
        del ax
        plt.close()
        plt.clf()

        if scree_only == True:
            return

    if gene_labels == False:

        #Prep data_scaled by adding labels from info to column (samples are rows)
        labels = pd.Series(info[1].values,index=info[0]).to_dict()

        scaled['label'] = scaled.index.to_series().map(labels)

        #Plot PCa & PCb
        pc_list = []
        for p in principle_components:
            pc_list.append('PC' + str(p))

        pc_list.append('label')

        df_pca = scaled[pc_list] #Prepare pca data

        if _3d_pca == False:
            df_pca.columns = ['PCa', 'PCb', 'label']
            unique_labels = df_pca['label'].unique() #Gather unique labels
            pca_plot = sns.scatterplot(df_pca.PCa, df_pca.PCb, hue=df_pca['label'], palette=palette)

            #Add confidence intervals
            for x in unique_labels:

                #slice df into label specific datasets
                df_slice = df_pca[df_pca['label'] == x]

                #make numpy array from label-specific dfs
                x_slice = df_slice.PCa.values
                y_slice = df_slice.PCb.values

                #pca maths
                cov = np.cov(x_slice, y_slice)
                lambda_, v = np.linalg.eig(cov)
                theta = np.degrees(np.arctan2(*v[:,0][::-1]))
                lambda_ = np.sqrt(lambda_)

                #plot
                pca_plot.add_patch(patches.Ellipse(xy=(np.mean(x_slice), np.mean(y_slice)),
                                  width=lambda_[0]*ci*2, height=lambda_[1]*ci*2,
                                  angle=theta,
                                  alpha=0.3, facecolor=palette[x], edgecolor='black', linewidth=1, linestyle='solid')
                                  )

            # Put the legend out of the figure
            handles,labels = pca_plot.get_legend_handles_labels()

            if order_legend != None:
                if type(order_legend) is list:
                    plt.legend([handles[idx] for idx in order_legend],[labels[idx] for idx in order_legend], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                else:
                    plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
                    print('order_legend datatype is invalid -- plotting samples in default order...')
            else:
                plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

            plt.xlabel('PC' + str(principle_components[0]) + ' (' + str(round(scree[(principle_components[0] - 1)],2)) + '%)')
            plt.ylabel('PC' + str(principle_components[1]) + ' (' + str(round(scree[(principle_components[1] - 1)],2)) + '%)')

            if grid == False:
                plt.grid(False)

            if save_fig != None:

                #Save plot
                plt.title(str(title))
                plt.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=bbox_to_anchor)

        elif _3d_pca == True:

            if plotly_login != None:
                import plotly
                import plotly.plotly as py
                import plotly.graph_objs as go
                plotly.tools.set_credentials_file(username=str(plotly_login[0]), api_key=str(plotly_login[1]))

                df_pca.columns = ['PCa', 'PCb', 'PCc', 'label']
                unique_labels = df_pca['label'].unique() #Gather unique labels

                pca0 = df_pca.loc[df_pca['label'] == unique_labels[0]]
                pca1 = df_pca.loc[df_pca['label'] == unique_labels[1]]
                pca2 = df_pca.loc[df_pca['label'] == unique_labels[2]]

                x0 = pca0.PCa.values
                y0 = pca0.PCb.values
                z0 = pca0.PCc.values
                trace0 = go.Scatter3d(
                    x=x0,
                    y=y0,
                    z=z0,
                    name=unique_labels[0],
                    mode='markers',
                    marker=dict(
                        size=size,
                        color=palette[unique_labels[0]],
                        opacity=0.8
                    )
                )

                x1 = pca1.PCa.values
                y1 = pca1.PCb.values
                z1 = pca1.PCc.values
                trace1 = go.Scatter3d(
                    x=x1,
                    y=y1,
                    z=z1,
                    name=unique_labels[1],
                    mode='markers',
                    marker=dict(
                        size=size,
                        color=palette[unique_labels[1]],
                        opacity=0.8
                    )
                )

                x2 = pca2.PCa.values
                y2 = pca2.PCb.values
                z2 = pca2.PCc.values
                trace2 = go.Scatter3d(
                    x=x2,
                    y=y2,
                    z=z2,
                    name=unique_labels[2],
                    mode='markers',
                    marker=dict(
                        size=size,
                        color=palette[unique_labels[2]],
                        opacity=0.8
                    )
                )

                data = [trace0, trace1, trace2]
                layout = go.Layout(
                    margin=dict(
                        l=0,
                        r=0,
                        b=0,
                        t=0
                    )
                )
                fig = go.Figure(data=data, layout=layout)

                if title != None:
                    py.iplot(fig, filename=title)
                else:
                    py.iplot(fig, filename='3d_pca')

            else:

                from mpl_toolkits.mplot3d import Axes3D

                fig = plt.figure()
                ax = Axes3D(fig)

                df_pca.columns = ['PCa', 'PCb', 'PCc', 'label']
                unique_labels = df_pca['label'].unique() #Gather unique labels

                pca0 = df_pca.loc[df_pca['label'] == unique_labels[0]]
                pca1 = df_pca.loc[df_pca['label'] == unique_labels[1]]
                pca2 = df_pca.loc[df_pca['label'] == unique_labels[2]]

                x0 = pca0.PCa.values
                y0 = pca0.PCb.values
                z0 = pca0.PCc.values
                ax.scatter(x0, y0, z0, c=palette[unique_labels[0]], label=str(unique_labels[0]))

                x1 = pca1.PCa.values
                y1 = pca1.PCb.values
                z1 = pca1.PCc.values
                ax.scatter(x1, y1, z1, c=palette[unique_labels[1]], label=str(unique_labels[1]))

                x2 = pca2.PCa.values
                y2 = pca2.PCb.values
                z2 = pca2.PCc.values
                ax.scatter(x2, y2, z2, c=palette[unique_labels[2]], label=str(unique_labels[2]))

                ax.set_xlabel(str(pc_list[0]))
                ax.set_ylabel(str(pc_list[1]))
                ax.set_zlabel(str(pc_list[2]))
                ax.legend()

                plt.show()

                if save_fig != None:
                    plt.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=bbox_to_anchor)

        else:
            return

    else:
        print('This feature has not been implemented yet')

    if return_pca_dataframe is True:
        return df_pca

    #Reset aesthetics
    sns.set_style("darkgrid")

"""
DESCRIPTION: Plot boxplot overlaid with swarmplot of each sample type's gene expression for the given gene

VARIABLES:
data= Dataframe (can be sample-normalized, prep_data() scaled, or log_scale() scaled)
info= MICARtools formatted sample info dataframe
gene_name= Name of gene to plot
palette= Dictionary of matplotlib compatible colors for samples
order= List of samples in order to plot
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_to_anchor= Format saved figure (often useful for making sure no text is cut off)
grid= Control plot gridlines (default: False)

ASSUMPTIONS:
data and info dataframes are properly formatted for MICARtools and any appropriate sample/gene normalizations have been performed
"""
def gene_overview(data, info, gene_name, palette, order=None, save_fig=None, dpi=600, bbox_to_anchor='tight', grid=False, whitegrid= False):

    if whitegrid == True:
        sns.set_style("whitegrid")

    data_copy = data.copy()
    data_copy = data_copy.dropna(axis=0)

    #Prep data_scaled by adding labels from info to column (samples are rows)
    labels = pd.Series(info[1].values,index=info[0]).to_dict()

    data_copy = data_copy.T
    data_copy['label'] = data_copy.index.to_series().map(labels)

    gene_df = data_copy[[str(gene_name), 'label']]
    gene_df[[str(gene_name)]] = gene_df[[str(gene_name)]].apply(pd.to_numeric, errors='coerce')

    ax = sns.catplot(x='label', y=str(gene_name), data=gene_df, color='black', order=order, kind='swarm') #Swarm plot
    ax = sns.boxplot(x='label', y=str(gene_name), data=gene_df, width =0.3, fliersize=0, order=order, palette=palette) #Boxplot, fliersize=0 removes outlier diamonds from sns
    plt.setp(ax.collections, sizes=[12]) #Resize markers for catplot

    if grid == False:
        ax.grid(False)

    fig = ax.get_figure()

    if save_fig != None:
        fig.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=bbox_to_anchor)

    #Reset plot aesthetics
    sns.set_style("darkgrid")

"""
DESCRIPTION: Calculates r, r^2 values, and p-values for every gene against target gene for given dataset

VARIABLES:
data= Dataframe of expression data (samples are columns, genes are rows)
gene_name= Gene to run correlations against
save_file= File path and name to save correlations matrix with
delimiter= Delimiter to use for save_file

ASSUMPTIONS:
data dataframes is properly formatted for MICARtools, sample scaling required but gene scaling is not necessary
data should not contain labels
"""
def linreg(data, gene_name, save_file, delimiter=','):

    df_corr = data.copy()
    df_corr = df_corr.dropna(axis=0)

    #Gene expression array for gene of interest
    interest = df_corr.loc[str(gene_name)].values.tolist()
    interest = np.array(interest).astype(np.float)
    interest = np.ndarray.tolist(interest)

    #Get expression arrays for all genes and run linear regression, save values to matrix
    lm_interest=[]

    for row in df_corr.iterrows():
        index, data = row
        gene = df_corr.loc[index].values.tolist()
        gene = np.array(gene).astype(np.float)
        gene = np.ndarray.tolist(gene)

        if len(gene) is not len(interest):
            continue
        else:
            slope, intercept, r_value, p_value, std_err = linregress(gene, interest)
            lm_interest.append([index, slope, intercept, r_value, (r_value**2), p_value, std_err])

    df_lm_interest = pd.DataFrame(lm_interest, columns=['gene', 'slope', 'intercept', 'r_value', 'r_squared', 'p_value', 'std_err'])

    # Save linear modeling metrics to .csv
    df_lm_interest.to_csv(str(save_file), sep=delimiter)

"""
If add_linreg used, title variable is void
"""
def scatter(data, info, gene1, gene2, palette, add_linreg=False, order_legend=None, title=None, save_fig=None, dpi=600, bbox_to_anchor='tight', grid=False, whitegrid=False, alpha=1):

    if whitegrid == True:
        sns.set_style("whitegrid")

    data_c = data.copy()
    data_c = data_c.dropna(axis=0)

    #Prep data_scaled by adding labels from info
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_c.loc['label'] = data_c.columns.map(labels.get)
    data_c = data_c.T

    ax = sns.scatterplot(data_c[str(gene1)], data_c[str(gene2)], hue=data_c['label'], palette=palette, alpha=alpha)

    if add_linreg == True:

        gene_a = data_c[str(gene1)].values.tolist()
        gene_a = np.array(gene_a).astype(np.float)
        gene_a = np.ndarray.tolist(gene_a)

        gene_b = data_c[str(gene2)].values.tolist()
        gene_b = np.array(gene_b).astype(np.float)
        gene_b = np.ndarray.tolist(gene_b)

        slope, intercept, r_value, p_value, std_err = linregress(gene_a, gene_b)

        title = 'r = ' + "%.2f" % round(r_value,4)

        x = np.linspace(data_c[str(gene1)].min(), data_c[str(gene1)].max(), 100)
        y = (slope * x) + intercept
        ax.plot(x, y, '-k')



    # Put the legend out of the figure
    handles,labels = ax.get_legend_handles_labels()

    if order_legend != None:
        if type(order_legend) is list:
            plt.legend([handles[idx] for idx in order_legend],[labels[idx] for idx in order_legend], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        else:
            plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            print('order_legend datatype is invalid -- plotting samples in default order...')
    else:
        plt.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    plt.xlabel(str(gene1))
    plt.ylabel(str(gene2))

    if title != None:
        plt.title(str(title))

    if grid == False:
        plt.grid(False)

    if save_fig != None:
        plt.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=bbox_to_anchor)

    #Revert to default styles
    sns.set_style('darkgrid')

"""
DESCRIPTION:

VARIABLES:
data= Sample-normalized, MICARtools-formatted data -- should NOT be gene-normalized
info= MICARtools formatted sample info dataframe
label_comp= Sample label to compare against base
label_base= Sample label to use as base for comparison
highlight_genes= If provided with a list, or a file path and name to a .csv-type series of gene names, will highlight genes of interest in a different color (Gene names are case-sensitive)
highlight_color= Color to use for highlighted genes
y_threshold= -log10(P-Value) threshold to use to identify significant hits. Will create a dotted line on the plot and use to pull significant hits if return_threshold_hits is not None
x_threshold= log2(Fold Change) threshold to use to identify significant hits. Will take the positive and negative of the number. Will create a dotted line on the plot and use to pull significant hits if return_threshold_hits is not None
return_threshold_hits= Will return a .csv-type matrix of significant genes/hits to the file path and name specified
return_threshold_hits_delimiter= Delimiter to use in exporting return_threshold_hits
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_to_anchor= Format saved figure (often useful for making sure no text is cut off)

LIMITATIONS:
Will only perform comparison between 2 sample types

ASSUMPTIONS:
data should ONLY be sample normalized. If using a previous function that returned a modified original
y_threshold must be a postive integer or float
"""
def volcano(data, info, label_comp, label_base, highlight_genes=None, highlight_color='DarkRed', alpha=1, alpha_highlights=1,
            y_threshold=10, x_threshold=1, return_threshold_hits=False, export_threshold_hits=None, export_threshold_hits_delimiter=',',
            save_fig=None, dpi=600, bbox_to_anchor='tight', whitegrid=False):

    if whitegrid == True:
        sns.set_style("whitegrid")

    #Prep dataframes
    data_c = data.copy()
    data_c = data_c.dropna(axis=0)
    info_c = info.copy()

    #Add labels to column name
    if 'label' in data_c.index:
        data_c = data_c.drop(labels='label',axis=0)

    info_c["id"] = info_c[0] + '_' + info_c[1]
    labels = pd.Series(info_c['id'].values,index=info_c[0]).to_dict()
    data_c = data_c.rename(labels, axis='columns')

    # Average every by cell line
    data_c['log2 Fold Change'] = np.log2((data_c.filter(regex=str(label_comp)).mean(axis=1)) / \
                                      (data_c.filter(regex=str(label_base)).mean(axis=1)))
    data_c['-log10 P-Value'] = ''

    # Calculate p-value using 1-way ANOVA with replicates and append to df_oxsm_volc
    for row in data_c.iterrows():
        index, data = row
        comp_row = data_c.loc[index].filter(regex=str(label_comp)).values.tolist()
        ground_row = data_c.loc[index].filter(regex=str(label_base)).values.tolist()

        # Append p_value to df_oxsm_volc
        statistic, p_value = stats.f_oneway(comp_row, ground_row)
        data_c.loc[index,'-log10 P-Value'] = float(-1 * (np.log10(p_value)))

    #Plot all genes
    ax = sns.scatterplot(x='log2 Fold Change', y='-log10 P-Value', data=data_c, color='Black', alpha=alpha)

    #Plot selected genes if user-specified
    if highlight_genes != None:
        if type(highlight_genes) is list:
            df_genes = data_c.reindex(labels=highlight_genes, axis=0)
        elif type(highlight_genes) is str:
            genes = custom_list(str(highlight_genes))
            df_genes = data_c.reindex(labels=genes, axis=0)
        else:
            print('Invalid highlight_genes option provided.')
            return

        ax = sns.scatterplot(x='log2 Fold Change', y='-log10 P-Value', data=df_genes, color=str(highlight_color), alpha=alpha_highlights)

    #Plot thresholds
    if y_threshold != None and type(y_threshold) is int or type(y_threshold) is float:
        if y_threshold > 0:
            ax.axhline(y_threshold, ls='--', color='b')
        else:
            print('Invalid y_threshold provided, must be a positive integer or float (Y-axis threshold will not be plotted)')
    else:
        print('Invalid y_threshold provided, must be a positive integer or float (Y-axis threshold will not be plotted)')

    if x_threshold != None and type(x_threshold) is int or type(x_threshold) is float:
        ax.axvline(-x_threshold, ls='--', color='b')
        ax.axvline(x_threshold, ls='--', color='b')
    else:
        print('Invalid x_threshold provided, must be an integer or float (X-axis threshold will not be plotted)')

    #Set labels and other plotting aesthetics
    ax.set_ylabel('-log$_1$$_0$(P-Value)')
    ax.set_xlabel('-log$_2$(Fold Change)')
    plt.grid(False)

    #Save plot if user-specified
    if save_fig != None:
        plt.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=bbox_to_anchor)

    #Save hits if user-specified
    if export_threshold_hits != None or return_threshold_hits == True:
        df_c = data_c[['log2 Fold Change', '-log10 P-Value']].copy()
        df_up = df_c.loc[(df_c['log2 Fold Change'] > x_threshold) & (df_c['-log10 P-Value'] > y_threshold)] #get upregulated hits
        df_down = df_c.loc[(df_c['log2 Fold Change'] < -x_threshold) & (df_c['-log10 P-Value'] > y_threshold)] #get downregulated hits
        thresh_hits = df_up.append(df_down) #append hits tables
        if export_threshold_hits != None:
            thresh_hits.to_csv(str(return_threshold_hits), sep=export_threshold_hits_delimiter) #export table for user
        else:
            return thresh_hits

    #Revert to default styles
    sns.set_style('darkgrid')

"""

"""
def jointplot(data, info, gene1, gene2, kind='reg', palette=None, order=None, save_fig=None, dpi=600, bbox_to_anchor='tight', whitegrid=False, grid=False):

    if whitegrid == True:
        sns.set_style("whitegrid")

    data_c = data.copy()
    data_c = data_c.dropna(axis=0)

    #Prep data_scaled by adding labels from info
    if 'label' not in data_c.index:
        labels = pd.Series(info[1].values,index=info[0]).to_dict()
        data_c.loc['label'] = data_c.columns.map(labels.get)

    #Get r
    data_c = data_c.T
    gene_a = data_c[str(gene1)].values.tolist()
    gene_a = np.array(gene_a).astype(np.float)
    gene_a = np.ndarray.tolist(gene_a)

    gene_b = data_c[str(gene2)].values.tolist()
    gene_b = np.array(gene_b).astype(np.float)
    gene_b = np.ndarray.tolist(gene_b)

    r_value = stats.pearsonr(gene_a, gene_b)[0]

    #Plot
    ax = sns.jointplot(x=gene_a, y=gene_b, kind='reg')
    ax.ax_joint.collections[0].set_visible(False)
    ax = sns.scatterplot(x=str(gene1), y=str(gene2), data=data_c, hue='label', palette=palette, hue_order=order)
    ax.set_title('r: ' + str(round(r_value,2)), y=0.92, x=0.09)

    if r_value > 0:
        plt.legend(loc='lower right')

    if grid == False:
        plt.grid(False)

    fig = ax.get_figure()
    plt.show()

    if save_fig != None:
        fig.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=bbox_to_anchor)

    #Revert to default styles
    sns.set_style('darkgrid')

"""

"""
def violin():

    print('')
