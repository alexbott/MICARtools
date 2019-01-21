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

"""
DESCRIPTION: Plots heatmap for prep_data formatted data

METHODS: We use centroid method as default as it looks at the gene set for each sample holistically to determine a 'centroid' then used for clustering. Other valid options are available (see VARIABLES)

VARIABLES:
data_scaled= Scaled dataframe as created with the MICARtools prep_data function
data_labeled= Unscaled dataframe with sample labels as created with the MICARtools prep_data function
col_colors= Dictionary of labels and colors for plotting, or valid seaborns.clustermap col_colors option
gene_list= List of genes (either as list or as .csv file path and name with list of genes) to plot (IMPORTANT: Gene names are case-sensitive)
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
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
def heatmap(data_scaled, data_labeled, col_colors=None, gene_list=None, save_fig=None, dpi=600, bbox_inches='tight', font_scale=.8,
    cmap=jakes_cmap, center=0, metric='euclidean', method='centroid', xticklabels=True, linewidths=.03, linecolor='#DCDCDC', col_cluster=True,
    row_cluster=False, figsize=(16,6.5)):

    #Set colors for plotting
    if type(col_colors) is dict:
        labels = data_labeled.loc['label']
        color_map = labels.map(col_colors)


    #Custom panel heatmap
    if gene_list != None:

        #Check file formats
        if type(gene_list) is list:
            custom_data = data_scaled.loc[gene_list]
        elif type(gene_list) is str:
            genes = custom_list(str(gene_list))
            custom_data = data_scaled.loc[genes]
        else:
            print('Incorrect gene_list type provided')
            return

        #Prep data for plotting
        plot_data = custom_data.dropna(axis=0)

    else:
        plot_data = data_scaled.dropna(axis=0)


    #Generate clustermap
    sns.set(font_scale=float(font_scale))
    if col_colors is None:
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
        plt.savefig(str(save_fig), dpi=dpi, bbox_inches=str(bbox_inches))

"""
DESCRIPTION: Create violin plots of subset of gene expressions or all gene expression by sample

VARIABLES:
data_labeled= Unscaled dataframe with sample labels as created with the MICARtools prep_data function
gene_list= List of genes (either as list or as .csv file path and name with list of genes) to plot (IMPORTANT: Gene names are case-sensitive)
order= List of samples in order to plot
palette= List of colors for samples in order of plotting
save_fig= If not None, provide full file path, name, and extension to save the file as
dpi= Set dpi of saved figure
bbox_inches= Format saved figure (often useful for making sure no text is cut off)
title= Provide title for figure and saved file if save_fig option used

ASSUMPTIONS:
Data has been scaled and labeled with the MICARtools prep_data function
"""
def sample_overview(data_scaled, info , gene_list=None, order=None, palette=None, save_fig=None, dpi=600, bbox_inches='tight', title=None, grid=False):

    #For some reason, function is returning data_scaled with labels, even though not specified to return anything
    #This solves the issue, but still throws the error
    scaled = data_scaled.copy()

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
            plt.savefig('./MICARtools_violin_plot.pdf', dpi=dpi, bbox_inches=bbox_inches)
        else:
            ax.set_title(str(title))
            plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

"""
DESCRIPTION:

METHODS: Confidence interval number to plot (i.e. 1 == CI1 == 68%, 2 == CI2 == 95%, 3 == CI3 == 99%,)

VARIABLES:
palette= Dictionary of categories and their corresponding plotting colors
n_components= Number of components to plot on the scree plot
order_legend= List of integers in which to reorder samples in legend
plotly_login= ['userid','api key']
USAGE:

ASSUMPTIONS:
Plotly and api key generation, addition

FEATURES TO ADD:
Allow for compatibility with adding labels for gene classes for plotting
Option to order legend
"""
def pca(data_scaled, info, palette, grouping='samples', principle_components=[1,2], _3d_pca=False, gene_labels=False, ci=2,
    gene_list=None, save_fig=None, scree_only=False, save_scree=None, n_components=10, dpi=600, bbox_inches='tight', title=None,
    return_pca_dataframe=False, order_legend=None, grid=False, plotly_login=None, fig_size=(10,10)):

    #Initial variable checks
    if len(principle_components) != 2 and _3d_pca == False:
        print('Incompatible options provided for principle_components and _3d_pca')
        return

    if plotly_login != None and _3d_pca == False:
        print('Only provide plotly login for _3d_pca')
        return

    if _3d_pca == True:
        fig_size=15
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

    if gene_list != None:

        #Check file formats and get dataframe with genes of interest
        if type(gene_list) is list:
            scaled = scaled.reindex[gene_list]

        elif type(gene_list) is str:
            genes = custom_list(str(gene_list))
            scaled = scaled.reindex[genes]

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
        plt.savefig(str(save_fig[:-4]) + '_scree.pdf', dpi=dpi, bbox_inches=bbox_inches)

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
                plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)

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
                    mode='markers',
                    marker=dict(
                        size=fig_size,
                        line=dict(
                            color=palette[unique_labels[0]],
                            width=0.25
                        ),
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
                    mode='markers',
                    marker=dict(
                        size=fig_size,
                        line=dict(
                            color=palette[unique_labels[1]],
                            width=0.25
                        ),
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
                    mode='markers',
                    marker=dict(
                        size=fig_size,
                        line=dict(
                            color=palette[unique_labels[2]],
                            width=0.25
                        ),
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
                ax.scatter(x0, y0, z0, c=palette[unique_labels[0]])

                x1 = pca1.PCa.values
                y1 = pca1.PCb.values
                z1 = pca1.PCc.values
                ax.scatter(x1, y1, z1, c=palette[unique_labels[1]])

                x2 = pca2.PCa.values
                y2 = pca2.PCb.values
                z2 = pca2.PCc.values
                ax.scatter(x2, y2, z2, c=palette[unique_labels[2]])

                ax.set_xlabel(str(pc_list[0]))
                ax.set_ylabel(str(pc_list[1]))
                ax.set_zlabel(str(pc_list[2]))
                plt.show()

                if save_fig != None:
                    plt.savefig(str(save_fig), dpi=dpi, bbox_to_anchor=bbox_to_anchor)

        else:
            return

    else:
        print('This feature has not been implemented yet')

    if return_pca_dataframe is True:
        return df_pca

"""

"""
def gene_overview():

    print('')

"""

"""
def linreg():

    print('')
