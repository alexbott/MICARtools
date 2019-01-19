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
        col_colors = labels.map(col_colors)

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
        plot_data = custom_data

    else:
        plot_data = data_scaled

    #Generate clustermap
    sns.set(font_scale=float(font_scale))
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
                    col_colors=col_colors,
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

ASSUMPTIONS:
Data has been scaled and labeled with the MICARtools prep_data function
"""
def sample_overview(data_scaled, info , gene_list=None, order=None, palette=None, save_fig=None, dpi=600, bbox_inches='tight', title=None):

    #Prep data_scaled by adding labels from info
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data_scaled.loc['label'] = data_scaled.columns.map(labels.get)

    #Output collapsed dataframe
    newIndex = ['label'] + [ind for ind in data_scaled.index if ind != 'label']
    data_scaled = data_scaled.reindex(index=newIndex)

    data_scaled = data_scaled.T
    data_scaled = data_scaled.set_index('label', drop=True)
    data_scaled = data_scaled.T

    data_unstacked = data_scaled.unstack().reset_index()
    data_unstacked = data_unstacked.rename(columns={data_unstacked.columns[0]: 'type', data_unstacked.columns[1]: 'gene', data_unstacked.columns[2]: 'expr'})

    #Custom panel heatmap
    if gene_list != None:

        #Check file formats
        if type(gene_list) is list:
            genes = gene_list

        elif type(gene_list) is str:
            genes = custom_list(str(gene_list))

        else:
            print('Incorrect gene_list type provided')
            return

    swarm = data_unstacked.loc[(data_unstacked.ix[:,1].isin(genes))]

    if order is None or palette is None:
        ax = sns.violinplot(x="type", y="expr", data=swarm)
    elif type(order) is list or type(palette) is list:
        ax = sns.violinplot(x="type", y="expr",
                        data=swarm, order=order,
                       palette=palette)
    else:
        return

    if save_fig is not None:
        if title is not None:
            ax.set_xlabel('')
            ax.set_ylabel('Expression')
            plt.savefig('./MICARtools_violin_plot.pdf', dpi=dpi, bbox_inches=bbox_inches)
        else:
            ax.set_xlabel('')
            ax.set_ylabel('Expression')
            ax.set_title(str(title))
            plt.savefig(str(save_fig), dpi=dpi, bbox_inches=bbox_inches)
