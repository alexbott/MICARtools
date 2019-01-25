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
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from multiprocessing import cpu_count, Pool
import scipy.stats as stats
import numpy as np
import pandas as pd
from functools import partial

"""
DESCRIPTION: Axis-agnostic list-reader from dataframe
VARIABLES:
USAGE:
ASSUMPTIONS:
"""
def custom_list(file, delimiter=','):

    gene_df = pd.read_csv(str(file), sep=delimiter, index_col=False, header=None)
    genes = gene_df.values.tolist()
    gene_list = [val for sublist in genes for val in sublist]

    return gene_list

"""
DESCRIPTION:
"""
#Truncate 45 nt
def truncate(gtf):

    gtf['plus'] = gtf[[2,3,4,6,8]].apply(lambda x:
        (x[3] + 45) if x[2] == "exon" and x[3] + 45 <= x[4] and x[6] == "+" and "exon_number \"1\"" in x[8] else (
        "delete_this" if x[2] == "exon" and x[3] + 45 > x[4] and x[6] == "+" and "exon_number \"1\"" in x[8] else x[3]),axis=1)

    gtf['minus'] = gtf[[2,3,4,6,8]].apply(lambda x:
        (x[4] - 45) if x[2] == "exon" and x[3] <= x[4] - 45 and x[6] == "-" and "exon_number \"1\"" in x[8] else (
        "delete_this" if x[2] == "exon" and x[3] > x[4] - 45 and x[6] == "-" and "exon_number \"1\"" in x[8] else x[4]),axis=1)

    #remove exon1s that are too short
    gtf = gtf[~gtf['plus'].isin(['delete_this'])]
    gtf = gtf[~gtf['minus'].isin(['delete_this'])]

    #copy new coordinates back to original columns
    gtf[3] = gtf['plus']
    gtf[4] = gtf['minus']

    #remove placeholder columns
    gtf = gtf.drop(columns=['plus','minus'])
    return gtf

"""
DESCRIPTION
"""
def calculate_fc(data, label_comp, label_base):

    # Average every by cell line
    data['log2 Fold Change'] = np.log2((data.filter(regex=str(label_comp)).mean(axis=1)) / \
                                      (data.filter(regex=str(label_base)).mean(axis=1)))
    data['-log10 P-Value'] = ''

    return data

"""
DESCRIPTION
"""
def calculate_p(data, label_comp, label_base, drop_index):

    # Calculate p-value using 1-way ANOVA with replicates and append to df_oxsm_volc
    for row in data.iterrows():
        index, row_data = row
        comp_row = data.loc[index].filter(regex=str(label_comp)).values.tolist()
        base_row = data.loc[index].filter(regex=str(label_base)).values.tolist()

        # Append p_value to df_oxsm_volc
        try:
            statistic, p_value = stats.ttest_ind(comp_row, base_row)
            data.loc[index,'-log10 P-Value'] = float(-1 * (np.log10(p_value)))
        except:
            drop_index.append(index)

    data = data.drop(labels=drop_index, axis=0)

    return data

"""
DESCRIPTION: Parallelize function on a chunk of a dataframe
"""
def parallelize(func, *args):

    cores = cpu_count() #Number of CPU cores on your system
    partitions = cpu_count() #Define as many partitions as you want

    data_split = np.array_split(args[0], partitions)
    pool = Pool(cores)

    if len(args) == 4:
        func = partial(calculate_p, label_comp=args[1], label_base=args[2], drop_index=args[3])
    elif len(args) == 3:
        func = partial(calculate_fc, label_comp=args[1], label_base=args[2])
    else:
        return

    data = pd.concat(pool.map(func, data_split))

    pool.close()
    pool.join()

    return data

"""
DESCRIPTION: Reset plotting object to avoid bleed through
"""
def reset_plot(whitegrid, ax=False):

    if ax == True:
        del ax

    plt.close()
    plt.clf()

    if whitegrid == True:
        sns.set_style("whitegrid")
    else:
        sns.set_style("darkgrid")
