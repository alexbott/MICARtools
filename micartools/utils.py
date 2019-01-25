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

"""
DESCRIPTION:
VARIABLES:
USAGE:
ASSUMPTIONS:
"""
def custom_list(file):
    gene_list = []

    with open(file, 'r') as f:
        for line in f:
            gene_list = line.split(",")

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
DESCRIPTION: Parallelize function on a chunk of a dataframe
"""
def parallelize(data, func):

    cores = cpu_count() #Number of CPU cores on your system
    partitions = cpu_count() #Define as many partitions as you want

    data_split = np.array_split(data, partitions)
    pool = Pool(cores)
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
