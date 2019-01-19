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
import pandas as pd
from sklearn import preprocessing
from .utils import custom_list

"""
DESCRIPTION: Create dictionary of probes and their gene names for the microarray probe collapser

VARIABLES:
reference= Full path and file name for GPL reference file, accessed from NCBI (should be a .txt file)
gene_list= Full path and file name to .csv file listing gene names to get probes for. Resulting function dictionary will only contain these genes and their probes (default: None)
no_multimappers= Do not allow ambiguous probes in the probe-gene dictionary (default: True)

USAGE:
import micartools as mat
gene_dict = mat.prep_collapser("~/Desktop/GPL570.txt")

ASSUMPTIONS:
If using gene_list option, file must be a .csv
Assumes GPL .txt file from NCBI is tab delimited
"""
def prep_collapser(reference, gene_list=None, no_multimappers=True):

    #Get probe/gene_name reference file and import into a dataframe
    df = pd.read_csv(str(reference), sep="\t", low_memory=False, comment='#') #no index
    df = df[['ID','Gene Symbol']]
    df = df.dropna()

    #If a custom gene list is given so that only certain probes are looked at and collapsed downstream
    #Take the gene list file, join genes for regex search
    #Create dictionary for just the probes corresponding with the genes of interest
    if gene_list != None:
        if type(gene_list) == list:

            gene_list = custom_list(gene_list)

            search = '|'.join(gene_list)
            search = search.upper()

            df = df[df['Gene Symbol'].str.contains(search)] #only df where probes of interest are
            df_dict = df.set_index('ID')['Gene Symbol'].to_dict()
    #Take full probe set and create dictionary
    else:
        df_dict = df.set_index('ID')['Gene Symbol'].to_dict()

    #If allowing for multimappers, return dictionary as is
    if no_multimappers == False:
        return df_dict
    #Remove any probes that map to several genes
    else:
        df_dict_mod = {}
        for key, value in df_dict.items():
            if '///' in value:
                continue
            else:
                df_dict_mod[key] = value

        return df_dict_mod

"""
DESCRIPTION: Collapse probes of microarray dataset using previously prepared probe collapser dictionary (see prep_collapser function))

METHODS: Works by taking all probe expression data of a particular gene and determines a mean expression value for each sample for the given gene

VARIABLES:
data= Dataframe of microarray probe data to be collapsed
dict= Probe collapser dictionary created in the prep_collapser function

USAGE:
import micartools as mat
df_collapsed = mat.probe_collapse(df, collapser_dict)

ASSUMPTIONS:
A probe collapser dictionary has been previously prepared using the prep_collapser function
"""
def probe_collapse(data, dict):

    print('Presumably, a SettingwithCopyWarning error will now appear. Ignore this warning, the collapser is working properly.')

    #Get list of probes to find in df
    dict_list = list(dict.keys())

    #only keep df rows where probes of interest are
    #(in cases where only looking at certain genes, default: all genes)
    #column 'name' header for probes in these files
    try:
        data = data[data['name'].isin(dict_list)]
    except:
        data['name'] = data.index
        data = data[data['name'].isin(dict_list)]

    #Map gene names in place of probes
    #Will set off SettingwithCopyWarning -- should be fine as we are replacing values in same column and it appears to change as expected
    data['name'] = data['name'].map(dict)

    #Set gene names as indices (allows for multiple indices with same name)
    #Needed to remove strings from df to allow for next step
    data = data.set_index('name', drop=True)

    #force data to float
    data_numeric = data.apply(pd.to_numeric)
    #Reset indices to its own column to allow for sorting in next step
    data_numeric['name_sort'] =  data_numeric.index

    #groupby index (gene name) and collapse, taking the mean of rows with same name in 'name_sort'
    #Sets name_sort column names (post-collapse) as indices
    data_collapsed = data_numeric.groupby('name_sort').mean()

    #Remove double header for indices
    del data_collapsed.index.name

    return data_collapsed

"""
DESCRIPTION: Ties prep_collapser and probe_collapse functions together

VARIABLES:
data= Dataframe of microarray probe data to be collapsed
reference= Full path and file name for GPL reference file, accessed from NCBI (should be a .txt file)
gene_list= Full path and file name to .csv file listing gene names to get probes for. Resulting function dictionary will only contain these genes and their probes (default: None)
no_multimappers= Do not allow ambiguous probes in the probe-gene dictionary (default: True)

USAGE:
import micartools as mat
df_collapsed = mat.auto_collapse(df, "~/Desktop/GPL570.csv")

ASSUMPTIONS:
See assumptions for prep_collapser and probe_collapse functions
"""
def auto_collapse(data, reference, gene_list=None, no_multimappers=True):

    dict = prep_collapser(reference, gene_list=gene_list, no_multimappers=no_multimappers)
    data_collapsed = probe_collapse(data, dict)

    return data_collapsed

"""
DESCRIPTION: Normalize samples, prints sample axis means for verification

METHODS: For each sample axis, divide each cell by the sum the axis divided by the factor provided (default: 1e6)

VARIABLES:
data= Dataframe of microarray probe data
axis= Axis where samples are found in the dataframe
factor= Numeric value to scale samples by
print_means= Print appropriate means that were scaled for verification

USAGE:
import micartools as mat
df_norm = mat.sample_norm(df)
"""
def sample_norm(data, axis=1, factor=1e6, print_means=False):

    #Initialize axis variables based on user input
    if axis == 0:
        axis_2 = 1
    elif axis == 1:
        axis_2 = 0
    else:
        pass

    #Perform normalization
    data_norm = data.divide((data.sum(axis=axis_2) / float(factor)),axis=axis)

    if print_means == True:
        print(data_norm.mean(axis=axis_2))

    return data_norm

"""
DESCRIPTION: Prepare dataframes for analysis plotting functions found within analyze.py

METHODS:
Original dataframe is unformatted besides adding labels from info to the first row of the dataframe
Formatted dataframe is scaled if option provided and dataframe is converted to float

VARIABLES:
data= MICARtools formatted dataframe of expression values
info= MICARtools formatted sample info dataframe
gene_scale= Scale genes (rows) of data
print_means= Print appropriate means that were scaled for verification

USAGE:
import micartools as mat
df_scaled, df_collapsed = mat.prep_df(df_collapsed, info)

ASSUMPTIONS:
Requires properly formatted df and info dataframes for MICARtools usage
"""
def prep_df(data, info, gene_scale=True, print_means=False):

    #Convert data to float
    data_scaled = data.astype(dtype='float')

    #gene normalization
    if gene_scale == True:
        data_scaled[data_scaled.columns] = preprocessing.scale(data_scaled[data_scaled.columns], axis=1)

    if print_means == True:
        print(data_scaled.mean(axis=1))

    #Map labels to samples
    labels = pd.Series(info[1].values,index=info[0]).to_dict()
    data.loc['label'] = data.columns.map(labels.get)

    #Output collapsed dataframe
    newIndex = ['label'] + [ind for ind in data.index if ind != 'label']
    data = data.reindex(index=newIndex)

    return data_scaled, data
