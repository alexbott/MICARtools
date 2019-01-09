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
IMPORT DEPENDENCIES
"""
import sys
import pandas as pd
import GEOparse

"""
FUNCTIONS
"""

"""
DESCRIPTION: Get dataframe from user file
VARIABLES:
file_name= full path of file to import into pandas dataframe
delimiter= delimiter type for importing file, default: ','
low_memory= Specify memory limits for importing large files, default: False (allows for large imports)
gene_axis= Orientiation of the data, where categorical data is either column-wise, (default: 'col') or row-wise ('row'). Case insensitive
USAGE:
import micartools as mat
data = mat.get_df("~/Desktop/data.csv")
ASSUMPTIONS:
Dataset does not contain axis labels (i.e. a column header for 'gene names')
Dataset only has gene names and sample_ids as column headers and row indices. Orientation is flexible, but needs to be specified in options if genes are not rows
If orientation is not default, it is then specified or else function will not be able to properly format the dataframe for downstream application
"""
def get_df(file_name, delimiter=",", low_memory=False, gene_axis='row'):

    #Read in file
    df = pd.read_csv(str(file_name), sep=delimiter, index_col=0, header=0, low_memory=low_memory)

    #Check data orientation
    if str(gene_axis).lower() == 'row':
        df = df
    elif str(gene_axis).lower() == 'col':
        df = df.T
    else:
        print("Incorrect gene_axis option specified")

    return df

"""
DESCRIPTION: Get GEO dataframe
VARIABLES:
geo_id= GEO ID for dataset of interest, input is case insensitive (ex: GSE20716)
USAGE:
import micartools as mat
data = mat.get_geo("GSE20716")
"""
def get_geo(geo_id):

    #Import GSE dataset
    gse = GEOparse.get_GEO(geo=str(geo_id).upper())
    df = gse.pivot_samples('VALUE')
    del df.index.name
    return df

"""
DESCRIPTION: Get user file with pertinent sample information
VARIABLES:
file_name= full path of file to import into pandas dataframe
delimiter= delimiter type for importing file, default: ','
axis= Orientiation of the data, where categorical data is either column-wise, (default: 'col') or row-wise ('row'). Case insensitive
sample_ids= Column or row number where sample IDs are found (default: 0)
labels= Column or row number where categorical label data are found (default: 1)
USAGE:
import micartools as mat
sample_info = mat.get_info("~/Desktop/sample_info.csv")
ASSUMPTIONS:
Data categories are not labeled
If orientation is not default, it is then specified or else function will not be able to properly format the dataframe for downstream application
"""
def get_info(file_name, delimiter=",", axis="col", sample_ids=0, labels=1):

    #Read in file
    df = pd.read_csv(str(file_name), sep=delimiter, header=None)

    #Reorganize dataframe as necessary so that data is column-wise and
    #sample_ids are the first column, labels are the second column
    if str(axis).lower() == 'col':
        if sample_ids == 0 and labels == 1:
            df = df
        else:
            df = df[[sample_ids, labels]] #Reorder dataframe columns
            df.columns = [0, 1]
    elif str(axis).lower() == 'row':
        df = df.T #Rotate dataframe to make it column-wise
        if sample_ids == 0 and labels == 1:
            df = df
        else:
            df = df[[sample_ids, labels]] #Reorder dataframe columns
            df.columns = [0, 1]
    else:
        print("Incorrect axis option specified")

    return df

"""
DESCRIPTION: Drop samples by sample IDs -- pass in a list of names
VARIABLES:
data= Dataframe containing expression data
ids= List of sample IDs to remove from the dataframe
USAGE:
import micartools as mat
df = mat.drop_samples(df, sample_list)
ASSUMPTIONS:
Dataframe axes have been properly formatted (samples are columns, genes are rows)
"""
def drop_samples(data, ids):

    #Check file formats
    if type(ids) is not list:
        return

    #Drop samples in list
    else:
        df = data.drop(ids, axis=1)
        return df

"""
DESCRIPTION: Drop samples by label group name
VARIABLES:
data= Dataframe containing expression data
info= Dataframe containing sample information data
label= Name of sample type to drop (string)
USAGE:
import micartools as mat
df = mat.drop_label(df, sample_info, "WT")
ASSUMPTIONS:
Dataframe axes have been properly formatted (samples are columns, genes are rows)
Only one string is given to drop per call instance of function
"""
def drop_label(data, info, label):

    #Check file formats
    if type(label) is not str:
        return

    #Drop samples by name (will grab from info df)
    else:
        #Create list of sample_ids based on name provided
        drop_ids = info[info[1] == str(label)]
        drop_ids_list = list(drop_ids[0])

        #Remove these samples from data
        df = data.drop(drop_ids_list, axis=1)

        return df

"""
DESCRIPTION: Keep samples by list of label names
VARIABLES:
data= Dataframe containing expression data
info= Dataframe containing sample information data
labels= List of sample types to keep
USAGE:
import micartools as mat
df = mat.keep_labels(df, sample_info, ['normal','adenoma'])
ASSUMPTIONS:
Dataframe axes have been properly formatted (samples are columns, genes are rows)
Labels provided are in list format
"""
def keep_labels(data, info, label_list):

    #Check file formats
    if type(label_list) is not list:
        return

    #Keep samples by name (will grab from info df)
    else:
        #Create list of sample_ids based on what is not provided in keep list
        drop_ids = info[~info[1].isin(label_list)]
        drop_ids_list = list(drop_ids[0])

        #Drop samples not given in list to keep
        df = data.drop(drop_ids_list, axis=1)

        return df
