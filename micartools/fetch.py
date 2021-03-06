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
import os, sys
import pandas as pd
import GEOparse
from .quality import clean_df

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
    data = pd.read_csv(str(file_name), sep=delimiter, index_col=0, header=0, low_memory=low_memory)

    #Check data orientation
    if str(gene_axis).lower() == 'row':
        data = data
    elif str(gene_axis).lower() == 'col':
        data = data.T
    else:
        print("Incorrect gene_axis option specified")

    data = clean_df(data)

    return data

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
    data = gse.pivot_samples('VALUE')
    del data.index.name

    data = clean_df(data)
    
    return data

"""
DESCRIPTION:
VARIABLES:
USAGE:
ASSUMPTIONS:
"""
def get_geo_info(geo_id):
    gse = GEOparse.get_GEO(geo_id) #Import GSE dataset

    with open('test.txt', 'w+') as f: #Create a text file for all the information for reference
        for gsm_name, gsm in gse.gsms.items():
            f.write(gsm_name + '\n')
            for key, value in gsm.metadata.items():
                f.write(" - %s : %s" % (key, ", ".join(value)) + '\n')

    for gsm_name, gsm in gse.gsms.items(): #Immediately print the most relevant information
        print("Name: ", gsm_name)
        for key, value in gsm.metadata.items():
            if key == 'title' or key == 'source_name_ch1' or key == 'characteristics_ch1' \
            or key == 'treatment_protocol_ch1' or key == 'data_processing':
                print(" - %s : %s" % (key, ", ".join(value)))

    dict = {}
    for gsm_name, gsm in gse.gsms.items(): #Create a dictionary for gsm identifier : sample description 
        for key, value in gsm.metadata.items():
            if key == 'title':
                dict[gsm_name] = value
    return dict



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
    info = pd.read_csv(str(file_name), sep=delimiter, header=None)

    #Reorganize dataframe as necessary so that data is column-wise and
    #sample_ids are the first column, labels are the second column
    if str(axis).lower() == 'col':
        if sample_ids == 0 and labels == 1:
            info = info
        else:
            info = info[[sample_ids, labels]] #Reorder dataframe columns
            info.columns = [0, 1]
    elif str(axis).lower() == 'row':
        info = info.T #Rotate dataframe to make it column-wise
        if sample_ids == 0 and labels == 1:
            info = info
        else:
            info = info[[sample_ids, labels]] #Reorder dataframe columns
            info.columns = [0, 1]
    else:
        print("Incorrect axis option specified")

    return info

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
        data_dropped = data.drop(ids, axis=1)

        return data_dropped

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
        data_dropped = data.drop(drop_ids_list, axis=1)

        return data_dropped

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
def keep_labels(data, info, label_list=None):

    if label_list == None:
        label_list = list(set(info[1]))

        keep_ids = info[info[1].isin(label_list)]
        keep_ids_list = list(keep_ids[0])

        #Drop samples not given in list to keep
        data_dropped = data[keep_ids_list]

    else:
        #Check file formats
        if type(label_list) is not list:

            return

        #Keep samples by name (will grab from info df)
        else:

            #Create list of sample_ids based on what is not provided in keep list
            drop_ids = info[~info[1].isin(label_list)]
            drop_ids_list = list(drop_ids[0])

            #Drop samples not given in list to keep
            data_dropped = data.drop(drop_ids_list, axis=1)

    return data_dropped

"""
DESCRIPTION: Compiles expression counts from multiple files into one table

VARIABLES:
delimiter= Delimiter style for expression files, will also output files if saved in this same format

ASSUMPTIONS:
File length of each is the same and ordered the same (same genes in the same order)
Files to parse are expected to be header-less and column[0] should be gene identifiers and column[1] should be expression values
"""
def catenate_files(directory, file_suffix='txt', gene_dictionary=None, sample_dictionary=None, save_file=None, delimiter='\t', drop_rows=0):

    #Walk through raw data files within given directory
    if directory[:-1] != '/':
        directory = directory + '/'

    file_list = []
    for subdir, dirs, files in os.walk(directory):
        for f in files:
            if f.endswith(file_suffix): #ignore hidden files and other uninterested files (particular for some submodules)
                file_list.append(f)
            else:
                pass

    #Sort files in alphabetical order (helps in formatting the count tables correctly)
    file_list = sorted(file_list)

    #get gene list from the first, as well as length(#rows)
    with open(str(directory) + str(file_list[0])) as f:
        gene_names = pd.read_csv(f, header=None, usecols=[0], dtype=str, sep=delimiter)
        row_num = len(gene_names) - drop_rows
        gene_names = gene_names[:-drop_rows]

    #populate dataframe with expression values
    data = pd.DataFrame(index=range(row_num))
    data['gene_names'] = gene_names

    for x in file_list:

        with open(str(directory) + x) as f:
            reader = pd.read_csv(f, header=None, usecols=[1], sep=delimiter)
            data[x] = pd.Series(data.index)
            length = len(reader)
            data[x] = reader

    #Change row names
    if gene_dictionary != None:
        data_set = rename_rows(data, gene_dictionary, 'gene_names')

    #Change column names
    if sample_dictionary != None:
        data_set = rename_cols(data, sample_dictionary)

    #Remove gene_names label
    data_set = data_set.rename({'gene_names': ''}, axis='columns')

    if save_file != None:
        data_set.to_csv(str(save_file),sep=delimiter)

    return data_set

"""
DESCRIPTION: Rename column names using dictionary

VARIABLES:
data= Dataframe to rename column names
dictionary= Keys are old name and values are new names
"""
def rename_cols(data, dictionary):

    data_c = data.copy()
    data_set = data_c.rename(columns=dictionary, inplace=True)
    return data_set

"""
DESCRIPTION: Rename values in a column (selected by providing column name) with a dictionary of keys and sort_values

VARIABLES:
data= Dataframe to rename row values
dictionary= Keys are old name and values are new names
label= Name of column to convert names; if 'index' is provided, will rename the index of the dataframe
"""
def rename_rows(data, dictionary, label):

    data_c = data.copy()

    if label == 'index':
        data_c['index'] = data_c.index
        data_c[label] = data_c['index'].replace(dictionary)
        data_c = data_c.set_index('index')

    else:
        data_c[label] = data_c[label].replace(dictionary)

    return data_c
