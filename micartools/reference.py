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
import numpy as np
import csv
import re
pd.options.mode.chained_assignment = None

from .utils import truncate, parallelize

"""
DESCRIPTION: Create a GTF reference file with only protein coding genes and the first n nucleotides of each first exon

RETURNS: Return the truncated, coding only table when running the function

METHODS: Considers strandedness
Multiprocesses chunks of a dataframe on cores of computer

VARIABLES:
gtf_file= GTF reference file path and name
truncate= Number of nucleotides to truncate from the 5' end of each exon 1 (considers strandedness)
save_file= Save truncated coding only GTF reference to file path and name
save_coding= Save coding only GTF reference to file path and name

ASSUMPTIONS:
Input file is a properly formatted GTF file
"""
def truncate_gtf(gtf_file, truncate=45, save_file=None, save_coding=None):

    #Import gtf reference file to
    if str(args_dict['input']).endswith('.gtf'):
        gtf = pd.read_csv(str(gtf_file), sep="\t", header=None, comment='#', low_memory=False)
    else:
        print('Error: A GTF-formatted file was not provided')
        return

    #Get only protein_coding coordinates
    gtf_coding = gtf[gtf.iloc[:, 8].str.contains('protein_coding') == True]

    #Save to .gtf file (tsv)
    if save_coding != None:
        gtf_coding.to_csv(str(save_coding), sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)

    gtf_coding_c = gtf_coding.copy()

    print("Multiprocessing reference chunks -- this may take a while...")
    gtf_truncated = parallelize(truncate, gtf_coding_c)

    if save_file != None:
        gtf_truncated.to_csv(str(save_file), sep='\t', header=None, index=False, quoting=csv.QUOTE_NONE)

    return gtf_truncated

"""
DESCRIPTION: Create a dictionary of gene IDs and common names from a GTF file

RETURNS: Dictionary of gene IDs and names

METHODS: Can be customized to create dictionary from any dataframe

VARIABLES:
file= Reference file path and name (IMPORTANT: If providing a csv or tsv, must be header-less)
file_type= Default 'gtf', other options: 'csv', 'tsv'
save_file= Save reference dictionary as a dataframe to file path and name (Only used for file_type='gtf')
gene_id= Position or column where gene_id information is found (default set for gtf)
gene_name= Position or column where gene_name information is found (default set for gtf)

ASSUMPTIONS:
Input file is a properly formatted GTF file
"""
def gene_dictionary(file, file_type='gtf', save_file=None, gene_id=0, gene_name=2):

    #Create dictionary from gtf file
    if file_type == 'gtf':
        gtf = pd.read_csv(str(file), sep="\t", header=None, comment='#', low_memory=False)

        #Get gene record rows from gtf
        gtf_genes = gtf.loc[gtf[2] == 'gene']

        #From the gene info column, split out gene_id and gene_name
        gtf_genes['gene_id'] = gtf[8].str.split(';').str[gene_id]
        gtf_genes['gene_name'] = gtf[8].str.split(';').str[gene_name]
        gtf_genes['length'] = abs(gtf[4] - gtf[3])
        gtf_genes['gene_id'] = gtf_genes['gene_id'].map(lambda x: x.lstrip('gene_id \"').rstrip('\"').rstrip(' '))
        gtf_genes['gene_name'] = gtf_genes['gene_name'].map(lambda x: x.lstrip('gene_name \"').rstrip('\"').rstrip(' '))

        #Create dictionary
        dict_df = gtf_genes[['gene_id','gene_name']].copy()
        gene_dict = pd.Series(gtf_genes.gene_id.values, index=gtf_genes.gene_name).to_dict()

        if save_file != None:
            dict_df.to_csv(str(save_file), sep='\t', index=False)

        return gene_dict

    #Create dictionary from dataframe where columns of file are gene_id and gene_name
    elif file_type == 'csv' or file_type == 'tsv':
        if file_type == 'csv':
            delimiter = ','
            data = pd.read_csv(str(file), sep=delimiter, header=None, comment='#', low_memory=False)
        elif file_type == 'tsv':
            delimiter = '\t'
            data = pd.read_csv(str(file), sep=delimiter, header=None, comment='#', low_memory=False)
        else:
            print('Error: Invalid file_type input')

        #Make sure file doesn't contain any spaces in names or IDs
        data[gene_id] = data[gene_id].map(lambda x: x.rstrip(' '))
        data[gene_name] = data[gene_name].map(lambda x: x.rstrip(' '))

        #Create dictionary
        dict_df = data[[gene_id,gene_name]].copy()
        dict_df.columns = ['gene_id','gene_name']
        gene_dict = pd.Series(data.gene_id.values, index=data.gene_name).to_dict()

        if save_file != None:
            dict_df.to_csv(str(save_file), sep=delimiter, index=False)

        return gene_dict

    else:
        print('Error: Invalid input provided')
        return

"""
DESCRIPTION: Rename .... maybe this is redundant. Needs to be tested anyways

VARIABLES:
gene_dictionary= MICARtools formatted gene_dictionary returned from gene_dictionary()
data= Dataframe where gene_id or gene_name is in a column
gene_id= Column number (int) or header name (str) in data of gene_id if converting gene_id to gene_name
gene_name= Column number (int) or header name (str) in data of gene_name if converting gene_name to gene_id

IMPORTANT:
gene_id and gene_name should not be used together (one should be None)

ASSUMPTIONS:
No column is provided as index
gene_dictionary keys are IDs and values are
"""
def apply_dictionary(gene_dictionary, data, gene_id=None, gene_name=None, save_file=None, delimiter=','):

    data_c = data.copy()

    if gene_id == None and gene_name != None:
        data_c[gene_name] = data_c[gene_name].map(lambda x: x.rstrip(' '))
        data_c['gene_id'] = data_c[gene_name].map(gene_dictionary)
    elif gene_name == None and gene_id != None:
        data_c[gene_id] = data_c[gene_id].map(lambda x: x.rstrip(' '))
        data_c['gene_name'] = data_c[gene_id].map(gene_dictionary)
    else:
        print('Error: Invalid gene_id or gene_name input provided')
        return

    if save_file != None:
        data_c.to_csv(str(save_file), sep=delimiter, index=False)

    return data_c
