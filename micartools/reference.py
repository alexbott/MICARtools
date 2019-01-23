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
from multiprocessing import cpu_count, Pool

from .utils import truncate, parallelize

"""
DESCRIPTION: Create a GTF reference file with only protein coding genes and the first n nucleotides of each first exon

METHODS: Considers strandedness
Multiprocesses chunks of a dataframe on cores of computer

VARIABLES:

ASSUMPTIONS:
"""
def truncate_gtf(gtf_file, truncate=45, save_file=None, save_coding_only=None, delimiter='\t', return_dataframe=False):

    #Import gtf reference file to
    if str(args_dict['input']).endswith('.gtf'):
        gtf = pd.read_csv(str(gtf_file), sep="\t", header=None, comment='#', low_memory=False)
    else:
        print('Error: A GTF-formatted file was not provided')
        return

    #Get only protein_coding coordinates
    gtf_coding = gtf[gtf.iloc[:, 8].str.contains('protein_coding') == True]

    #Save to .gtf file (tsv)
    if save_coding_only != None:
        gtf_coding.to_csv(str(save_coding_only), sep=delimiter, header=None, index=False, quoting=csv.QUOTE_NONE)

    gtf_coding_c = gtf_coding.copy()

    print("Multiprocessing reference chunks -- this may take a while...")
    gtf_truncated = parallelize(gtf_coding_c, truncate)

    if return_dataframe == True:
        return gtf_truncated

    if save_file != None:
        gtf_truncated.to_csv(str(save_file), sep=delimiter, header=None, index=False, quoting=csv.QUOTE_NONE)
