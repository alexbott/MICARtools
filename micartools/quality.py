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
import matplotlib as plt
import seaborn as sns

"""
DESCRIPTION: Check sample means and medians
METHODS: Output density plot for dataframe, barplot for each sample
VARIABLES:
USAGE:
ASSUMPTIONS:
Dataframe has been properly formatted so that probes or genes are rows and samples are columns
"""
def check_samples(data):

    wid = len(list(data))
    ax = data.boxplot(column=list(data), figsize=(wid,wid/3))
    ax.set_xlabel('Samples')
    ax.set_ylabel('Expression')

"""
DESCRIPTION: Cleans axis of NULL values
VARIABLES:
USAGE:
ASSUMPTIONS:
If dataframe has been properly formatted previously and genes are in rows, the default parameters will remove along the gene axis
"""
def clean_df(data, axis=0):

    data = data.dropna(axis=axis)
    return data
