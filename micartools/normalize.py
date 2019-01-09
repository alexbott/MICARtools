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




def prep_data(data, info, gene_axis='row'):

    #Make sure data is properly formatted
    if gene_axis == 'row':
        data = data
    elif gene_axis == 'col':
        data = data.T
    else:
        #Exit out of function
        print('Inappropriate gene_axis input provided')
        return

    #Import data and drop NAs
    orig_df = read_df(file)
    orig_df = orig_df.T
    orig_df = orig_df.dropna(axis=1)

    #Save gene list
    feat_cols = list(orig_df.columns[2:].values)

    #Drop labels
    if label == True and grade == False:
        df_scale1 = orig_df.loc[:, orig_df.columns != 'label']
    if label == False and grade == True:
        df_scale1 = orig_df.loc[:, orig_df.columns != 'grade']
    elif label == True and grade == True:
        df_scale = orig_df.loc[:, orig_df.columns != 'label']
        df_scale1 = df_scale.loc[:, df_scale.columns != 'grade']
    else:
        pass

    #Convert data to float
    df_scale2 = df_scale1.astype(dtype='float')

    #gene normalization (column)
    df_scale2[df_scale2.columns] = preprocessing.scale(df_scale2[df_scale2.columns])

    #Format dataframe for downstream
    df_scale2 = df_scale2.T

    return df_scale2, orig_df, feat_cols
