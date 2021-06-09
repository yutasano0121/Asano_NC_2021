import pandas as pd
import re
import os

workDir = os.getcwd()  # Or enter a path to your directory of choice.
inputDir = os.path.join(workDir, 'input')


# load the Immgen data. Make the column of the gene symbols the index.
df = pd.read_csv(
    os.path.join(inputDir, 'immgen_GSE15907_data.csv'),
    index_col=1
)
df = df.drop('ProbeSetID', axis=1)

# load the cell type annotation, which corresponds to columns.
anno_cellType = pd.read_csv(
    os.path.join(inputDir, 'immgen_annotation.csv')
)

# aggregate rows with a same gene symbols
df = df.groupby(df.index).agg('sum')

# fetch B cells in periphery
B = [col for col in df.columns.values if re.match(r"B.*Sp|B.*PC", col)]
df_B = df.loc[:, B]
df_B.to_csv(os.path.join(inputDir, 'immgen_Bcells.csv'))

# average values by cell types
cellTypes = [cellType for cellType in anno_cellType['cellType'].values if re.match(r"B.*Sp|B.*PC", cellType)]
df_B_mean = df_B.groupby(cellTypes, axis=1).agg('mean')

# filter out unexpressed genes
genefilter = df_B_mean.apply(max, axis=1) >= 100
df_out = df_B_mean.loc[genefilter, :]


# calculate correlation coefficient with AHNAK
corr_AHNAK = pd.DataFrame(df_out.corrwith(df_out.loc["Ahnak", ], axis=1))
corr_AHNAK.columns = ["correlationCoefficient"]

# save the results
df_out.to_csv(os.path.join(inputDir, 'immgen_Bcells_mean.csv'))
corr_AHNAK.to_csv(
    os.path.join(inputDir, 'immgen_Bcells_correlationWithAhnak.csv')
)
print("Result saved!")
