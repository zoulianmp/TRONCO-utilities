# Load Tronco files
source('../visualization/oncoprint.R')
source('../correctness/is.compliant.R')
source('../loading/import.genotypes.R')
source('../clustering/subtypes.split.R')

# Example with TCGA data for colorectal cancer

# Load all dataset (boolean matrix, gene names, clustering indexes, sample ids + stages)
alldata = read.table('./gene_indiv_mat')
names = read.table('./gene_id_symbol')
idxclust = read.table('./idx.k4.txt')
samples = read.table('./sample_id')
stages = read.table('sample_stage')

# assign column and row names
rownames(alldata) = unlist(samples)
colnames(alldata) = unlist(names)
rownames(idxclust) = unlist(samples)

# Events (randomly sampled)
n = c('KRAS', 'TP53', 'PTEN', 'ARID1A', 'BRAF', 'APC')

alldata = alldata[, n]

tcga = import.genotypes(alldata, stage.annot=stages)
is.compliant(tcga)

# Extract only cluster 1
x = subtypes.split(tcga, idxclust, idx=1)
is.compliant(x)

# Work with x..