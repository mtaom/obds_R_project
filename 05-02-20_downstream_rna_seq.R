## We are now in the git-tracked repository
## all paths are relative to it

##gitignore - things you don't want git to track
## is a text file -> write the file names you don't git to track
## .Rproj --> will have project options, packrat things, etc.

## rhs panel -- Git tab --> can choose which files you wan to commit, write a 
# message.etc --> committed locally on your PC --> still need to push to your
# GitHub

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install('DESeq2')
BiocManager::install('rhdf5')
## having troubles as packagfes was not available as binary for Windows
## installed from source
BiocManager::install('GenomeInfoDbData', type = 'source')


library('tximport')
library('DESeq2')
library('rhdf5')
library('tidyverse')
library('biomaRt')

## loading in information about the sample
samples <- read_tsv('data/obds_sampletable.csv')
samples

## stealing these from our abundace file, a list of transcripts we are
# interested in
transcripts_of_interest <- read_tsv('data/pseudoaligned/ERR1755082/abundance.tsv')
transcripts_of_interest <- transcripts_of_interest$target_id

## need to create pair of IDs and genes
## taking these from the EnsDb package we installed
# mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# # then pulling out specific columns that we are interested in #
# annot_genes <- getBM(attributes = c("ensembl_transcript_id_version", "ensembl_gene_id", "mgi_symbol"),
#                      values = transcripts_of_interest, mart = mouse)
### run the first bit and then write it into a file that we read in below
## biomaRt crushes a lot, so saves a lot of trouble!!
# write_tsv(annot_genes, 'data/annot-genes.tsv')

## reading in our transcript annotation
annot_genes <- read_tsv('data/annot-genes.tsv')
## making new annotation dataframe (only has 2 columns now!)
annot_genes_filtered <- annot_genes[1:2]

## steps to get to importing gene-level abundancies in the end
files <- file.path('data/pseudoaligned', samples$Sample_accession, "abundance.h5")
## renaming files
names(files) <- samples$Sample_accession
## here importing all of our files, they are kallisto type; txOut -- make sure also to get the gene info, not transcript!
## tx2gene -- we are giving it the file we created for it (needs two columns transcript ID and gene ID)
txi.kallisto <- tximport(files, type = "kallisto", txOut = FALSE, tx2gene = annot_genes_filtered)
head(txi.kallisto$counts)


## $ makes a vector; [] make a dataframe
### making a sample table - which is a dataframe containing genotype of mice
## and then adding our row names (which are our sample names -- Sample_accession)
###### also converting it into a factor! wa want the column sample_title to 
# be a factor - will be easier for further downstream analysis
## not factor anymore - need to get rid of the _123 (replicate) part first
## going to use gsub
sampleTable <- data.frame(condition = samples$sample_title)
rownames(sampleTable) <- samples$Sample_accession

######Converting to a DESEQ2 object ########

## for Deseq2 need to specify that the values in the column are factors
## gsub - issue in our condition column (would make 12 factors because had  _123)
## now we removed them AND made it into factor
## last command to check what the factors are
sampleTable$condition <- gsub('_[123]', '',sampleTable$condition)
## need to also gsub 2/3 -- deseq2 doesnt like it!
## gsub(paramter, replacement, what on/data/x)
sampleTable$condition <- gsub('/', '_',sampleTable$condition)
## have 4 conditions thus will end up with 4 levels/factors
sampleTable$condition <- as.factor(sampleTable$condition)
levels(sampleTable$condition)

### Makign a Deseqdataset (dds) ####
## needed to be able to run Deseq2
dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)


#### DESEQ2 ###### 
# normalisation  - estimateSizeFactors
# variance estimation - tries to normalise across samples
# identify DEGs nbinomWaldTest 
 ### all three are wrapped in one fucntion called DESeq function

##minimal pre-filtering, keep only rows that have more than 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)
#res <- results(dds) ## does a contrast of its own - alphabetical order

## specify the groups you want to do DGE on 
res <- results(dds, contrast=c("condition","Egr2_Kin_CD4","Egr2_Kin_CD8"))


## volcano plot & pheatmap
## x axis should have log2 and -log10(padj) on y axis

## need to export the results into a table so can read it in and plot
## res - a big object! access next level with @
## ensure it keeps row names/geneID
res_df <- data.frame(res, row.names = res@rownames)

## removing NAs - keep all the rows where the padj column is not NA
res_df <- res_df[!is.na(res_df$padj),]


### losing my row names when using mutate :(
#res_df <- res_df %>% mutate(neg_log10 = -log10(padj))
## making a new column called neg_log10 and assigning what it should be
res_df$neg_log10 <- -log10(res_df$padj)
## ordering by increasing pvalue
resOrdered <- res_df[order(res_df$pvalue),]
summary(resOrdered)

## keep only results that are significant(ie padj < 0.05)
res_significant <- subset(resOrdered, padj < 0.05) 
## think about also changing the value of p where it is equal to zero!!
#negative log10 of values will be infinity and thus plotted like crazy
 
## volcano plot ## plot all genes
##logical vector, if true, output TRUE - will be used to colour by
## can ask to output whatever! e.g. instead of TRUE -  significantly DEGs
## this is our threshold for significant genes (significatnly DEGs)
# absolute value of log2fold change - to have negative fold change positive
threshold <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, TRUE, FALSE)

ggplot(res_df, aes(x=log2FoldChange, y=neg_log10, colour = threshold)) +
  geom_point() +
  scale_y_continuous(limits = c(0,500))

##what are the dots up there?? INFINITY!
highs <- subset(res_df, neg_log10 > 350)
## really significantly DEGs - consider changing their p value from
## zero to 0.00000000000000000000001 - yet for plotting it moves it all
## down to the same place and it looks a bit silly - edgeR - plots better?

## heatmap


