#setwd
setwd('/Users/josoga2/Documents/wale_docs/phd/HB/bash/module_10/') #adjust to your own pc

#libraries
library(DESeq2)
library(pheatmap)

#count file
c_e_count <- read.delim('counts.txt', header = T)
c_e_meta <- read.delim('metadata.tsv', header = T, stringsAsFactors = T)

#preview
head(c_e_count)
head(c_e_meta)

#keep to the important columns
raw_counts <- c_e_count[, c_e_meta$sample]
head(raw_counts)

#add rownames
rownames(raw_counts) <- c_e_count$Geneid

head(raw_counts)

#create desqdataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = c_e_meta,
                              design = ~gender)

#preview
dds
dds$sample
dds$gender

#Perform the differential expression analysis
dds <- DESeq(dds)

final_res <- results(dds)

#look at your result
head(final_res)

#we have a truncated data, let's see the distro of p-values
plot(density(x = na.omit(final_res$pvalue)))

#ok let's look at our differentially expressed genes
plot(x = final_res$log2FoldChange, 
     y = -log10(final_res$padj),
     cex = 0.25,
     pch = 19, 
     col = 'grey',
     ylim = c(0,20),
     ylab = 'Adjusted P-Value',
     xlab = 'Log2 FC')

abline(v = c(-2, 2), h = -log10(0.05), lwd = 0.5, lty = 2)
#where are the upregulated
upregulated <- subset(final_res, padj < 0.05 & log2FoldChange > 2)
points(upregulated$log2FoldChange,
       y = -log10(upregulated$padj), 
       cex = 0.35,
       pch = 19,
       col = 'salmon')

#where are the downregulated
downregulated <- subset(final_res, padj < 0.05 & log2FoldChange < -2)
points(downregulated$log2FoldChange,
       y = -log10(downregulated$padj), 
       cex = 0.35,
       pch = 19,
       col = 'lightblue')

mtext('A simple volcano')
#we can merge the two to do a clean and less memory efficient heatmap
degs <- rbind(raw_counts[rownames(upregulated),], 
              raw_counts[rownames(downregulated),])
pheatmap(degs, 
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = F,
         scale = 'row',
         show_colnames = T)

#what are the genes that are upregulated
rownames(upregulated)
rownames(downregulated)

#exporting the files
write.csv(upregulated, 'upregulated.csv')
write.csv(downregulated, 'downregulated.csv')
write.csv(raw_counts, 'raw_counts.csv')

#Functional Enrichment Analysis
# Visit https://bioinformatics.sdstate.edu/go/
