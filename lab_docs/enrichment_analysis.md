# Functional Enrichment Analysis within RERconverge
---

This section describes how to run pathway enrichment analysis using RERconverge output and functions included with the RERconverge package. Enrichment analysis detects groups of genes that are evolving faster or slower with a phenotype of interest. In the RERconverge package, the enrichment function is implemented as a Wilcoxon Rank-Sum Test on a list of genes ranked based on their correlation statistics. It detects distribution shifts in groups of genes compared to all genes, and it thereby bypasses the need to set a foreground significance cutoff like some other enrichment methods.

**Input** to the enrichment function is the output from RERconverge correlation functions and pathways of interest with gene symbols (gene names).

**Output** is enrichment statistics for each pathway, including genes in the pathways and their ranks.

## Import Pathway Annotations

Now that we have our gene statistics, we need pathway annotations. Download all curated gene sets, gene symbols (c2.all.v6.2.symbols.gmt) from [GSEA-MSigDB](http://software.broadinstitute.org/gsea/downloads.jsp) as gmtfile.gmt. You must register for an account prior to downloading. The rest of the vignette expects "gmtfile.gmt" to be in the current working directory.

With the file in the appropriate place, simply use the *read.gmt* function to read in the annotation data.
**Remember** in your own studies, that these could be many different forms of annotation and even custom annotations made for your project from functional genomics data.

```
# https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#C7
aC5bp <- read.gmt("annotations/c5.go.bp.v2023.2.Hs.symbols.gmt")
aC5hpo <- read.gmt("annotations/c5.hpo.v2023.2.Hs.symbols.gmt")
```

## Format Pathway Annotations

RERconverge enrichment functions expect pathways to be in named pathway-group lists contained within a list (see diagram above). A pathway-group is a group of similar pathways stored together (for example KEGG pathways might be one pathway-group and MGI pathways might be another pathway-group). Each pathway-group list contains a named list called *genesets* with each element of the list a character vector of gene names in a particular pathway, and the names of the elements the names of the pathways. The names of the genesets are also contained as a character vector in the second element of the pathway-group list named *geneset.names*. Each element of the *genesets* list is a character vector of gene names corresponding to a particular pathway, and the names of the elements in *genesets* are the names of the pathways. The *geneset.names* vector contains the names of the elements in *genesets*.

To convert our gmt file to the correct format, we simply need to put it inside another list; the *read.gmt* function automatically reads in gmt files in the format described above.

```
annotlist <- list( aC5bp, aC5hpo )
names(annotlist) <- c( "GOBP", "HumanPheno" )
```

---

## Calculate Enrichment Using *fastwilcoxGMTall*

We can now use *annotlist* and *stats* to calculate pathway enrichment. We will use the function *fastwilcoxGMTAll*, which calculates the estimated Wilcoxon Rank-Sum statistics based on the ranks of the genes in *stats*. The test essentially measures the distribution shift in ranks of genes of interest (each pathway) compared to the background rank distribution (ranks of all other genes included in the pathway-group annotations that are not part of the pathway of interest). We set *outputGeneVals* to true so the function returns names of the genes in the pathway and their ranks. *num.g* specifies the minimum number of genes required to be present to run the enrichment (10 by default).

```
# Calculate enrichments
enrichment=fastwilcoxGMTall( stats, annotlist, outputGeneVals=T, num.g=10 )

# To aid in sorting, we will add a new column 'signedLogP' of -log(P-value) with sign of the direction of enrichment (up or down).
enrichment$GOBP$signedLogP = -log10(enrichment$GOBP$pval) * sign(enrichment$GOBP$stat)
enrichment$HumanPheno$signedLogP = -log10(enrichment$HumanPheno$pval) * sign(enrichment$HumanPheno$stat)

# Inspect enrichments as such, clicking on column 'signedLogP' to sort ascending or descending.
View(enrichment$HumanPheno)
View(enrichment$GOBP)

# or when not using RStudio
head( enrichment$HumanPheno[ order( enrichment$HumanPheno$signedLogP, decreasing=TRUE ) , ] , n=20)
head( enrichment$GOBP[ order( enrichment$GOBP$signedLogP, decreasing=TRUE ) , ] , n=20)
```

Examine the resulting enrichments in both categories. Take time to sort in both positive and negative directions for the signedLogP column at far right.
- What is AMELOGENESIS_IMPERFECTA? Is this important?
- Why could it be that some other annotations not directly related to teeth appear "more significant"?

---
# Your future studies
Exploring these functional associations could also reveal new pathways and functions involved in your phenotype of interest.
