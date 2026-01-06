# Continuous Trait Analysis
## Association of branch RERs with Lifespan

---
  
RERconverge can also calculate correlations between rates of evolution of genes and the change in a *continuous* trait. To perform a continuous trait analysis, start with a named vector in R. Vector names must match the names in the trees read in previously, in this case the genome version names (*e.g.*, "hg38").

An interesting continuous trait, lifespan, is available for these species in `mammal108phenotypes`
```{r}
View(mammal108phenotypes)
# or
head( mammal108phenotypes )
```

We must convert the lifespan column into a named vector for use as a continuous trait in RERconverge, as below.
```
lifespan = mammal108phenotypes$Longevity.yrs
names(lifespan) = rownames(mammal108phenotypes)
```

We now convert the trait vector to paths comparable to the paths in the RER matrix using the function 'char2Paths' as shown here:

```
charpaths = char2Paths( lifespan, treesObj )
```

We are using `metric = diff`, which means that branch lengths are assigned to the trait tree based on the difference in trait values on the nodes connected to that branch. Other choices are available.

The `char2Paths` function creates a paths vector with length equal to the number of columns in the RER matrix. The phylogenetic relationships represented in the "char2Paths" output are the same as those represented in the RER matrix.

Finally, we use `correlateWithContinuousPhenotype` to find correlations between the rate of evolution of genes (in RER matrix `rers`) and the rate of change of a phenotype (encoded in `charpaths`). `correlateWithContinuousPhenotype` is a wrapper function for continuous trait analysis using the `getAllCor` function. By default it performs Pearson correlation. To perfom a rank-based correlation (Spearman) use the `getAllCor` function with `method="s"`. For Spearman, Winsorizing extreme values is not necessary for rank-based tests.\

This function uses the following input variables (all the options set here are also the defaults):
  
- `min.sp`: the minimum number of species in the gene tree for that gene to be included in the analysis. The default is 10, but you may wish to modify it depending upon the number of species in your master tree.
- `winsorizeRER`/`winsorizetrait`: pulls the most extreme N values (default N=3) in both the positive and negative tails to the value of the N+1 most extreme value. This process mitigates the effect of extreme outliers before calculating correlations.

```{r}
corLifespan = correlateWithContinuousPhenotype(rers, charpaths, min.sp = 10, winsorizeRER = 2, winsorizetrait = 2)

# As with binary trait, we set a signed log P-value for sorting
corLifespan$stat = -log10(corLifespan$P) * sign(corLifespan$Rho)

View(corLifespan)
# or
head( corLifespan[order(corLifespan$stat),] )
```

In these results, Rho is the standard statistic for a Pearson correlation, N is the number of branches included in the analysis, and P is the uncorrected correlation p-value. Since huge numbers of statistical tests are being performed in these analyses, it is essential to correct p-values using a method such as the Benjamini-Hochberg correction (p.adj).

---
# Enrichment Analysis for Lifespan Genes

Which biological functions are under additional constraint (negative Rho) or relaxed constraint (positive Rho) in long-lived species?

### If you have not yet run the enrichment lab activity, use it to load the annotation list!
[RERconverge's gene set enrichment analysis](https://github.com/nclark-lab/PhyloGenomicsLab/tree/main/lab_docs/enrichment_analysis.md)

```
# Calculate correlation statistic for genes versus lifespan
statLifespan = getStat(corLifespan)

# Calculate enrichments
enrichment=fastwilcoxGMTall( statLifespan, annotlist, outputGeneVals=T, num.g=10 )

# To aid in sorting, add a new column 'signedLogP' of -log(P-value) with sign of the direction of enrichment (up or down).
enrichment$GOBP$signedLogP = -log10(enrichment$GOBP$pval) * sign(enrichment$GOBP$stat)
enrichment$HumanPheno$signedLogP = -log10(enrichment$HumanPheno$pval) * sign(enrichment$HumanPheno$stat)

# Inspect enrichments as such, clicking on column 'signedLogP' to sort ascending or descending.
View(enrichment$GOBP)

# or when not using RStudio
head( enrichment$GOBP[ order( enrichment$GOBP$signedLogP, decreasing=TRUE ) , ] , n=20)
```


---
# Validating Genes with Plots
One important consideration of the results is the impact of one or a few species on the overall correlation. To assess this risk, we can examine individual correlation plots as follows:
  
```{r}
x=charpaths
y=rers['TET1',]
pathnames=namePathsWSpecies(treesObj$masterTree) 
names(y)=pathnames
plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Lifespan Change", 
     ylab="Evolutionary Rate", main="Gene TET1 Pearson Correlation",
     pch=19, cex=1)
text(x,y, labels=names(y), pos=4)
abline(lm(y~x), col='red',lwd=3)
```

In this case, we see that the positive correlation is driven by all species and not just a single clade. Note that labelled points are terminal nodes in the phylogeny and unlabelled points are internal nodes.
---

# Questions:

1. Which processes did you expect to be correlated with long lifespan *a priori*?
2. Which processes appeared in this comparative analysis?
3. List 2 other continuous traits you would like to examine for any set of species.
