---
# Convergent Evolution of Tooth Enamel Loss in Mammals
## An RERconverge Demonstration
### Nathan Clark
### University of Pittsburgh
### 2026 January 21
---

## Overview

Convergent evolution provides a powerful natural experiment for identifying the genetic basis of complex traits. Selective pressures are constantly shifting and their action on individual genes changes over evolutionary time as species encounter new conditions.
As a particular gene becomes more important for fitness, its evolutionary rate will slow due to increased constraint and rejection of nearly neutral mutations. Conversely, if a gene is less important its rate will increase as some deleterious mutations fix in the species.
Finally, if positively selected changes to a gene accompany changing conditions in a species, the gene's evolutionary rate will similarly increase.
Thus, studying shifts in relative evolutionary rates (RERs) between species can reveal genes which are responding to selective pressures.
However, many genes will change rates in every species, so our ability to infer such connections relies on repeated evolutionary events, convergent evolution, so that we can associate a gene's shift with a specific selective pressure, phenotype, or environmental change.

One convergent evolutionary example in mammals is the **independent loss of tooth enamel**, which has occurred in multiple lineages including **pangolins**, **anteaters**, **armadillos** and **baleen whales**.
In this workshop, I demonstrate how to use **`RERconverge`** to identify genes whose evolutionary rates show **convergent shifts** associated with enamel loss. In other published contexts, we have similarly identified genes responding to convergent transitions of mammals to an aquatic and subterranean life, long lifespan, loss of hair, and high altitude.

---

## Biological Background
Tooth enamel is the hardest tissue in vertebrates and requires a coordinated developmental program involving genes such as:
- *ENAM*
- *AMELX*
- *AMBN*
- *AMTN*
- *ANKRD11*
- *DSPP*
- *MMP20*

Lineages that no longer form mineralized teeth often exhibit **relaxed selection** or **pseudogenization** of these genes. RERconverge allows us to detect these patterns genome-wide by searching for those genes that consistently accelerate their RERs in toothless or enameless species.

---
## Data Requirements

This study requires:
1. A **species phylogeny** covering all taxa of interest
2. A set of **gene trees** (one per gene, same taxa)
3. A **binary phenotype** indicating enamel presence or absence

All trees are assumed to be in *Newick format*.

---
### Load or Install [R](https://www.r-project.org/) and possibly [R Studio](https://posit.co/download/rstudio-desktop/)
`module load r/4.5.0`

### Install packages (run once)
```{install r libraries}
install.packages("devtools")
install.packages("ape")
install.packages("phytools")
install.packages("phangorn")

install.packages("BiocManager")
BiocManager::install("ggtree")
BiocManager::install("impute")

devtools::install_github("nclark-lab/RERconverge")
```

---
### Load Required Libraries

```{r libraries}
# Load RERconverge package. Dependencies will load automatically.
library(RERconverge)
```

---
## 108 mammalian species in the study
More than 400 mammalian species genomes are available today, the workshop uses a carefully chosen set of ~100 to improve runtime and to contain specific phenotypes for this study.
Read in the phenotypes table and examine the species using `read.csv()` and `View()` (only available in RStudio).

```
mammal108phenotypes <- read.table("mammal108phenotypes.tsv", header=TRUE, sep="\t")
View(mammal108phenotypes)
```

## Species Phylogeny

**Mammalian species phylogeny used for RERconverge analysis.**  
Load and visualize the 108 mammal species tree we will use for the workshop.
The **`ape`** package provides the `read.tree()` function to read the tree from the "speciesTree108.nwk" file into a phylo object.
Then, `plot()` will automaticaly detect the phylo object and plot the phylogram.
This tree represents the shared evolutionary history of the mammalian taxa included in the analysis. Branch lengths reflect neutral evolutionary divergence and serve as the baseline against which gene-specific rate shifts are measured.

```{r species-tree}
# Load the species tree
speciesTree <- read.tree("speciesTree108.nwk")

# Plot tree basic, with genome version names
par(mar=c(1,1,2,1)+0.1) # Change margin to make title fit.
plot(speciesTree, cex=0.5, edge.width = 0.5, main="Mammalian Species Tree", align.tip.label=F)

# Plot tree, with common names
speciesTreeCommon <- speciesTree
speciesTreeCommon$tip.label <- mammal108phenotypes[speciesTreeCommon$tip.label,"common_name"]
plot(speciesTreeCommon, cex=0.5, edge.width = 0.5, main="Mammalian Species Tree", align.tip.label=F)

```
![A phylogeny of 108 mammal species](images/speciesTree.jpg)

---
# Preparing gene trees
Next, the study requires multiple sequence alignments (MSAs) of many orthologous gene sets from the common set of species  in the species phylogeny.
Not all species are required for each gene, RERconverge can handle sparse data, since orthologous gene sets always experience gains and losses naturally.

There are multiple potential sources for orthologous gene MSAs:
1. Download MSAs from projects such as [Zoonomia](https://zoonomiaproject.org/), the [Vertebrate Genomes Project](https://vertebrategenomesproject.org/), EMSEMBL, etc...
2. Extract genes from a whole-genome alignment of many species in HAL or MAF format.
3. Group and Align annotated genes from your species. Suggested to use OrthoFinder followed by OrthoSnap.
4. Download [Clark lab gene trees](https://github.com/nclark-lab/ComparativeData/wiki). (Alignments also available).

We will begin using a demonstration set of MSAs available in the **`alignments`** directory.

# Examining a multiple sequence alignment with ggmsa
## Lens Instrinsic Membrane protein 2 (LIM2) encodes a protein important in lens function and hence vision.
Which species have the most amino acid changes? Why these species?

[alignments.zip](https://pitt-my.sharepoint.com/:f:/g/personal/nclark_pitt_edu/IgBYR3CPWrWjQZ5mfrloQ8hpARMYAm84WkMXLu786CbPy_g?e=5y5NbH)

[SeaView](https://doua.prabi.fr/software/seaview) alignment program

# BiocManager::install("Biostrings")
# library(Biostrings)
# BiocManager::install("DECIPHER")
# library(DECIPHER)
```{r LIM2 alignment}
alnLIM2 <- readAAStringSet("alignments/LIM2.mfa")
colors <- c(`-`="#000000", `A`="#BDB1E8", `R`="#EFA2C5", `N`="#F6602F",
+             `D`="#FD5559", `C`="#12C7FE", `Q`="#DDACB4", `E`="#FEA097", `G`="#F46802",
+             `H`="#FCA708", `I`="#369BD9", `L`="#2E95EC", `K`="#CF7690", `M`="#4B8EFE",
+             `F`="#76997D", `P`="#FD2AE3", `S`="#A08A9A", `T`="#9A84D5", `W`="#74C80D",
+             `Y`="#9BB896", `V`="#89B9F9")
BrowseSeqs(alnLIM2,colors=colors,patterns=names(colors))
```


---
# Calculating branch lengths
RERconverge includes tree-building functions that perform maximum likelihood branch length estimation given a fixed tree topology and alignments for each sequence of interest. These functions are built directly on [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html) functions `pml` and `optim.pml`, including arguments for parameters passed directly to those functions. For more details on those functions, refer to [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html) documentation.

## Input file specification
Tree building functions require two inputs: a master tree topology and alignments from which to estimate branch lengths.

```{r results='hide', message = FALSE, warning = FALSE}

estimatePhangornTreeAll(alndir=alignmentfn, treefile=mastertreefn, output.file=outputfn)

```


---

# Phenotype Definition: Tooth Enamel Loss

**Figure 2. Independent losses of tooth enamel across mammals.**  
Tooth enamel has been lost independently in several mammalian lineages. These convergent events form the biological basis of this analysis, allowing us to search for genes showing repeated evolutionary rate shifts.

We encode enamel loss as a binary trait:

- **1** = enamel absent
- **0** = enamel present

```{r phenotype}
# Binary phenotype vector
enamelTrait <- c(
  Homo_sapiens       = 0,
  Mus_musculus       = 0,
  Canis_lupus        = 0,
  Bos_taurus         = 0,
  Loxodonta_africana = 0,
  Manis_javanica     = 1,  # pangolin
  Myrmecophaga_tridactyla = 1,  # anteater
  Balaenoptera_musculus   = 1   # baleen whale
  Armadillo = 1 # armadillo
)

# Ensure correct ordering
enamelTrait <- enamelTrait[speciesTree$tip.label]

knitr::kable(data.frame(Species = names(enamelTrait), EnamelAbsent = enamelTrait))
```

---
## Reading in gene trees with `readTrees`

To run RERconverge, you will first need to supply a file containing **gene trees** for all genes to be included in your analysis. This is a tab delimited file with the following information on each line:

Gene_name Newick_tree

Now in R, read in the gene trees. The `useSpecies` input variable can be provided to most RERconverge functions. Excluding one or more species from this vector will exclude them from the analyses. We leave `useSpecies` null here, but you can feed in a list of species you want to limit the trees to if you desire:

```{r, cache = TRUE}
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), useSpecies=NULL)
```

First, the code tells us that there are 500 items, or gene trees, in the file. Then it says that the maximum number of tips in the gene trees is 62 and, later, it reports that it will use the 71 genes in this set that have data for all 62 species to estimate a **master tree**. The master tree will be used for subsequent analyses.

RERconverge is intended to be used on genome-scale datasets, containing a large number of gene trees with data present for all species. It thus has a minimum number of such gene trees required for `readTrees` to use to estimate a **master tree**; this is set with the `minTreesAll` option and is 20 by default. If your dataset is smaller, you may adjust this or supply your own master tree using the option `masterTree` (this should be a `phylo` object generated using ape's `read.tree`); however, we recommend interpreting results with caution in this case. If you want to read in less trees than your whole dataset for time purposes, set the `max.read` argument to however many gene trees you want it to read in.

---
## Estimating relative evolutionary rates (RER) with `getAllResiduals`

The next step is to estimate **relative evolutionary rates**, or RERs, for all branches in the tree for each gene. Intuitively, a gene's RER for a given branch represents how quickly or slowly the gene is evolving on that branch relative to its overall rate of evolution throughout the tree.

Briefly, RERs are calculated by normalizing branch lengths across all trees by the master branch lengths. Branch lengths are then corrected for the heteroskedastic relationship between average branch length and variance using weighted regression. For a more detailed description of how RERs are computed, see [@Chikina2016] and [@Partha2017].

We will use the `getAllResiduals` function to calculate RERs. This uses the following input variables (all the options set here are also the defaults):

-   `useSpecies`: a vector that can be used to specify a subset of species to use in the analysis. Here we will use the species in our AdultWeightLog vector that will be used for continuous trait analysis. Note that these are also the same species used for binary trait analysis. These species should be a subset of species included in `toyTrees$masterTree$tip.label`.
-   `transform`: the method used to transform the raw data. By transforming the raw data, we reduce the heteroscedasticity (relationship between mean and variance) and the influence of outliers. Here we will use a square-root transform ("sqrt"), which has performed the best at reducing heteroskedasticity in our datasets. Also available are "none" (no transformation) and "log" (natural logarithm transformation) and "asinh" (inverse hyperbolic sine transformation).
-   `n.pcs`: Number of principal components to normalize by (default: 0, mean normalization).
-   `use.weights`: whether to use a weighted regression to estimate RER. Weighting allows further correction for the relationship between mean and variance, which can be directly estimated from the data.
-   `weights`: Manual weights for a weighted regression. If these aren't provided, it uses ones generated in `readTrees`.
-   `norm`: A character string specifying the normalization method. Options include "scale," "zscore," or "quantile" (default is "scale").
-   The documentation provides other specific parameters for those who want to fine tune their analysis.

Here is the basic method, with the recommended settings:

```{r, cache = TRUE}
data("logAdultWeightcm")
mamRERw = getAllResiduals(toyTrees,useSpecies=names(logAdultWeightcm), 
                          transform = "sqrt", n.pcs = 0, use.weights = T,
                          weights=NULL,norm="scale")
```

Part of the output of this function tells you that the cutoff is set to 0.05. Any branches shorter than this will be excluded from the analysis. It then  calculates relative evolutionary rates for sets of gene trees.

The plots generated by this function show the log variance of the RERs resulting from the original method (on the left) and the variance after transformation and weighted regression (on the right). Notice the heteroscedasticity (positive trend between values and their variance) that is present before the new transformation method is applied, is now gone. The x-axis displays bins of branch lengths on the tree, and the y-axis is the (log-scaled) variance in these branch lengths across trees. As you can see by comparing the right plot to the left plot, transforming and performing a weighted regression reduces the relationship between the mean branch length (x-axis) and the variance in branch length (y-axis). You can alter values for `transform`, `n.pcs`, `use.weights`, and `norm` to attempt to optimize heteroskedasticity correction.

If you wish to save this RER object for later, you can use R's `saveRDS` function. This will allow you to load it later with `readRDS`, using a different name, if you wish.

```{r, message = FALSE, cache = TRUE}
saveRDS(mamRERw, file="mamRERw.rds") 
newmamRERw = readRDS("mamRERw.rds")
```

---

# Computing Relative Evolutionary Rates (RERs)

**Figure 3. Schematic of Relative Evolutionary Rate calculation.**  
For each gene tree, branch lengths are regressed against the species tree. Residuals from this regression represent relative evolutionary rates, normalized across the phylogeny.

RERconverge calculates **branch-length residuals** for each gene relative to the species tree, producing normalized evolutionary rates.

```{r rers}
rerResults <- getAllResiduals(
  geneTrees,
  speciesTree,
  min.sp = 5
)
```

---

# Association with Enamel Loss

**Figure 4. Conceptual framework for RER–phenotype association.**  
Branches corresponding to enamel-less lineages are designated as foreground. Genes whose RERs are consistently elevated or reduced on these branches are inferred to be associated with enamel loss.

We now test whether evolutionary rates are **systematically correlated** with enamel loss across the phylogeny.

```{r correlation}
traitResults <- correlateWithBinaryPhenotype(
  rerResults,
  enamelTrait,
  method = "spearman"
)

# Multiple testing correction
traitResults$FDR <- p.adjust(traitResults$p.value, method = "fdr")
```

---

# Significant Genes

Genes with **FDR < 0.05** are considered significantly associated with enamel loss.

```{r significant-genes}
significantGenes <- subset(traitResults, FDR < 0.05)

significantGenes <- significantGenes[
  order(abs(significantGenes$Rho), decreasing = TRUE),
]

head(significantGenes)
```

---

# Interpretation

- **Positive correlation (Rho > 0)**  
  → Accelerated evolution in enamel-less mammals (relaxed constraint)

- **Negative correlation (Rho < 0)**  
  → Increased constraint or compensatory selection

Expected findings include:

- Strong signals in enamel matrix genes (*ENAM*, *AMELX*)
- Enrichment for tooth development and epithelial pathways
- Evidence of repeated pseudogenization events

---

# Visualization (Optional)

```{r visualization, eval=FALSE}
# Plot RERs for a candidate gene
plotRers(
  rerResults[["ENAM"]],
  speciesTree,
  enamelTrait,
  main = "Relative Evolutionary Rates of ENAM"
)
```

---

# Summary and Learning Objectives

By completing this handout, students should be able to:

- Explain how convergent evolution can be used to infer gene–trait associations
- Describe how RERconverge computes relative evolutionary rates
- Interpret correlations between evolutionary rates and binary traits
- Critically evaluate biological conclusions from comparative genomics analyses

---

# Conclusions

This analysis illustrates how **RERconverge** leverages convergent phenotypic evolution to identify genes underlying complex traits. Tooth enamel loss provides a clear, biologically interpretable example where evolutionary rate shifts reflect repeated relaxation of developmental constraints.

---

# References

- Kowalczyk et al. *Genome Biology* (2019)
- Chikina et al. *PLoS Genetics* (2016)
