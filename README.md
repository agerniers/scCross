# scCross

R and Scala implementation of scCross [1].

### References

[1] ...


## RScala

To run scCross in R, the user must install the `rscala` package (https://cran.r-project.org/web/packages/rscala/index.html), that can call the Scala implementation of scCross from within R. Once installed, he can check wether there is a version of Scala on his computer using:

``` R
rscala::scalaConfig()
```
This function tries to find Scala and Java on the user's computer and, if needed, downloads and installs Scala and Java in the user's `~/.rscala` directory.


### Initializing scCross

To run scCross, the user must first set his working directory (`setwd`) to the location of the `scCross` directory, and then load the functions contained in the `sc_cross.R` script.
``` R
source("./src/main/R/sc_cross.R")
```
Then, he must initialise an instance of RScala containing the scCross source code.
``` R
scX.rs = initRscala()
```
This `mcc.rs` object needs to be passed as the first argument to any R function using Scala code.


## Preparing the data

scCross assumes the data is represented in a matrix (or data frame) containing positive values for presence of expression, and negative values in case of absence of (or negligible) expression. Such data can be obtained from (normalize) count data using an appropriate transformation, such as $\log_{10}(x + 0.1)$.

**Note:** Applying such a transformation means the data can no longer be stored in a sparse structure, which can be problematic when dealing with large structures of data. Therefore, all given `R` functions *will assume by default the given data are count values (or normalized counts)* and will apply the transformation internally. By default, they apply a $\log_{10}(x + 0.1)$ transformation. This behavior can be changed by setting the `data.tfo` parameter accordingly :

* `"log10"` : $\log_{10}(x + 0.1)$ transformation (default)
* `"log2"` : $\log_{2}(x + 0.5)$ transformation
* a custom `R` function with one parameter (corresponding to a single numerical value) that transforms the value passed as argument
* `"none"` if the data is already in the appropriate format and no transformation should be applied

Let's assume the (normalized) count data is contained in a matrix/dataframe called `mat` (assuming the rows represent genes and the columns cells). scCross requires genes expressed in too many cells of the dataset to be removed.
``` R
gene.expr = rowSums(mat > 0) # Can be useful to plot a histogram to define a threshold. Typically keeping genes expressed in less that 25% keeps the vast majority of the genes
mat.filt = mat[gene.expr >= 2 & gene.expr <= 0.25 * ncol(mat), ]
```

**Note:** Instead of removing these genes, one could instead opt to remove the median value (to search for overexpression) or invert the values (to look for absence of expression within some cells).

MicroCellClust needs, for each sample, a vector containing the sum of positive expressions of each gene. This can be computed internally, but as it might be time-consuming for large instances, we recommend to pre-compute it once and store it in memory. Assuming `samples` contains the sample annotation of each cell:
```R
gs = geneSumPerSample(mat.filt, samples) # If needed, don't forget to set the `data.tfo` argument
```
**Remark:** By default, the R implementation represents the cells by the columns, and the genes by the rows. In the inverse case, the user can indicate this by giving the parameter `cellsOnCol = FALSE` to the function. *Note that in the Scala implementation, the samples (here cells) are represented by the rows, which is traditionally the case in machine learning.*


### Rareness score (optional)
To operate on large instances, scCross uses a rareness score as heuristics. One should compute a vector `rs` that contains a value for each cell. Existing methods include :
* FiRE [Jindal *et al.*, 2018] https://github.com/princethewinner/FiRE
* DoRC [Chen *et al.*, 2019] https://github.com/chenxofhit/DoRC


## Running scCross

scCross is run using the following function (the `gene.sum` and `rareness.score` parameters are optinal):
``` R
scX.res = runScCross(scX.rs, mat.filt, samples, gene.sum = gs, rareness.score = rs) # If needed, don't forget to set the `data.tfo` and `cellesOnCol` arguments
```
As this function calls the Scala solver, the user must provide the `scX.rs ` object as first argument. `runScCross` returns a list containing the indices and names of the cells (`scX.res$cells.idx` and `scX.res$cells.names`) and genes (`scX.res$genes.idx` and `scX.res$genes.names`) composing the identified bicluster, as well as the latter's objective value (`scX.res$obj.value`). `scX.res$info` contains additional information, including intermediate results.

### Optimization problem parameters

Arguments `kappa` and `nNeg` correspond to the $\kappa \ge 0$ and $\mu \in [0, 1]$ metaparameters of the optimisation problem [1], controlling respectively the out-of-cluster expression and the maximum proportion of negative values inside the bicluster. When omitted, default values of `kappa = 100 / ncol(mat)` and `nNeg = 0.1` are taken. 

To tune these values, one could re-run several times the `runScCross ` function until obtaining the desired result. However, as it operates in two steps (beam search + local search), one can perform this tuning by only re-running the second step (the local search) using:
``` R
runScCross.quick(scX.rs, mat.filt, samples, scX.res$info$res1.cells.idx, scX.res$info$res1.genes.idx, gene.sum = gs) # scX.res$info$res1.___ contains the result of the first step (beam search), which will be used as initial solution
```

### Looking for other solutions

Once a first cell/gene bicluster has been found, one can re-run `runScCross` while excluding the cells in the `scX.rs` solution to search for another bicluster. This can be done using:
``` R
scX.res.2 = runScCross(scX.rs, mat.filt, samples, gene.sum = gs, rareness.score = rs, cells.excl = scX.res$cells.idx) # If needed, don't forget to set the `data.tfo` and `cellesOnCol` arguments
```
Alternatively, one could also/instead exclude the genes from the first solution using `genes.excl = scX.res$genes.idx`. Excluding previously identified genes rather than cells could be useful to identify hierachical structures, as cells could be classiffied multiple times.


### Other parameters

* By default, `runScCross ` uses the FiRE threshold $Q3 + 1.5 \cdot IQR$ on the rareness scores to select the cells for the beam search. An other value than $1.5$ can be specified using the `iqr` parameter. 
* By default, `runScCross ` uses a minimum of 1000 cells (assuming the data is larger than that) for the beam search. This value can be changed using the `min.cells.beam` parameter.
* One can set the `seed` parameter for the random events in the local search.

