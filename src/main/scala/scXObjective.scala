import scala.collection.mutable.{ListBuffer, Map => MutMap}
import scala.math.{ceil, pow}

object scXObjective {

    /**
      * Get the assignment of genes given the assignment of cells
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param samplesIDs for each sample, the list of cells it contains
      * @param cells an assignment of cells
      * @param expr map containing for each cell the expressed genes
      * @param geneSumPerSample the sum of positive expression of each gene for each sample
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return - the assignment of genes respecting the maximal percentage of negative values using a greedy approach
      *         (genes ranked per number of negative values, then per obj value)
      *         - the corresponding objective value
      */
    def getGenes(m: Array[Array[Double]], samplesIDs: Array[List[Int]], cells: List[Int], expr: Map[Int, List[Int]],
                 geneSumPerSample: Array[Array[Double]], nNeg: Double = 0.2, kappa: Double = 1): (List[Int], Double) = {
        val geneExpr = cells.flatMap(expr(_))

        val candidateGenes = geneExpr.groupBy(identity)
            .map(j => (j._1, cells.length - j._2.size, getGeneObj(m, samplesIDs, cells, j._1, geneSumPerSample, kappa)))
            .toList.filter(j => j._3 >= 0 && j._2 <= ceil(cells.length * nNeg))
            .sortBy(j => (j._2, -j._3))

        return selectCandGenes(m, candidateGenes, cells.size, nNeg)
    }

    /**
      * Gives the contribution of a gene to the objective under the given assignment of cells
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param samplesIDs for each sample, the list of cells it contains
      * @param cells an assignment of cells
      * @param gene a gene
      * @param geneSumPerSample the sum of positive expression of each gene for each sample
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return the objective value
      */
    def getGeneObj(m: Array[Array[Double]], samplesIDs: Array[List[Int]], cells: List[Int], gene: Int,
                   geneSumPerSample: Array[Array[Double]], kappa: Double = 1): Double = {
        val nbSam = samplesIDs.length.toDouble
        var samSum = 0.0
        var samProdPos = 1.0
        var samSumPos = 0.0

        for (s <- samplesIDs.indices) {
            val samObj = getGeneMccObj(m, cells intersect samplesIDs(s), gene, geneSumPerSample(s), kappa)
            samSum += samObj
            samProdPos *= samObj max 0.1
            samSumPos += samObj max 0.1
        }

        return (nbSam * pow(samProdPos, 1.0 / nbSam) / samSumPos) * samSum
    }

    /**
      * Gives the contribution of a gene to the objective under the given assignment of cells for a given sample
      * (i.e. MicroCellClust objective)
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param cells an assignment of cells
      * @param gene a gene
      * @param geneSum the cumulative positive expression of the genes
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return the objective value
      */
    def getGeneMccObj(m: Array[Array[Double]], cells: List[Int], gene: Int, geneSum: Array[Double], kappa: Double): Double = {
        var obj = -kappa * geneSum(gene)

        for (i <- cells) {
            if (m(i)(gene) >= 0) obj += (1 + kappa) * m(i)(gene)
            else obj += m(i)(gene)
        }

        return obj
    }

    /**
      * Get the assignment of markers given the assignment of samples with MicroCellClust objective
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param cells an assignment of cells
      * @param expr a map containing for each cell the expressed genes
      * @param geneSum the cumulative positive expression of the genes
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return - the assignment of genes respecting the maximal percentage of negative values using a greedy approach
      *         (genes ranked per number of negative values, then per obj value)
      *         - the corresponding objective value
      */
    def getGenesMcc(m: Array[Array[Double]], cells: List[Int], expr: Map[Int, List[Int]], geneSum: Array[Double],
                    nNeg: Double = 0.1, kappa: Double = 1): (List[Int], Double) = {
        val geneExpr = cells.flatMap(expr(_))

        val candidateGenes = geneExpr.groupBy(identity)
            .map(j => (j._1, cells.length - j._2.size, getGeneMccObj(m, cells, j._1, geneSum, kappa)))
            .toList.filter(j => j._3 >= 0 && j._2 <= ceil(cells.length * nNeg))
            .sortBy(j => (j._2, -j._3))

        return selectCandGenes(m, candidateGenes, cells.size, nNeg)
    }

    /**
      * Selects the assignment of markers from a list of candidates for a certain assignment of samples
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param candidateGenes a list of genes of format (ID_gene, nb_negatives, obj_val)
      * @param nbCells the number of cells in the assignment
      * @param nNeg the maximum percentage of negatives allowed inside the cluster
      * @return a list of genes and corresponding objective value
      */
    def selectCandGenes(m: Array[Array[Double]], candidateGenes: List[(Int, Int, Double)], nbCells: Int, nNeg: Double): (List[Int], Double) = {
        val genes = new ListBuffer[Int]
        var obj = 0.0
        var finished = false
        var n = 0.0
        var j = 0
        while (!finished && j < candidateGenes.length) {
            val g = candidateGenes(j)
            if ((n + g._2) / (nbCells * (genes.length + 1)) <= nNeg) {
                genes += g._1
                n += g._2
                obj += g._3
                j += 1
            } else finished = true
        }
        return (genes.toList, obj)
    }

    /**
      * For an assignment of cells, get the list of objective contributions for each gene
      * @param m an expression matrix
      * @param nbSam the number of samples
      * @param samIdCells an array with the sample number of each cell
      * @param cells an assignment of cells
      * @param geneSumPerSample the sum of positive expression of each gene for each sample
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return a list of genes that can still satisfy the constraint with one cell added, containing:
      *         - gene number
      *         - the objective contributions in each sample
      *         - the number of negatives in each sample
      *         - the number of cells in each sample
      */
    def getObjAllGenes(m: Array[Array[Double]], nbSam: Int, samIdCells: Array[Int], cells: List[Int], geneSumPerSample: Array[Array[Double]],
                       nNeg: Double = 0.2, kappa: Double = 1): List[(Int, Array[Double], Array[Int], Array[Int])] = {

        val expSums = ListBuffer[(Int, Array[Double], Array[Int], Array[Int])]()
        val maxNegNext = ceil((cells.length + 1) * nNeg)

        for (j <- m(0).indices) {

            val sp = Array.fill(nbSam)(0.0)
            val sn = Array.fill(nbSam)(0.0)
            val nn = Array.fill(nbSam)(0)
            val nbCells = Array.fill(nbSam)(0)
            for (i <- cells) {
                nbCells(samIdCells(i)) += 1
                if (m(i)(j) >= 0) {
                    sp(samIdCells(i)) += m(i)(j)
                } else {
                    sn(samIdCells(i)) += m(i)(j)
                    nn(samIdCells(i)) += 1
                }
            }

            if (nn.sum <= maxNegNext) {
                val o = (0 until nbSam).toArray.map(p => -kappa * geneSumPerSample(p)(j) + (1 + kappa) * sp(p) + sn(p))
                expSums += ((j, o, nn, nbCells))
            }

        }

        return expSums.toList
    }

    /**
      * Get the assignment of genes given the assignment of cells
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param nbSam the number of samples
      * @param nbCells the number of cells
      * @param newCell the cell to be added
      * @param sampleId the sample number of that cell
      * @param prevObjs the list of objective values of the previously selected cells
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return - the assignment of genes respecting the maximal percentage of negative values using a greedy approach
      *         (genes ranked per number of negative values, then per obj value)
      *         - the corresponding objective value
      */
    def getGenesFromPrev(m: Array[Array[Double]], nbSam: Int, nbCells: Int, newCell: Int, sampleId: Int,
                         prevObjs: List[(Int, Array[Double], Array[Int], Array[Int])], nNeg: Double = 0.2, kappa: Double = 1): (List[Int], Double) = {
        val candidateGenes = prevObjs.map(x => getGeneObjFromPrev(m, nbSam, newCell, sampleId, x._1, x._3.clone, x._2.clone, x._4.clone, kappa))
            .filter(j => j._3 >= 0 && j._2 <= ceil(nbCells * nNeg))
            .sortBy(j => (j._2, -j._3))

        return selectCandGenes(m, candidateGenes, nbCells, nNeg)
    }

    /**
      * Gives the contribution of a gene to the objective under the given assignment of cells
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param nbSam the number of samples
      * @param newCell the cell to be added
      * @param sampleId the sample number of that cell
      * @param gene a gene
      * @param negNb the number of in-cluster negative values of this gene in each sample
      * @param objs the objective contributions of this gene for each sample
      * @param nbCells the number of selected cells in each sample
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return the gene id, the number of negative values, the objective value
      */
    def getGeneObjFromPrev(m: Array[Array[Double]], nbSam: Int, newCell: Int, sampleId: Int, gene: Int, negNb: Array[Int],
                           objs: Array[Double], nbCells: Array[Int], kappa: Double = 1): (Int, Int, Double) = {
        if (m(newCell)(gene) >= 0) {
            objs(sampleId) += (1 + kappa) * m(newCell)(gene)
        } else {
            objs(sampleId) += m(newCell)(gene)
            negNb(sampleId) += 1
        }
        nbCells(sampleId) += 1

        val posObjs = objs.map(x => if (x > 0.1) x else 0.1)
        val obj = (nbSam * pow(posObjs.product, 1.0 / nbSam) / posObjs.sum) * objs.sum

        return (gene, negNb.sum, obj)
    }

    /**
      * Get gene assignment using MicroCellClust objective when adding a pair of cells
      * @param m an expression matrix
      * @param nbCells the number of cells in the assignment
      * @param newPair the cells added to the previous assignment
      * @param prevObjs a list of genes with their objective values and nb. of in-cluster negatives before adding the pair
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return - the assignment of genes respecting the maximal percentage of negative values using a greedy approach
      *         (genes ranked per number of negative values, then per obj value)
      *         - the corresponding objective value
      */
    def getGenesMccFromPrevPair(m: Array[Array[Double]], nbCells: Int, newPair: List[Int], prevObjs: List[(Int, Double, Int)],
                                nNeg: Double = 0.2, kappa: Double = 1): (List[Int], Double) = {
        val maxNeg = ceil(nbCells * nNeg)

        val candidateMarks = prevObjs.map { case (j, obj, nn) => getGeneMccObjFromPrevPair(m, newPair, j, obj, nn, kappa) }
            .filter(j => j._3 >= 0 && j._2 <= maxNeg)
            .sortBy(j => (j._2, -j._3))

        return selectCandGenes(m, candidateMarks, nbCells, nNeg)
    }

    /**
      * Gives the contribution of a gene to the objective under the given assignment of cells
      * @param m an expression matrix
      * @param newPair a pair of cells
      * @param gene a gene
      * @param prevObj the objective value before adding the pair
      * @param prevNegNb the number of in-cluster negatives before adding the pair
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return the gene number, the number of in-cluster negatives, objective value
      */
    def getGeneMccObjFromPrevPair(m: Array[Array[Double]], newPair: List[Int], gene: Int,
                                  prevObj: Double, prevNegNb: Int, kappa: Double): (Int, Int, Double) = {
        var obj = prevObj
        var nn = prevNegNb

        for (newCell <- newPair) {
            if (m(newCell)(gene) >= 0) {
                obj += (1 + kappa) * m(newCell)(gene)
            } else {
                obj += m(newCell)(gene)
                nn += 1
            }
        }
        return (gene, nn, obj)
    }

    /**
      * Computes the sum of expression over the given genes for each cell
      * @param m an expression matrix
      * @param genes a list of genes
      * @return an array with the sum of maker expression for each sample
      */
    def getCellsObj(m: Array[Array[Double]], genes: List[Int]): Array[Double] = {
        return m.indices.toArray.map(i => genes.map(j => m(i)(j)).sum)
    }

    /**
      * Get the assignment of genes given the assignment of cells
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param cells an assignment of cells
      * @param objsNew for each gene, the objective values for each sample
      * @param nnNew for each gene, the number of in-cluster negative values
      * @param samplesIDs for each sample, the list of cells it contains
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @return - the assignment of genes respecting the maximal percentage of negative values using a greedy approach
      *         (genes ranked per number of negative values, then per obj value)
      *         - the corresponding objective value
      */
    def getGenesLS(m: Array[Array[Double]], cells: List[Int], objsNew: Array[Array[Double]], nnNew: Array[Int],
                   samplesIDs: Array[List[Int]], nNeg: Double): (List[Int], Double) ={

        val candGenes = nnNew.zipWithIndex.toList.filter { case (n, j) => n <= ceil(cells.size * nNeg) }
            .map { case (n, j) => (j, n, getGeneObjLS(samplesIDs, objsNew(j))) }
            .filter(j => j._3 >= 0).sortBy(j => (j._2, -j._3))

        return selectCandGenes(m, candGenes, cells.size, nNeg)
    }

    /**
      * Gives the contribution of a gene to the objective under the given assignment of cells
      * @param samplesIDs for each sample, the list of cells it contains
      * @param geneObjs the objective values of the gene for each sample
      * @return the objective value
      */
    def getGeneObjLS(samplesIDs: Array[List[Int]], geneObjs: Array[Double]): Double = {
        val nbSam = samplesIDs.length.toDouble
        var samSum = 0.0
        var samProdPos = 1.0
        var samSumPos = 0.0

        for (p <- samplesIDs.indices) {
            samSum += geneObjs(p)
            samProdPos *= geneObjs(p) max 0.1
            samSumPos += geneObjs(p) max 0.1
        }

        return (nbSam * Math.pow(samProdPos, 1.0 / nbSam) / samSumPos) * samSum
    }

    /**
      * Get the cumulative positive expression of each gene
      * @param m an expression matrix
      * @return an array with the sum of positive expressions of each gene
      */
    def getGeneSum(m: Array[Array[Double]]): Array[Double] = {
        val geneSum = Array.fill(m(0).length)(0.0)

        for (j <- m(0).indices) {
            for (i <- m.indices) {
                if (m(i)(j) >= 0) geneSum(j) += m(i)(j)
            }
        }

        return geneSum
    }

    /**
      * Get the cumulative positive expression of each gene for each sample
      * @param m an expression matrix
      * @param samplesIDs for each sample, the list of cells it contains
      * @return an array with, for each sample, an array with the sum of positive expressions of each gene
      */
    def getGeneSumPerSample(m: Array[Array[Double]], samplesIDs: Array[List[Int]]): Array[Array[Double]] = {
        return samplesIDs.map(s => getGeneSum(s.toArray.map(m(_))))
    }

    /**
      * Builds a mapping containing, for each cell, the genes it expresses
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @return a cell -> list of genes map
      */
    def buildExprMap(m: Array[Array[Double]]): Map[Int, List[Int]] = {
        val expr = MutMap[Int, List[Int]]()

        for (i <- m.indices) {
            expr += (i -> m(i).zipWithIndex.filter(_._1 >= 0).map(_._2).toList)
        }

        return expr.toMap
    }

    /**
      * Computes the percentage of the given genes expressed by the cell
      * @param cell a cell
      * @param genes a list of genes
      * @param expr a cell -> gene map
      * @return the percentage
      */
    def getSupport(cell: Int, genes: List[Int], expr: Map[Int, List[Int]]): Double = {
        val inter = genes intersect expr(cell)
        return inter.length.toDouble / genes.length
    }
}
