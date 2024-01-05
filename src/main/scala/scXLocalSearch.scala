import scXObjective._

import scala.collection.immutable.{Map, Set}
import scala.math.{ceil, exp, floor, pow}
import scala.util.Random

object scXLocalSearch {

    /**
      * Local search to improve initial solution
      * @param m an expression matrix with samples cells on the rows and markers genes on the columns
      * @param samIdCells an array with the sample number of each cell
      * @param initCells an initial assignment of samples
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param nbIter the number of neighbours to visit
      * @param preGeneSumPerSample the sum of positive values for each marker (required if m contains only part of the samples)
      * @param seed the seed of the pseudo-random number generator (0 == no seed)
      * @param verbose enable/disable printing
      * @return - an assignment of cells
      *         - an assignment of genes
      *         - the corresponding objective value
      */
    def localSearch(m: Array[Array[Double]], samIdCells: Array[Int], initCells: Set[Int], nNeg: Double = 0.2, kappa: Double = 1,
                    nbIter: Int = 1000, nbRestart: Int = 20, preGeneSumPerSample: Array[Array[Double]] = Array(),
                    seed: Int = 0, verbose: Boolean = true): (List[Int], List[Int], Double) = {
        val expr = buildExprMap(m)
        val samplesIDs = (0 to samIdCells.max).toArray.map(i => samIdCells.zipWithIndex.filter(_._1 == i).map(_._2).toList)
        val geneSumPerSample = if (preGeneSumPerSample.length == samplesIDs.length && preGeneSumPerSample(0).length == m(0).length) preGeneSumPerSample
                               else getGeneSumPerSample(m, samplesIDs)

        val (initGenes, initObj) = getGenes(m, samplesIDs, initCells.toList, expr, geneSumPerSample, nNeg, kappa)
        var best = (initCells, initGenes, initObj)

        if (verbose) {
            val i = ""
            val nb = "New best?"
            val c = "Nb. cells"
            val g = "Nb. genes"
            val o = "Obj. value"
            val t = "Exe. time [s]"
            println(f"| $i%-25s | $nb%-12s | $c%-12s | $g%-12s | $o%-12s | $t%-15s |")
            val i2 = "Initial solution"
            val nb2 = "    ****"
            val c2 = initCells.size
            val g2 = initGenes.length
            val o2 = (initObj * 100).round / 100.0
            println(f"| $i2%-25s | $nb2%-12s | $c2%-12s | $g2%-12s | $o2%-12s | $i%-15s |")
        }

        for (nbr <- 1 to nbRestart) {
            if (verbose) {
                val i = "Local search restart " + nbr
                print(f"| $i%-25s |")
            }
            var bestRestart = (Set[Int](), List[Int](), 0.0)
            val t0 = System.currentTimeMillis

            var objs = Array.fill(m(0).length)(Array.fill(samplesIDs.length)(0.0))
            var nn = Array.fill(m(0).length)(0)
            for (j <- m(0).indices) {
                objs(j) = geneSumPerSample.map(-kappa * _ (j))
                for (i <- initCells) {
                    if (m(i)(j) >= 0) {
                        objs(j)(samIdCells(i)) += (1 + kappa) * m(i)(j)
                    } else {
                        objs(j)(samIdCells(i)) += m(i)(j)
                        nn(j) += 1
                    }
                }
            }
            var current = (initCells, initGenes, initObj)

            val heurOrd = getCellsObj(m, best._2).zipWithIndex.sortBy(-_._1).map(_._2)

            val rand = if (seed == 0) new Random() else new Random(seed * nbr)
            var temp = temperature(initObj, 0)

            for (k <- 0 until nbIter) {

                val (curOrdered, oocOrdered) = heurOrd.partition(current._1 contains _)

                val (cells, objsNew, nnNew) = neighbour(m, samIdCells, curOrdered, oocOrdered, objs, nn, kappa, rand)

                val (genes, obj) = getGenesLS(m, cells.toList, objsNew, nnNew, samplesIDs, nNeg)

                if (obj > bestRestart._3) {
                    bestRestart = (cells, genes, obj)

                    if (obj > best._3) {
                        best = (cells, genes, obj)
                    }
                }

                if (acceptance(current._3, obj, temp) >= rand.nextDouble()) {
                    current = (cells, genes, obj)
                    objs = objsNew
                    nn = nnNew
                }

                temp = temperature(initObj, k)
            }

            val (cells, genes, obj) = greedyImprove(m, samplesIDs, samIdCells, bestRestart, expr, geneSumPerSample, nNeg, kappa)

            if (obj > bestRestart._3) {
                bestRestart = (cells, genes, obj)

                if (obj > best._3) {
                    best = (cells, genes, obj)
                }
            }

            val t1 = System.currentTimeMillis

            if (verbose) {
                val nb = if (bestRestart._3 >= best._3) "    ****" else ""
                val c = bestRestart._1.size
                val g = bestRestart._2.length
                val o = (bestRestart._3 * 100).round / 100.0
                val t = (t1 - t0).toDouble / 1000
                println(f" $nb%-12s | $c%-12s | $g%-12s | $o%-12s | $t%-15s |")
            }
        }

        return (best._1.toList.sorted, best._2, best._3)
    }

    /**
      * Stochastically generates a new neighbour
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param samIdCells an array with the sample number of each cell
      * @param curOrdered the cells currently selected, ordered decreasingly by sum of expression over initial genes
      * @param oocOrdered the cells currently not selected, ordered decreasingly by sum of expression over initial genes
      * @param objs for all genes, the objective values in each sample
      * @param nn number of negative values in the current cell selection for each gene
      * @param rand a pseudo-random number generator
      * @return a new assignment of cells, with updated àbjs and nn
      */
    def neighbour(m: Array[Array[Double]], samIdCells: Array[Int], curOrdered: Array[Int], oocOrdered: Array[Int],
                  objs: Array[Array[Double]], nn: Array[Int], kappa: Double, rand: Random):
                 (Set[Int], Array[Array[Double]], Array[Int]) = {

        if (rand.nextInt(m.length) < oocOrdered.length) { // Add sample to current
            val idx = floor(oocOrdered.length * pow(rand.nextDouble(), 2)).toInt
            val cell = oocOrdered(idx)

            val objsNew = objs.indices.toArray.map(j => {
                val a = if (m(cell)(j) >= 0) (1 + kappa) * m(cell)(j) else m(cell)(j)
                val b = objs(j).clone()
                b(samIdCells(cell)) += a
                b
            })
            val nnNew = nn.indices.toArray.map(j => {
                val a = if (m(cell)(j) < 0) 1 else 0
                nn(j) + a
            })

            return (curOrdered.toSet + cell, objsNew, nnNew)

        } else { // Remove sample from current
            val idx = (curOrdered.length - 1) - floor(curOrdered.length * pow(rand.nextDouble(), 2)).toInt
            val cell = curOrdered(idx)

            val objsNew = objs.indices.toArray.map(j => {
                val a = if (m(cell)(j) >= 0) (1 + kappa) * m(cell)(j) else m(cell)(j)
                val b = objs(j).clone()
                b(samIdCells(cell)) -= a
                b
            })
            val nnNew = nn.indices.toArray.map(j => {
                val a = if (m(cell)(j) < 0) 1 else 0
                nn(j) - a
            })

            return (curOrdered.toSet - cell, objsNew, nnNew)
        }
    }

    /**
      * Temperature for the simulated annealing
      * @param t0 the initial temperature
      * @param k the iteration number
      * @return the current temperature
      */
    def temperature(t0: Double, k: Int): Double = {
        t0 * pow(0.99, k)
    }

    /**
      * Computes the probability to accept the neighbour
      * @param obj the current objective value
      * @param objNew the objective value of the neighbour
      * @param temp the current temperature
      * @return a probability
      */
    def acceptance(obj: Double, objNew: Double, temp: Double): Double = {
        if (objNew >= obj) {
            return 1.0
        } else {
            return exp(-(obj - objNew) / temp) max 0.01
        }
    }

    /**
      * Greedily searches if a better solution can be found by adding out-of-cluster cells, in order of their sum
      * of expression over the selected genes
      * @param m an expression matrix with samples (cells) on the rows and markers (genes) on the columns
      * @param samplesIDs for each sample, the list of cells it contains
      * @param samIdCells an array with the sample number of each cell
      * @param init the current solution
      * @param expr map containing for each cell the expressed genes
      * @param geneSumPerSample the sum of positive expression of each gene for each sample
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @return a (possibly improved) solution
      */
    def greedyImprove(m: Array[Array[Double]], samplesIDs: Array[List[Int]], samIdCells: Array[Int], init: (Set[Int], List[Int], Double),
                      expr: Map[Int, List[Int]], geneSumPerSample: Array[Array[Double]], nNeg: Double, kappa: Double): (Set[Int], List[Int], Double) = {
        var current = init._1
        var best = init

        val cellsObj = getCellsObj(m, init._2).zipWithIndex
        val (inclSamUnsorted, oocSamUnsorted) = cellsObj.partition(init._1 contains _._2)
        val inclSam = inclSamUnsorted.sortBy(_._1)
        val oocSam = oocSamUnsorted.sortBy(-_._1)

        var objs = Array.fill(m(0).length)(Array.fill(samplesIDs.length)(0.0))
        var nn = Array.fill(m(0).length)(0)
        for (j <- m(0).indices) {
            objs(j) = geneSumPerSample.map(-kappa * _ (j))
            for (i <- init._1) {
                if (m(i)(j) >= 0) {
                    objs(j)(samIdCells(i)) += (1 + kappa) * m(i)(j)
                } else {
                    objs(j)(samIdCells(i)) += m(i)(j)
                    nn(j) += 1
                }
            }
        }

        for (i <- 0 until (oocSam.length min inclSam.length)) {
            if (oocSam(i)._1 > inclSam(i)._1) {
                current = current + oocSam(i)._2 - inclSam(i)._2

                objs = objs.indices.toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) >= 0) (1 + kappa) * m(oocSam(i)._2)(j) else m(oocSam(i)._2)(j)
                    val r = if (m(inclSam(i)._2)(j) >= 0) (1 + kappa) * m(inclSam(i)._2)(j) else m(inclSam(i)._2)(j)
                    objs(j)(samIdCells(oocSam(i)._2)) += a
                    objs(j)(samIdCells(inclSam(i)._2)) -= r
                    objs(j)
                })
                nn = nn.indices.toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) < 0) 1 else 0
                    val r = if (m(inclSam(i)._2)(j) < 0) 1 else 0
                    nn(j) + a - r
                })
            } else {
                current = current + oocSam(i)._2

                objs = objs.indices.toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) >= 0) (1 + kappa) * m(oocSam(i)._2)(j) else m(oocSam(i)._2)(j)
                    objs(j)(samIdCells(oocSam(i)._2)) += a
                    objs(j)
                })
                nn = nn.indices.toArray.map(j => {
                    val a = if (m(oocSam(i)._2)(j) < 0) 1 else 0
                    nn(j) + a
                })
            }

            val (genes, obj) = getGenesLS(m, current.toList, objs, nn, samplesIDs, nNeg)

            if (obj > best._3) {
                best = (current, genes, obj)
            }
        }

        return best
    }

}
