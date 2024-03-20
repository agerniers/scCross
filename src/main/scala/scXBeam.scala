import scXObjective._

import scala.collection.mutable.{PriorityQueue => MutPriorityQueue, Set => MutSet}

object scXBeam {

    /**
      * Find the bicluster maximising the objective using a beam search
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param samIdCells  an array with the sample number of each cell
      * @param preGeneSumPerSample the sum of positive values for each marker (required if m contains only part of the samples)
      * @param rarenessScore a list of rareness scores for each cell
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param kAdapt whether to adapt kappa if too many/few genes are selected
      * @param nHeurPair the number of pairs to form for each cell at level 2
      * @param nHeurKeep the number of top-solutions to consider for expansion at the next level
      * @param nhAdapt whether to decrease nHeurKeep during the search
      * @param maxNbCells the maximum number of cells in the bicluster
      * @param stopNoImprove stop search when no better solution is found after x levels
      * @param verbose enable/disable printing
      * @return - an assignment of cells
      *         - an assignment of genes
      *         - the corresponding objective value
      */
    def findCluster(m: Array[Array[Double]], samIdCells: Array[Int], preGeneSumPerSample: Array[Array[Double]] = Array(),
                    rarenessScore: Array[Double] = Array(), nNeg: Double = 0.2, kappa: Double = 1, kAdapt: Boolean = true,
                    nHeurPair: Int = 100, nHeurKeep: Int = 100, nhAdapt: Boolean = true, maxNbCells: Int = Int.MaxValue,
                    stopNoImprove: Int = 25, verbose: Boolean = true): (List[Int], List[Int], Double, Double) = {
        val expr = buildExprMap(m)
        val samplesIDs = (0 to samIdCells.max).toArray.map(i => samIdCells.zipWithIndex.filter(_._1 == i).map(_._2).toList)
        val geneSumPerSample = if (preGeneSumPerSample.length == samplesIDs.length && preGeneSumPerSample(0).length == m(0).length) preGeneSumPerSample else getGeneSumPerSample(m, samplesIDs)
        val nbCells = m.length
        val hOrd = if (rarenessScore.length != nbCells) Array[Array[(Double, Int)]]()
                   else samplesIDs.map(l => rarenessScore.zipWithIndex.filter(x => l contains x._2).sortBy(-_._1))
        val nbSam = samplesIDs.length

        var nH = nHeurKeep
        val geneSum = geneSumPerSample.transpose.map(_.sum)

        val t0 = System.currentTimeMillis
        var prevLvlNBest = start(m, samplesIDs, samIdCells, geneSum, hOrd, expr, nNeg, kappa, nHeurPair, nHeurKeep, verbose)
        val t1 = System.currentTimeMillis
        val (genes, obj) = getGenes(m, samplesIDs, prevLvlNBest.head.toList, expr, geneSumPerSample, nNeg, kappa)
        var best = (prevLvlNBest.head, genes, obj)

        if (verbose) {
            val l = "Level"
            val k = "(kappa)"
            val n = "(nb. heur.)"
            val p = "(nb. cand.)"
            val nb = "New best?"
            val g = "Nb. genes"
            val o = "Obj. value"
            val m = "(mean sim.)"
            val t = "Exe. time [s]"
            println(f"| $l%-6s | $k%-8s | $n%-12s | $p%-12s | $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |")
        }

        if (verbose) {
            val lvl = "start"
            val nb = "    ****"
            val g = best._2.length
            val o = (best._3 * 100).round / 100.0
            val m = ""
            val t = (t1 - t0).toDouble / 1000
            println(f"| $lvl%-6s | $kappa%-8s | $nH%-12s | $nbCells%-12s | $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |") // + best._1.toList.sorted.mkString(" "))
        }

        var kA = kappa

        var lvl = best._1.size + 1
        val startLvl = lvl
        val maxLvl = maxNbCells min nbCells
        var noImprove = 0
        var finished = false
        var cellFilt = false
        var cellPool = (0 until nbCells).toSet

        //
        // Use queue for each patient
        //

        while (!finished && lvl <= maxLvl.min(startLvl + 10)) {
            if (verbose) {
                val npl = prevLvlNBest.size
                val p = cellPool.size
                print(f"| $lvl%-6s | $kA%-8s | $npl%-12s | $p%-12s |")
            }
            val t0 = System.currentTimeMillis
            val nBestQueues = Array.fill(nbSam)(MutPriorityQueue[(Set[Int], List[Int], Double)]()(Ordering[Double].on(x => -x._3)))

            val hist = MutSet[Set[Int]]()

            for (prev <- prevLvlNBest) {
                val prevObjs = getObjAllGenes(m, nbSam, samIdCells, prev.toList, geneSumPerSample, nNeg, kA)
                for (i <- cellPool if !prev.contains(i)) {
                    val cells = prev + i
                    val l = hist.size
                    hist += cells
                    if (hist.size - l >= 1) {
                        val (genes, obj) = getGenesFromPrev(m, nbSam, cells.size, i, samIdCells(i), prevObjs, nNeg = nNeg, kappa = kA)
                        for (sam <- cells.map(samIdCells(_))) {
                            nBestQueues(sam) += ((cells, genes, obj))
                            if (nBestQueues(sam).size > nHeurKeep) nBestQueues(sam).dequeue()
                        }
                    }
                }
            }

            val t1 = System.currentTimeMillis

            if (nBestQueues.exists(_.nonEmpty)) {
                val nBest = samplesIDs.indices.flatMap(sam => nBestQueues(sam).dequeueAll).toList.distinct.sortBy(-_._3)
                val cBest = nBest.head

                val msim = if (nH <= 1) 1
                else nBest.drop(1).map(x => (x._2 intersect cBest._2).size.toDouble / cBest._2.size).sum / (nBest.length - 1)

                // If similarity of genes between solution is high: decrease nHeurKeep
                if (nhAdapt) {
                    if (msim >= 0.98 && lvl >= 200) nH = 1
                    else if (msim >= 0.98) nH = (nHeurKeep / 8) max 1
                    else if (msim >= 0.9) nH = nHeurKeep / 4
                    else if (msim >= 0.7) nH = nHeurKeep / 2
                    else nH = nHeurKeep
                }

                // If similarity of genes between solution is high: keep only samples that express these markers
                val t2 = System.currentTimeMillis()

                if (!cellFilt && msim >= 0.9) {
                    cellFilt = true
                    val geneUnion = nBest.flatMap(_._2).distinct
                    cellPool = (0 until nbCells).toSet.filter(i => scXObjective.getSupport(i, geneUnion, expr) >= 0.3)
                } else if (cellFilt && lvl % 10 == 0) {
                    val geneUnion = nBest.flatMap(_._2).distinct
                    cellPool = (0 until nbCells).toSet.filter(i => scXObjective.getSupport(i, geneUnion, expr) >= 0.3)
                }
                val t3 = System.currentTimeMillis()

                if (cBest._3 > best._3) {
                    best = (cBest._1, cBest._2, cBest._3)
                    noImprove = 0
                    if (verbose) {
                        val nb = "    ****"
                        val g = cBest._2.length
                        val o = (cBest._3 * 100).round / 100.0
                        val m = (msim * 10000).round / 100.0
                        val t = (t1 - t0 /* + t3 - t2 */).toDouble / 1000
                        println(f" $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |")
                    }

                    if (kAdapt) {
                        if (cBest._2.length <= 200 && kA <= 1000) {
                            if (kA < kappa) { // If kappa was temporarily decreased, check if we can go back to normal
                                val ng = getGenes(m, samplesIDs, cBest._1.toList, expr, geneSumPerSample, nNeg, kappa)._1.length

                                if (ng >= 20) {
                                    kA = kappa
                                    best = (Set(), List(), 0.0)
                                } else {
                                    prevLvlNBest = nBest.map(x => x._1)
                                    lvl += 1
                                }
                            } else {
                                prevLvlNBest = nBest.map(x => x._1)
                                lvl += 1
                            }
                        } else { // If the solution contains too many markers, increase kappa
                            var nTry = 5
                            var kTest = kA * 2
                            while (nTry > 0) {
                                val ng = getGenes(m, samplesIDs, cBest._1.toList, expr, geneSumPerSample, nNeg, kappa)._1.length
                                if (ng > 100) {
                                    nTry -= 1
                                    kTest = kTest * 2
                                } else if (ng < 50) {
                                    nTry -= 1
                                    kTest = kTest * 0.75
                                } else {
                                    nTry = 0
                                }
                            }
                            kA = kTest
                            best = (Set(), List(), 0.0)
                        }
                    } else {
                        prevLvlNBest = nBest.map(x => x._1)
                        lvl += 1
                    }
                } else {
                    noImprove += 1
                    if (verbose) {
                        val nb = ""
                        val g = cBest._2.length
                        val o = (cBest._3 * 100).round / 100.0
                        val m = (msim * 10000).round / 100.0
                        val t = (t1 - t0 + t3 - t2).toDouble / 1000
                        println(f" $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |")
                    }
                    if (noImprove >= stopNoImprove) {
                        finished = true
                        println()
                    } else if (cBest._3 <= 0.0) {
                        finished = true
                        println()
                    }
                    prevLvlNBest = nBest.map(x =>  x._1)
                    lvl += 1
                }
            } else {
                finished = true
                println()
            }
        }

        //
        // Use one queue
        //

        while (!finished && lvl <= maxLvl) {
            if (verbose) {
                val npl = prevLvlNBest.size
                val p = cellPool.size
                print(f"| $lvl%-6s | $kA%-8s | $npl%-12s | $p%-12s |")
            }
            val t0 = System.currentTimeMillis
            val nBestQueue = MutPriorityQueue[(Set[Int], List[Int], Double)]()(Ordering[Double].on(x => -x._3))

            val hist = MutSet[Set[Int]]()

            for (prev <- prevLvlNBest) {
                val prevObjs = getObjAllGenes(m, nbSam, samIdCells, prev.toList, geneSumPerSample, nNeg, kA)
                for (i <- cellPool if !prev.contains(i)) {
                    val cells = prev + i
                    val l = hist.size
                    hist += cells
                    if (hist.size - l >= 1) {
                        val (genes, obj) = getGenesFromPrev(m, nbSam, cells.size, i, samIdCells(i), prevObjs, nNeg = nNeg, kappa = kA)
                        nBestQueue += ((cells, genes, obj))
                        if (nBestQueue.size > nH) nBestQueue.dequeue()
                    }
                }
            }

            val t1 = System.currentTimeMillis

            if (nBestQueue.nonEmpty) {
                val nBest = nBestQueue.clone.dequeueAll.toList.reverse
                val cBest = nBest.head

                val msim = if (nH <= 1) 1
                else nBest.drop(1).map(x => (x._2 intersect cBest._2).size.toDouble / cBest._2.size).sum / (nBest.length - 1)

                // If similarity of markers between solution is high: decrease nHeurKeep
                if (nhAdapt) {
                    if (msim >= 0.98 && lvl >= 200) nH = 1
                    else if (msim >= 0.98) nH = (nHeurKeep / 8) max 1
                    else if (msim >= 0.9) nH = nHeurKeep / 4
                    else if (msim >= 0.7) nH = nHeurKeep / 2
                    else nH = nHeurKeep
                }

                // If similarity of markers between solution is high: keep only samples that express these markers
                val t2 = System.currentTimeMillis()

                if (!cellFilt && msim >= 0.9) {
                    cellFilt = true
                    val markUnion = nBest.flatMap(_._2).distinct
                    cellPool = (0 until nbCells).toSet.filter(i => scXObjective.getSupport(i, markUnion, expr) >= 0.3)
                } else if (cellFilt && lvl % 10 == 0) {
                    val markUnion = nBest.flatMap(_._2).distinct
                    cellPool = (0 until nbCells).toSet.filter(i => scXObjective.getSupport(i, markUnion, expr) >= 0.3)
                }
                val t3 = System.currentTimeMillis()

                if (cBest._3 > best._3) {
                    best = (cBest._1, cBest._2, cBest._3)
                    noImprove = 0
                    if (verbose) {
                        val nb = "    ****"
                        val g = cBest._2.length
                        val o = (cBest._3 * 100).round / 100.0
                        val m = (msim * 10000).round / 100.0
                        val t = (t1 -t0 /* + t3 - t2 */).toDouble / 1000
                        println(f" $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |") // + best._1.toList.sorted.mkString(" "))
                    }

                    if (kAdapt) {
                        if (cBest._2.length <= 200 && kA <= 1000) {
                            if (kA < kappa) { // If kappa was temporarily decreased, check if we can go back to normal
                                val ng = getGenes(m, samplesIDs, cBest._1.toList, expr, geneSumPerSample, nNeg, kappa)._1.length

                                if (ng >= 20) {
                                    kA = kappa
                                    best = (Set(), List(), 0.0)
                                } else {
                                    prevLvlNBest = nBest.map(x => x._1)
                                    lvl += 1
                                }
                            } else {
                                prevLvlNBest = nBest.map(x => x._1)
                                lvl += 1
                            }
                        } else { // If the solution contains too many markers, increase kappa
                            var nTry = 5
                            var kTest = kA * 2
                            while (nTry > 0) {
                                val ng = getGenes(m, samplesIDs, cBest._1.toList, expr, geneSumPerSample, nNeg, kappa)._1.length
                                if (ng > 100) {
                                    nTry -= 1
                                    kTest = kTest * 2
                                } else if (ng < 50) {
                                    nTry -= 1
                                    kTest = kTest * 0.75
                                } else {
                                    nTry = 0
                                }
                            }
                            kA = kTest
                            best = (Set(), List(), 0.0)
                        }
                    } else {
                        prevLvlNBest = nBest.map(x => x._1)
                        lvl += 1
                    }
                } else {
                    noImprove += 1
                    if (verbose) {
                        val nb = ""
                        val g = cBest._2.length
                        val o = (cBest._3 * 100).round / 100.0
                        val m = (msim * 10000).round / 100.0
                        val t = (t1 - t0 + t3 - t2).toDouble / 1000
                        println(f" $nb%-12s | $g%-12s | $o%-12s | $m%-12s | $t%-15s |") // + cBest._1.toList.sorted.mkString(" "))
                    }
                    if (noImprove >= stopNoImprove) {
                        finished = true
                        println()
                    } else if (cBest._3 <= 0.0) {
                        finished = true
                        println()
                    }
                    prevLvlNBest = nBest.map(x =>  x._1)
                    lvl += 1
                }
            } else {
                finished = true
                println()
            }
        }

        return (best._1.toList.sorted, best._2, best._3, kA)
    }

    /**
      * Search for a set of initial solutions by first defining pairs of cells among the same sample, and then combining
      * combining pairs coming from different samples using a beam search
      * @param m an expression matrix with cells on the rows and genes on the columns
      * @param samplesIDs for each sample, the list of cells it contains
      * @param samIdCells an array with the sample number of each cell
      * @param geneSum the sum of positive expression of each gene (across all samples)
      * @param hOrd the ordering of cells in decreasing order of rareness score; if empty, generates all pairs of cells
      * @param expr a map containing for each cell the expressed genes
      * @param nNeg the maximum percentage of negative values allowed inside the cluster
      * @param kappa a weighting constant for the out-of-cluster expression
      * @param nHeurPair the number of pairs to form for each sample at level 2
      * @param nHeurKeep the number of top-solutions to consider for expansion at the next level
      * @param verbose enable/disable printing
      * @return the nHeurKeep best solutions of level with highest objective (not necessarily involving all samples)
      * TODO : merge best of different levels ?
      * TODO if start fails (i.e. big subpop in each patient ==> no matching pairs at lvl2), redo without these cells ?
      */
    def start(m: Array[Array[Double]], samplesIDs: Array[List[Int]], samIdCells: Array[Int], geneSum: Array[Double],
              hOrd: Array[Array[(Double, Int)]], expr: Map[Int, List[Int]], nNeg: Double = 0.2, kappa: Double = 1,
              nHeurPair: Int = 100, nHeurKeep: Int = 100, verbose: Boolean = true): List[Set[Int]] = {
        val nbSam = samplesIDs.length

        val nBestPairsQueues = Array.fill(nbSam)(MutPriorityQueue[(Set[Int], List[Int], Double)]()(Ordering[Double].on(x => -x._3)))

        if (verbose) {
            println("*** Starting phase ***")
            print("Evaluating pairs for sample: ")
        }

        if (hOrd.isEmpty) { // All pairs
            for (samId <- samplesIDs.indices) {
                if (verbose) print(samId + " ")
                for (pair <- samplesIDs(samId).combinations(2)) {
                    val (genes, obj) = getGenesMcc(m, pair, expr, geneSum, nNeg, kappa)
                    nBestPairsQueues(samId) += ((pair.toSet, genes, obj))
                    if (nBestPairsQueues(samId).size > nHeurKeep) nBestPairsQueues(samId).dequeue()
                }
            }
        } else { // Use rareness score
            for (samId <- samplesIDs.indices) {
                if (verbose) print(+ samId + " ")
                for (a <- samplesIDs(samId).indices; b <- (a + 1) until (samplesIDs(samId).size min (a + nHeurPair))) {
                    val pair = List(hOrd(samId)(a)._2, hOrd(samId)(b)._2)
                    val (genes, obj) = getGenesMcc(m, pair, expr, geneSum, nNeg, kappa)
                    nBestPairsQueues(samId) += ((pair.toSet, genes, obj))
                    if (nBestPairsQueues(samId).size > nHeurKeep) nBestPairsQueues(samId).dequeue()
                }
            }
        }

        println()

        val nBestQueues = Array.fill(nbSam)(MutPriorityQueue[(Set[Int], List[Int], Double, Set[Int])]()(Ordering[Double].on(x => -x._3)))

        if (verbose) print("Evaluating combinations: ")

        for (samComb <- samplesIDs.indices.combinations(2)) {
            if (verbose) print(samComb.mkString("/") + " ")
            for (pairA <- nBestPairsQueues(samComb(0))) {
                val prevObjs = getObjAllGenes(m, 1, Array.fill(m.length)(0), pairA._1.toList, Array(geneSum), nNeg, kappa).map(x => (x._1, x._2.head, x._3.head))
                for (pairB <- nBestPairsQueues(samComb(1))) {
                    val cells = pairA._1 union pairB._1
                    val (genes, obj) = getGenesMccFromPrevPair(m, cells.size, pairB._1.toList, prevObjs, nNeg, kappa)
                    for (sam <- samComb) {
                        nBestQueues(sam) += ((cells, genes, obj, samComb.toSet))
                        if (nBestQueues(sam).size > nHeurKeep) nBestQueues(sam).dequeue()
                    }
                }
            }
        }

        var nBest = samplesIDs.indices.flatMap(sam => nBestQueues(sam).dequeueAll).toList.distinct.sortBy(-_._3)
        var prevLvlNBest = nBest.map(x => (x._1, x._4))

        if (verbose) println("\n2 samples: " + nBest.size + " solutions, new best: " + nBest.head._3)

        var lvl = 3
        var finished = false
        while (!finished && lvl <= nbSam) {
            if (verbose) print(lvl + " samples: ")
            for ((prevCells, prevSamples) <- prevLvlNBest) {
                val prevObjs = getObjAllGenes(m, 1, Array.fill(m.length)(0), prevCells.toList, Array(geneSum), nNeg, kappa).map(x => (x._1, x._2.head, x._3.head))
                for (newSam <- samplesIDs.indices if !prevSamples.contains(newSam)) {
                    for (pair <- nBestPairsQueues(newSam)) {
                        val cells = prevCells union pair._1
                        val (genes, obj) = getGenesMccFromPrevPair(m, cells.size, pair._1.toList, prevObjs, nNeg, kappa)
                        for (sam <- prevSamples + newSam) {
                            nBestQueues(sam) += ((cells, genes, obj, prevSamples + newSam))
                            if (nBestQueues(sam).size > nHeurKeep) nBestQueues(sam).dequeue()
                        }
                    }
                }
            }

            if (nBestQueues.exists(_.nonEmpty)) {
                val lvlBest = samplesIDs.indices.flatMap(sam => nBestQueues(sam).dequeueAll).toList.distinct.sortBy(-_._3)
                if (verbose) print(lvlBest.size + " solutions")
                if (lvlBest.head._3 >= nBest.head._3) {
                    nBest = lvlBest
                    if (verbose) print(", new best: " + nBest.head._3)
                }
                prevLvlNBest = lvlBest.map(x => (x._1, x._4))
                if (verbose) println()

                lvl += 1
            } else {
                if (verbose) print("no solutions")
                finished = true
            }
        }

        if (verbose) println("*** Beam search ***")

        return nBest.map(_._1)
    }

}
