Example of usage: ¡®aa_kde_em_clustering.m¡¯

Main function: Idx = kde_em_clustering(files,Par)

Input:
-¡®files¡¯ is a cell containing the directories to data files. Each data file contains data from one measurement. The number of data files can be 2 or more.
-¡®Par¡¯ is a structure containing the parameters for the clustering algorithm. There are 6 parameters:
-Par.numcluster: Number of clusters, needs to be specified.
-Par.normalize:  Optional, 2(default) - normalize positive numbers to [0.5,1] and normalize negative numbers to [0,0.5], 1 - normalize all numbers to [0,1], 0 - do not normalize. Normalizing data before clustering can avoid one measurement dominating others. 
-Par.anchor:     Optional, number of anchor points, default: 100. Increasing this number can improve accuracy of clustering result, but can increase running time.
-Par.maxit:      Optional, maximum iterations of EM algorithm, default: 400. Increasing this number can improve accuracy of clustering result, but can increase running time.
-Par.Leps:       Optional, termination criteria of EM algorithm, default: 1. Decreasing this number can improve accuracy of clustering result, but can increase running time.
-Par.plot:       Optional, 1(default) - plot clustering result. The algorithm plots the clustering results in 2D figures (PCA dimension reduction).

Output:
-¡®Idx¡¯ is a vector indicating the cluster label for each plant.