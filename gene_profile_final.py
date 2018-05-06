from gene_expression_comparison_functions import *
import numpy as np

## Define cDNA microarray data
group = np.array([0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1]) # 0= KO mice; 1= control mice
time = np.array([3, 3, 3, 3, 3, 3, 9, 9, 9, 9, 9, 9, 12, 12, 12, 12, 12, 12]) # Time points treatments in weeks (3, 9 and 12)
data_gene,data_array,labels,dic, mylist = load_affy_data('Affymetrix_dataframe.xls',group,time)

## Statistical analysis can be choose bettween oneway or twoway ANOVA (OPTIONAL: you can skip the FDR curves by adding the argument fast0=True)
p_oneway, fdr_oneway, labels_oneway = oneway_ANOVA(dic, labels, mylist,fast0=True)
#p_twoway, fdr_twoway, labels_groups, labels_interaction = twoway_ANOVA(dic, np.array(mylist), data_array, labels,fast0=True)


## Hierarchical Clustering analysis using Pearson Correlation as distance metric to cluster
clusters,matrix_data = Hierarchical_Clustering(dic,time, group, labels_oneway, fdr_oneway)


## Visualization Hierarchical Clustering (projecting the first two principal components)
covariance = (1/len(matrix_data)) * np.matmul(np.transpose(matrix_data), matrix_data)
U,S,V = np.linalg.svd(covariance)
matrix_reduce = np.matmul(matrix_data, U[:,0:2])

k=2 #k= number of clusters
myclus=  fcluster(clusters, k, criterion = 'maxclust')
plt.figure(16)
for kk in range(1, k+1):
    color_list = ['red', 'blue', 'orange', 'green', 'yellow', 'black', 'pink', 'purple']
    plt.scatter(matrix_reduce[myclus==kk, 0], matrix_reduce[myclus==kk, 1], c= color_list[kk-1])    
    cluster_t= labels_oneway[2, 0:fdr_oneway[2]]
    dic['cluster_'+str(kk)]= pd.DataFrame(labels[cluster_t[(myclus==kk)], :] )
    dic['cluster_'+str(kk)].to_csv('Gene_id_cluster-'+str(kk)+'.csv')


plt.show()
