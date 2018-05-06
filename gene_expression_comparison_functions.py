
from scipy import stats
from statsmodels.stats.multitest import local_fdr
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)

def load_affy_data(file_name,group,time,columns_label = ['group', 'time', 'gene'],sheet='sheet1'):
    # Data is  loaded in standard affymetrix array format where:
        # rows correspond to the signalintensity of studied genes 
        # columns correspond to: probe id, symbol gene, description of gene and subsequently number of samples analyzed
    # This outputs data_gene, the matrix with all the relevant data, data_array which also contains the respective groups and times as the first 
    #two rows, labels, with the gene labels, and a dictionary with various quantities of interest that can be used for 

    N=group.shape[0];
    myfile = pd.ExcelFile(file_name) 
    data = myfile.parse(sheet)
    data = data.values
    labels = np.array(data[:,[1, 2]])
    
    data_gene = np.transpose(np.array(data[:, 3:3+N], dtype='f'))
    data1 =  np.append(np.transpose([group]), data_gene, axis=1)
    data_array =  np.append(np.transpose([time]), data1, axis=1) 
    
    dic = {}
    
    timelist=list(sorted(set(time)))
    
    data_gene=np.transpose(data_gene)
    
    for i in timelist:
        dic['mean_ctrl_' + str(i)] = np.mean(data_gene[:, (group==1) * (time==i)], axis=1)
        dic['KO_' + str(i)] = data_gene[:,(group==0) * (time==i)]
        dic['ctrl_' + str(i)] = data_gene[:,(group==1) * (time==i)]
        dic['ratio_' + str(i)] = np.transpose(np.transpose(dic['KO_' + str(i)]) / dic['mean_ctrl_' + str (i)])

    return np.transpose(data_gene), data_array , labels, dic,  timelist
    

def gene_comparison(p_values, q_list, i, q0=0.27,fast=False):
    # Identification of differentially expressed genes after accounting for multiple testing measuring False Discovery Rate. 
    # Given the p_values of the genes, we measure the number of genes under the set FDR (q0).

    N = len(p_values)
    x = np.arange(1, N+1)
    p_values = np.sort(p_values)
    fdr = sum(p_values < (q0 * x)/N)
    
    if fast==False:
        for j in q_list:
            if j == 0 :
                fdr_total = sum(p_values < (j*x)/N)
            else:
                fdr = sum(p_values < (j*x)/N)
                fdr_total = np.append(fdr_total, [[fdr]])
        
        plt.figure(i)
        plt.plot(fdr_total,q_list, 'g.')
        plt.title(str(i)+' weeks')  
        plt.xlabel('number of genes')
        plt.ylabel('False Discovery Rate')
    
        plt.figure(i+1)
        
        plt.plot(x, p_values, 'b.', label='p(x)')
        plt.plot(x, (q0*x)/N, 'r', label='q(x)/N')
        plt.title(str(i)+' weeks')   
        plt.xlabel('number of genes')
        plt.ylabel('p values')
        plt.legend()
    
    return fdr
    
    
    
def oneway_ANOVA(dic, labels, timelist, fast0 ,q_list=np.arange(0,0.51,0.001)):
    # Given our dictionary, the labels, and the list of times, it returns:
        # p_table using one_way anova: the list of p_values for each time point 
        # fdr_pair: the number of genes under the set FDR foreach time point
        # A_table: list of gene sorted by their p_value (each row correspond to the gene location and each column correspond to a different time point)   
    
    # Gene comparison between KO and control of same time point (pairwise comparison)   
    fdr_pair = np.zeros(len(timelist), dtype='int')
    j= 0
    
    for i in timelist:
        dic['F_' + str(i)], dic['p_' + str(i)] = stats.f_oneway(np.transpose(dic['ctrl_' + str(i)]),np.transpose(dic['KO_' + str(i)]))
        p_values = dic['p_' + str(i)]
        
        
        fdr_pair[j] = gene_comparison(p_values, q_list, i,fast=fast0)
        
        A = np.argsort(dic['p_' + str(i)])
        gene_selec_pair = pd.DataFrame(np.append(labels[A[0:fdr_pair[j]], :], np.transpose([p_values[A[0:fdr_pair[j]]]]), axis=1), columns= ['Gene ID', 'Gene Description', 'p'])
        gene_selec_pair.to_csv('One-Way-ANOVA_'+str(i)+'weeks.csv')
        j = j+1
        if i == timelist[0]:
            p_table = [p_values]
            A_table = [A]
        else:
            p_table = np.append(p_table, [p_values], axis=0)
            A_table = np.append(A_table, [A], axis=0)
    
    FDR_pairwise = pd.DataFrame(fdr_pair, index=['FDR (3weeks)', 'FDR (9weeks)', 'FDR (12weeks)'], columns=['Gene selected'])
    FDR_pairwise.to_csv('FDR-pairwise_table.csv')
    
    # Validation p_values selected by Volcano Plot 
    for i in  timelist:
        Xi= np.transpose([np.mean(dic['ratio_'+str(i)],  axis=1)])
        if i == timelist[0]:
            matrix = Xi
        else:
            matrix = np.append(matrix, Xi,axis=1)
        
    Xt = np.log2(matrix)
    Y = - np.log(np.transpose(p_table))
    lab_oneway = np.transpose(A_table)
    
    for k in range (0, 3):
        lab=lab_oneway[0:fdr_pair[k],k]
        plt.figure(k)
        plt.scatter(Xt[:,k], Y[:,k], c='r')
        plt.scatter(Xt[lab,k], Y[lab,k], c='b')   
        plt.title('Volcano Plot of '+str(timelist[k])+' weeks')
        plt.xlabel('Log2(ratio)')
        plt.ylabel('-Log10(p_value)')   
    
    return p_table, fdr_pair, A_table
        

    
def twoway_ANOVA(dic, timelist, data_array, labels,fast0, q_list=np.arange(0,0.51,0.001)):
    # Given our dictionary, the labels, and the list of times, it returns:
        # p_table using two_way anova: the list of p_values: between group, versus time or the interaction of both factors 
        # fdr: the number of genes under the set FDR under the following conditions (between group, versus time or the interaction of both factors)
        # A_groups: list of gene sorted by their p_value resulting  from the comparison between groups (each row correspond to the gene location)   
        # A_interaction: list of gene sorted by their p_value resulting from the interaction between groups and time (each row correspond to the gene location)  

    N = data_array.shape[1]
    
    for i in range(2, N):
        selec = [0, 1, i]
        columns_label = ['group', 'time', 'gene']
        data_frame1 = pd.DataFrame(data=data_array[:, selec], columns=columns_label, dtype='f')
        formula = 'gene ~ C(group) + C(time) + C(group):C(time)'
        model = ols(formula, data_frame1).fit()
        aov_table = anova_lm(model, typ=2)     
        p_values = np.array(aov_table.iloc[0:3, 3])
        if i == 2:
            p_table = [p_values]
        else:
            p_table = np.append(p_table, [p_values], axis=0)
     
    #Comparison over group
    q_list = np.arange(0, 0.51, 0.001)
    fdr = np.zeros(3, dtype= 'int')
    
    for i in range(0,3):
        fdr[i]= gene_comparison(p_table[:,i], q_list, timelist[i],fast=fast0) # Column order: p values by groups, p values by time and p values by interaction groups-time
    
    fdr_table = pd.DataFrame(fdr, index=['FDR (groups)', 'FDR (time)', 'FDR (interaction)'], columns=['Gene selected'])
    fdr_table.to_csv('FDR_table.csv')
    
    A_groups = np.argsort(p_table[:, 0])
    gene_selec_groups = pd.DataFrame(np.append(labels[A_groups[0:fdr[0]], :], p_table[A_groups[0:fdr[0]], :], axis=1), columns= ['Gene ID', 'Gene Description', 'p(group)', 'p(time)', 'p(interaction)'])
    gene_selec_groups.to_csv('Gene_selec_groups.csv')
    
    A_interaction = np.argsort(p_table[:, 2])
    gene_selec_interaction = pd.DataFrame(np.append(labels[A_interaction[0:fdr[2]], :], p_table[A_interaction[0:fdr[2]], :], axis=1), columns= ['Gene ID', 'Gene Description', 'p(group)', 'p(time)', 'p(interaction)'])
    gene_selec_interaction.to_csv('Gene_selec_interaction.csv')    
    
    return p_table, fdr, A_groups, A_interaction


def Hierarchical_Clustering(dic, time, group, labels_oneway, fdr_oneway):
    # Hierarchical clustering using the Pearson correlation as the distance metric to cluster the selected genes from one_way ANOVA. 
    # Returns the linkage matrix Z and plots the respective clustered heat maps. 
    
    
    timelist=list(sorted(set(time)))
    for i in timelist:
        dic['selec_KO_'+str(i)] = dic['KO_'+str(i)][labels_oneway[2, 0:fdr_oneway[2]], :]
        dic['selec_ctrl_'+str(i)] = dic['ctrl_'+str(i)][labels_oneway[2, 0:fdr_oneway[2]], :]    
        dic['selec_KO_mean'+str(i)] = np.mean(dic['selec_KO_'+str(i)],  axis=1)
        dic['selec_ctrl_mean'+str(i)] = np.mean(dic['selec_ctrl_'+str(i)], axis=1)
        if i ==timelist[0]:
            Xc = np.append(dic['selec_ctrl_'+str(i)], dic['selec_KO_'+str(i)], axis=1)
        else:
            Xc0 = np.append(Xc, dic['selec_ctrl_'+str(i)], axis=1)
            Xc = np.append(Xc0, dic['selec_KO_'+str(i)],axis=1)

    
    def pearson_dist(a,b):
       return 1 - stats.pearsonr(a,b)[0]
    
    Z = linkage(Xc, 'average', metric=pearson_dist)
    
    c, coph_dists = cophenet(Z, pdist(Xc))
    c
    plt.figure(figsize=(25, 10))
    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(Z, leaf_rotation=90.,leaf_font_size=8.)
#    dendrogram(Z,leaf_rotation=90.,  # rotates the x axis labels
#        leaf_font_size=8.,  # font size for the x axis labels
#    )
    
    #def greentoredpalette():
    #    numpoints=254;
    #    numblack=-2;
    #    half=int(numpoints/2)
    #    redColorMap = np.append(np.zeros([1,half+numblack]), np.linspace(0,255,num= half-numblack).astype(int));
    #    greenColorMap = np.append(np.linspace(255,0,num= half-numblack).astype(int),np.zeros([1,half+numblack]));
    #    mymap=np.append(np.array([redColorMap,greenColorMap]),np.zeros([1,numpoints]),axis=0);
    #    mymap=np.transpose(mymap.astype(int))
    #    clist=[]
    #    for i in range(0,mymap.shape[0]):
    #    	clist.append("#{:02x}{:02x}{:02x}".format(mymap[i,0],mymap[i,1],mymap[i,2]))
    #    
    #    return clist
    #
    #mypal=sns.color_palette(greentoredpalette())
    mypal = sns.light_palette("green")
    
    group_label=['KO','Control'];
    
    xprelabels=[];
    for i in range(0,time.shape[0]):
        xprelabels.append(group_label[group[i]]+' w'+str(time[i]))
    xprelabels=np.array(xprelabels);  

#  g = sns.clustermap(Xc, metric= pearson_dist,vmin= 3, vmax= 15, cmap=mypal)
    g = sns.clustermap(Xc, metric= pearson_dist,vmin= 3, vmax= 15, cmap=mypal, yticklabels=False)

    
    g.ax_heatmap.set_xticklabels(xprelabels[leaves_list(g.dendrogram_col.linkage)])
    return Z, Xc