AD = read.csv("AD_scaled.csv", header = TRUE)
row.names(AD)=AD[,1]
AD = AD[,-1]
AD_transpose = t(AD)

library(flexclust)

for(k in 2:100){
  for(b in 1:20){
    kDot = kcca(AD_transpose, k, family = kccaFamily("kmeans"))
    #calculate S, report k, group number and S in the results
    group = kDot@cluster
    AD_transpose_group=as.data.frame(cbind(AD_transpose,group))
    S_total=0
    for(a in 1:k){
      data_matrix = AD_transpose_group[AD_transpose_group$group==a,]
      data_matrix = as.matrix(data_matrix [,-length(data_matrix)])
      if(nrow(data_matrix)==0){
        S_total=Inf
      }else if(nrow(data_matrix) >1){
        S_matrix = as.matrix(dist(data_matrix, method = "euclidean"))
        sum = 0
        count = 0
        for(i in 1:ncol(S_matrix)){
          for(j in 1:ncol(S_matrix)){
            if(i<j){
              sum=sum + S_matrix[i,j]
              count=count+1}
          }
        }
        S_subgroup=sum/count
        print("subgroup")
        print(a)
        print(S_subgroup)
        S_total = S_total + S_subgroup
        print("S_total")
        print(S_total)
      }else{
        S_subgroup = 0
      }
    }
    S = S_total/k
    print(S)
    results = paste (k, "\t", S)
    write(results, "euclidean_k_S.txt", append=TRUE)
    
    center = as.matrix(kDot@centers)
    if( nrow(kDot@centers) == k){
      center_matrix = as.matrix(dist(center, method = "euclidean"))
      sum = 0
      count = 0
      for(i in 1:ncol(center_matrix)){
        for(j in 1:ncol(center_matrix)){
          if(i<j){
            sum=sum + center_matrix[i,j]
            count=count+1}
        }
      }
      D = sum/count
      results = paste (k, "\t",D)
      write(results,"eucliean_k_D.txt", append=TRUE)
    }
  }
}