AD = read.csv("AD_scaled.csv", header = TRUE)
row.names(AD)=AD[,1]
AD = AD[,-1]
AD_transpose = t(AD)

distDot <- function(x, centers)
{
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  for(k in 1:nrow(centers)){
    for(i in 1:nrow(x)){
      z[i,k] <- (x[i,] %*% centers[k,]) / sqrt(sum(x[i,]*x[i,]))/ sqrt(sum(centers[k,]*centers[k,]))
      z[i,k] <- acos(ifelse(z[i,k] > 1, 1, z[i,k])) / pi * 180
    }
  }
  z
}

library(flexclust)
library(lsa)
for(k in 2:100){
  for(b in 1:10){
    kDot = kcca(AD_transpose, k,control = list(iter.max = 8000,tolerance=1e-08, min.size = 1), family = kccaFamily(dist = distDot, cent = colMeans))
    #calculate S, report k, group number and S in the results
    group = kDot@cluster
    AD_transpose_group=as.data.frame(cbind(AD_transpose,group))
    S_total=0
    for(a in 1:k){
      data_matrix = AD_transpose_group[AD_transpose_group$group==a,]
      data_matrix = as.matrix(data_matrix [,-length(data_matrix)])
      if(nrow(data_matrix)==0){
        S_total=Inf
      }else if(nrow(data_matrix)>1){
        S_matrix = cosine(t(data_matrix))
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
        S_subgroup = 1
      }
      }
    S = S_total/k
    print(S)
    results = paste (k, "\t", S)
    write(results, "dot_product_k_S.txt", append=TRUE)
    
    center = as.matrix(kDot@centers)
    if( nrow(kDot@centers) == k){
      center_matrix = cosine(t(center))
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
      write(results,"dot_product_k_D.txt", append=TRUE)
      }
    }
  }

