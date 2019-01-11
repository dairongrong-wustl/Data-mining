euclidean_k_S= read.delim("euclidean_k_S.txt", header = FALSE)
colnames(euclidean_k_S)=c("k","S")
euclidean_k_D= read.delim("eucliean_k_D.txt", header = FALSE)
colnames(euclidean_k_D)=c("k","D")

dot_product_k_S= read.delim("dot_product_k_S.txt", header = FALSE)
colnames(dot_product_k_S)=c("k","S")
dot_product_k_D= read.delim("dot_product_k_D.txt", header = FALSE)
colnames(dot_product_k_D)=c("k","D")

euclidean_k = cbind(euclidean_k_S, euclidean_k_D[2])
euclidean_k$S_OVER_D = euclidean_k$S/euclidean_k$D
dot_product_k = cbind(dot_product_k_S, dot_product_k_D[2])
dot_product_k$S_OVER_D = dot_product_k$S/dot_product_k$D

library("dplyr")
euclidean_results = euclidean_k %>% grouped_df ("k") %>% arrange(S_OVER_D, .by_group = TRUE) %>%
  filter(row_number()==1)

dot_product_results = dot_product_k %>% grouped_df ("k") %>% arrange(S_OVER_D, .by_group = TRUE) %>%
  filter(row_number()==1)
write.csv(euclidean_results, "euclidean_results.csv")
write.csv(dot_product_results, "dot_product_results.csv")