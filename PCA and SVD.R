data = read.csv("training_data_complete_scaled.csv", header = TRUE)
row.names(data)=data[,1]
data = data[,-1]

rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}
gene_pca = data
gene_centralized= as.matrix(gene_pca - rep.row(apply(gene_pca, 2, mean),nrow(gene_pca))) 


cov_gene_centra = 1/((nrow(gene_centralized))-1) * t(gene_centralized) %*% gene_centralized
results = eigen(cov_gene_centra)
pca_eigen_value = results$values
sum_information = sum(results$values)
accumulative = 0
for (i in 1:length(results$values)){
  accumulative = accumulative + results$values[i]
  accu_write = paste(i,"\t",accumulative/sum_information)
  write(accu_write,"accumulate_information.txt", append = TRUE)
}

accumulative_information = read.delim("accumulate_information.txt", header = FALSE)
plot(accumulative_information, xlab ="k",ylab="Accumulative information ")
write.csv(results$vectors[,1:3], "eigen vectors.csv")


#########################################################################
top_3_eigv <- results$vectors[, 1:3]
top_3_coord <- as.matrix(gene_centralized) %*% top_3_eigv
color_code <- rep("red", nrow(top_3_coord))
#change control to "green"
color_code[grep("CON", row.names(top_3_coord))] = "green"

# 3d plot

car::scatter3d(x = top_3_coord[, 1], y = top_3_coord[, 2], z = top_3_coord[, 3],
               point.col = color_code,sphere.size=2,
               grid = FALSE, surface = FALSE)

# 2d plot
top_2_coord <- as.matrix(gene_centralized) %*% results$vectors[, 1:2]
#########################################################################

library("plot3D")
scatter3D(results$vectors[,1], results$vectors[,2], results$vectors[,3], bty = "g", type = "h", 
          pch = 16, cex = 1.5, xlab = "eigen vector 1", ylab = "eigen vector 2",
          zlab = "eigen vector 3", main= "AD cases and normal controls",theta = 20, ticktype = "detailed")


##############explore begin#######################
try_data_x <-1/((nrow(gene_centralized))-1) * t(gene_centralized) %*% gene_centralized
try_pc <- (princomp(try_data_x))
plot(try_pc)
try_pc$loadings
eigen(t(try_data_x) %*% try_data_x)
##############explore end   ######################

##############explore begin#######################
test = matrix(c(1:15), ncol = 3)
r = c(0,0,0,1,1)
test = cbind(r, test)
scatter3D(test[,2], test[,3], test[,4], 50, colvar= test[,1],'filled', ticktype = "detailed")

##############explore end   ######################

##############explore begin#######################
color_code2 = rep("red", nrow(top_3_coord))
color_code2[grep("CON", row.names(top_3_coord))] = "green"
car::scatter3d(top_3_coord[,1],top_3_coord[,2],top_3_coord[,3],point.col = color_code2, 
               sphere.size =2,  grid = FALSE, surface = FALSE)
plot(top_3_coord[,1],top_3_coord[,2],type = "p", col = color_code2, pch=19)

##############explore end   ######################


######################svD: head##########################
data_matrix <- as.matrix(data)
data_svd <- svd(data_matrix)
svd_eigen_value <- data_svd$d
plot(svd_eigen_value)
left_singular <- data_svd$u
right_singular <- data_svd$v
######################SVD: tail##########################
