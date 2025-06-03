y <- matrix(c(4.45, 4.61, 5.27, 5.00, 5.82, 5.79), byrow=FALSE, nrow=6)
#y
X <- matrix(c(1,1,1,0,0,0, 0, 0,0, 1, 1,1), byrow=FALSE, nrow=6)
# X
Z <- matrix(c(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,1,0,0, 0,0,0,1, 0,0,1,0),byrow=TRUE, nrow=6)
#Z


A <- matrix(c(2,1,11/16,7/8, 1,2,43/32,27/16, 11/16, 43/32,2,91/64, 7/8, 27/16, 91/64,2), nrow=4, byrow=TRUE)
#A
R <- matrix(c(1/18,0,0,0,0,0, 
              0,1/18,0,0,0,0, 
              0,0,1/18,0,0,0,
              0,0,0,1/9,0,0, 
              0,0,0,0,1/9,0, 
              0,0,0,0,0,1/9), nrow=6, byrow=T)
a11 <- t(X) %*% solve(R) %*% X
a12 <- t(X) %*% solve(R) %*% Z
a21 <- t(Z) %*% solve(R) %*% X
a22 <- t(Z) %*% solve(R) %*% Z + solve(A) * 5
lhs <- rbind(cbind(a11, a12), cbind(a21, a22))
lhs


b1 <- t(X) %*% solve(R) %*% y
b2 <- t(Z) %*% solve(R) %*% y
rhs <- rbind(b1, b2)


solve(lhs) %*% rhs
out <- solve(lhs) %*% rhs


library(data.table)
d = fread("/Users/subhashmahamkali/Downloads/pop_info.txt")
# Keep only the second column (V2)
land <- subset(d, V2 == "Landrace")[, .(V1)]
wild <- subset(d, V2 %in% c("Wild"))[, .(V1)]
imp  <- subset(d, V2 == "Improved")[, .(V1)]


fwrite(land, "landrace.txt", sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(wild, "wild.txt",     sep = "\t", quote = FALSE, col.names = FALSE)
fwrite(imp,  "improved.txt", sep = "\t", quote = FALSE, col.names = FALSE)

#bcftools view -S improved.txt -m2 -M2 --min-ac 1 -v snps -Oz -o improved_clean.vcf.gz SorGSD.289snp.miss05.vcf.gz
