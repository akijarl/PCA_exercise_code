library(reshape2)
#library(tidyverse)


whale.pca <- prcomp(whale[,-1], center = TRUE,scale. = TRUE)

lapply(whale, FUN=class)
by(whale, INDICES=whale$Location, FUN=function(z){
  # Use numeric fields for the PCA
  pca <- prcomp(z[,unlist(lapply(z, FUN=class))=="integer"])
  plot(pca$x[,1:2], pch=16, main=z[1,"Location"]) # 2 first principal components
  z
})

biplot(whale.pca, pc.biplot = T)
biplot(pca)

whale %>%
  gather("Type", "Value",-Species) %>%
  ggplot(aes(Species, Value, fill = Type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw()

whale %>%
  gather("Type", "Value",-Location) %>%
  ggplot(aes(Location, Value, fill = Type)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw()

pca <- prcomp(whale[,unlist(lapply(whale, FUN=class))=="integer"])
plot(pca$x[,1:2], pch=16, col=as.numeric(whale[,"Location"]), main="Blue Whale Populations") # 2 first principal components
legend("bottomlef", pch=16, col=unique(as.numeric(whale[,"Location"])), legend=unique(whale[,"Location"]))






str(whale.pca)
plot(whale.pca)
biplot(whale.pca)




colnames(dat) <- snp$type
rownames(dat) <- whale[,1]

whale.pca <- prcomp(dat,
                  center = TRUE,
                  scale. = TRUE)


summary(whale.pca)
plot(whale.pca)
biplot(whale.pca)
biplot(cows.pca)

prcomp(mtcars[,c(1:7,10,11)], center = TRUE,scale. = TRUE)

mat <- scale(matrix(rnorm(20), 4, 5))
dimnames(mat) <- list(paste("Sample", 1:4), paste("Var", 1:5))

# Perform PCA
myPCA <- prcomp(mat, scale. = F, center = F)
myPCA$rotation # loadings
myPCA$x

require(adegenet)
data(microbov)
x.cows <- tab(microbov, freq=TRUE, NA.method="mean")

pca.cows <- dudi.pca(x.cows, center=TRUE, scale=FALSE)
s.label(pca.cows$li)

plot(whale.pca)
# log transform 
log.p <- log(snp$p)
snp.type <- snp$type

Eigenvalues <- eigen(cov(us.bis.yield))$values
Eigenvectors <- eigen(cov(us.bis.yield))$vectors

AllCallAnnot <- snpgdsOpen("All_calls_annot_NRecip_TSL.gds")

pcaC <- snpgdsPCA(AllCallAnnot,autosome.only=F)

pc.percent <- pcaC$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pcaC$sample.id,
                  EV1 = pcaC$eigenvect[,1],
                  EV2 = pcaC$eigenvect[,2],
                  EV3 = pcaC$eigenvect[,3],
                  EV4 = pcaC$eigenvect[,4],
                  EV5 = pcaC$eigenvect[,5],
                  EV6 = pcaC$eigenvect[,6],
                  stringsAsFactors = FALSE)
head(tab)
