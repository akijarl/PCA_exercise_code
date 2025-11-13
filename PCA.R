#test code for PCA concepts

x <- rep(1:10)
y <- rep(1,10)

plot(x,y)

mean(x)
x-mean(x)
(x-mean(x))^2
sum((x-mean(x))^2)
sum((x-mean(x))^2) * 1/(length(x)-1)

cov(x,x)
var(x)

var(y)

cov(x,y)

z <- x+2

plot(x,z)

var(x)
var(z)

cov(x,z)


T0 <- prcomp(cbind(x,x),center = TRUE)
T1 <- prcomp(cbind(x,z),center = TRUE)

biplot(T0)
biplot(T1)

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
