library(onion) # 四元数
library(geometry) # 3次元ベクトル・クロス積関数
library(Matrix) # 疎行列
library(rgl) # 3Dプロット
library(knitr) # 文書作成・画像埋め込み
library(igraph) # グラフ理論
library(gplots)
library(sphereplot) # mercator
library(geosphere)
library(MASS)
library(foreach)
library(doSNOW)

# cumpute Gaussian curvature of j_th vertex
# input: index of vertices of quanterion vertices object
gaussian_curvature <- function(j, vertices, faces.v){
  #for(j in seq(n_vertices)){
    plane <- t(faces.v[, apply(faces.v==j, 2, any)])
    # shuffle the order of summation of triangle areas
    p <- p1 <- 1 # start
    s <- c(j, setdiff(plane[p,], j)[1]) # 2nd vertex
    p1 <- c(p1, setdiff(which(apply(plane == tail(s, 1), 1, any)), p))
    while(p < nrow(plane) - 1){
      s <- c(s, setdiff(plane[tail(p1, 1),], s))
      p1 <- c(p1, which(sapply(apply(plane, 1, setdiff, s), length) == 1)[2])
      p <- p + 1
    }
    plane <- plane[p1,] # ordered sequentially

    circum_mat <- matrix(1, 3, 3)
    diag(circum_mat) <- -1
    theta_j <- rep(0, nrow(plane))
    circumcenter <- matrix(0, nrow(plane), 3)
    for(p in seq(nrow(plane))){
      hogemat <- t(as.matrix(vertices[plane[p,]])[-1,])
      rname <- plane[p,]
      e12 <- sweep(hogemat[!(rname %in% j),], 2, hogemat[ (rname %in% j),], "-")
      theta_j[p] <- acos( e12[1,] %*% e12[2,]/prod(sqrt(rowSums(e12^2))) )
      if(theta_j[p] < pi/2){ # non-obtuse angle compute circumcenter point
        abc <- mapply(function(k) sum(apply(hogemat[c(k%%3+1, (k+1)%%3+1),], 2, diff)^2), seq(3))
        # abc <- rowSums(rbind(apply(e12, 2, diff), e12)^2)
        abc2 <- circum_mat %*% abc * abc
        circumcenter[p,] <- colSums(diag(c(abc2)) %*% hogemat)/sum(abc2)
        # circumcenter <- (abc2[1]*hogemat[1,] + abc2[2]*hogemat[2,] + abc2[3]*hogemat[3,])/sum(abc2)
      } else { # obtuse angle
        circumcenter[p,] <- colSums(hogemat[!(rname %in% j),])/2 # center point
      }
    }
    mat <- sweep(circumcenter, 2, as.numeric(vertices[j])[-1], "-") # vector from vertex j_th
    n <- nrow(circumcenter)
    Amixed <- mapply(function(k) sqrt(prod(rowSums(mat[c(k, k%%n+1),]^2)) - (mat[k,]%*%mat[k%%n+1,])^2) / 2, seq(n))
    Gcurvature <- (2*pi - sum(theta_j))/sum(Amixed, na.rm=TRUE)
  #}
  return(Gcurvature)
}

b <- mapply(function(j) gaussian_curvature(j, vertices, faces.v), seq(length(vertices)))
cut0 <- cut(b, quantile(b, c(0, 0.02, 0.98, 1)))
cols <- bluered(length(levels(cut0)))[cut0]

plot3d(mesh.tri, type="dits", axes=FALSE, xlab="", ylab="", zlab="")
shade3d(mesh.tri, col="grey")
wire3d(mesh.tri)
spheres3d(t(mesh.tri$vb), col=cols, radius=0.02)
