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


# cpp で計算したやつ
w <- "~/Desktop/par/"
files <- paste(w, list.files(w, pattern="gcf"), sep="")
ord <- order(as.numeric(gsub("out|_gcf.txt", "", gsub(w, "", files))))
files <- files[ord]
gcfs <- do.call(cbind.data.frame, mapply(read.table, files))
matplot(t(gcfs), type="l", xlab="# of step", ylab="Gaussian curvature", lty=1, lwd=2, cex.lab=1.3)


w <- "~/Desktop/out/"
files <- paste(w, list.files(w, pattern="out"), sep="")
ord <- order(as.numeric(gsub("out|.obj", "", gsub(w, "", files))))
files <- files[ord]
v0 <- mapply(function(x) read.table(x), files, SIMPLIFY=FALSE)
v1 <- mapply(function(x) lapply(split(x[,2:4], x$V1), as.matrix), v0, SIMPLIFY=FALSE)

d <- read.table("/spintransformation/spinxFaringFast/example/bunny.obj")
d0 <- split(d[, 2:4], d$V1)
d0 <- lapply(d0, as.matrix)
bun.f <- d0$f
# v1 <- c(list(d0), v1)


w <- "~/Desktop/par/"
files <- paste(w, list.files(w, pattern="lambda"), sep="")
ord <- order(as.numeric(gsub("out|_lambda.txt", "", gsub(w, "", files))))
files <- files[ord]
l0 <- mapply(function(x) read.table(x), files, SIMPLIFY=FALSE)


theta <- mapply(function(x) 2*acos(x[, 1]), l0) # lambda から求まる回転角
# d = rad * 180 / pi

yl <- c(179/180*pi, pi)
matplot(t(theta), type="l", xlab="# of step", ylab="Radian", ylim=yl, lty=1)
abline(h=178:180/180*pi, lty=3)


w <- "~/Desktop/par/"
# mean curvature of vartice
files <- paste(w, list.files(w, pattern="mcn.txt"), sep="")
ord <- order(as.numeric(gsub("out|_mcn.txt", "", gsub(w, "", files))))
files <- files[ord]
mcn <- mapply(function(x) read.table(x), files, SIMPLIFY=FALSE)
# mcn <- do.call(cbind.data.frame, mcn)
mcn <- mapply(function(z) sqrt(rowMeans(z[, 2:4]^2))/2, mcn)

# mean curvature of face (triangle)
files <- paste(w, list.files(w, pattern="rho.txt"), sep="")
ord <- order(as.numeric(gsub("out|_rho.txt", "", gsub(w, "", files))))
files <- files[ord]
rho <- mapply(function(x) read.table(x), files, SIMPLIFY=FALSE)
rho <- do.call(cbind.data.frame, rho)

w <- "~/Desktop/par/"
files <- paste(w, list.files(w, pattern="rhov"), sep="")
ord <- order(as.numeric(gsub("out|_rhov.txt", "", gsub(w, "", files))))
files <- files[ord]
rhov <- mapply(function(x) read.table(x), files, SIMPLIFY=FALSE)
rhov <- mapply(function(z) sqrt(rowMeans(z[, 2:4]^2))/2, rhov)

w <- "~/Desktop/par/"
files <- paste(w, list.files(w, pattern="vol"), sep="")
ord <- order(as.numeric(gsub("out|_vol.txt", "", gsub(w, "", files))))
files <- files[ord]
vol <- sapply(mapply(function(x) read.table(x), files, SIMPLIFY=FALSE), sum)
plot(vol)

w <- "~/Desktop/par/"
files <- paste(w, list.files(w, pattern="area"), sep="")
ord <- order(as.numeric(gsub("out|_area.txt", "", gsub(w, "", files))))
files <- files[ord]
area0 <- do.call(cbind.data.frame, mapply(function(x) read.table(x), files, SIMPLIFY=FALSE))
area <- colSums(area0)
plot(area)



u <- gcfs
bun.f <- v1[[1]]$f
for(i in seq(v1)){
  bun.v <- v1[[i]]$v
  p <- 10
  cut0 <- c(-Inf, seq(-p, p, length=300), Inf)
  cut1 <- cut(u[,i], cut0, include.lowest=TRUE)
  cols <- colorpanel(length(cut0), "green", grey(0.4), "red")[cut1]
  mesh.tri <- tmesh3d(t(bun.v), t(bun.f),homogeneous=FALSE)
  plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
  shade3d(mesh.tri, col=rep(cols, each=3))
  #rgl.viewpoint(theta=200, phi=0, zoom=0.65)
  #rgl.snapshot(paste("mesh", i, ".png", sep=""))
}

# 正球のガウス曲率か平均曲率に近い色で塗り分け
u <- rho

#conv <- 1/ra^2 # gauss 1/ra^2 mean 2/ra
conv <- median(u[, ncol(u)]) # gauss 1/ra^2 mean 2/ra

delta <- 0 # 収束値より ± 微小区間を設定する場合
for(i in seq(v1)){
  bun.v <- v1[[i]]$v
  cut0 <- c(-Inf, conv-delta, conv+delta, Inf)
  cut1 <- cut(u[,i], cut0, include.lowest=TRUE)
  cols <- c("green", grey(0.4), "red")[cut1]
  cols <- ifelse(u[,i]<conv, "green", "red")
  mesh.tri <- tmesh3d(t(bun.v), t(bun.f),homogeneous=FALSE)
  plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
  shade3d(mesh.tri, col=rep(cols, each=3))
  #rgl.viewpoint(theta=200, phi=0, zoom=0.65)
  #rgl.snapshot(paste("mesh", i, ".png", sep=""))
}


# 球面に色を付けた時、同じ色が付いている閉じた領域をチェックする
K0 <- NULL
for(k0 in 1:nrow(v1[[i]]$f)){
  if( !any(mapply(function(z) k0 %in% z, K0)) ){ # いままでのK に kが含まれていなかったら新規のtriangle
    print(paste("Processing", k0))
    k <- k0
    foo0 <- which(mapply(function(x) any(v1[[i]]$f[x,] %in% v1[[i]]$f[k,]), seq(nrow(v1[[i]]$f))))
    foo1 <- foo0[cols[foo0] == cols[k[1]]] # k の三角形と同じ色をもって隣り合っている三角形たち
    #plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
    #shade3d(mesh.tri, col=rep(cols, each=3))
    #rgl.viewpoint(10, 30, zoom=0.7)
    #m2 <- tmesh3d(t(v1[[i]]$v), t(v1[[i]]$f[foo1,]), homogeneous=FALSE)
    #wire3d(m2)
    K <- list(k, foo1)
    while( !all(K[[length(K)]] %in% K[[length(K)-1]]) ){
      k <- setdiff(K[[length(K)]], K[[length(K)-1]])
      foo0 <- which(mapply(function(x) any(v1[[i]]$f[x,] %in% v1[[i]]$f[k,]), seq(nrow(v1[[i]]$f))))
      foo1 <- foo0[cols[foo0] == cols[k[1]]]
      K <- c(K, list(foo1))
      #open3d()
      #plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
      #shade3d(mesh.tri, col=rep(cols, each=3))
      #m <- v1[[i]]$f[setdiff(K[[length(K)]], K[[length(K)-1]]),]
      #m2 <- tmesh3d(t(v1[[i]]$v), t(v1[[i]]$f[foo1,]), homogeneous=FALSE)
      #points3d(matrix(v1[[i]]$v[m,], nc=3))
      #texts3d(matrix(v1[[i]]$v[m,], nc=3), texts=unique(c(v1[[i]]$f[k,])))
      #wire3d(m2, col=4)
    }
    K1 <- sort(unique(unlist(K)))
    K0 <- c(K0, list(K1))
  }
}

# K0 は同じ色で塗られるtriangle のひとつづきの領域
# 隣り合う領域を探す
V0 <- mapply(function(z) unique(c(v1[[i]]$f[z,])), K0) # 含まれる頂点たち
adj_mat <- mapply(function(z2) mapply(function(z1) any(z2 %in% z1)+0, V0), V0)
g0 <- graph.adjacency(adj_mat, mode="undirected", diag=FALSE)
V(g0)$size <- mapply(function(z) sum(area0[z, i]), K0) * 10
V(g0)$color <- mapply(function(z) cols[z[1]], K0)
plot(g0)


# 変化差分の大きいところをプロットする
# extract large change portion
udif <- t(tail(t(u), -1)-head(t(u), -1))

cut0 <- cut(udif[,i], seq(-max(abs(udif[,i])), max(abs(udif[,i])), length=101), include.lowest=TRUE) # 単純に分ける
cut0 <- lapply(split(u[,i], u[,i] >= 0), quantile, seq(0, 1, by=0.1))
cut0 <- c(cut0[[1]], 0, cut0[[2]]); cut0 <- cut(u[,i], cut0)

cols <- ifelse(udif[,i] > 0, "red", "green")
cols <- bluered(length(levels(cut0)))[cut0]
plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
shade3d(mesh.tri, col=rep(cols, each=3))
rgl.viewpoint(-50, 50, zoom=0.7)


# faring object にプロット
# i to i+1 step の変化率なので、i 番目(スクリプト的にはi-1)
i <- i+1
udif <- t(tail(t(u), -1)-head(t(u), -1))
idx <- which(xor(apply(udif>0, 1, all), apply(udif<0, 1, all))) # 差分の符号変化がないところ
idx <- tail(order(udif[,i]), 100)
idf <- which(apply(apply(v1[[i]]$f, 1, "%in%", idx), 2, any)) # vertex を含むtriangle


cut0 <- cut(udif[,i], seq(-max(abs(udif[,i])), max(abs(udif[,i])), length=101), include.lowest=TRUE) # 単純に分ける
cut0 <- lapply(split(u[,i], u[,i] >= 0), quantile, seq(0, 1, by=0.1))
cut0 <- c(cut0[[1]], 0, cut0[[2]]); cut0 <- cut(u[,i], cut0)

cols <- replace(rep("white", nrow(u)), idx, "blue")
# vertex なら
plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
shade3d(tmesh3d(t(bun.v), t(bun.f), homogeneous=FALSE), col=rep(replace(rep("white", nrow(v1[[i]]$f)), idf, "pink"), each=3))
spheres3d(bun.v[idx,], col=cols[idx], radius=0.002)
rgl.viewpoint(-95.5,59.35,zoom=0.001)

# 耳のダイナミックな変化を追跡する
bun.v <- ratio[i] * v1[[i]]$v
plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
shade3d(tmesh3d(t(bun.v), t(bun.f), homogeneous=FALSE), col=rep(replace(rep("white", nrow(v1[[i]]$f)), idf, "pink"), each=3))
spheres3d(bun.v[idx,], col=cols[idx], radius=0.002)



cols <- ifelse(udif[,i] > 0, "red", "green")
cols <- bluered(length(levels(cut0)))[cut0]
bun.v <- ratio[i-1] * v1[[i-1]]$v
plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
shade3d(tmesh3d(t(bun.v), t(bun.f), homogeneous=FALSE), col=rep(cols, each=3))
rgl.viewpoint(-50, 50, zoom=0.7)




### 等積変換する isovolumic transformation
# 体積比の3乗根がなると信じてやってみる
# (1/3) of volume ratio
i <- 1
ratio <- (vol[1]/vol)^(1/3)
ra <- ((3*vol[1])/4/pi)^(1/3) # 基準となる体積の球の半径
Ara <- 4*pi*ra^2 # 基準となる体積の球の表面積

par(mfrow=c(2, 1), mar=c(2, 4.5, 2, 2), cex.lab=1.5)
plot(vol*ratio^3, xlab="", ylab="Volume", type="o", pch=16)
par(mar=c(5, 4.5, 1, 2))
plot(area*ratio^2, xlab="Step", ylab="Surface area", type="o", pch=16, ylim=c(0, max(area*ratio^2)))
abline(h=Ara, lty=3, lwd=3, col=grey(0.6))


new_v <- list(v = ratio[i] * v1[[i]]$v, f=v1[[i]]$f)


open3d()
mesh.tri <- tmesh3d(t(new_v$v), t(bun.f),homogeneous=FALSE)
plot3d(bun.v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
shade3d(mesh.tri, col=rep(cols, each=3))


# メルカトル図法を試す
# Mercator projection
s <- car2sph(x=new_v$v[,1], y=new_v$v[,2], z=new_v$v[,3], deg=FALSE)
merc2d <- mercator(s[,1:2], r=ra)
plot(merc2d, pch=16, cex=0.03)

# 頂点の経度・緯度移動
# transision of latitude and longitude
for(i in seq(ratio)){
  new_v <- list(v = ratio[i] * v1[[i]]$v, f=v1[[i]]$f)
  s <- car2sph(x=new_v$v[,1], y=new_v$v[,2], z=new_v$v[,3], deg=FALSE)
  merc2d <- mercator(s[,1:2], r=ra)
  plot(merc2d, pch=16, cex=0.03)
}


new_v <- list(v = ratio[i] * v1[[i]]$v, f=v1[[i]]$f)
s <- car2sph(x=new_v$v[,1], y=new_v$v[,2], z=new_v$v[,3], deg=FALSE)
sdeg <- car2sph(x=new_v$v[,1], y=new_v$v[,2], z=new_v$v[,3], deg=TRUE)
merc2d <- mercator(s[,1:2], r=ra)
plot(merc2d, pch=16, cex=0.03)

cut0 <- mapply(function(i) cut(sdeg[,i], seq(-180/i, 180/i, length=101), include.lowest=TRUE), 1:2, SIMPLIFY=FALSE)
cut1 <- cut0[[1]]
cols <- colorpanel(length(levels(cut1)), "green", grey(0.4), "red")[cut1]

mer01 <- mapply(function(i) seq(-180/i, 180/i, by=5), 1:2)

plot(merc2d, pch=16, cex=0.03, col=cols)
plot(sdeg[,1:2], pch=16, cex=0.03, col=cols)
abline(v=mer01[[1]], h=mer01[[2]], lty=3)


# 距離の分布
di <- c(mapply(function(x) c(dist(merc2d[ new_v$f[x,], ])), seq(nrow(new_v$f))))
dq <- quantile(di, 0.992)

# 三角形の面積の分布(参考)
ai <- mapply(function(x){abs(det(merc2d[ new_v$f[x,][1:2],] - merc2d[ new_v$f[x,][3],]))/2}, seq(nrow(new_v$f)))

# coloring mercator projection
plot(merc2d, pch=16, cex=0.03)
for(j in 1:nrow(new_v$f)){
  tmp <- merc2d[ new_v$f[j, c(1:3,1)],]
  if( all(c(dist(tmp[1:3,])) < dq) ) polygon(tmp, col=cols[j], border=1)
}

l <- 2
xyz <- mapply(function(i){
  new_v <- list(v = ratio[i] * v1[[i]]$v, f=v1[[i]]$f)
  s <- car2sph(x=new_v$v[,1], y=new_v$v[,2], z=new_v$v[,3], deg=FALSE)
  xyz <- sph2car(s[l,"long"], s[l, "lat"], r=seq(0, 0.9, length=100), deg=FALSE)
  xyz
}, seq(length(v1)-1), SIMPLIFY=FALSE)
for(n in seq(xyz)) points3d(xyz[[n]], col=rainbow(length(xyz))[n])


mesh.tri <- tmesh3d(t(new_v$v), t(new_v$f),homogeneous=FALSE)
plot3d(new_v$v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
shade3d(mesh.tri, col=rep(cols, each=3))
for(m in seq_along(xyz)) points3d(xyz[[m]], col=rainbow(length(xyz))[m])


i <- 1
new_v <- list(v = ratio[i] * v1[[i]]$v, f=v1[[i]]$f)
s <- car2sph(x=new_v$v[,1], y=new_v$v[,2], z=new_v$v[,3], deg=FALSE)
idx <- mapply(function(l){
    norm_sph <- c(cos(s[l,1])*cos(s[l,2]), sin(s[l,1])*cos(s[l,2]), sin(s[l,2]))
    abc <- mapply(function(i){
        hoge_tri <- v1[[1]]$v[new_v$f[i,],]
        Rad <- 1/sum(base::solve(hoge_tri, rep(-1, 3)) * norm_sph)
        abc <- -base::solve(t(hoge_tri), Rad*norm_sph)
        if( all(abc > 0) & colSums(abc*hoge_tri)%*%norm_sph > 0) colSums(abc*hoge_tri)
      }, seq(nrow(new_v$f)))
    triangle_points <- do.call(rbind.data.frame, abc)
    colnames(triangle_points) <- NULL
    triangle_points
  }, 1:5, SIMPLIFY=FALSE)

m <- which(sapply(abc, length)>0)
i <- 1
new_v <- list(v = ratio[i] * v1[[i]]$v, f=v1[[i]]$f)
mesh.tri <- tmesh3d(t(new_v$v), t(new_v$f),homogeneous=FALSE)
plot3d(new_v$v, type="n", axes=FALSE, xlab="", ylab="", zlab="")
wire3d(mesh.tri, col=rep(cols, each=3))
n <- 30
points3d(xyz[[n]], col=rainbow(length(xyz))[n])
triangles3d(v1[[1]]$v[ v1[[1]]$f[m[3],], ], col=rep(cols[m[3]], each=3))
points3d(v1[[1]]$v[ v1[[1]]$f[m[3],], ])


tri_xyz <- mapply(function(i){
    new_v <- list(v = ratio[i] * v1[[i]]$v, f=v1[[i]]$f)
    s <- car2sph(x=new_v$v[,1], y=new_v$v[,2], z=new_v$v[,3], deg=FALSE)
    norm_sph <- c(cos(s[l,1])*cos(s[l,2]), sin(s[l,1])*cos(s[l,2]), sin(s[l,2]))
    abc <- mapply(function(f){
        hoge_tri <- v1[[1]]$v[new_v$f[f,],]
        Rad <- 1/sum(base::solve(hoge_tri, rep(-1, 3)) * norm_sph)
        abc <- -base::solve(t(hoge_tri), Rad*norm_sph)
        if( all(abc > 0) & colSums(abc*hoge_tri)%*%norm_sph > 0) colSums(abc*hoge_tri)
      }, seq(nrow(new_v$f)))
    triangle_points <- do.call(rbind.data.frame, abc)
    colnames(triangle_points) <- NULL
    triangle_points
  }, seq(length(v1) - 1), SIMPLIFY=FALSE)

col1 <- rainbow(length(tri_xyz))
for(m in 2:length(tri_xyz)) points3d(tri_xyz[[m]], col=col1[m])





