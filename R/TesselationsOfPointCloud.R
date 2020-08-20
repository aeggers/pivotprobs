###########################################################################################
# code for Andy Eggers, written 8/16/20 by John Nolan,  jpnolan@american.edu
simplicesFromVertices <- function( V, epsilon=1.0e-15 ) {
# compute a tessellation of a collection of vertices in R^n
# V should be an (nV x n) matrix, with each column S[,i], i=1,...,nV giving a point on the unit simplex in R^n
# (note that R and package mvmesh generally uses rows for vectors, so this is nonstandard)
# The code tests that the points are on the unit simplex by testing that column sums are
# (within epsilon) of 1

# test that all points are on the unit simplex
diff <-  max( abs( colSums(V) - 1.0 ))
stopifnot( diff < epsilon )

# project down to (n-1) dimensional space & transpose matrix for call to mvmeshFromVertices( )
n <- nrow(V)
V1 <- t( V[1:(n-1),] )
SVI <- mvmeshFromVertices( V1 )$SVI  # indices of the simplices
S <- array( 0.0, c(n,n,ncol(SVI)) )
for (i in 1:ncol(SVI)) {
  for (j in 1:nrow(SVI)) {
    S[,j,i] <- V[ ,SVI[j,i] ]
  }
}
return(S) }
###########################################################################################

test <- F

if(test){

  ## AE checking what it does compared to my functions
  V <- cbind(c(1,0,0), c(.5, .5,0), c(.5, 0, .5), c(1/3, 1/3, 1/3))
  simplicesFromVertices(V)

  # does this do the same as
  geometry::convhulln(t(V)[,-3])
  mvmeshFromVertices(t(V)[,-3])$SVI
  # TODO: check whether my method is producing redundant simplices. does it get the right result, even if slow?


  # 3d plotting library and mvmesh library (both on CRAN)
  library(rgl)
  library(mvmesh)
  # vertices
  V <- cbind( c(1,0,1e-20), c(0,1,0), c(0,0,1), c(.25,.5,.25), c(.1,.4,.5) )
  V
  S <- simplicesFromVertices( V )
  S

  # plot vertices and simplices
  open3d()
  points3d( t(V), size=10,col="red" )
  for (i in 1:dim(S)[3]) {
    wrap <- cbind( S[,,i], S[,1,i] )
    lines3d(t(wrap), width=3, col='blue')
  }

}
