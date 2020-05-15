# R function to compute LU factorisation of Lower Hessenberg matrix (W.L.O.G.) in optimal time complexity with 
# principles of numerical and conventional linear algebra
LUHessenberg <- function(input){ 
  n <- ncol(input)
  for (i in 1:(n-1)){
    for (j in (i+1):n){ 
      input[j,i] <- input[j,i]/input[i,i] #store the coefficients for L part
      input[j,i+1] <- input[j,i+1] - input[j,i]*input[i,i+1] #update the upper triangular, non diagonal values (diagonal values already 0)
    }
  }
  return(input) #output is both L and U matrices composed into one dense matrix, per space complexity of numerical linear algebra
}

# function to create test input (lower) Hessenberg matrix. This step can be can be processed trivially using 
# the tril() and diag() methods in matlab
makemat <- function(v){ 
  matv <- matrix(v, nrow = 1, ncol = length(v)) #make vector v into a row matrix
  col1 <- t(matv) #transpose
  part1<-NULL
  for (i in 1:length(v)){
    part1 <- cbind(part1, col1 + v[i])
  }
  part1[upper.tri(part1)] <- 0 # set upper tri to 0
  part2 <- diag(-v[1:length(v)-1])
  part2 <- cbind(rep(0, length(v)), rbind(part2, rep(0, length(v)-1)))
  return(part1 + part2)
}
#test round(LUHessenberg(makemat(c(1,2,3,4,5,6,7))), digits = 4)
#test2 round(LUHessenberg(makemat(c(1,2,10,1,1,4,7))), digits = 4)
