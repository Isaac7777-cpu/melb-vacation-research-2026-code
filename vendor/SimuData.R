# theta0 = c(range, other parameters,  nugget effect, psill)


SimuData <- function (theta0, rate, loc, cov_model, nug, Srep)
{
  # n: sample size
  n = loc$n
  D = matrix(0, nrow = n, ncol = n)            
  for(i in 1:n) {
    D[i, ] <- sqrt((loc$x[i] - loc$x)^2 + (loc$y[i] - loc$y)^2)
  }

  q = length(theta0)
  eta0 = theta0[-q]
  if (nug == TRUE){
    eta0[q-1] = theta0[q]/(theta0[q-1] + theta0[q])
  } 
  cormat = cor.mat(D, eta0, cov_model, nug)
  
  if (nug == TRUE){
    Sigma0 = (theta0[q-1] + theta0[q])*cormat
  } else{
    Sigma0 = theta0[q]*cormat
  }

  Ymat=rmvnorm(Srep,sigma=Sigma0)
  #ymean = apply(Ymat, 1, mean)
  #Ymat = Ymat - matrix(ymean, nrow = S, ncol = n, byrow = FALSE)
  object <- list (Ymat=Ymat, loc=loc, D=D)
  invisible(object)
}


