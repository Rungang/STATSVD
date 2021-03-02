# Source code for STATSVD and two partial algorithm version.
# We also implement several baselines including HOSVD, HOOI, SHOSVD, SHOOI.


require(rTensor)
require(ssvd)
require(MASS)

#' Calculate the Square of l2 norm of a vector.
#' @param x Input vector.
#' @return Square of l2 norm of a vector.
vector_sum_square <- function(x){
  return(sum(x^2))
}


#' Calculate the sine-theta distance between two subspaces
#' @param U First subspace
#' @param hatU Second subspace
#' @param q Hilbert-Schmidt norm
#' @return q-sine-theta distances
sine.theta <- function(U, hatU, q){ # Sine-theta distance between two singular subspaces.
  # U and hatU should be of the same dimensions
  try(if(missing("q")) q = 2)
  try(if(nrow(U)!=nrow(hatU) | ncol(U)!=ncol(hatU)) stop("Matrix does not match") )

  r = ncol(U)
  v = 1 - (svd(t(U) %*% hatU)$d)^2
  if(is.infinite(q))
    return(max(v))
  else
    return((sum(v^q))^(1/q))
}



#' Estimate Gaussian noise level via MAD.
#'
#' @param Y Input tensor, with dimension p1*p2*...*pd.
#' @return Estimate of standard deviation of Gaussian noise.
get.sigma <- function(Y){
  if(missing("Y")) stop("missing argument: Y is required as the tensor type.")
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  X = k_unfold(Y,1)@data
  V = as.vector(X)
  return(1.4826*median(abs(V)))
}



#' Estimate Tucker rank of a given tensor
#' @param Y Input tensor, with dimension p1*p2*...*pd.
#' @param sigma (Estimated) Gaussian noise level
#' @param sparse_mode A d-dimensional logicical vector that indicates whether each mode is sparse.
#' @return Estimated rank as a d-dimensinal integer tuple
get.rank <- function(Y, sigma, sparse_mode){
  if(missing("Y")) stop("missing argument: Y is required as the tensor type.")
  if(missing("sigma")) stop("missing argument: sigma is required as the non-negative double")
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  try(if(sigma <= 0) stop("invalid input: sigma must be positive."))
  p = dim(Y)
  d = length(p)
  r = rep(0,d)
  try(if(missing("sparse_mode")) sparse_mode=rep(TRUE, d))
  p.prod = prod(p)
  p.minus.k = prod(p) / p
  p.ast = max(p)

  I_0 = list()
  for(i in 1:d){ # Initialization step 1: find significant index set I_k
    Y_i = k_unfold(Y, i)@data
    Y_i_row_norm = apply(Y_i, 1, vector_sum_square)
    #consider sparse model!
    if(isTRUE(sparse_mode[i])){
      I_0 = c(I_0, list((apply(Y_i, 1, vector_sum_square) > sigma^2*(p.minus.k[i]+2*sqrt(p.minus.k[i]*log(p.prod))+2*log(p.prod)))
                          |   (apply(abs(Y_i), 1, max) > 2*sigma*sqrt(log(p.prod)))))
    }
    else
      I_0 = c(I_0, list(rep(TRUE,p[i])))
    # select the significant indices
  }
  tilde_Y = Y
  s = rep(0,d)
  for (i in 1:d){ # Initialization step 2: construct tilde_Y
    tilde_Y = ttm(tilde_Y, diag(I_0[[i]]*1), i)
    s[i] = sum(I_0[[i]]*1)
  }

  for(i in 1:d){
    MY = k_unfold(tilde_Y,d)@data
    p.k = dim(MY)[1]
    p.minus.k = dim(MY)[2]
    s.value = svd(MY)$d
    ci = s[i]
    cj = prod(s)/s[i]
    delta = sqrt(ci) + sqrt(cj) + sqrt(2*ci*(1+log(p.k/ci)) + 2*cj*(1+log(p.minus.k/cj))
                                           + 4*log(p.k))
    r[i] = length(which(s.value >= sigma*delta))

  }
  return(r)
}



#' Sparse low-rank tensor decomposition using High-order SVD
#' @param Y Input tensor, with dimension p1*p2*...*pd.
#' @param r (Estimated) Tucker ranks.
#' @param sparse_mode A d-dimensional logicical vector that indicates whether each mode is sparse.
#' @return Estimated subspace in each mode.
HOSVD <- function(Y, r, sparse_mode){
  try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))

  p = dim(Y)
  d = length(p)
  try(if(missing("sparse_mode")) {
    sparse_mode = rep(FALSE,3)
    })

  if(is.atomic(r) && length(r)==1){
    r = rep(r, d)
  }


  U_0 = list();
  for (i in 1:d){
    MY = k_unfold(Y, i)@data
    if(!isTRUE(sparse_mode[i])){
      U_0 = c(U_0, list(svd(MY)$u[,1:r[i]]))
    }
    else{
      # pc = nsprcomp(t(MY), ncomp = r[i], k = s[i], center = FALSE)
      # temp = pc$rotation
      temp = ssvd(MY, method="method",r=r[i])$u
      colnames(temp) = NULL
      U_0 = c(U_0, list(temp))
    }
  }

  return(U_0)
}


#' Sparse low-rank tensor decomposition using High-order orthogonal iteration
#' @param Y Input tensor, with dimension p1*p2*...*pd.
#' @param r (Estimated) Tucker ranks.
#' @param sparse_mode A d-dimensional logicical vector that indicates whether each mode is sparse.
#' @param tmax Maximal iteration times.
#' @param vartol Covergence telorance.
#' @return Estimated subspace in each mode.
HOOI <- function(Y, r, sparse_mode , tmax, vartol){
  try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  try(if(missing("tmax")) tmax = 5)
  try(if(missing("vartol")) vartol = 1e-6)

  p = dim(Y)
  d = length(p)

  try(if(missing("sparse_mode")) {
    sparse_mode = rep(FALSE,3)
  })

  if(is.atomic(r) && length(r)==1){
    r = rep(r, d)
  }


  U_t = list();


  for (i in 1:d){#Initialization
    MY = k_unfold(Y, i)@data
    if(!isTRUE(sparse_mode[i])){
      U_t = c(U_t, list(t(svd(MY)$u[,1:r[i]])))
    }
    else{
      temp = ssvd(MY, method="method",r=r[i])$u
      colnames(temp) = NULL
      U_t = c(U_t, list(t(temp)))
    }
  }

  t = 1;
  approx = -1

  while(t<tmax){ # Stop criterion: convergence or maximum number of iteration reached
    for(i in 1:d){
      A = ttl(Y, U_t[-i], (1:d)[-i])
      A_matrix = k_unfold(A, i)@data
      if(!isTRUE(sparse_mode[i])){
        svd.result = svd(A_matrix)
        U_t[[i]] = t(svd.result$u[,1:r[i]])
        svector = svd.result$d[1:r[i]]
      }
      else{
        temp = ssvd(A_matrix, method="method",r=r[i])
        U_t[[i]] = t(temp$u)
        svector = temp$d
      }
    }
    if (abs(sum(svector^2) - approx) > vartol){
      #print(t)
      t = t+1
      approx = sum(svector^2)
    }
    else {
      break
    }
  }
  for(i in 1:d){
    U_t[[i]] = t(U_t[[i]])
  }
  return(U_t)
}


#' Sparse low-rank tensor decomposition using STATSVD.
#' See Zhang and Han (2019) for reference.
#' @param Y Input tensor, with dimension p1*p2*...*pd.
#' @param r (Estimated) Tucker ranks.
#' @param sigma (Estimated) noise level.
#' @param tmax Maximal iteration times.
#' @param sparse_mode A d-dimensional logicical vector that indicates whether each mode is sparse.
#' @param vartol Covergence telorance.
#' @return Estimated subspace in each mode.
STATSVD <- function(Y, r, sigma, tmax, sparse_mode, vartol){
  # Sparse tensor PCA algorithm for general order-d tensors

  try(if(missing("Y")) stop("missing argument: Y is required as the tensor type."))
  try(if(missing("r")) stop("missing argument: r is required as a scalar or vector."))
  try(if(class(Y) != "Tensor") stop("invalid input: Y should be of tensor type."))
  try(if(missing("tmax")) tmax = 10)
  try(if(missing("vartol")) vartol = 1e-6)
  p = dim(Y)
  d = length(p)

  try(if(missing("sparse_mode")) sparse_mode=rep(TRUE, d))

  p.prod = prod(p)
  p.minus.k = prod(p) / p
  if(is.atomic(r) && length(r)==1){
    r = rep(r, d)
  }
  r.minus.k = prod(r) / r # Introduce r_{-k}

  try(if(d != length(r)) stop("invalid input: r and the order of Y is incompatible."))
  tag.warning = FALSE
  U_t = list(); I_0 = list()
  for(i in 1:d){ # Initialization step 1: find significant index set I_k
    Y_i = k_unfold(Y, i)@data
    Y_i_row_norm = apply(Y_i, 1, vector_sum_square)
    #sparse mode
    if(isTRUE(sparse_mode[i])){
      I_0 = c(I_0, list((apply(Y_i, 1, vector_sum_square) > (sigma^2*(p.minus.k[i]+2*sqrt(p.minus.k[i]*log(p.prod))+2*log(p.prod))))
                        | (apply(abs(Y_i), 1, max) > 2*sigma*sqrt(log(p.prod)))))
    }
    else
      I_0 = c(I_0, list(rep(TRUE,p[i])))
    # select the significant indices
  }
  tilde_Y = Y
  for (i in 1:d){ # Initialization step 2: construct tilde_Y
    tilde_Y = ttm(tilde_Y, diag(I_0[[i]]*1), i)
  }
  for (i in 1:d){ # Initialization step 3: find loading U_0k
    if(isTRUE(sparse_mode[i])){ # some optimization
      Ui = matrix(0, nrow = p[i], ncol = r[i])
      datai = k_unfold(tilde_Y, i)@data
      selec = (1:p[i])[I_0[[i]]]
      if(length(selec)<r[i]){  # in case that sk < rk
        tag.warning = TRUE
        newdata = matrix(0, nrow=nrow(datai), ncol = ncol(datai))
        newdata[selec,] = datai[selec,]
        Ui = (svd(newdata)$u[,1:r[i]])
      }
      else{
        Ui[I_0[[i]],] = svd(datai[I_0[[i]],])$u[,1:r[i]]
      }
      U_t = c(U_t, list(t(Ui)))
    }
    else{
      U_t = c(U_t, list(t(svd(k_unfold(tilde_Y, i)@data)$u[,1:r[i]])))
    }

  }


  t = 1; approx = -1;
  svector = 0;
  while(t<tmax){ # Stop criterion: convergence or maximum number of iteration reached
    #print(t)
    for(i in 1:d){
      A = ttl(Y, U_t[-i], (1:d)[-i])
      A_matrix = k_unfold(A, i)@data
      if(!isTRUE(sparse_mode[i])){
        svd.result = svd(A_matrix)
        #print(svd.result)
        U_t[[i]] = t(svd.result$u[,1:r[i]])
        svector = svd.result$d[1:r[i]]
      }
      else{
        A_k_row_norm = apply(A_matrix, 1, vector_sum_square)
        I_k = A_k_row_norm > sigma^2 * (r.minus.k[i] + 2*(sqrt(r.minus.k[i]*log(p.prod))+log(p.prod)))
        selec = (1:p[i])[I_k]
        if(length(selec) < r[i]){
          #tag.warning = TRUE
          next
        }



        B_matrix = A_matrix[I_k,]

        hat.V = svd(B_matrix)$v[,1:r[i]]


        bar.A = A_matrix %*% hat.V
        bar.I_k = apply(bar.A, 1, vector_sum_square) > sigma^2 * (r[i] + 2*(sqrt(r[i]*log(p.prod))+log(p.prod)))
        This.U = matrix(0, nrow(bar.A), r[i])

        if(length((1:p[i])[bar.I_k]) < r[i]){
          #tag.warning = TRUE
          next
        }


        C = bar.A[bar.I_k,]
        svd.result = svd(C)
        This.U[bar.I_k,] = svd.result$u[,1:r[i]]
        U_t[[i]] = t(This.U)
        svector = svd.result$d[1:r[i]]
      }
    }
    if (abs(sum(svector^2) - approx) > vartol & t<tmax){
      t = t+1
      approx = sum(svector^2)
    }
    else {
      break
    }
  }

  for(i in 1:d){
    U_t[[i]] = t(U_t[[i]])
  }

  if(tag.warning) warning("Input noise level is too large and the estimation may be not reliable")

  return(U_t)
}


