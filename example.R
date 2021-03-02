<<<<<<< HEAD
library(STATSVD)
=======
>>>>>>> a215e08dc581744950ca249aee3bf2e4922f578f
set.seed(100)
# Parameter setting
p = c(50, 50, 50) # dimension.
r = c(5, 5, 5) # rank.
s = c(10, 10, 50) # support cardinality. Here the 3rd mode is dense by construction.
lambda = 60 # Least singular value of signal tensor.
sigma = 0.5 # Noise level.

# Generate a low-rank tensor.
S.array = array(rnorm(prod(r)), dim=r)
S = as.tensor(S.array)
d1 = svd(k_unfold(S,1)@data)$d
d2 = svd(k_unfold(S,2)@data)$d
d3 = svd(k_unfold(S,3)@data)$d
lambda0 = min(c(d1[length(d1)], d2[length(d2)], d3[length(d3)]))
S = S*lambda/lambda0 # Core tensor

U1 = matrix(0, p[1], r[1]); U2 = matrix(0, p[2], r[2]); U3 = matrix(0, p[3], r[3])
U1[sample(p[1], s[1]),] = svd(matrix(rnorm(s[1]*r[1]), nrow=s[1], ncol=r[1]))$u
U2[sample(p[2], s[2]),] = svd(matrix(rnorm(s[2]*r[2]), nrow=s[2], ncol=r[2]))$u
U3[sample(p[3], s[3]),] = svd(matrix(rnorm(s[3]*r[3]), nrow=s[3], ncol=r[3]))$u # loadings
X = ttm(ttm(ttm(S, U1, 1), U2, 2), U3, 3)
Z = sigma * as.tensor(array(rnorm(prod(p)), dim=p))
Y = X + Z # observation

# estimate sigma.
hat.sigma = get.sigma(Y);

# estimate rank.
hat.rank = get.rank(Y, hat.sigma, c(TRUE, TRUE, FALSE))

# STATSVD.
U.hat = STATSVD(Y, hat.rank, hat.sigma, sparse_mode = c(TRUE, TRUE, FALSE))
P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
X.hat = ttm(ttm(ttm(Y, P.U1, 1), P.U2, 2), P.U3, 3)
cat(sprintf("Estimation error from STATSVD is %f", fnorm(X.hat-X)/fnorm(X)))

# Compare with High-order Orthogonal Iteration (HOOI)
U.hat = HOOI(Y, hat.rank)
P.U1 = U.hat[[1]]%*%t(U.hat[[1]])
P.U2 = U.hat[[2]]%*%t(U.hat[[2]])
P.U3 = U.hat[[3]]%*%t(U.hat[[3]])
X.hat = ttm(ttm(ttm(Y, P.U1, 1), P.U2, 2), P.U3, 3)
cat(sprintf("Estimation error from HOOI is %f", fnorm(X.hat-X)/fnorm(X)))


