# STAT-SVD
R code for Sparse Tensor Alternating Thresholding Singular Value Decomposition (STAT-SVD). 

# Citation
@article{zhang2019optimal,
  title={Optimal Sparse Singular Value Decomposition for High-Dimensional High-Order Data},
  author={Zhang, Anru and Han, Rungang},
  journal={Journal of the American Statistical Association},
  volume={114},
  number={528},
  pages={1708--1725},
  year={2019},
  publisher={Taylor and Francis Ltd}
}

# Following packages are prerequisited: 'rTensor', 'ssvd', 'MASS'.

# Function Descriptions:

sine.theta(U1, U2); 
% U1, U2 are two subspaces with the same dimension, return the sine theta distance between U1 and U2.

get.sigma(T);
% T is a tensor-type variable, return the estimated standard deviation of noise under the regular Gaussian model.
% You are highly recommended to use empirical estimator when the noise is far from Gaussian.

get.rank(T);
% T is a tensor-type variable, return the estimated Tucker rank as an array.
% You are highly recommended to use empirical estimator when the noise is far from Gaussian.

HOSVD(T, r, sparse_mode);
% T is a noisy tensor, r is the Tucker rank array (can be automatically estimated when missing). sparse_mode is a logic array with the same length as r.
% Return the estimated singular space in each mode calculated by HOSVD.
% e.g.  >>>> U = HOSVD(T, c(3,4,5), sparse_mode = c(TRUE, TRUE, FALSE)); # U[[1]] is of size (p1*3); U[[1]] and U[[2]] would be sparse while U[[3]] is nonsparse.

