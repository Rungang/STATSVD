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

# Instructions
Algorithm code and auxiliary functions are all contained in STATSVD.R.  

STATSVD requires the following packages for full functionality: 'MASS', 'rTensor', 'ssvd'. 

The algorithm requires the tensor rank, sparsity mode indexes and noise level as prior. It returns the estimated orthogonal loading (or subspace) for each mode. An estimation of the tensor can be then obtained via projection. See example.R for illustration.

Adaptive estimations for rank and noise level are provided, but may only be valid under Gaussian additive model (i.e. Y = X + Z), see the paper for more details. For Non-Gaussian data, you are highly recommended to provide an empirical estimation for these two tuning parameters. Over-estimating noise level may lead to some warning and inaccurate estimation.

# Auxiliary Functions
```R
sine.theta(U1, U2); % Calculate the sine theta distance between subspaces U1 and U2.
sigma.hat = get.sigma(Y); % Estimate standard deviation of noise under the regular Gaussian additive model given observed tensor Y.
r.hat = get.rank(Y); % Estimate Tucker rank given observed tensor Y.
```

Please contact rhan32@stat.wisc.edu with any questions or comments.

