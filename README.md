# tmvn

This repository houses the `R` package `tmvn`, which provides efficient methods
for sampling from a truncated multivariate normal and Student t distribution 
subject to convex polytope restrictions. The method implemented in this package 
leverages a rejection sampling method for the truncated univariate normal 
distribution which has superior acceptance rates compared to popular existing methods. 

## Installation
1. Install and load the `devtools` package.
2. Install `tmvn` from GitHub using `install_github("suchitm/tmvn")` 

## References
Li, Y., & Ghosh, S. K. (2015). [Efficient sampling methods for 
  truncated multivariate normal and student-t distributions subject to 
  linear inequality constraints](http://www.stat.ncsu.edu/information/library/papers/mimeo2649_Li.pdf). Journal of Statistical Theory and 
  Practice, 9(4), 712-732.

