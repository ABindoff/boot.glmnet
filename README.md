README
================
Bindoff, A. D.

2021-03-25

## Background

Penalized regression methods are used to select variables for prediction
models (or diagnostic models). In some fields the data may be very
high-dimensional, resulting in very large model matrixes. Frequently,
users must fit many boostrap replications to obtain useful statistics
such as non-parametric confidence intervals. This can quickly exceed
memory and computation constraints, even using sparse matrix
representations.

The `boot.glmnet` R package aims to extend the popular `glmnet` R
package \[1\] by providing functions to fit bootstrap replications in
blocks which can be computed sequentially to be saved and combined, or
computed in parallel using the `mc_boot_glmnet` function which takes
advantage of the `future.apply` R package \[2\] to distribute computing
in all supported environments (including Windows).

Additionally, `boot.glmnet` will calculate the proportion of replicates
for which a predictor term is selected via LASSO (or elastic-net), and
return an approximate confidence interval around this.

Methods are provided for conveniently summarising results above
significance thresholds and plotting results.

## Example

TBA

### References

\[1\] Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
Regularization Paths for Generalized Linear Models via Coordinate
Descent. Journal of Statistical Software, 33(1), 1-22. URL
<https://www.jstatsoft.org/v33/i01/>.

\[2\] H. Bengtsson, A Unifying Framework for Parallel and Distributed
Processing in R using Futures, arXiv:2008.00553, 2020
