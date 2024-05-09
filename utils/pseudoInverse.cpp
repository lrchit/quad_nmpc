#include "pseudoInverse.h"
/*!
 * Compute the pseudo inverse of a matrix
 * @param matrix : input matrix
 * @param sigmaThreshold : threshold for singular values being zero
 * @param invMatrix : output matrix
 */

// inline 为了防止mutiple definition
void pseudoInverse(DMat const& matrix, double sigmaThreshold,
                   DMat& invMatrix) {
  if ((1 == matrix.rows()) && (1 == matrix.cols())) {
    invMatrix.resize(1, 1);
    if (matrix.coeff(0, 0) > sigmaThreshold) {
      invMatrix.coeffRef(0, 0) = 1.0 / matrix.coeff(0, 0);
    } else {
      invMatrix.coeffRef(0, 0) = 0.0;
    }
    return;
  }

  Eigen::JacobiSVD<DMat> svd(matrix,
                                Eigen::ComputeThinU | Eigen::ComputeThinV);
  // not sure if we need to svd.sort()... probably not
  int const nrows(svd.singularValues().rows());
  DMat invS;
  invS = DMat::Zero(nrows, nrows);
  for (int ii(0); ii < nrows; ++ii) {
    if (svd.singularValues().coeff(ii) > sigmaThreshold) {
      invS.coeffRef(ii, ii) = 1.0 / svd.singularValues().coeff(ii);
    } else {
      // invS.coeffRef(ii, ii) = 1.0/ sigmaThreshold;
      // printf("sigular value is too small: %f\n",
      // svd.singularValues().coeff(ii));
    }
  }
  invMatrix = svd.matrixV() * invS * svd.matrixU().transpose();
}