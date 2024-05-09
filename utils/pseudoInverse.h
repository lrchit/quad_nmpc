#ifndef PSEUDO_INVERSE_H
#define PSEUDO_INVERSE_H

#include <stdio.h>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/SVD>
#include "cppTypes.h"
using namespace std;

/*!
 * Compute the pseudo inverse of a matrix
 * @param matrix : input matrix
 * @param sigmaThreshold : threshold for singular values being zero
 * @param invMatrix : output matrix
 */

void pseudoInverse(DMat const& matrix, double sigmaThreshold, DMat& invMatrix);
#endif