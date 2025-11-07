#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Logistic sigmoid function used in the Rasch model
inline double logistic(double x)
{
  if (x < 0)
  {
    double z = std::exp(x);
    return z / (1.0 + z);
  }
  else
  {
    double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
}

/*
 * Analytical First-Order Bias Correction for JML Estimators in the Rasch Model.
 *
 * The function corrects the first-order incidental parameter bias from the person parameters (theta)
 * and item parameters (beta) estimated by Joint Maximum Likelihood (JML).
 *
 * The correction is based on the score and information terms of the Rasch model log-likelihood:
 * For each person n and item i, define:
 *   p_ni = logistic(theta_n - beta_i)
 *   u_ni = X_ni - p_ni            (score function component)
 *   v_ni = p_ni * (1 - p_ni)      (information component, variance of Bernoulli)
 *
 * Then the bias terms B_theta and B_beta are approximated as ratios of sums of v_ni and u_ni terms:
 * B_theta[n] = sum_i v_ni * u_ni / sum_i v_ni^2
 * B_beta[i]  = sum_n v_ni * (p_ni - X_ni) / sum_n v_ni^2
 *
 * Finally, corrected parameter estimates are:
 * theta_corrected = theta - B_theta / I
 * beta_corrected = beta - B_beta / I
 * where I is the number of items.
 *
 * Missing responses (NA) are skipped.
 *
 * @param theta NumericVector, person parameter estimates to correct (modified in place)
 * @param beta NumericVector, item parameter estimates to correct (modified in place)
 * @param X NumericMatrix, response data (persons x items) coded 0/1 with possible NAs
 * @param N Integer, number of persons (rows in X)
 * @param I Integer, number of items (columns in X)
 */

// [[Rcpp::export]]
Rcpp::List biasCorrectionJMLE(NumericVector &theta, NumericVector &beta, const NumericMatrix &X,
                              int N, int I)
{

  NumericVector biasTheta(N);
  NumericVector biasBeta(I);

  // Compute bias term for each person parameter theta_n
  for (int n = 0; n < N; ++n)
  {
    double numerator = 0.0;
    double denominator = 0.0;
    for (int i = 0; i < I; ++i)
    {
      if (NumericVector::is_na(X(n, i)))
        continue; // skip missing

      double eta = theta[n] - beta[i];
      double p = logistic(eta);

      double u = X(n, i) - p;   // score component
      double v = p * (1.0 - p); // information component

      numerator += v * u;
      denominator += v * v;
    }
    biasTheta[n] = (denominator > 0.0) ? numerator / denominator : 0.0;
  }

  // Compute bias term for each item parameter beta_i
  for (int i = 0; i < I; ++i)
  {
    double numerator = 0.0;
    double denominator = 0.0;
    for (int n = 0; n < N; ++n)
    {
      if (NumericVector::is_na(X(n, i)))
        continue;

      double eta = theta[n] - beta[i];
      double p = logistic(eta);

      double u = p - X(n, i); // score component flipped for items
      double v = p * (1.0 - p);

      numerator += v * u;
      denominator += v * v;
    }
    biasBeta[i] = (denominator > 0.0) ? numerator / denominator : 0.0;
  }

  // Apply bias correction scaled by 1/I (number of items)
  for (int n = 0; n < N; ++n)
  {
    theta[n] -= biasTheta[n] / static_cast<double>(I);
  }
  for (int i = 0; i < I; ++i)
  {
    beta[i] -= biasBeta[i] / static_cast<double>(I);
  }

  return List::create(
      Named("theta") = theta,
      Named("beta") = beta);
}
