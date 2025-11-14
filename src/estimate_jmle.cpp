#include <Rcpp.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;

/* -----------------------------------------------------------------------
Helper: numerically stable logistic (sigmoid)
------------------------------------------------------------------------- */
inline double logistic(double x)
{
  if (x >= 0)
  {
    double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  else
  {
    double z = std::exp(x);
    return z / (1.0 + z);
  }
}

/* -----------------------------------------------------------------------
Helper: validate X is a (logical/integer/double) matrix, coerce to NumericMatrix
------------------------------------------------------------------------- */
inline Rcpp::NumericMatrix require_real_matrix(SEXP X_)
{
  if (!Rf_isMatrix(X_))
    Rcpp::stop("X must be a matrix.");
  int t = TYPEOF(X_);
  if (!(t == REALSXP || t == INTSXP || t == LGLSXP))
    Rcpp::stop("X must be a numeric matrix.");
  Rcpp::NumericMatrix X = Rcpp::as<Rcpp::NumericMatrix>(X_);
  return X;
}

/* -----------------------------------------------------------------------
1) WLE estimation (unchanged from original)
------------------------------------------------------------------------- */
Rcpp::List estimate_wle(NumericMatrix X, NumericVector beta,
                        int max_iter = 1000, double tol = 1e-10,
                        double lower_ext = NA_REAL, double upper_ext = NA_REAL,
                        double wle_adj = 1e-8)
{
  const int N = X.nrow();
  const int K = X.ncol();
  if (K == 0)
    stop("Number of items must be > 0.");
  if (beta.size() != K)
    stop("Length of 'beta' must equal number of columns in 'X'.");

  // Sanity: only 0/1/NA
  for (int i = 0; i < N; ++i)
  {
    for (int k = 0; k < K; ++k)
    {
      double v = X(i, k);
      if (NumericMatrix::is_na(v))
        continue;
      if (!(v == 0.0 || v == 1.0))
        Rcpp::stop("X must contain only 0, 1, or NA (found %.8g at row %d, col %d).",
                   v, i + 1, k + 1);
    }
  }

  NumericVector raw_score(N), wle(N), se(N);
  IntegerVector conv(N), n_iter(N);

  // Main loop over persons calculating raw scores and WLE
  for (int i = 0; i < N; i++)
  {
    int score = 0, J = 0;
    for (int k = 0; k < K; k++)
    {
      double x = X(i, k);
      if (!NumericMatrix::is_na(x))
      {
        score += (int)x;
        J += 1;
      }
    }
    raw_score[i] = score; // store raw score

    if (J == 0) // no valid items
    {
      wle[i] = NA_REAL;
      se[i] = NA_REAL;
      conv[i] = 0;
      n_iter[i] = 0;
      continue;
    }

    // initial theta estimate
    bool is_extreme = (score == 0 || score == J);
    double score_d = static_cast<double>(score);                      // observed score
    double J_d = static_cast<double>(J);                              // number of valid items
    double theta = std::log((score_d + 0.5) / (J_d - score_d + 0.5)); // initial theta
    double xstar_const = NA_REAL;                                     // for extreme scores

    // handle extreme scores
    if (is_extreme)
    {
      double a = std::min(std::max(wle_adj, 1e-8), J_d - 1e-8); // guard for extremes
      double target_sum = (score == 0) ? a : (J_d - a);         // adjusted target sum
      xstar_const = target_sum / J_d;                           // effective score for all items
    }

    int converged = 0;
    double fi = 0.0, dfi = 0.0, delta = 0.0;
    int iter = 0;

    // Newton-Raphson iteration
    for (iter = 0; iter < max_iter; iter++)
    {
      fi = 0.0;
      dfi = 0.0;
      double wle_bias_sum = 0.0;

      for (int k = 0; k < K; k++)
      {
        double x = X(i, k);
        if (NumericMatrix::is_na(x))
          continue;
        double z = theta - beta[k];                 // theta - b_k
        double p = logistic(z);                     // P(X=1|theta)
        double var = p * (1.0 - p);                 // Variance
        double xeff = is_extreme ? xstar_const : x; // effective score
        fi += (xeff - p);                           // Score function
        dfi += var;                                 // Fisher information
        wle_bias_sum += var * (1.0 - 2.0 * p);      // correction term for WLE
      }

      // Check for zero information
      if (std::abs(dfi) < 1e-12)
        break;

      // WLE bias correction
      double bias = 0.5 * wle_bias_sum / (-dfi); // bias term
      double fi_wle = fi - bias;                 // adjusted score function
      delta = fi_wle / (-dfi);                   // Newton-Raphson update
      if (!R_finite(delta))                      // guard against numerical issues
        break;
      if (std::abs(delta) > 5.0)          // limit step size
        delta = (delta > 0 ? 5.0 : -5.0); // limit step size

      double theta_new = theta - delta; // update theta

      // Check convergence
      if (std::abs(delta) < tol)
      {
        theta = theta_new;
        converged = 1;
        break;
      }
      theta = theta_new;
    }
    // Compute standard error if converged
    double I_final = 0.0;
    if (converged)
    {
      for (int k = 0; k < K; k++)
      {
        double x = X(i, k);
        if (NumericMatrix::is_na(x))
          continue;
        double z = theta - beta[k]; // theta - b_k
        double p = logistic(z);     // P(X=1|theta)
        I_final += p * (1.0 - p);   // Fisher information
      }
    }

    wle[i] = theta;
    se[i] = converged ? std::sqrt(1.0 / std::max(I_final, 1e-12)) : NA_REAL;
    conv[i] = converged;
    n_iter[i] = converged ? (iter + 1) : max_iter;
  }
  // Return results
  return List::create(
      Named("raw_score") = raw_score,
      Named("wle") = wle,
      Named("standard_error") = se,
      Named("conv") = conv,
      Named("iterations") = n_iter);
}

/* -----------------------------------------------------------------------
2) OPTIMIZED JML estimation - all fixes incorporated
------------------------------------------------------------------------- */
Rcpp::List estimate_jmle(NumericMatrix X,
                         int max_iter = 1000,
                         double conv = 1e-6,
                         double eps = 0.0,
                         bool bias_correction = false,
                         std::string center = "items",
                         double max_update = 1.5,
                         bool verbose = false,
                         bool estimatewle = false,
                         double wle_adj = 1e-8)
{
  const int N = X.nrow();
  const int I = X.ncol();
  if (I == 0)
    stop("Number of items must be > 0.");

  // Sanity check: only 0/1/NA
  for (int i = 0; i < N; ++i)
  {
    for (int k = 0; k < I; ++k)
    {
      double v = X(i, k);
      if (NumericMatrix::is_na(v))
        continue;
      if (!(v == 0.0 || v == 1.0))
        Rcpp::stop("X must contain only 0, 1, or NA (found %.8g at row %d, col %d).",
                   v, i + 1, k + 1);
    }
  }

  // SPEED: Pre-compute row/col statistics once
  IntegerVector row_obs(N, 0), col_obs(I, 0);
  NumericVector row_sum(N, 0.0), col_sum(I, 0.0);
  LogicalVector is_extreme(N, false);

  // Compute observed counts and sums
  for (int p = 0; p < N; ++p)
  {
    for (int i = 0; i < I; ++i)
    {
      if (!NumericMatrix::is_na(X(p, i)))
      {
        row_obs[p]++;
        col_obs[i]++;
        double val = X(p, i);
        row_sum[p] += val;
        col_sum[i] += val;
      }
    }

    // Mark extreme persons (for deletion in standard JML when eps=0)
    if (row_obs[p] > 0 && eps <= 0.0)
    {
      is_extreme[p] = (row_sum[p] <= 0.0 || row_sum[p] >= row_obs[p]);
    }
  }

  // Initialize theta and beta with logit of proportions
  NumericVector theta(N, 0.0), beta(I, 0.0);

  // Initial person parameters
  for (int p = 0; p < N; ++p)
  {
    if (row_obs[p] > 0)
    {
      double prop = (row_sum[p] + 0.5) / (row_obs[p] + 1.0);
      prop = std::min(std::max(prop, 1e-12), 1.0 - 1e-12);
      theta[p] = std::log(prop / (1.0 - prop));
    }
  }
  // Initial item parameters
  for (int i = 0; i < I; ++i)
  {
    if (col_obs[i] > 0)
    {
      double prop = (col_sum[i] + 0.5) / (col_obs[i] + 1.0);
      prop = std::min(std::max(prop, 1e-12), 1.0 - 1e-12);
      beta[i] = -std::log(prop / (1.0 - prop)); // Note: negative for difficulty
    }
  }
  // Centering functions if items are required
  auto center_items = [&]()
  {
    double mean_beta = 0.0;
    for (int i = 0; i < I; ++i)
      mean_beta += beta[i];
    mean_beta /= (double)I;
    for (int i = 0; i < I; ++i)
      beta[i] -= mean_beta;
    for (int p = 0; p < N; ++p)
      theta[p] += mean_beta;
  };
  // Centering functions if persons are required
  auto center_persons = [&]()
  {
    double mean_theta = 0.0;
    int cnt = 0;
    for (int p = 0; p < N; ++p)
    {
      if (!is_extreme[p])
      {
        mean_theta += theta[p];
        cnt++;
      }
    }
    if (cnt > 0)
    {
      mean_theta /= (double)cnt;
      for (int p = 0; p < N; ++p)
        theta[p] -= mean_theta;
      for (int i = 0; i < I; ++i)
        beta[i] += mean_theta;
    }
  };

  // Main iteration loop
  int iter = 0;
  String converged = "FALSE";
  const double denom_guard = 1e-12; // to avoid division by zero

  // SPEED: Allocate once, reuse
  NumericVector theta_grad(N), theta_info(N);
  NumericVector beta_grad(I), beta_info(I);

  // Initialize with centered values BEFORE the loop
  if (center == "items")
    center_items();
  else if (center == "persons")
    center_persons();

  // Initialize old values ONCE
  NumericVector theta_old = clone(theta);
  NumericVector beta_old = clone(beta);

  // Iterative updates
  while (iter < max_iter)
  {
    // ========== STEP 1: Update PERSONS (items fixed) ==========
    std::fill(theta_grad.begin(), theta_grad.end(), 0.0);
    std::fill(theta_info.begin(), theta_info.end(), 0.0);

    // Compute gradients and information for persons
    for (int p = 0; p < N; ++p)
    {
      if (is_extreme[p])
        continue;

      for (int i = 0; i < I; ++i)
      {
        double x = X(p, i);
        if (NumericMatrix::is_na(x))
          continue;

        double z = theta[p] - beta[i]; // theta - b_i
        double P = logistic(z);        // P(X=1|theta,beta)
        double W = P * (1.0 - P);      // Variance

        theta_grad[p] += x - P; // Score function
        theta_info[p] += W;     // Fisher information
      }
    }

    // Apply person parameter updates
    for (int p = 0; p < N; ++p)
    {
      // Skip extreme persons
      if (is_extreme[p])
        continue;
      // Guard against zero information
      if (theta_info[p] < 1e-12)
        theta_info[p] = 1.0;
      // Compute step size with max_update limit
      double step = theta_grad[p] / theta_info[p];
      if (step > max_update)
        step = max_update;
      else if (step < -max_update)
        step = -max_update;

      theta[p] += step;
    }

    // ========== STEP 2: Update ITEMS (persons fixed) ==========
    std::fill(beta_grad.begin(), beta_grad.end(), 0.0);
    std::fill(beta_info.begin(), beta_info.end(), 0.0);

    // Compute gradients and information for items
    for (int i = 0; i < I; ++i)
    {
      for (int p = 0; p < N; ++p)
      {
        // Skip extreme persons
        if (is_extreme[p])
          continue;

        double x = X(p, i);
        if (NumericMatrix::is_na(x))
          continue;

        double z = theta[p] - beta[i]; // theta - b_i
        double P = logistic(z);        // P(X=1|theta,beta)
        double W = P * (1.0 - P);      // Variance

        beta_grad[i] += P - x; // Score function
        beta_info[i] += W;     // Fisher information
      }
    }

    // Apply item parameter updates
    for (int i = 0; i < I; ++i)
    {
      // Guard against zero information
      if (beta_info[i] < 1e-12)
        beta_info[i] = 1.0;
      // Compute step size with max_update limit
      double step = beta_grad[i] / beta_info[i];
      if (step > max_update)
        step = max_update;
      else if (step < -max_update)
        step = -max_update;

      beta[i] += step;
    }

    // ========== STEP 3: Center after both updates ==========
    if (center == "items")
      center_items();
    else if (center == "persons")
      center_persons();

    // ========== STEP 4: Check convergence ==========
    double max_change_items = 0.0;
    for (int i = 0; i < I; ++i)
    {
      // Compute maximum change for items
      double change = std::fabs(beta[i] - beta_old[i]);
      if (change > max_change_items)
        max_change_items = change;
    }

    // Compute maximum change for persons
    double max_change_persons = 0.0;
    for (int p = 0; p < N; ++p)
    {
      if (is_extreme[p])
        continue;
      double change = std::fabs(theta[p] - theta_old[p]);
      if (change > max_change_persons)
        max_change_persons = change;
    }

    // Overall maximum change for convergence
    double max_change = std::max(max_change_items, max_change_persons);

    if (verbose && (iter < 5 || iter % 10 == 0))
    {
      Rprintf("Iteration %4d: max_items=%.8f, max_persons=%.8f\n",
              iter, max_change_items, max_change_persons);
      R_FlushConsole();
    }

    if (max_change < conv && iter >= 5)
    {
      converged = "TRUE";
      if (verbose)
        Rprintf("Converged: max_items=%.8f, max_persons=%.8f\n",
                max_change_items, max_change_persons);
      break;
    }

    theta_old = clone(theta);
    beta_old = clone(beta);

    iter++;
  }

  // Optional bias correction (post-hoc scaling)
  if (bias_correction && I > 1)
  {
    double factor = (double)(I - 1) / (double)I;
    for (int i = 0; i < I; ++i)
      beta[i] *= factor;

    // Re-center after bias correction
    if (center == "items")
      center_items();
    else if (center == "persons")
      center_persons();
  }

  // Optional WLE estimation using final beta
  NumericVector theta_wle(N, NA_REAL);
  if (estimatewle)
  {
    List w = estimate_wle(X, beta, max_iter, 1e-10, NA_REAL, NA_REAL, wle_adj);
    theta_wle = w["wle"];
  }

  // Assign ±Inf for extreme persons (only if eps == 0)
  if (eps <= 0.0)
  {
    for (int p = 0; p < N; ++p)
    {
      if (row_obs[p] > 0)
      {
        if (row_sum[p] <= 0.0)
          theta[p] = R_NegInf;
        else if (row_sum[p] >= row_obs[p])
          theta[p] = R_PosInf;
      }
    }
  }

  return List::create(
      _["data"] = X,
      _["theta"] = theta,
      _["beta"] = beta,
      _["iterations"] = iter,
      _["converged"] = converged,
      _["bias_correction"] = bias_correction,
      _["center"] = center,
      _["wle_estimate"] = theta_wle);
}

/* -----------------------------------------------------------------------
Exported wrappers (unchanged)
------------------------------------------------------------------------- */
// [[Rcpp::export]]
Rcpp::List estimate_wle(SEXP X_, Rcpp::NumericVector beta,
                        int max_iter = 100, double tol = 1e-10,
                        double lower_ext = NA_REAL, double upper_ext = NA_REAL,
                        double wle_adj = 1e-8)
{
  try
  {
    Rcpp::NumericMatrix X = require_real_matrix(X_);
    return estimate_wle(X, beta, max_iter, tol, lower_ext, upper_ext, wle_adj);
  }
  catch (std::exception &e)
  {
    Rcpp::stop(e.what());
  }
  catch (...)
  {
    Rcpp::stop("estimate_wle: unknown error.");
  }
}

// [[Rcpp::export]]
Rcpp::List estimate_jmle(SEXP X_,
                         int max_iter = 1000,
                         double conv = 1e-6,
                         double eps = 0.0,
                         bool bias_correction = false,
                         std::string center = "items",
                         double max_update = 1.5,
                         bool verbose = false,
                         bool estimatewle = false,
                         double wle_adj = 1e-8)
{
  try
  {
    Rcpp::NumericMatrix X = require_real_matrix(X_);
    return estimate_jmle(X, max_iter, conv, eps, bias_correction,
                         center, max_update, verbose, estimatewle, wle_adj);
  }
  catch (std::exception &e)
  {
    Rcpp::stop(e.what());
  }
  catch (...)
  {
    Rcpp::stop("estimate_jmle: unknown error.");
  }
}
