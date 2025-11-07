#include <Rcpp.h>
#include <cmath>
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
1) WLE estimation (function name unchanged)
- Preserves API; adds numeric guards and extreme stabilization.
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
  const double eps01 = 0.0;
  for (int i = 0; i < N; ++i)
  {
    for (int k = 0; k < K; ++k)
    {
      double v = X(i, k);
      if (NumericMatrix::is_na(v))
        continue;
      bool is_zero = (eps01 == 0.0) ? (v == 0.0) : (std::fabs(v - 0.0) <= eps01);
      bool is_one = (eps01 == 0.0) ? (v == 1.0) : (std::fabs(v - 1.0) <= eps01);
      if (!(is_zero || is_one))
        Rcpp::stop("X must contain only 0, 1, or NA (found %.8g at row %d, col %d).",
                   v, i + 1, k + 1);
    }
  }

  NumericVector raw_score(N), wle(N), se(N);
  IntegerVector conv(N), n_iter(N);

  // Track min/max beta as in original (compatibility)
  if (K > 0)
  {
    double min_b = beta[0];
    double max_b = beta[0];
    for (int k = 1; k < K; k++)
    {
      if (beta[k] < min_b)
        min_b = beta[k];
      if (beta[k] > max_b)
        max_b = beta[k];
    }
    (void)min_b;
    (void)max_b;
  }

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
    raw_score[i] = score;

    if (J == 0)
    {
      wle[i] = NA_REAL;
      se[i] = NA_REAL;
      conv[i] = 0;
      n_iter[i] = 0;
      continue;
    }

    bool is_extreme = (score == 0 || score == J);
    double score_d = static_cast<double>(score);
    double J_d = static_cast<double>(J);

    // Warm start from logit with 0.5 continuity correction
    double theta = std::log((score_d + 0.5) / (J_d - score_d + 0.5));
    double xstar_const = NA_REAL;
    if (is_extreme)
    {
      // target sum to avoid 0/1 fixation in iteration; bounded
      double a = std::min(std::max(wle_adj, 1e-8), J_d - 1e-8);
      double target_sum = (score == 0) ? a : (J_d - a);
      xstar_const = target_sum / J_d;
    }

    int converged = 0;
    double fi = 0.0, dfi = 0.0, delta = 0.0;
    int iter = 0;

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
        double z = theta - beta[k];
        double p = logistic(z);
        double var = p * (1.0 - p);
        double xeff = is_extreme ? xstar_const : x;
        fi += (xeff - p);
        dfi += var;
        wle_bias_sum += var * (1.0 - 2.0 * p);
      }

      if (std::abs(dfi) < 1e-12)
        break;

      double bias = 0.5 * wle_bias_sum / (-dfi);
      double fi_wle = fi - bias;
      delta = fi_wle / (-dfi);
      if (!R_finite(delta))
        break;
      if (std::abs(delta) > 5.0)
        delta = (delta > 0 ? 5.0 : -5.0);

      double theta_new = theta - delta;
      if (std::abs(delta) < tol)
      {
        theta = theta_new;
        converged = 1;
        break;
      }
      theta = theta_new;
    }

    double I_final = 0.0;
    if (converged)
    {
      for (int k = 0; k < K; k++)
      {
        double x = X(i, k);
        if (NumericMatrix::is_na(x))
          continue;
        double z = theta - beta[k];
        double p = logistic(z);
        I_final += p * (1.0 - p);
      }
    }

    wle[i] = theta;
    se[i] = converged ? std::sqrt(1.0 / std::max(I_final, 1e-12)) : NA_REAL;
    conv[i] = converged;
    n_iter[i] = converged ? (iter + 1) : max_iter;
  }

  return List::create(
      Named("raw_score") = raw_score,
      Named("wle") = wle,
      Named("standard_error") = se,
      Named("conv") = conv,
      Named("iterations") = n_iter);
}

/* -----------------------------------------------------------------------
2) JML estimation
- Robust JML with symmetric x*-use, no ±Inf mid-iteration, stable centering.
------------------------------------------------------------------------- */
Rcpp::List estimate_jmle(NumericMatrix X,
                         int max_iter = 1000,
                         double conv = 1e-6,
                         double eps = 0.0,
                         bool bias_correction = false,
                         std::string center = "items",
                         double max_update = 1.0,
                         bool verbose = false,
                         bool estimatewle = false,
                         double wle_adj = 1e-8)
{

  const int N = X.nrow();
  const int I = X.ncol();
  if (I == 0)
    stop("Number of items must be > 0.");

  // Sanity: only 0/1/NA
  const double eps01 = 0.0;
  for (int i = 0; i < N; ++i)
  {
    for (int k = 0; k < I; ++k)
    {
      double v = X(i, k);
      if (NumericMatrix::is_na(v))
        continue;
      bool is_zero = (eps01 == 0.0) ? (v == 0.0) : (std::fabs(v - 0.0) <= eps01);
      bool is_one = (eps01 == 0.0) ? (v == 1.0) : (std::fabs(v - 1.0) <= eps01);
      if (!(is_zero || is_one))
        Rcpp::stop("X must contain only 0, 1, or NA (found %.8g at row %d, col %d).",
                   v, i + 1, k + 1);
    }
  }

  // Row/col counts and sums
  IntegerVector row_obs(N, 0), col_obs(I, 0);
  NumericVector row_sum(N, 0.0), col_sum(I, 0.0);
  for (int p = 0; p < N; ++p)
  {
    for (int i = 0; i < I; ++i)
    {
      if (!NumericMatrix::is_na(X(p, i)))
      {
        row_obs[p]++;
        col_obs[i]++;
        if (X(p, i) != 0.0)
        {
          row_sum[p] += 1.0;
          col_sum[i] += 1.0;
        }
      }
    }
  }

  // Initialize theta, beta via logit with 0.5 continuity
  NumericVector theta(N, 0.0), beta(I, 0.0);
  for (int p = 0; p < N; ++p)
  {
    int n = std::max(1, row_obs[p]);
    double s = row_sum[p];
    double pr = (s + 0.5) / (n + 1.0);
    pr = std::min(std::max(pr, 1e-12), 1.0 - 1e-12);
    theta[p] = std::log(pr / (1.0 - pr));
  }
  for (int i = 0; i < I; ++i)
  {
    int m = std::max(1, col_obs[i]);
    double s = col_sum[i];
    double pr = (s + 0.5) / (m + 1.0);
    pr = std::min(std::max(pr, 1e-12), 1.0 - 1e-12);
    beta[i] = -std::log(pr / (1.0 - pr));
  }

  // Early centering to harmonize first steps
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
  auto center_persons = [&]()
  {
    double mean_theta = 0.0;
    int cnt = 0;
    for (int p = 0; p < N; ++p)
    {
      mean_theta += theta[p];
      cnt++;
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
  if (center == "items")
    center_items();
  else if (center == "persons")
    center_persons();

  int iter = 0;
  double max_diff = R_PosInf;
  const double denom_guard = 1e-6; // slightly larger for stability on tiny data

  // x* helper, symmetric for both sums; only affects extremes when eps>0
  auto xstar = [&](int p, int i) -> double
  {
    double x = X(p, i);
    if (NumericMatrix::is_na(x))
      return NA_REAL;
    int n = row_obs[p];
    double s = row_sum[p];
    if (eps > 0.0 && n > 0)
    {
      if (s <= 0.0)
        return eps / (double)n;
      if (s >= n)
        return 1.0 - (eps / (double)n);
    }
    return x;
  };

  // Mark extremes (persons) for eps==0 case; we won’t set ±Inf until the end
  LogicalVector is_extreme(N, false);

  while (iter < max_iter)
  {

    // Refresh extreme flags for eps==0 (informational; no mid-iteration ±Inf)
    for (int p = 0; p < N; ++p)
    {
      int n = row_obs[p];
      double s = row_sum[p];
      is_extreme[p] = (n > 0 && (s <= 0.0 || s >= n)) && (eps <= 0.0);
    }

    NumericVector theta_grad(N, 0.0), theta_info(N, 0.0);
    NumericVector beta_grad(I, 0.0), beta_info(I, 0.0);

    // Person side
    // Compute person gradients and information using x* when eps>0
    for (int p = 0; p < N; ++p)
    {
      if (is_extreme[p] && eps <= 0.0)
      {
        // fully inert: no contribution and no update
        theta_grad[p] = 0.0;
        theta_info[p] = 1.0;
        continue;
      }
      double sumP = 0.0, infoP = 0.0, sumX = 0.0;
      for (int i = 0; i < I; ++i)
      {
        double x = X(p, i);
        if (NumericMatrix::is_na(x))
          continue;
        // if person p is extreme and eps==0, we already continued; else:
        double P = logistic(theta[p] - beta[i]);
        double xs = xstar(p, i);
        sumX += xs;
        sumP += P;
        infoP += P * (1.0 - P);
      }
      theta_grad[p] = sumX - sumP;
      theta_info[p] = std::max(infoP, denom_guard);
    }
    // Item side (consistent x*)
    // Compute item gradients and information consistently using x*,
    // and drop extreme persons when eps==0 (pure JML deletion)
    for (int i = 0; i < I; ++i)
    {
      double sumP = 0.0, infoP = 0.0, sumX = 0.0;
      for (int p = 0; p < N; ++p)
      {
        if (is_extreme[p] && eps <= 0.0)
          continue; // deletion for pure JML
        double x = X(p, i);
        if (NumericMatrix::is_na(x))
          continue;
        double P = logistic(theta[p] - beta[i]);
        double xs = xstar(p, i);
        sumX += xs;
        sumP += P;
        infoP += P * (1.0 - P);
      }
      beta_grad[i] = sumP - sumX;
      beta_info[i] = std::max(infoP, denom_guard);
    }

    // Updates with clipping and smooth damping
    max_diff = 0.0;
    auto damp = [](double step)
    {
      // smooth damping: step *= 1/(1+|step|/0.5)
      double c = 0.5;
      double f = 1.0 / (1.0 + std::fabs(step) / c);
      return step * f;
    };

    // Separate caps for persons/items derived from max_update
    double cap_theta = max_update;
    double cap_beta = 0.5 * max_update;

    for (int p = 0; p < N; ++p)
    {
      if (is_extreme[p] && eps <= 0.0)
        continue;
      double step = theta_grad[p] / theta_info[p];
      if (step > cap_theta)
        step = cap_theta;
      if (step < -cap_theta)
        step = -cap_theta;
      step = damp(step);
      theta[p] += step;
      if (std::fabs(step) > max_diff)
        max_diff = std::fabs(step);
    }

    for (int i = 0; i < I; ++i)
    {
      double step = beta_grad[i] / beta_info[i];
      if (step > cap_beta)
        step = cap_beta;
      if (step < -cap_beta)
        step = -cap_beta;
      step = damp(step);
      beta[i] += step;
      if (std::fabs(step) > max_diff)
        max_diff = std::fabs(step);
    }

    /// Centering on each iteration
    if (center == "items")
    {
      double mean_beta = 0.0;
      for (int i = 0; i < I; ++i)
        mean_beta += beta[i];
      mean_beta /= (double)I;
      for (int i = 0; i < I; ++i)
        beta[i] -= mean_beta;
      // shift ALL persons equally (including extremes, which remain finite here)
      for (int p = 0; p < N; ++p)
        theta[p] += mean_beta;
    }
    else if (center == "persons")
    {
      double mean_theta = 0.0;
      for (int p = 0; p < N; ++p)
        mean_theta += theta[p];
      mean_theta /= (double)N;
      for (int p = 0; p < N; ++p)
        theta[p] -= mean_theta;
      for (int i = 0; i < I; ++i)
        beta[i] += mean_theta;
    }
    if (verbose && (iter % 10 == 0))
      Rprintf("Iteration %d, max parameter change = %g\n", iter, max_diff);

    iter++;
    if (max_diff <= conv)
      break;
  }

  // Optional post-hoc bias correction (global scaling), then re-center
  if (bias_correction && I > 0)
  {
    double factor = (I - 1.0) / (double)I;
    for (int i = 0; i < I; ++i)
      beta[i] *= factor;
    if (center == "items")
    {
      double mean_beta = 0.0;
      for (int i = 0; i < I; ++i)
        mean_beta += beta[i];
      mean_beta /= (double)I;
      for (int i = 0; i < I; ++i)
        beta[i] -= mean_beta;
      for (int p = 0; p < N; ++p)
        theta[p] += mean_beta;
    }
    else if (center == "persons")
    {
      double mean_theta = 0.0;
      for (int p = 0; p < N; ++p)
        mean_theta += theta[p];
      mean_theta /= (double)N;
      for (int p = 0; p < N; ++p)
        theta[p] -= mean_theta;
      for (int i = 0; i < I; ++i)
        beta[i] += mean_theta;
    }
  }

  // Optional WLE using final beta
  NumericVector theta_wle(N, NA_REAL);
  if (estimatewle)
  {
    List w = estimate_wle(X, beta, max_iter, conv, NA_REAL, NA_REAL, wle_adj);
    theta_wle = w["wle"];
  }

// Assign ±Inf for extremes only at the end and only if eps == 0
  if (eps <= 0.0)
  {
    for (int p = 0; p < N; ++p)
    {
      if ((row_obs[p] > 0) && (row_sum[p] <= 0.0))
        theta[p] = R_NegInf;
      else if ((row_obs[p] > 0) && (row_sum[p] >= row_obs[p]))
        theta[p] = R_PosInf;
    }
  }
  return List::create(
      _["data"] = X,
      _["theta"] = theta,
      _["beta"] = beta,
      _["iterations"] = iter,
      _["converged"] = (max_diff <= conv),
      _["bias_correction"] = bias_correction,
      _["center"] = center,
      _["wle_estimate"] = theta_wle);
}

/* -----------------------------------------------------------------------
Exported wrappers (names/signatures unchanged)
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
