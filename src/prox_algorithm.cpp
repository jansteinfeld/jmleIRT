#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// Safe logit transform
inline double safe_logit(double p)
{
  constexpr double lower = 1e-6;
  constexpr double upper = 1 - 1e-6;
  if (p < lower)
    return std::log(lower / (1.0 - lower));
  if (p > upper)
    return std::log(upper / (1.0 - upper));
  return std::log(p / (1.0 - p));
}

// Weighted mean ignoring NAs
double weighted_mean(const NumericVector &x, const NumericVector &w)
{
  double sum_wx = 0.0, sum_w = 0.0;
  int n = x.size();
  for (int i = 0; i < n; i++)
  {
    if (!NumericVector::is_na(x[i]) && w[i] > 0)
    {
      sum_wx += w[i] * x[i];
      sum_w += w[i];
    }
  }
  return (sum_w > 0) ? sum_wx / sum_w : NA_REAL;
}

// Weighted variance ignoring NAs
double weighted_var(const NumericVector &x, const NumericVector &w, double mean)
{
  double sum_w = 0.0, sum_w_diff2 = 0.0;
  int n = x.size();
  for (int i = 0; i < n; i++)
  {
    if (!NumericVector::is_na(x[i]) && w[i] > 0)
    {
      double diff = x[i] - mean;
      sum_w_diff2 += w[i] * diff * diff;
      sum_w += w[i];
    }
  }
  return (sum_w > 1) ? sum_w_diff2 / (sum_w - 1) : NA_REAL;
}

// [[Rcpp::export]]
List prox_rasch(NumericMatrix dat,
                    NumericMatrix dat_resp,
                    NumericVector freq,
                    double conv = 0.001,
                    int maxiter = 30)
{

  int IP = dat.nrow();
  int I = dat.ncol();
  if (dat_resp.size() == 0)
  {
    dat_resp = clone(dat);
    for (int i = 0; i < dat.nrow(); ++i)
    {
      for (int j = 0; j < dat.ncol(); ++j)
      {
        dat_resp(i, j) = !NumericVector::is_na(dat(i, j)) ? 1.0 : 0.0;
      }
    }
  }
  if (freq.size() == 0)
  {
    freq = NumericVector(dat.nrow(), 1.0);
  }
  

  NumericVector s_i(I, 0.0), n_i(I, 0.0);
  NumericVector r_n(IP, 0.0), n_n(IP, 0.0);
  NumericVector logit_freq_i(I);
  NumericVector logit_freq_n(IP);

  for (int j = 0; j < I; j++)
  {
    double sum_s = 0.0, sum_n = 0.0;
    for (int i = 0; i < IP; i++)
    {
      double resp = dat_resp(i, j);
      double val = !NumericVector::is_na(dat(i, j)) ? dat(i, j) : 0.0;
      sum_s += val * resp * freq[i];
      sum_n += resp * freq[i];
    }
    s_i[j] = sum_s;
    n_i[j] = sum_n;
    double prop = (sum_n > 0) ? (sum_s / sum_n) : 0.5;
    logit_freq_i[j] = safe_logit((prop + 0.01) / 1.02);
  }
  for (int i = 0; i < IP; i++)
  {
    double sum_r = 0.0, sum_n = 0.0;
    for (int j = 0; j < I; j++)
    {
      double resp = dat_resp(i, j);
      double val = !NumericVector::is_na(dat(i, j)) ? dat(i, j) : 0.0;
      sum_r += val * resp * freq[i];
      sum_n += resp * freq[i];
    }
    r_n[i] = sum_r;
    n_n[i] = sum_n;
    double prop = (sum_n > 0) ? (sum_r / sum_n) : 0.5;
    logit_freq_n[i] = safe_logit((prop + 0.01) / 1.02);
  }

  NumericVector d_i(I, 0.0), mu_i(I, 0.0), sigma_i(I, 1.0);
  NumericVector b_n(IP, 0.0), mu_n(IP, 0.0), sigma_n(IP, 1.0);

  NumericVector N_i = clone(n_i);
  NumericVector N_n(IP, 0.0);
  for (int i = 0; i < IP; i++)
  {
    double count = 0.0;
    for (int j = 0; j < I; j++)
    {
      count += dat_resp(i, j);
    }
    N_n[i] = count;
  }

  double par_change = 1.0;
  int iter = 0;

  while (par_change > conv && iter < maxiter)
  {
    NumericVector d_i_old = clone(d_i);

    // Update item difficulty
    for (int j = 0; j < I; j++)
    {
      d_i[j] = mu_i[j] - std::sqrt(1 + sigma_i[j] * sigma_i[j] / 2.9) * logit_freq_i[j];
    }

    // Update ability estimate for pattern
    for (int i = 0; i < IP; i++)
    {
      b_n[i] = mu_n[i] + std::sqrt(1 + sigma_n[i] * sigma_n[i] / 2.9) * logit_freq_n[i];
    }

    // Center persons
    double weighted_mean_b = 0.0, total_freq = 0.0;
    for (int i = 0; i < IP; i++)
    {
      weighted_mean_b += b_n[i] * freq[i];
      total_freq += freq[i];
    }
    weighted_mean_b /= total_freq;
    for (int i = 0; i < IP; i++)
    {
      b_n[i] -= weighted_mean_b;
    }

    // Update mu_i and sigma_i (mean and SD of logit abilities of persons who answered item i)
    for (int j = 0; j < I; j++)
    {
      NumericVector b_for_item(IP);
      NumericVector weights(IP);
      for (int i = 0; i < IP; i++)
      {
        b_for_item[i] = b_n[i];
        weights[i] = dat_resp(i, j) * freq[i];
      }
      double mean_b = weighted_mean(b_for_item, weights);
      mu_i[j] = mean_b;
      double var_b = weighted_var(b_for_item, weights, mean_b);
      sigma_i[j] = (var_b < 0) ? 0.0001 : std::sqrt(var_b);
    }

    // Update mu_n and sigma_n (mean and SD of logit difficulties encountered by person n)
    for (int i = 0; i < IP; i++)
    {
      NumericVector d_for_person(I);
      NumericVector weights(I);
      for (int j = 0; j < I; j++)
      {
        d_for_person[j] = d_i[j];
        weights[j] = dat_resp(i, j);
      }
      double mean_d = weighted_mean(d_for_person, weights);
      mu_n[i] = mean_d;
      double var_d = weighted_var(d_for_person, weights, mean_d);
      sigma_n[i] = (var_d < 0) ? 0.0001 : std::sqrt(var_d);
    }

    iter++;
    par_change = 0.0;
    for (int j = 0; j < I; j++)
    {
      double diff = std::abs(d_i[j] - d_i_old[j]);
      if (diff > par_change)
        par_change = diff;
    }
  }

  return List::create(
      Named("b") = d_i,
      Named("theta") = b_n,
      Named("iter") = iter,
      Named("sigma.i") = sigma_i,
      Named("sigma.n") = sigma_n);
}
