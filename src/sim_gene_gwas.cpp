#include <RcppEnsmallen.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::depends(RcppEnsmallen)]]

using namespace Rcpp;

// Define a differentiable objective function by implementing both Evaluate()
// and Gradient() separately.
class LogisticRegressionFunction
{
public:
  // Construct the object with the given the design
  // matrix and responses.
  LogisticRegressionFunction(const arma::mat& X,
                             const arma::vec& y) :
  X(X), y(y) { }

  // Return the objective function for model parameters beta.
  double Evaluate(const arma::mat& beta)
  {

    // Negative log likelihood
    //   sum(log(1 + exp(X * beta))) - y' * X * beta

    arma::vec xbeta = X * beta;
    double yxbeta = dot(y, xbeta);
    // X * beta => exp(X * beta)
    arma::vec exp_xbeta_1 = exp(xbeta) + 1.0;
    return sum(log(exp_xbeta_1)) - yxbeta;
  }

  // Compute the gradient for model parameters beta
  void Gradient(const arma::mat& beta, arma::mat& g)
  {

    // Gradient
    //   X' * (p - y), p = exp(X * beta) / (1 + exp(X * beta))

    // exp(X * beta) => p
    arma::vec exp_xbeta = exp(X * beta);
    exp_xbeta /= (exp_xbeta + 1.0);
    g = X.t() * (exp_xbeta - y);
  }

private:
  // The design matrix.
  const arma::mat& X;
  // The responses to each data point.
  const arma::vec& y;
};


// [[Rcpp::export]]
arma::mat logit_reg_lbfgs(const arma::mat& X, const arma::vec& y) {

  // Construct the first objective function.
  LogisticRegressionFunction lrf(X, y);

  // Create the L_BFGS optimizer with default parameters.
  // The ens::L_BFGS type can be replaced with any ensmallen optimizer that can
  // handle differentiable functions.
  ens::L_BFGS lbfgs;

  lbfgs.MaxIterations() = 10;

  // Create a starting point for our optimization randomly.
  // The model has p parameters, so the shape is p x 1.
  arma::mat beta(X.n_cols, 1, arma::fill::zeros);

  // Run the optimization
  lbfgs.Optimize(lrf, beta);

  // Compute the Hessian matrix with the final beta values:
  // exp(X * beta) => p
  arma::vec exp_xbeta = exp(X * beta);
  arma::vec pred_p = exp_xbeta / (1.0  + exp_xbeta);
  arma::mat diag_w = diagmat(pred_p % (1.0 - pred_p));
  arma::mat info_mat = inv(X.t() * diag_w * X);
  // Get the standard errors:
  arma::vec beta_se = sqrt(info_mat.diag());
  // Compute the z-scores:
  arma::vec z_scores = beta / beta_se;
  // Compute the two-sided p-values:
  arma::vec pvals = 2 * normcdf(-1 * abs(z_scores));

  arma::mat final_results = join_rows(beta, beta_se,
                                      z_scores, pvals);

  return final_results;
}


// Resample a genotype matrix the provided number of times
//' @rdname resample_genotype_data
//' @export
// [[Rcpp::export]]
arma::mat resample_genotype_data(arma::mat genotype_data,
                                 int n_resamples) {

  int n_subjects = genotype_data.n_rows;

  arma::mat sim_genotype_data;

  arma::uvec row_indices = arma::linspace<arma::uvec>(0, n_subjects - 1, n_subjects);
  arma::uvec resample_rows_i = Rcpp::RcppArmadillo::sample(row_indices, n_resamples, true);
  sim_genotype_data = genotype_data.rows(resample_rows_i);

  return sim_genotype_data;
}


// Given a genotype matrix generate simulated GWAS data with provided settings
//' @rdname simulate_gene_gwas_data
//' @export
// [[Rcpp::export]]
arma::mat simulate_gene_gwas_data(arma::mat genotype_data,
                                  bool is_non_null,
                                  arma::uvec causal_snp_i,
                                  double causal_or,
                                  double case_rate) {
  int n_subjects = genotype_data.n_rows;
  arma::vec beta = arma::zeros(genotype_data.n_cols);
  // Initialize the intercept based on the provided case-rate,
  double prob_intercept = log(case_rate / (1.0 - case_rate));

  arma::vec geno_effect = arma::zeros(n_subjects);
  double double_n_subjects = genotype_data.n_rows;
  int n_snps = genotype_data.n_cols;

  // Next determine the effect sizes and predicted probabilities based on the
  // type of gene - if null we don't need to do anything else, but if non-null
  if (is_non_null) {
    // 1 - Get the genotype row sums using only the causal SNPs:
    arma::vec causal_snp_geno_sums = arma::vectorise(arma::sum(genotype_data.cols(causal_snp_i), 1));
    // 2 - What are the different possible unique counts:
    arma::vec unique_geno_counts = arma::unique(causal_snp_geno_sums);
    // 3 - For each unique count, see how many subjects have it:
    arma::vec geno_count_freq(unique_geno_counts.n_elem);
    arma::uvec geno_count_i;
    double geno_i_val;
    for (int i = 0; i < unique_geno_counts.n_elem; i++) {
      geno_i_val = unique_geno_counts[i];
      geno_count_i = find(causal_snp_geno_sums == geno_i_val);
      geno_count_freq[i] = geno_count_i.n_elem;
    }
    // Now update the intercept:
    prob_intercept -= (log(causal_or) * sum(unique_geno_counts % geno_count_freq)) / double_n_subjects;

    arma::vec causal_beta(causal_snp_i.n_elem);
    causal_beta.fill(log(causal_or));

    geno_effect = causal_snp_geno_sums * log(causal_or);

  }

  // Next the probability of being a case:
  arma::vec exp_xbeta = exp(geno_effect + prob_intercept);
  arma::vec case_prob = exp_xbeta / (1.0  + exp_xbeta);
  // Create numeric versions of this to generate the vector of coin tosses for
  // each subject:
  // Generate the vector of coin toss results:
  arma::vec case_status(n_subjects);
  for (int subject_i = 0; subject_i < n_subjects; subject_i++) {
    case_status[subject_i] = R::rbinom(1, case_prob[subject_i]);
  }

  // Now create the matrix of GWAS results:
  arma::mat sim_gwas(n_snps, 4); // since the logit returns four columns
  arma::mat design_mat = arma::ones(n_subjects, 1);
  arma::mat snp_results;
  for (int j = 0; j < n_snps; j++) {
    snp_results = logit_reg_lbfgs(arma::join_horiz(design_mat, genotype_data.col(j)), case_status);
    sim_gwas.row(j) = snp_results.row(1); // Ignore's the intercept row
  }

  return sim_gwas;

}



// Compute correlation matrix with Armadillo
//' @rdname compute_cor_matrix
//' @export
// [[Rcpp::export]]
arma::mat compute_cor_matrix(arma::mat X) {
  int ncols = X.n_cols;

  // Initialize the correlation matrix
  arma::mat cor_matrix(ncols, ncols);

  // Calculate it:
  cor_matrix = arma::cor(X, 0);
  return cor_matrix;
}

// Compute cholesky matrix with Armadillo
//' @rdname compute_chol_matrix
//' @export
// [[Rcpp::export]]
arma::mat compute_chol_matrix(arma::mat X) {
  arma::mat chol_matrix = arma::chol(X, "lower");
  return chol_matrix;
}

// Check if the correlation matrix is positive definite
//' @rdname check_is_pd_cor_matrix
//' @export
// [[Rcpp::export]]
bool check_is_pd_cor_matrix(arma::vec eig_val, double eps) {

  bool result;
  int n_ev = eig_val.n_elem;

  // Compute the tolerance and delta for modifying the eigenvalues:
  double max_ev = max(abs(eig_val));
  double tol = n_ev * max_ev * eps;

  // How many of the eigenvalues are greater than the tolerance:
  arma::uvec eig_val_tol_idx = find(eig_val > tol);
  int pos_eig_val = eig_val_tol_idx.n_elem;
  if (pos_eig_val == n_ev) {
    result = true;
  } else {
    result = false;
  }
  return result;
}

// Make the correlation matrix be positive definite
//' @rdname make_pd_cor_matrix
//' @export
// [[Rcpp::export]]
arma::mat make_pd_cor_matrix(arma::mat cor_matrix,
                             arma::vec eig_val, // vector of eigenvalues
                             arma::mat eig_vec, // matrix of eigenvectors
                             double eps) {
  int ncols = cor_matrix.n_cols;

  // Compute the tolerance and delta for modifying the eigenvalues:
  double max_ev = max(abs(eig_val));
  double tol = ncols * max_ev * eps;
  double delta = 2 * tol;
  double new_ev_i;

  // Loop through and modify the eigen-values:
  for (int i = 0; i < ncols; i++) {
    new_ev_i = delta - eig_val(i);
    if (new_ev_i < 0) {
      new_ev_i = 0;
    }
    eig_val(i) = new_ev_i;
  }

  // Create a diagonal matrix of these new eigen values:
  arma::mat d_eig_val = arma::diagmat(eig_val);

  // Next get the transformation matrix:
  arma::mat trans_matrix = eig_vec * d_eig_val * eig_vec.t();

  arma::mat new_cor_matrix = cor_matrix + trans_matrix;

  return new_cor_matrix;
}

// Force the correlation matrix be positive definite
//' @rdname force_pd_cor_matrix
//' @export
// [[Rcpp::export]]
arma::mat force_pd_cor_matrix(arma::mat cor_matrix, double eps) {
  // Just following this approach: https://vegas2.qimrberghofer.edu.au/vegas2v2
  arma::mat new_cor_matrix(cor_matrix.n_cols, cor_matrix.n_cols);

  int max_loops = 3;
  bool is_pd = false;

  // Loop through max_loops times to modify diagonal if need be:
  for (int i = 0; i < max_loops; i++) {
    new_cor_matrix = cor_matrix + (arma::eye(cor_matrix.n_cols, cor_matrix.n_cols) * (1 + .0001 * pow(10, i)));
    arma::vec new_eig_val = arma::eig_sym(new_cor_matrix);
    is_pd = check_is_pd_cor_matrix(new_eig_val, eps);
    if (is_pd) {
      break;
    }
  }

  return new_cor_matrix;
}


// Compute the cholesky of the correlation matrix
//' @rdname compute_chol_cor_matrix
//' @export
// [[Rcpp::export]]
arma::mat compute_chol_cor_matrix(arma::mat cor_matrix, double eps) {
  int ncols = cor_matrix.n_cols;

  // Initialize the result matrix to eventually return:
  arma::mat final_chol_matrix(ncols, ncols);

  arma::vec eig_val; // vector of eigenvalues
  arma::mat eig_vec; // matrix of eigenvectors
  arma::eig_sym(eig_val, eig_vec, cor_matrix);

  // First check if Armadillo can compute the cholesky
  bool is_pd = check_is_pd_cor_matrix(eig_val, eps);

  if (is_pd == false) {
    // Make it positive definite with the make_pd_cor_matrix function:
    arma::mat new_cor_matrix = make_pd_cor_matrix(cor_matrix, eig_val, eig_vec, eps);
    arma::vec new_eig_val = eig_sym(new_cor_matrix);
    bool now_pd = check_is_pd_cor_matrix(new_eig_val, eps);

    if (now_pd == false) {
      new_cor_matrix = force_pd_cor_matrix(new_cor_matrix, eps);
    }
    final_chol_matrix = compute_chol_matrix(new_cor_matrix);
  } else {
    final_chol_matrix = compute_chol_matrix(cor_matrix);
  }

  return final_chol_matrix;
}



// Compute Brown's covariance approximation sum given correlation matrix
//' @rdname compute_brown_cov_sum
//' @export
// [[Rcpp::export]]
double compute_brown_cov_sum(arma::mat genotype_cor_mat) {

  // Number of SNPs in gene:
  int n_snps = genotype_cor_mat.n_cols;
  double snps_rho;
  double snps_cov;
  double cov_sum = 0;

  for (int i = 0; i < n_snps; i++) {
    for (int j = i+1 ; j < n_snps; j++) {
      snps_rho = genotype_cor_mat(i,j);
      if (snps_rho >= 0) {
        snps_cov = snps_rho * (3.25 + 0.75 * snps_rho);
      } else {
        snps_cov = snps_rho * (3.27 + 0.71 * snps_rho);
      }
      cov_sum += snps_cov;
    }
  }
  cov_sum *= 2;

  return cov_sum;
}

// Compute the two-sided update to Brown's covariance approximation sum given correlation matrix
//' @rdname compute_two_sided_cov_sum
//' @export
// [[Rcpp::export]]
double compute_two_sided_cov_sum(arma::mat genotype_cor_mat) {

  // Number of SNPs in gene:
  int n_snps = genotype_cor_mat.n_cols;
  double cov_sum = 0;

  for (int i = 0; i < n_snps; i++) {
    for (int j = i+1 ; j < n_snps; j++) {
      cov_sum += 4 * pow(genotype_cor_mat(i,j), 2);
    }
  }
  cov_sum *= 2;

  return cov_sum;
}


// Compute the MAGMA version of Brown's combination test
//' @rdname compute_magma
//' @export
// [[Rcpp::export]]
arma::mat compute_magma(arma::vec snp_pvals, arma::mat genotype_cor_mat) {

  // Number of SNPs in gene:
  int n_snps = snp_pvals.n_elem;

  // Expected value for chi-squared:
  int test_exp_value = 2 * n_snps;

  // Get the test statistic:
  double test_stat = -2 * sum(log(snp_pvals));

  // Initialize the matrix to return
  arma::mat result_data(1, 4);
  result_data(0, 0) = test_stat;

  // Initialize the chi-squared parameters:
  double scale;
  double df;
  double test_variance;
  double init_pval;
  double magma_pval;
  // And the MAGMA fudging:
  double chisq_asym_fudge = 1.025;
  double fudge_correction;

  // Square the correlation matrix first:
  arma::mat square_genotype_cor_mat  = square(genotype_cor_mat);

  // First compute the covariance sum using Brown's approximation using the
  // compute_brown_cov_sum function:
  double brown_cov_sum = compute_brown_cov_sum(square_genotype_cor_mat);

  // Compute the test stat variance:
  test_variance = 4 * n_snps + brown_cov_sum;

  // Next the scaling:
  scale = test_variance / (2 * test_exp_value);
  // Degrees of freedom:
  df = (2 * pow(test_exp_value, 2)) / test_variance;
  // Next compute the initial p-value:
  init_pval = R::pchisq(test_stat / scale, df, 0, 0);
  fudge_correction = pow(chisq_asym_fudge, log10(init_pval));
  if (fudge_correction < 0.5) {
    fudge_correction = 0.5;
  }
  // Compute the MAGMA p-value:
  magma_pval = pow(init_pval, fudge_correction);

  result_data(0, 1) = scale;
  result_data(0, 2) = df;
  result_data(0, 3) = magma_pval;

  return result_data;
}

// Compute Brown's combination test given type of p-values
//' @rdname compute_browns
//' @export
// [[Rcpp::export]]
arma::mat compute_browns(arma::vec snp_pvals, arma::mat genotype_cor_mat,
                         bool two_sided) {

  // Number of SNPs in gene:
  int n_snps = snp_pvals.n_elem;

  // Expected value for chi-squared:
  int test_exp_value = 2 * n_snps;

  // Get the test statistic:
  double test_stat = -2 * sum(log(snp_pvals));

  // Initialize the matrix to return
  arma::mat result_data(1, 4);
  result_data(0, 0) = test_stat;

  // Initialize the chi-squared parameters:
  double scale;
  double df;
  double test_variance;
  double browns_pval;
  double brown_cov_sum;

  // Compute the covariance approximation depending on type of p-value:
  if (two_sided) {
    brown_cov_sum = compute_two_sided_cov_sum(genotype_cor_mat);
  } else {
    brown_cov_sum = compute_brown_cov_sum(genotype_cor_mat);
  }

  // Compute the test stat variance:
  test_variance = 4 * n_snps + brown_cov_sum;

  // Next the scaling:
  scale = test_variance / (2 * test_exp_value);
  // Degrees of freedom:
  df = (2 * pow(test_exp_value, 2)) / test_variance;
  // Next compute the p-value:
  browns_pval = R::pchisq(test_stat / scale, df, 0, 0);

  result_data(0, 1) = scale;
  result_data(0, 2) = df;
  result_data(0, 3) = browns_pval;

  return result_data;
}

// Compute gene level statistic using Armadillo
//' @rdname compute_sim_pvals
//' @export
// [[Rcpp::export]]
arma::mat compute_sim_pvals(arma::vec snp_pvals, arma::vec snp_zstats,
                            arma::mat genotype_cor_mat,
                            int n_blocks, int block_size, double eps) {

  double fisher_pval;
  double vegas_pval;
  double fisher_test_stat = -2 * sum(log(snp_pvals));
  double vegas_test_stat = sum(square(snp_zstats));
  int n_snps = genotype_cor_mat.n_cols;
  //int pval_direction = 2;
  //int test_sign = -1;
  arma::mat result_data(1, 4);

  // Initialize the matrices to hold the VEGAS and Fisher test statistics
  // for the given number of simulations and blocks
  arma::mat vegas_data(n_blocks * block_size, 1);
  arma::mat fisher_data(n_blocks * block_size, 1);

  // Compute the cholesky of the correlation matrix:
  arma::mat genotype_cor_chol = compute_chol_cor_matrix(genotype_cor_mat, eps);

  // Now loop through the number of blocks to generate the simulations for the
  // test statistics:
  int i;
  for (i = 0; i < n_blocks; i++) {
    // Generate the uncorrelated Gaussian data:
    arma::mat Y = arma::randn(block_size, n_snps);
    // Generate the z-statistics:
    arma::mat sim_data = trans(genotype_cor_chol * Y.t());

    // Compute and assign the squared versions of the z-statistics:
    vegas_data.rows(i * block_size,
                    i * block_size + (block_size - 1)) = arma::sum(square(sim_data), 1);

    // Compute the Fisher test statistics - just the log transformation of the
    // two-sided p-value for now:
    // Need to first take the negative of each z-statistic:
    arma::mat trans_data = -1 * abs(sim_data);
    // Now compute the fisher statistics:
    arma::mat pval_data = 2 * normcdf(trans_data);
    fisher_data.rows(i * block_size,
                     i * block_size + (block_size - 1)) = -2 * arma::sum(log(pval_data), 1);

  }

  // Compute the p-values given the simulated data, first for Fisher's
  arma::uvec fisher_idx = find(fisher_data <= fisher_test_stat);
  arma::uvec vegas_idx = find(vegas_data <= vegas_test_stat);
  double n_sims = fisher_data.n_rows;
  double fisher_count = fisher_idx.n_elem;
  double vegas_count = vegas_idx.n_elem;
  fisher_pval = 1.0 - (fisher_count / n_sims);
  vegas_pval = 1.0 - (vegas_count / n_sims);

  result_data(0, 0) = fisher_test_stat;
  result_data(0, 1) = fisher_pval;
  result_data(0, 2) = vegas_test_stat;
  result_data(0, 3) = vegas_pval;

  //arma::mat final_data = join_rows(vegas_data, fisher_data);

  return result_data;
}


// Given matrix of p-values and z-statistics, and reference genotype data - compute
// the various types of gene level statistic using Armadillo
//' @rdname compute_gene_level_test
//' @export
// [[Rcpp::export]]
arma::mat compute_gene_level_test(arma::vec snp_pvals, arma::vec snp_zstats,
                                  arma::mat genotype_data,
                                  int sim_n_blocks, int sim_block_size,
                                  double eps) {

  // Need to add checks to see that number of rows for snp_stats is the same
  // as the number of columns for genotype_data

  // First compute the correlation matrix of the genotype data using the wrapper
  // for the Armadillo correlation function:
  arma::mat genotype_cor_matrix = compute_cor_matrix(genotype_data);

  // Proceed to generate the results for the different gene-level tests, where
  // each function should return a matrix of one row and the necessary columns
  // corresponding to the calculations by the gene test

  // First the MAGMA test with fudging
  arma::mat magma_fudge_gene_test = compute_magma(snp_pvals, genotype_cor_matrix);
  // Next the MAGMA test without fudging - just Brown's but with squared correlation values
  arma::mat magma_no_fudge_gene_test = compute_browns(snp_pvals, square(genotype_cor_matrix),
                                                      false);
  // Use the correct covariance approximation
  arma::mat browns_two_sided_gene_test = compute_browns(snp_pvals, genotype_cor_matrix,
                                                        true);
  // Use the incorrect approximation
  arma::mat browns_gene_test = compute_browns(snp_pvals, genotype_cor_matrix, false);
  // The simulated truth statistics - four columns of data, VEGAS and Fisher sims
  arma::mat sim_gene_tests = compute_sim_pvals(snp_pvals, snp_zstats,
                                               genotype_cor_matrix,
                                               sim_n_blocks, sim_block_size,
                                               eps);

  // Join them together:
  arma::mat brown_tests = join_rows(magma_fudge_gene_test, magma_no_fudge_gene_test,
                                    browns_two_sided_gene_test, browns_gene_test);
  arma::mat gene_tests = join_rows(brown_tests, sim_gene_tests);

  // Return
  return gene_tests;
}

// Given a genotype matrix generate simulated GWAS data with provided settings
//' @rdname simulate_gene_level_stats
//' @export
// [[Rcpp::export]]
arma::mat simulate_gene_level_stats(arma::mat init_genotype_data,
                                  bool resample,
                                  int n_resamples,
                                  bool is_non_null,
                                  arma::uvec causal_snp_i,
                                  double causal_or,
                                  double case_rate,
                                  int sim_n_blocks, int sim_block_size,
                                  double eps) {

  arma::mat genotype_data;
  arma::mat sim_gwas_results;
  if (resample) {
    genotype_data = resample_genotype_data(init_genotype_data,n_resamples);
  } else {
    genotype_data = init_genotype_data;
  }
  sim_gwas_results =  simulate_gene_gwas_data(genotype_data,
                                              is_non_null,
                                              causal_snp_i,
                                              causal_or,
                                              case_rate);

  arma::vec sim_pvals = sim_gwas_results.col(3);
  arma::vec sim_zstats = sim_gwas_results.col(2);

  arma::mat gene_stats = compute_gene_level_test(sim_pvals, sim_zstats,
                                                 genotype_data,
                                                 sim_n_blocks, sim_block_size,
                                                 eps);

  return gene_stats;

}


// Generate simulated truth data using Armadillo
//' @rdname generate_sim_truth_data
//' @export
// [[Rcpp::export]]
arma::mat generate_sim_truth_data(arma::mat genotype_cor_matrix,
                                  int n_blocks,
                                  int block_size,
                                  double eps) {

  int n_snps = genotype_cor_matrix.n_cols;

  // Initialize the matrix to hold the VEGAS and Fisher test statistics
  // for the given number of simulations and blocks:
  arma::mat vegas_data(n_blocks * block_size, 1);
  arma::mat fisher_data(n_blocks * block_size, 1);

  // Compute the cholesky of the correlation matrix:
  arma::mat genotype_cor_chol = compute_chol_cor_matrix(genotype_cor_matrix, eps);

  // Now loop through the number of blocks to generate the simulations for the
  // test statistics:
  for (int i = 0; i < n_blocks; i++) {
    // Generate the uncorrelated Gaussian data:
    arma::mat Y = arma::randn(block_size, n_snps);
    // Generate the z-statistics:
    arma::mat sim_data = trans(genotype_cor_chol * Y.t());

    // Compute and assign the squared versions of the z-statistics:
    vegas_data.rows(i * block_size,
                    i * block_size + (block_size - 1)) = arma::sum(square(sim_data), 1);

    // Compute the Fisher test statistics - just the log transformation of the
    // two-sided p-value for now:
    // Need to first take the negative of each z-statistic:
    arma::mat trans_data = -1 * abs(sim_data);
    // Now compute the fisher statistics:
    arma::mat pval_data = 2 * normcdf(trans_data);
    fisher_data.rows(i * block_size,
                     i * block_size + (block_size - 1)) = -2 * arma::sum(log(pval_data), 1);

  }

  arma::mat sim_truth_test_data = arma::join_horiz(vegas_data, fisher_data);

  return sim_truth_test_data;
}

// Compute gene level statistic using Armadillo
//' @rdname compute_fixed_sim_pvals
//' @export
// [[Rcpp::export]]
arma::mat compute_fixed_sim_pvals(arma::vec snp_pvals, arma::vec snp_zstats,
                                  arma::mat sim_truth_matrix) {

  double fisher_pval;
  double vegas_pval;
  double fisher_test_stat = -2 * sum(log(snp_pvals));
  double vegas_test_stat = sum(square(snp_zstats));

  arma::mat result_data(1, 4);
  arma::vec vegas_data = sim_truth_matrix.col(0);
  arma::vec fisher_data = sim_truth_matrix.col(1);


  // Compute the p-values given the simulated data, first for Fisher's
  arma::uvec fisher_idx = find(fisher_data <= fisher_test_stat);
  arma::uvec vegas_idx = find(vegas_data <= vegas_test_stat);
  double n_sims = sim_truth_matrix.n_rows;
  double fisher_count = fisher_idx.n_elem;
  double vegas_count = vegas_idx.n_elem;
  fisher_pval = 1.0 - (fisher_count / n_sims);
  vegas_pval = 1.0 - (vegas_count / n_sims);

  result_data(0, 0) = fisher_test_stat;
  result_data(0, 1) = fisher_pval;
  result_data(0, 2) = vegas_test_stat;
  result_data(0, 3) = vegas_pval;

  return result_data;
}



// Given simulated GWAS results with simulated truth test statistics according
// to the correlation matrix, compute the gene level test statistics
//' @rdname compute_fixed_gene_level_test
//' @export
// [[Rcpp::export]]
arma::mat compute_fixed_gene_level_test(arma::mat sim_gwas_data,
                                        arma::mat sim_truth_matrix,
                                        arma::mat genotype_cor_matrix) {

  // Intialize the vector of pvalues and z-statistics given column structure
  arma::vec snp_pvals = sim_gwas_data.col(3);
  arma::vec snp_zstats = sim_gwas_data.col(2);

  // Next compute each of the different tests
  // First the MAGMA test with fudging
  arma::mat magma_fudge_gene_test = compute_magma(snp_pvals, genotype_cor_matrix);
  // Next the MAGMA test without fudging - just Brown's but with squared correlation values
  arma::mat magma_no_fudge_gene_test = compute_browns(snp_pvals, square(genotype_cor_matrix),
                                                      false);
  // Use the correct covariance approximation
  arma::mat browns_two_sided_gene_test = compute_browns(snp_pvals, genotype_cor_matrix,
                                                        true);
  // Use the incorrect approximation
  arma::mat browns_gene_test = compute_browns(snp_pvals, genotype_cor_matrix, false);
  // The simulated truth statistics - four columns of data, VEGAS and Fisher sims
  arma::mat sim_gene_tests = compute_fixed_sim_pvals(snp_pvals, snp_zstats,
                                                     sim_truth_matrix);

  // Join them together:
  arma::mat brown_tests = join_rows(magma_fudge_gene_test, magma_no_fudge_gene_test,
                                    browns_two_sided_gene_test, browns_gene_test);
  arma::mat gene_tests = join_rows(brown_tests, sim_gene_tests);

  // Return
  return gene_tests;
}



// Given a genotype matrix generate simulated GWAS data with provided settings
//' @rdname null_fixed_cor_gene_test_sim
//' @export
// [[Rcpp::export]]
arma::mat fixed_cor_gene_test_sim(arma::mat genotype_data,
                                  int n_gene_sims,
                                  bool is_non_null,
                                  arma::uvec causal_snp_i,
                                  double causal_or,
                                  double case_rate,
                                  int truth_sim_n_blocks,
                                  int truth_sim_block_size,
                                  double eps) {

  // First compute the correlation matrix of the genotype data using the wrapper
  // for the Armadillo correlation function:
  arma::mat genotype_cor_matrix = compute_cor_matrix(genotype_data);

  // Generate the matrix of z-stats and p-values to use for computing the
  // simulated truth p-values:
  arma::mat sim_truth_matrix = generate_sim_truth_data(genotype_cor_matrix,
                                                       truth_sim_n_blocks,
                                                       truth_sim_block_size,
                                                       eps);
  // Initialize the sim truth vectors to pass into the function for computing
  // the simulation based p-values:
  arma::vec truth_sim_pvals = sim_truth_matrix.col(0);
  arma::vec truth_sim_zstats = sim_truth_matrix.col(1);

  // Intialize the simulation results matrix:
  arma::mat sim_gene_tests(n_gene_sims, 20); // TEMPORARY STORE 20 COLS
  for (int i = 0; i < n_gene_sims; i++) {

    // Generate the GWAS results:
    arma::mat sim_gwas_data =  simulate_gene_gwas_data(genotype_data,
                                                       is_non_null, causal_snp_i,
                                                       causal_or, case_rate);
    sim_gene_tests.row(i) = compute_fixed_gene_level_test(sim_gwas_data,
                                                          sim_truth_matrix,
                                                          genotype_cor_matrix);
  }
  // Return the final matrix of results
  return sim_gene_tests;

}





