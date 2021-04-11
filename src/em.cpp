#include "utils.h"

/* Helper Functions */
double calculate_sd(double *data, int size_data);
double function_dnorm(double x, double mean, double sigma);
double * function_forward(double * f_line, double para_mu_pos, double para_mu_neg, double para_sd0, double para_tao1, double para_tao2, double para_rho, double * para_tran_p, double o_b, int dist_b);
double * function_backward(double * b_line, double para_mu_pos, double para_mu_neg, double para_sd0, double para_tao1, double para_tao2, double para_rho, double * para_tran_p, double o_b, int dist_b);
double function_update_p(double hk_init, int * obs_dist, double * P_kl_all, double rho, int k, int Bins);
double function_update_rho(double rho_init, int * obs_dist, double * P_kl_all, double * tran_p_new, int Bins);
static double Brent_fmin(double ax, double bx, double tol, double f(double x, int * obs_dist, double * P_kl_all, double rho, int hmmDMR_k, int Bins), double *x, int * obs_dist, double * P_kl_all, double rho, int hmmDMR_k, int Bins);

/* hummingbirdEMinternal: Expectation-Maximization Algorithm for Fitting the Hidden Markov Model.
 * This function reads in methylated and unmethylated read count data,
 * transforms it into logarithm bin-wise data, sets up initial values
 * and implements the EM algorithm to estimate HMM parameters and find
 * the best sequence of hidden states based on model fitting.
 */
// [[Rcpp::export]]
SEXP hummingbirdEMinternal(SEXP normM, SEXP normUM, SEXP abnormM, SEXP abnormUM, SEXP pos, SEXP binSize)
{
	Rprintf("Reading input...\n");

	// Bin size...
	int bin_size = 0;

	if(TYPEOF(binSize) == NILSXP){
		bin_size = 40; //default
	}
	else{
		bin_size = asReal(binSize);
	}

	Rprintf("Bin size: %d.\n", bin_size);

	// Reading matrices...
	SEXP S_dim_matrices = getAttrib(normM, R_DimSymbol);
	int total_lines = INTEGER(S_dim_matrices)[0];
	int total_repl = INTEGER(S_dim_matrices)[1];
	
	SEXP S_dim_matrices_ab = getAttrib(abnormM, R_DimSymbol);
	int total_repl_ab = INTEGER(S_dim_matrices_ab)[1];

	Rprintf("Total lines: %d, total replicates in normal group: %d and in abnormal group: %d\n", total_lines, total_repl, total_repl_ab);

	double *matrix_normM, *matrix_normUM, *matrix_ABnormM, *matrix_ABnormUM, *matrix_normC, *matrix_ABnormC;
	int *matrix_pos;

	matrix_normM = (double *)malloc(sizeof(double)*total_lines*total_repl);
	matrix_normUM = (double *)malloc(sizeof(double)*total_lines*total_repl);
	matrix_ABnormM = (double *)malloc(sizeof(double)*total_lines*total_repl_ab);
	matrix_ABnormUM = (double *)malloc(sizeof(double)*total_lines*total_repl_ab);
	matrix_normC = (double *)malloc(sizeof(double)*total_lines*total_repl);
	matrix_ABnormC = (double *)malloc(sizeof(double)*total_lines*total_repl_ab);
	matrix_pos = (int *)malloc(sizeof(int)*total_lines);

	int i, j, j_ab;
	for(i=0;i<total_lines;i++){
		matrix_pos[i] = INTEGER(pos)[i];

		for(j=0;j<total_repl;j++){
			matrix_normM[idx2c(i,j,total_lines)] = INTEGER(normM)[idx2c(i,j,total_lines)] + 0.5;
			matrix_normUM[idx2c(i,j,total_lines)] = INTEGER(normUM)[idx2c(i,j,total_lines)] + 0.5;

		}
		
		for(j_ab=0;j_ab<total_repl_ab;j_ab++){
			matrix_ABnormM[idx2c(i,j_ab,total_lines)] = INTEGER(abnormM)[idx2c(i,j_ab,total_lines)] + 0.5;
			matrix_ABnormUM[idx2c(i,j_ab,total_lines)] = INTEGER(abnormUM)[idx2c(i,j_ab,total_lines)] + 0.5;
		}
		
	}

	for(i=0;i<total_lines;i++){
		for(j=0;j<total_repl;j++){
			matrix_normC[idx2c(i,j,total_lines)] = matrix_normM[idx2c(i,j,total_lines)] + matrix_normUM[idx2c(i,j,total_lines)];
		}
		
		for(j_ab=0;j_ab<total_repl_ab;j_ab++){
			matrix_ABnormC[idx2c(i,j_ab,total_lines)] = matrix_ABnormM[idx2c(i,j_ab,total_lines)] + matrix_ABnormUM[idx2c(i,j_ab,total_lines)];
		}
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/******************************* Read Process ********************************/
	/*****************************************************************************/
	/*****************************************************************************/

	Rprintf("Processing input...\n");

	int dist_flag = 1;
	int start_pos = matrix_pos[0];
	int end_pos = matrix_pos[total_lines-1];
	int numofBins = ceil((end_pos - start_pos)/bin_size)+1;
	int start = 0;
	int end = 0;
	int k = 0;
	int kcol, kcol_ab;
	int found;
	int temp_dist = 0;

	double sum_normM = 0.0;
	double sum_normC = 0.0;
	double sum_abnormM = 0.0;
	double sum_abnormC = 0.0;

	int full_bins = 0;

	int *obs_dist = NULL, *obs_start = NULL, *obs_end = NULL;
	double *obs_norm_p = NULL, *obs_abnorm_p = NULL, *obs_o = NULL;

	for(i = 0; i < numofBins; i++)
	{
		start = start_pos + i*bin_size;
		end = start_pos + (i+1)*bin_size-1;

		found = 0;
		sum_normM = 0.0;
		sum_normC = 0.0;
		sum_abnormM = 0.0;
		sum_abnormC = 0.0;
		while(matrix_pos[k] >= start && matrix_pos[k] <= end)
		{
			found = 1;

			for(kcol = 0; kcol < total_repl; kcol++){
				sum_normM = sum_normM + matrix_normM[idx2c(k,kcol,total_lines)];
				sum_normC = sum_normC + matrix_normC[idx2c(k,kcol,total_lines)];
			}
			
			for(kcol_ab = 0; kcol_ab < total_repl_ab; kcol_ab++){
				sum_abnormM = sum_abnormM + matrix_ABnormM[idx2c(k,kcol_ab,total_lines)];
				sum_abnormC = sum_abnormC + matrix_ABnormC[idx2c(k,kcol_ab,total_lines)];
			}

			k++;
		}

		if(found == 1)
		{
			full_bins++;

			obs_norm_p = (double*) realloc (obs_norm_p, full_bins * sizeof(double));
			if(obs_norm_p == NULL){
				free(obs_norm_p);
				error("Error allocating memory.\n");
			}
			obs_abnorm_p = (double*) realloc (obs_abnorm_p, full_bins * sizeof(double));
			if(obs_abnorm_p == NULL){
				free(obs_abnorm_p);
				error("Error allocating memory.\n");
			}
			obs_o = (double*) realloc (obs_o, full_bins * sizeof(double));
			if(obs_o == NULL){
				free(obs_o);
				error("Error allocating memory.\n");
			}
			obs_dist = (int*) realloc (obs_dist, full_bins * sizeof(int));
			if(obs_dist == NULL){
				free(obs_dist);
				error("Error allocating memory.\n");
			}
			obs_start = (int*) realloc (obs_start, full_bins * sizeof(int));
			if(obs_start == NULL){
				free(obs_start);
				error("Error allocating memory.\n");
			}
			obs_end = (int*) realloc (obs_end, full_bins * sizeof(int));
			if(obs_end == NULL){
				free(obs_end);
				error("Error allocating memory.\n");
			}

			obs_norm_p[full_bins-1] = sum_normM/sum_normC;
			obs_abnorm_p[full_bins-1] = sum_abnormM/sum_abnormC;

			obs_o[full_bins-1] = log(obs_abnorm_p[full_bins-1]/(1-obs_abnorm_p[full_bins-1])) - log(obs_norm_p[full_bins-1]/(1-obs_norm_p[full_bins-1]));

			obs_start[full_bins-1] = start;
			obs_end[full_bins-1] = end;

			if(dist_flag == 0)
			{
				obs_dist[full_bins-1] = (temp_dist + 1)*bin_size;
				dist_flag = 1;
			}
			else if(dist_flag == 1)
			{
				obs_dist[full_bins-1] = (1)*bin_size;
				dist_flag = 1;
			}
		}
		else if(found == 0 && dist_flag == 0)
		{
			temp_dist = temp_dist + 1;
			dist_flag = 0;
		}
		else if(found == 0 && dist_flag == 1)
		{
			temp_dist = 1;
			dist_flag = 0;
		}

	}

	Rprintf("Processing input completed...\n");

	// Free memory
	free(matrix_normM);
	free(matrix_normUM);
	free(matrix_ABnormM);
	free(matrix_ABnormUM);
	free(matrix_normC);
	free(matrix_ABnormC);

	/*****************************************************************************/
	/*****************************************************************************/
	/******************************* Initial Value *******************************/
	/*****************************************************************************/
	/*****************************************************************************/

	Rprintf("Calculation of the initial value...\n");

	double sd;
	sd = calculate_sd(obs_o, full_bins);

	int sum_initial_state_1 = 0, sum_initial_state_2 = 0;
	double sum_obs_1 = 0.0, sum_obs_2 = 0.0;

	for(i = 0; i < full_bins; i++)
	{
		if(obs_o[i] > sd)
		{
			sum_initial_state_1++;

			sum_obs_1 = sum_obs_1 + obs_o[i];
		}
		else if(obs_o[i] < -sd)
		{
			sum_initial_state_2++;

			sum_obs_2 = sum_obs_2 + obs_o[i];
		}
	}

	// Emission parameters
	double tmp_mu_pos, tmp_mu_neg;

	if((sum_initial_state_1 != 0) && (sum_initial_state_1 != 1))
	{
		tmp_mu_pos = sum_obs_1/sum_initial_state_1;
	}
	else
	{
		tmp_mu_pos = (3/2)*sd;
	}

	if((sum_initial_state_2 != 0) && (sum_initial_state_2 != 1))
	{
		tmp_mu_neg = sum_obs_2/sum_initial_state_2;
	}
	else
	{
		tmp_mu_neg = -(3/2)*sd;
	}

	double mu_pos, mu_neg, sd0, tao1, tao2;

	mu_pos = (tmp_mu_pos+fabs(tmp_mu_neg))/2;
	mu_neg = (-fabs(tmp_mu_pos)+tmp_mu_neg)/2;
	sd0 = sd;
	tao1 = sd;
	tao2 = sd;

	// Transition parameters
	double p0 = (double)1/3, p1 = (double)1/3, p2 = (double)1/3, rho = 0.05;
	double * tran_p;

	tran_p = (double *)malloc(sizeof(double)*3*3);

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_p[idx2c(i,j,3)] = 0.1;
		}
	}

	Rprintf("Initial Value calculated...\n");

	/*****************************************************************************/
	/*****************************************************************************/
	/************************************* EM ************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	Rprintf("EM begins...\n");

	int iter = 1;

	// Final states
	int *direction;
	direction = (int *)malloc(sizeof(int)*full_bins);

	int count_states = 0;

	while(1){
	Rprintf("Iteration: %d\n", iter);
	if(iter % 20 == 0){
		R_CheckUserInterrupt();
	}

	int B = full_bins, bin_cnt = 0;

	// Forward algorithm
	double *f;
	f = (double *)malloc(sizeof(double)*B*3);

	f[idx2c(bin_cnt,0,B)] = p0*function_dnorm(obs_o[bin_cnt], 0.0, sd0);
	f[idx2c(bin_cnt,1,B)] = p1*function_dnorm(obs_o[bin_cnt], mu_pos, sd0);
	f[idx2c(bin_cnt,2,B)] = p2*function_dnorm(obs_o[bin_cnt], mu_neg, sd0);

	double *f_line;
	f_line = (double *)malloc(sizeof(double)*3);

	while(bin_cnt < B-1)
	{
		bin_cnt++;

		f_line[0] = f[idx2c(bin_cnt-1,0,B)];
		f_line[1] = f[idx2c(bin_cnt-1,1,B)];
		f_line[2] = f[idx2c(bin_cnt-1,2,B)];

		f_line = function_forward(f_line, mu_pos, mu_neg, sd0, tao1, tao2, rho, tran_p, obs_o[bin_cnt], obs_dist[bin_cnt]);

		f[idx2c(bin_cnt,0,B)] = f_line[0];
		f[idx2c(bin_cnt,1,B)] = f_line[1];
		f[idx2c(bin_cnt,2,B)] = f_line[2];

	}

	// Free memory
	free(f_line);

	// Backward algorithm
	bin_cnt = B-1;

	double *b;
	b = (double *)malloc(sizeof(double)*B*3);

	double *b_line;
	b_line = (double *)malloc(sizeof(double)*3);

	b[idx2c(bin_cnt,0,B)] = 1.0;
	b[idx2c(bin_cnt,1,B)] = 1.0;
	b[idx2c(bin_cnt,2,B)] = 1.0;

	b_line[0] = b[idx2c(bin_cnt,0,B)];
	b_line[1] = b[idx2c(bin_cnt,1,B)];
	b_line[2] = b[idx2c(bin_cnt,2,B)];

	while(bin_cnt > 0)
	{
		bin_cnt--;

		b_line = function_backward(b_line, mu_pos, mu_neg, sd0, tao1, tao2, rho, tran_p, obs_o[bin_cnt+1], obs_dist[bin_cnt+1]);

		b[idx2c(bin_cnt,0,B)] = b_line[0];
		b[idx2c(bin_cnt,1,B)] = b_line[1];
		b[idx2c(bin_cnt,2,B)] = b_line[2];
	}

	// Free memory
	free(b_line);

	// E step continuing
	double *mult_fb, *sum_mult_fb, *P_k;

	mult_fb = (double *)malloc(sizeof(double)*B*3);
	sum_mult_fb = (double *)malloc(sizeof(double)*B);
	P_k = (double *)malloc(sizeof(double)*B*3);

	for(i=0;i<B;i++){
		for(j=0;j<3;j++){
			mult_fb[idx2c(i,j,B)] = f[idx2c(i,j,B)]*b[idx2c(i,j,B)];
		}
		sum_mult_fb[i] = mult_fb[idx2c(i,0,B)] + mult_fb[idx2c(i,1,B)] + mult_fb[idx2c(i,2,B)];
	}

	for(i=0;i<B;i++){
		for(j=0;j<3;j++){
			P_k[idx2c(i,j,B)] = mult_fb[idx2c(i,j,B)]/sum_mult_fb[i];
		}
	}

	double sum_P_kl_loga_kl, sum_P_k_logOb;
	sum_P_kl_loga_kl = 0.0;
	sum_P_k_logOb = 0.0;

	double *P_kl_all;
	P_kl_all = (double *)malloc(sizeof(double)*(B-1)*9);

	double * tran_a;
	tran_a = (double *)malloc(sizeof(double)*3*3);

	double * P_kl_temp;
	P_kl_temp = (double *)malloc(sizeof(double)*9);

	double dist_rate;

	double * temporary_tran_a;
	temporary_tran_a = (double *)malloc(sizeof(double)*9);

	double * Ob;
	Ob = (double *)malloc(sizeof(double)*3);

	for(bin_cnt=0; bin_cnt<=B-2; bin_cnt++){

		dist_rate = 1 - exp(-rho*obs_dist[bin_cnt+1]);

		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				tran_a[idx2c(i,j,3)] = tran_p[idx2c(i,j,3)]*dist_rate;
			}
		}

		tran_a[idx2c(0,0,3)] = 1 -(tran_a[idx2c(0,1,3)] + tran_a[idx2c(0,2,3)]);
		tran_a[idx2c(1,1,3)] = 1 -(tran_a[idx2c(1,0,3)] + tran_a[idx2c(1,2,3)]);
		tran_a[idx2c(2,2,3)] = 1 -(tran_a[idx2c(2,0,3)] + tran_a[idx2c(2,1,3)]);

		for(i=0;i<3;i++){
			P_kl_temp[i*3] = f[idx2c(bin_cnt,i,B)]*tran_a[idx2c(i,0,3)]*function_dnorm(obs_o[bin_cnt+1], 0.0, sd0)*b[idx2c(bin_cnt+1,0,B)];
			P_kl_temp[i*3+1] = f[idx2c(bin_cnt,i,B)]*tran_a[idx2c(i,1,3)]*function_dnorm(obs_o[bin_cnt+1], mu_pos, tao1)*b[idx2c(bin_cnt+1,1,B)];
			P_kl_temp[i*3+2] = f[idx2c(bin_cnt,i,B)]*tran_a[idx2c(i,2,3)]*function_dnorm(obs_o[bin_cnt+1], mu_neg, tao2)*b[idx2c(bin_cnt+1,2,B)];
		}

		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				if(tran_a[idx2c(i,j,3)] == 0.0){
					tran_a[idx2c(i,j,3)] = 0.0001;
				}
			}
		}

		double sum_P_kl_temp = 0.0;

		for(i=0;i<9;i++){
			sum_P_kl_temp = sum_P_kl_temp + P_kl_temp[i];
		}

		for(j=0;j<9;j++){
			P_kl_all[idx2c(bin_cnt,j,B-1)] = P_kl_temp[j]/sum_P_kl_temp;
		}

		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				temporary_tran_a[i*3+j] = log(tran_a[idx2c(i,j,3)]);
			}
		}

		for(j=0;j<9;j++){
			sum_P_kl_loga_kl = sum_P_kl_loga_kl + P_kl_all[idx2c(bin_cnt,j,B-1)]*temporary_tran_a[j];
		}

		Ob[0] = function_dnorm(obs_o[bin_cnt], 0.0, sd0);
		if(Ob[0] == 0.0){
			Ob[0] = 0.0001;
		}
		Ob[1] = function_dnorm(obs_o[bin_cnt], mu_pos, tao1);
		if(Ob[1] == 0.0){
			Ob[1] = 0.0001;
		}
		Ob[2] = function_dnorm(obs_o[bin_cnt], mu_neg, tao2);
		if(Ob[2] == 0.0){
			Ob[2] = 0.0001;
		}

		for(j=0;j<3;j++){
			sum_P_k_logOb = sum_P_k_logOb + P_k[idx2c(bin_cnt,j,B)]*log(Ob[j]);
		}

		if(bin_cnt == (B-2))
		{

			Ob[0] = function_dnorm(obs_o[B-1], 0.0, sd0);
			if(Ob[0] == 0.0){
				Ob[0] = 0.0001;
			}
			Ob[1] = function_dnorm(obs_o[B-1], mu_pos, tao1);
			if(Ob[1] == 0.0){
				Ob[1] = 0.0001;
			}
			Ob[2] = function_dnorm(obs_o[B-1], mu_neg, tao2);
			if(Ob[2] == 0.0){
				Ob[2] = 0.0001;
			}

			for(j=0;j<3;j++){
				sum_P_k_logOb = sum_P_k_logOb + P_k[idx2c(B-1,j,B)]*log(Ob[j]);
			}

		}

	}

	//Keep old ps
	double p0_current, p1_current, p2_current;
	p0_current = p0;
	p1_current = p1;
	p2_current = p2;

	if(p0 == 0.0){
		p0 = 0.0001;
	}
	if(p1 == 0.0){
		p1 = 0.0001;
	}
	if(p2 == 0.0){
		p2 = 0.0001;
	}

	double sum_fin_eloglik = 0.0;

	for(j=0;j<3;j++){
		if(j == 0)
			sum_fin_eloglik = sum_fin_eloglik + P_k[idx2c(0,j,B)]*log(p0);
		if(j == 1)
			sum_fin_eloglik = sum_fin_eloglik + P_k[idx2c(0,j,B)]*log(p1);
		if(j == 2)
			sum_fin_eloglik = sum_fin_eloglik + P_k[idx2c(0,j,B)]*log(p2);
	}

	// Free memory
	free(mult_fb);
	free(sum_mult_fb);
	free(tran_a);
	free(P_kl_temp);
	free(temporary_tran_a);

	// Keep values of E step, i.e. current
	double * current_tran_p;
	current_tran_p = (double *)malloc(sizeof(double)*3*3);

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			current_tran_p[idx2c(i,j,3)] = tran_p[idx2c(i,j,3)];
		}
	}

	double current_mu_pos, current_sd0, current_rho;
	current_mu_pos = mu_pos;
	current_sd0 = sd0;
	current_rho = rho;

	/*****************************************************************************/

	// Update p_kl where (l != k).
	double * hks;
	hks = (double *)malloc(sizeof(double)*3);

	int l1, l2;
	double * temporary_tran_p;
	temporary_tran_p = (double *)malloc(sizeof(double)*9);

	for(k=0;k<3;k++){

		if(k == 0){
			l1 = 1;
			l2 = 2;
		}
		if(k == 1){
			l1 = 3;
			l2 = 5;
		}
		if(k == 2){
			l1 = 6;
			l2 = 7;
		}

		for(i=0;i<3;i++){
			for(j=0;j<3;j++){
				temporary_tran_p[i*3+j] = tran_p[idx2c(i,j,3)];
			}
		}

		double sum_of_l1 = 0.0, sum_of_l2 = 0.0;

		for(i=0;i<B-1;i++){
			sum_of_l1 = sum_of_l1 + P_kl_all[idx2c(i,l1,B-1)];
			sum_of_l2 = sum_of_l2 + P_kl_all[idx2c(i,l2,B-1)];
		}

		double hk_init;
		hk_init = (sum_of_l1/temporary_tran_p[l1] + sum_of_l2/temporary_tran_p[l2])/2.0;

		double brent_ax, brent_bx, brent_tol;

		brent_tol = 0.01;

		brent_ax = sum_of_l1 + sum_of_l2;
		brent_bx = ((double)100)*brent_ax;

		Brent_fmin(brent_ax, brent_bx, brent_tol, function_update_p, &hk_init, obs_dist, P_kl_all, rho, k, B);

		hks[k] = hk_init;
	}

	double * tran_p_new;
	tran_p_new = (double *)malloc(sizeof(double)*3*3);

	double * sum_of_all;
	sum_of_all = (double *)malloc(sizeof(double)*9);
	for(j=0;j<9;j++){
		sum_of_all[j] = 0.0;
	}

	for(i=0;i<B-1;i++){
		for(j=0;j<9;j++){
			sum_of_all[j] = sum_of_all[j] + P_kl_all[idx2c(i,j,B-1)];
		}
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_p_new[idx2c(i,j,3)] = sum_of_all[i*3+j]/hks[i];
		}
	}

	// Using grid search to find a new value for rho.
	double o_upd_rho = 0.0;
	double rho_init, temp_rho_new, rho_new;
	int counter_grid = 0;

	rho_init = 0.002;

	while(rho_init <= 0.1){

		o_upd_rho = function_update_rho(rho_init, obs_dist, P_kl_all, tran_p_new, B);

		if(counter_grid == 0){
			temp_rho_new = o_upd_rho;
			rho_new = rho_init;
			counter_grid = 1;
		}
		else{
			if(o_upd_rho > temp_rho_new){
				temp_rho_new = o_upd_rho;
				rho_new = rho_init;
			}
		}

		rho_init = rho_init + 0.002;

	}

	// Update mu+, mu-, sd0, tao1, and tao2.
	double tmp_mu_pos_new, tmp_mu_neg_new, mu_pos_new, mu_neg_new;
	double sd0_new, tao1_new, tao2_new, tmp_sigma_new;
	double p0_new, p1_new, p2_new;

	tmp_mu_pos_new = 0.0;
	tmp_mu_neg_new = 0.0;

	double sum_pos_mult = 0.0, sum_pos = 0.0;
	double sum_neg_mult = 0.0, sum_neg = 0.0;
	double sum_sd_pow_calc = 0.0, sum_sd_calc = 0.0;
	double sum_tao_pos_calc = 0.0, sum_tao_neg_calc = 0.0;

	for(i = 0; i < B; i++)
	{
		sum_pos_mult = sum_pos_mult + obs_o[i]*P_k[idx2c(i,1,B)];
		sum_pos = sum_pos + P_k[idx2c(i,1,B)];

		sum_neg_mult = sum_neg_mult + obs_o[i]*P_k[idx2c(i,2,B)];
		sum_neg = sum_neg + P_k[idx2c(i,2,B)];

		sum_sd_pow_calc = sum_sd_pow_calc + pow(obs_o[i],2)*P_k[idx2c(i,0,B)];
		sum_sd_calc = sum_sd_calc + P_k[idx2c(i,0,B)];

	}

	tmp_mu_pos_new = sum_pos_mult/sum_pos;
	tmp_mu_neg_new = sum_neg_mult/sum_neg;

	mu_pos_new = (tmp_mu_pos_new + fabs(tmp_mu_neg_new))/2;
	mu_neg_new = (-fabs(tmp_mu_pos_new) + tmp_mu_neg_new)/2;

	sd0_new = sqrt(sum_sd_pow_calc/sum_sd_calc);

	for(i = 0; i < B; i++)
	{
		sum_tao_pos_calc = sum_tao_pos_calc + pow((obs_o[i]-mu_pos_new),2)*P_k[idx2c(i,1,B)];
		sum_tao_neg_calc = sum_tao_neg_calc + pow((obs_o[i]-mu_neg_new),2)*P_k[idx2c(i,2,B)];
	}

	tao1_new = sqrt(sum_tao_pos_calc/sum_pos);
	tao2_new = sqrt(sum_tao_neg_calc/sum_neg);

	tmp_sigma_new = (sd0_new + tao1_new + tao2_new)/3;

	sd0_new = tmp_sigma_new;
	tao1_new = tmp_sigma_new;
	tao2_new = tmp_sigma_new;

	p0_new = P_k[idx2c(0,0,B)];
	p1_new = P_k[idx2c(0,1,B)];
	p2_new = P_k[idx2c(0,2,B)];

	// updating values
	p0 = p0_new;
	p1 = p1_new;
	p2 = p2_new;
	mu_pos = mu_pos_new;
	mu_neg = mu_neg_new;
	sd0 = sd0_new;
	tao1 = tao1_new;
	tao2 = tao2_new;
	rho = rho_new;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_p[idx2c(i,j,3)] = tran_p_new[idx2c(i,j,3)];
		}
	}

	// EM continuing
	double tol = 0.01;

	double sum_curr_tranp = 0.0, sum_tranp = 0.0;
	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			sum_curr_tranp = sum_curr_tranp + current_tran_p[idx2c(i,j,3)];
			sum_tranp = sum_tranp + tran_p[idx2c(i,j,3)];
		}
	}

	count_states = 0;

	if(fabs(current_mu_pos - mu_pos) < tol && \
	fabs(current_sd0 - sd0) < tol && \
	fabs(current_rho - rho) < tol && \
	fabs(current_tran_p[idx2c(0,1,3)] - tran_p_new[idx2c(0,1,3)]) < tol && \
	fabs(current_tran_p[idx2c(0,2,3)] - tran_p_new[idx2c(0,2,3)]) < tol && \
	fabs(current_tran_p[idx2c(1,0,3)] - tran_p_new[idx2c(1,0,3)]) < tol && \
	fabs(current_tran_p[idx2c(1,2,3)] - tran_p_new[idx2c(1,2,3)]) < tol && \
	fabs(current_tran_p[idx2c(2,0,3)] - tran_p_new[idx2c(2,0,3)]) < tol && \
	fabs(current_tran_p[idx2c(2,1,3)] - tran_p_new[idx2c(2,1,3)]) < tol && \
	fabs(p0_current - p0_new) < tol && \
	fabs(p1_current - p1_new) < tol && \
	fabs(p2_current - p2_new) < tol ){
		// Find the sequence of best states
		Rprintf("Calculation of states...\n");

		for(i=0;i<B;i++){
			if((P_k[idx2c(i,0,B)] >= P_k[idx2c(i,1,B)]) && (P_k[idx2c(i,0,B)] >= P_k[idx2c(i,2,B)])){
				direction[i] = 0;
				count_states = count_states +1;
			}
			else if((P_k[idx2c(i,1,B)] >= P_k[idx2c(i,0,B)]) && (P_k[idx2c(i,1,B)] >= P_k[idx2c(i,2,B)])){
				direction[i] = 1;
				count_states = count_states +1;
			}
			else if((P_k[idx2c(i,2,B)] >= P_k[idx2c(i,0,B)]) && (P_k[idx2c(i,2,B)] >= P_k[idx2c(i,1,B)])){
				direction[i] = -1;
				count_states = count_states +1;
			}
		}

		// Free memory
		free(f);
		free(b);
		free(P_k);
		free(P_kl_all);
		free(temporary_tran_p);
		free(hks);
		free(tran_p_new);
		free(current_tran_p);

		// Convergence
		Rprintf("EM converged after %d iterations.\n", iter);
		break;
	}
	else{
		iter = iter + 1;
	}

	// Free memory
	free(f);
	free(b);
	free(P_k);
	free(P_kl_all);
	free(temporary_tran_p);
	free(hks);
	free(tran_p_new);
	free(current_tran_p);

	}

	Rprintf("Saving output...\n");

	SEXP S_em = PROTECT(allocMatrix(INTSXP, full_bins, 4));
	for(i=0;i<full_bins;i++){
		INTEGER(S_em)[idx2c(i,0,full_bins)] = direction[i];
		INTEGER(S_em)[idx2c(i,1,full_bins)] = obs_dist[i];
		INTEGER(S_em)[idx2c(i,2,full_bins)] = obs_start[i];
		INTEGER(S_em)[idx2c(i,3,full_bins)] = obs_end[i];
	}

	SEXP S_em2 = PROTECT(allocMatrix(REALSXP, full_bins, 2));
	for(i=0;i<full_bins;i++){
		REAL(S_em2)[idx2c(i,0,full_bins)] = obs_norm_p[i];
		REAL(S_em2)[idx2c(i,1,full_bins)] = obs_abnorm_p[i];
	}

	/**********/
	const char *names[] = {"obs", "normAbnorm", ""};
	SEXP S_res = PROTECT(mkNamed(VECSXP, names));
	SET_VECTOR_ELT(S_res, 0, S_em);
	SET_VECTOR_ELT(S_res, 1, S_em2);

	UNPROTECT(3);
	/**********/

	// Free memory
	free(direction);
	free(obs_dist);
	free(obs_start);
	free(obs_end);
	free(obs_norm_p);
	free(obs_abnorm_p);
	free(matrix_pos);

	free(obs_o);
	free(tran_p);

	Rprintf("****** Program ended. ******\n");

	return (S_res);
}


/* Helper Functions */

double calculate_sd(double *data, int size_data){
	double sum = 0.0, mean, total = 0.0, ret_val_sd = 0.0;
	int i;

    for(i = 0; i < size_data; ++i)
    {
        sum += data[i];
    }

    mean = sum/size_data;

	for(i = 0; i < size_data; ++i)
        total += pow(data[i] - mean, 2);

    ret_val_sd = sqrt(total/(size_data-1));

	return(ret_val_sd);
}

double function_dnorm(double x, double mean, double sigma)
{
	double dnorm_f;
	dnorm_f = exp(-(pow((x-mean),2))/(2*(pow(sigma,2))))/(sigma*sqrt(2*M_PI));

	return(dnorm_f);
}

double * function_forward(double * f_line, double para_mu_pos, double para_mu_neg, double para_sd0, double para_tao1, double para_tao2, double para_rho, double * para_tran_p, double o_b, int dist_b)
{
	double *ret_val_forward = NULL;
	ret_val_forward = (double *)malloc(sizeof(double)*3);

	double *f_bin_b = NULL;
	f_bin_b = (double *)malloc(sizeof(double)*3);

	double dist_rate = 1 - exp(-para_rho*dist_b);
	int i, j;
	double * tran_a, *f_bin_b_init;

	tran_a = (double *)malloc(sizeof(double)*3*3);
	f_bin_b_init = (double *)malloc(sizeof(double)*3*3);

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_a[idx2c(i,j,3)] = para_tran_p[idx2c(i,j,3)]*dist_rate;
		}
	}

	tran_a[idx2c(0,0,3)] = 1 -(tran_a[idx2c(0,1,3)] + tran_a[idx2c(0,2,3)]);
	tran_a[idx2c(1,1,3)] = 1 -(tran_a[idx2c(1,0,3)] + tran_a[idx2c(1,2,3)]);
	tran_a[idx2c(2,2,3)] = 1 -(tran_a[idx2c(2,0,3)] + tran_a[idx2c(2,1,3)]);

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			f_bin_b_init[idx2c(i,j,3)] = f_line[i]*tran_a[idx2c(i,j,3)];
		}
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			if(j == 0)
				f_bin_b_init[idx2c(i,j,3)] = f_bin_b_init[idx2c(i,j,3)]*function_dnorm(o_b, 0.0, para_sd0);
			if(j == 1)
				f_bin_b_init[idx2c(i,j,3)] = f_bin_b_init[idx2c(i,j,3)]*function_dnorm(o_b, para_mu_pos, para_tao1);
			if(j == 2)
				f_bin_b_init[idx2c(i,j,3)] = f_bin_b_init[idx2c(i,j,3)]*function_dnorm(o_b, para_mu_neg, para_tao2);
		}
	}

	f_bin_b[0] = 0.0;
	f_bin_b[1] = 0.0;
	f_bin_b[2] = 0.0;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			f_bin_b[j] = f_bin_b[j] + f_bin_b_init[idx2c(i,j,3)];
		}
	}

	double which_min = 0.0;
	int which_min_i = 0;

	which_min = (double)1/f_bin_b[0];
	for(i=0;i<3;i++){
		if(which_min > (double)1/f_bin_b[i]){
			which_min = (double)1/f_bin_b[i];
			which_min_i = i;
		}
	}

	for(i=0;i<3;i++){
		ret_val_forward[i] = f_bin_b[i]*((double)1/f_bin_b[which_min_i]);
	}

	//Free memory
	free(tran_a);
	free(f_bin_b_init);
	free(f_bin_b);

	return (ret_val_forward);
}

double * function_backward(double * b_line, double para_mu_pos, double para_mu_neg, double para_sd0, double para_tao1, double para_tao2, double para_rho, double * para_tran_p, double o_b, int dist_b)
{
	double *ret_val_backward = NULL;
	ret_val_backward = (double *)malloc(sizeof(double)*3);

	double *b_bin_b = NULL;
	b_bin_b = (double *)malloc(sizeof(double)*3);

	double dist_rate = 1 - exp(-para_rho*dist_b);

	int i, j;
	double * tran_a, *b_bin_b_init;

	tran_a = (double *)malloc(sizeof(double)*3*3);
	b_bin_b_init = (double *)malloc(sizeof(double)*3*3);

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			tran_a[idx2c(i,j,3)] = para_tran_p[idx2c(i,j,3)]*dist_rate;
		}
	}

	tran_a[idx2c(0,0,3)] = 1 -(tran_a[idx2c(0,1,3)] + tran_a[idx2c(0,2,3)]);
	tran_a[idx2c(1,1,3)] = 1 -(tran_a[idx2c(1,0,3)] + tran_a[idx2c(1,2,3)]);
	tran_a[idx2c(2,2,3)] = 1 -(tran_a[idx2c(2,0,3)] + tran_a[idx2c(2,1,3)]);

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			if(j == 0)
				b_bin_b_init[idx2c(i,j,3)] = tran_a[idx2c(i,j,3)]*function_dnorm(o_b, 0.0, para_sd0);
			if(j == 1)
				b_bin_b_init[idx2c(i,j,3)] = tran_a[idx2c(i,j,3)]*function_dnorm(o_b, para_mu_pos, para_tao1);
			if(j == 2)
				b_bin_b_init[idx2c(i,j,3)] = tran_a[idx2c(i,j,3)]*function_dnorm(o_b, para_mu_neg, para_tao2);
		}
	}

	b_bin_b[0] = 0.0;
	b_bin_b[1] = 0.0;
	b_bin_b[2] = 0.0;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			b_bin_b_init[idx2c(i,j,3)] = b_bin_b_init[idx2c(i,j,3)]*b_line[j];
		}
	}

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			b_bin_b[i] = b_bin_b[i] + b_bin_b_init[idx2c(i,j,3)];
		}
	}

	double which_min = 0.0;
	int which_min_i = 0;

	which_min = (double)1/b_bin_b[0];
	for(i=0;i<3;i++){
		if(which_min > (double)1/b_bin_b[i]){
			which_min = (double)1/b_bin_b[i];
			which_min_i = i;
		}
	}

	for(i=0;i<3;i++){
		ret_val_backward[i] = b_bin_b[i]*((double)1/b_bin_b[which_min_i]);
	}

	//Free memory
	free(tran_a);
	free(b_bin_b_init);
	free(b_bin_b);

	return (ret_val_backward);
}

double function_update_p(double hk_init, int * obs_dist, double * P_kl_all, double rho, int k, int Bins)
{
	double ret_val_update_p;
	ret_val_update_p = 0.0;
	int i = 0;

	double tmp1 = 0.0, term1 = 0.0, term2 = 0.0;
	double sum_pkl_a = 0.0, sum_pkl_b = 0.0;

	int c1, c2, c3;

	if(k == 0){
		c1 = 1;
		c2 = 2;
		c3 = 0;
	}
	if(k == 1){
		c1 = 3;
		c2 = 5;
		c3 = 4;
	}
	if(k == 2){
		c1 = 6;
		c2 = 7;
		c3 = 8;
	}

	for(i=0;i<Bins-1;i++){
		tmp1 = tmp1 + P_kl_all[idx2c(i,c1,Bins-1)] + P_kl_all[idx2c(i,c2,Bins-1)];
	}

	for(i=0;i<Bins-1;i++){
		term1 = term1 + P_kl_all[idx2c(i,c3,Bins-1)]*log(1-tmp1*(1-exp(-rho*obs_dist[i+1]))/(hk_init));
	}

	for(i=0;i<Bins-1;i++){
		sum_pkl_a = sum_pkl_a + P_kl_all[idx2c(i,c1,Bins-1)];
		sum_pkl_b = sum_pkl_b + P_kl_all[idx2c(i,c2,Bins-1)];
	}

	for(i=0;i<Bins-1;i++){
		term2 = term2 + P_kl_all[idx2c(i,c1,Bins-1)]*log(sum_pkl_a*(1-exp(-rho*obs_dist[i+1]))/(hk_init)) + P_kl_all[idx2c(i,c2,Bins-1)]*log(sum_pkl_b*(1-exp(-rho*obs_dist[i+1]))/(hk_init));
	}

	ret_val_update_p = term1 + term2;

	return(-ret_val_update_p);
}

double function_update_rho(double rho_init, int * obs_dist, double * P_kl_all, double * tran_p_new, int Bins)
{
	double ret_val_update_rho;
	ret_val_update_rho = 0.0;
	int i, j;
	int num_Bins = Bins-1;
	double exp_component = 0.0, log_component = 0.0;

	for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			if(tran_p_new[idx2c(i,j,3)] == 0.0){
				tran_p_new[idx2c(i,j,3)] = 0.0001;
			}
		}
	}

	double G2 = 0.0, G3 = 0.0;
	double sum_tran_p_new_a, sum_tran_p_new_b, sum_tran_p_new_c;

	sum_tran_p_new_a = tran_p_new[idx2c(0,1,3)] + tran_p_new[idx2c(0,2,3)];
	sum_tran_p_new_b = tran_p_new[idx2c(1,0,3)] + tran_p_new[idx2c(1,2,3)];
	sum_tran_p_new_c = tran_p_new[idx2c(2,0,3)] + tran_p_new[idx2c(2,1,3)];

	double tran_p_new_loop2_1, tran_p_new_loop2_2, tran_p_new_loop2_3, tran_p_new_loop2_4, tran_p_new_loop2_5, tran_p_new_loop2_6;
	tran_p_new_loop2_1 = log(tran_p_new[idx2c(0,1,3)]);
	tran_p_new_loop2_2 = log(tran_p_new[idx2c(0,2,3)]);
	tran_p_new_loop2_3 = log(tran_p_new[idx2c(1,0,3)]);
	tran_p_new_loop2_4 = log(tran_p_new[idx2c(1,2,3)]);
	tran_p_new_loop2_5 = log(tran_p_new[idx2c(2,0,3)]);
	tran_p_new_loop2_6 = log(tran_p_new[idx2c(2,1,3)]);

	for(i=0;i<num_Bins;i++){

		exp_component = (1-exp(-rho_init*obs_dist[i+1]));
		log_component = log(exp_component);

		G2 += P_kl_all[idx2c(i,0,num_Bins)]*log(1-sum_tran_p_new_a*exp_component) + \
		P_kl_all[idx2c(i,4,num_Bins)]*log(1-sum_tran_p_new_b*exp_component) + \
		P_kl_all[idx2c(i,8,num_Bins)]*log(1-sum_tran_p_new_c*exp_component);

		G3 += P_kl_all[idx2c(i,1,num_Bins)]*(tran_p_new_loop2_1+log_component) + \
		P_kl_all[idx2c(i,2,num_Bins)]*(tran_p_new_loop2_2+log_component) + \
		P_kl_all[idx2c(i,3,num_Bins)]*(tran_p_new_loop2_3+log_component) + \
		P_kl_all[idx2c(i,5,num_Bins)]*(tran_p_new_loop2_4+log_component) + \
		P_kl_all[idx2c(i,6,num_Bins)]*(tran_p_new_loop2_5+log_component) + \
		P_kl_all[idx2c(i,7,num_Bins)]*(tran_p_new_loop2_6+log_component);
	}

	ret_val_update_rho = G2 + G3;

	return(ret_val_update_rho);
}

/*
Comment:
 The function Brent_fmin is the function fmin copied straight from
 'optimize.c' in the 'stats' package of R's core installation,
 with a few modifications in the names of the input arguments.

 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 2003-2004  The R Foundation
 *  Copyright (C) 1998--2014  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/


 R's  optimize() :   function	fmin(ax,bx,f,tol)
   =    ==========		~~~~~~~~~~~~~~~~~

        an approximation  x  to the point where  f  attains a minimum  on
    the interval  (ax,bx)  is determined.

    INPUT..

    ax    left endpoint of initial interval
    bx    right endpoint of initial interval
    f     function which evaluates  f(x, info)  for any  x
          in the interval  (ax,bx)
    tol   desired length of the interval of uncertainty of the final
          result ( >= 0.)

    OUTPUT..

    fmin  abcissa approximating the point where  f  attains a minimum
*/
static double Brent_fmin(double ax, double bx, double tol, double f(double x, int * obs_dist, double * P_kl_all, double rho, int hmmDMR_k, int Bins), double *x, int * obs_dist, double * P_kl_all, double rho, int hmmDMR_k, int Bins)
{
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w; //, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    (*x) = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = f((*x), obs_dist, P_kl_all, rho, hmmDMR_k, Bins);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs((*x)) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs((*x) - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */

	    r = ((*x) - w) * (fx - fv);
	    q = ((*x) - v) * (fx - fw);
	    p = ((*x) - v) * q - ((*x) - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - (*x)) || p >= q * (b - (*x))) { /* a golden-section step */

	    if ((*x) < xm) e = b - (*x); else e = a - (*x);
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */

	    d = p / q;
	    u = (*x) + d;

	    /* f must not be evaluated too close to ax or bx */

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if ((*x) >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d) >= tol1)
	    u = (*x) + d;
	else if (d > 0.)
	    u = (*x) + tol1;
	else
	    u = (*x) - tol1;

	fu = f(u, obs_dist, P_kl_all, rho, hmmDMR_k, Bins);

	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < (*x)) b = (*x); else a = (*x);
	    v = w;    w = (*x);   (*x) = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < (*x)) a = u; else b = u;
	    if (fu <= fw || w == (*x)) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == (*x) || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    return (*x);
}


