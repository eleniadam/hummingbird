#include "utils.h"

/* hummingbirdPostAdjustment: Post Adjustment algorithm for the output of the EM.
 * This function adjusts HMM output such that each detected DMR
 * has a minimum length and a minimum number of CpGs in each DMR.
 */
// [[Rcpp::export]]
SEXP hummingbirdPostAdjustment(SEXP em, SEXP pos, SEXP minCpGs, SEXP minLength, SEXP maxGap)
{
	Rprintf("Reading input...\n");

	// User defined values
	int min_CpGs = 0, min_length = 0, max_gap = 0;

	// minimum number of CpGs
	if(TYPEOF(minCpGs) == NILSXP){
		min_CpGs = 10; //default
	}
	else{
		min_CpGs = asReal(minCpGs);
	}

	// minimum length
		if(TYPEOF(minLength) == NILSXP){
		min_length = 500; //default
	}
	else{
		min_length = asReal(minLength);
	}

	// maximum gap
	if(TYPEOF(maxGap) == NILSXP){
		max_gap = 300; //default
	}
	else{
		max_gap = asReal(maxGap);
	}

	Rprintf("Min CpGs: %d, Min Length: %d, Max gap: %d.\n", min_CpGs, min_length, max_gap);

	// Matrix EM
	SEXP S_dim_em = getAttrib(em, R_DimSymbol);
	int full_bins = INTEGER(S_dim_em)[0];

	int *direction = NULL, *obs_dist = NULL, *obs_start = NULL, *obs_end = NULL;
	direction = (int *)malloc(sizeof(int)*full_bins);
	obs_dist = (int *)malloc(sizeof(int)*full_bins);
	obs_start = (int *)malloc(sizeof(int)*full_bins);
	obs_end = (int *)malloc(sizeof(int)*full_bins);

	int i;
	for(i=0;i<full_bins;i++){
		direction[i] = INTEGER(em)[idx2c(i,0,full_bins)];
		obs_dist[i] = INTEGER(em)[idx2c(i,1,full_bins)];
		obs_start[i] = INTEGER(em)[idx2c(i,2,full_bins)];
		obs_end[i] = INTEGER(em)[idx2c(i,3,full_bins)];
	}

	// Matrix Pos
	SEXP S_dim_pos = getAttrib(pos, R_DimSymbol);
	int total_lines = INTEGER(S_dim_pos)[0];

	int *matrix_pos;
	matrix_pos = (int *)malloc(sizeof(int)*total_lines);

	for(i=0;i<total_lines;i++){
		matrix_pos[i] = INTEGER(pos)[i];
	}

	/****************************************************************/

	Rprintf("Post Adjustment begins...\n");

	int flag = 0, region_cnt = 0;
	int tmp_direction = 0, start_pos_adj = 0, end_pos_adj = 0;
	int start_bin = 0, end_bin = 0, bin_cnt = 0;
	int DMR_length = 0, num_CpGs = 0, region_length = 0;
	char region_state[10];

	int *mat_region = NULL, *mat_start_pos = NULL, *mat_end_pos = NULL, *mat_length = NULL, *mat_num_CpGs = NULL, *int_mat_state = NULL;

	for(bin_cnt=0; bin_cnt<full_bins; bin_cnt++){

		if((direction[bin_cnt] != 0) && (flag == 0)){
			flag = 1;
			tmp_direction = direction[bin_cnt];
			start_pos_adj = obs_start[bin_cnt];
			start_bin = bin_cnt;
		}
		else if ((direction[bin_cnt] != 0) && (flag == 1)){

			if((tmp_direction != direction[bin_cnt]) || (tmp_direction == direction[bin_cnt] && obs_dist[bin_cnt] > max_gap) ){

				end_pos_adj = obs_end[bin_cnt-1];
				end_bin = bin_cnt-1;
				DMR_length = end_pos_adj - start_pos_adj;
				num_CpGs = 0;
				for(i=0; i<total_lines; i++){
					if(matrix_pos[i] >= start_pos_adj && matrix_pos[i] <= end_pos_adj){
						num_CpGs = num_CpGs + 1;
					}
				}

				if(DMR_length < min_length || num_CpGs < min_CpGs){
					for(i=start_bin; i<=end_bin; i++){
						direction[i] = 0;
					}
				}
				else{
					region_cnt = region_cnt + 1;
					if(tmp_direction == 1){
						strcpy(region_state, "hyper");
					}
					if(tmp_direction == 2){
						strcpy(region_state, "hypo");
					}
					region_length = end_pos_adj - start_pos_adj;

					mat_region = (int*) realloc (mat_region, region_cnt * sizeof(int));
					mat_start_pos = (int*) realloc (mat_start_pos, region_cnt * sizeof(int));
					mat_end_pos = (int*) realloc (mat_end_pos, region_cnt * sizeof(int));
					int_mat_state = (int*) realloc (int_mat_state, region_cnt * sizeof(int));
					mat_length = (int*) realloc (mat_length, region_cnt * sizeof(int));
					mat_num_CpGs = (int*) realloc (mat_num_CpGs, region_cnt * sizeof(int));

					mat_start_pos[region_cnt-1] = start_pos_adj;
					mat_end_pos[region_cnt-1] = end_pos_adj;
					mat_length[region_cnt-1] = region_length + 1;
					mat_region[region_cnt-1] = region_cnt;
					mat_num_CpGs[region_cnt-1] = num_CpGs;
					int_mat_state[region_cnt-1] = tmp_direction;

				}
				start_pos_adj = obs_start[bin_cnt];
				start_bin = bin_cnt;
				tmp_direction = direction[bin_cnt];

			}
			if(tmp_direction == direction[bin_cnt] && obs_dist[bin_cnt] <= max_gap && bin_cnt == full_bins-1){

				end_pos_adj = obs_end[bin_cnt];
				end_bin = bin_cnt;
				DMR_length = end_pos_adj - start_pos_adj;
				num_CpGs = 0;
				for(i=0; i<total_lines; i++){
					if(matrix_pos[i] >= start_pos_adj && matrix_pos[i] <= end_pos_adj){
						num_CpGs = num_CpGs + 1;
					}
				}

				if(DMR_length < min_length || num_CpGs < min_CpGs){
					for(i=start_bin; i<=end_bin; i++){
						direction[i] = 0;
					}

				}
				else{
					region_cnt = region_cnt + 1;
					if(tmp_direction == 1){
						strcpy(region_state, "hyper");
					}
					if(tmp_direction == 2){
						strcpy(region_state, "hypo");
					}
					region_length = end_pos_adj - start_pos_adj;

					mat_region = (int*) realloc (mat_region, region_cnt * sizeof(int));
					mat_start_pos = (int*) realloc (mat_start_pos, region_cnt * sizeof(int));
					mat_end_pos = (int*) realloc (mat_end_pos, region_cnt * sizeof(int));
					int_mat_state = (int*) realloc (int_mat_state, region_cnt * sizeof(int));
					mat_length = (int*) realloc (mat_length, region_cnt * sizeof(int));
					mat_num_CpGs = (int*) realloc (mat_num_CpGs, region_cnt * sizeof(int));

					mat_start_pos[region_cnt-1] = start_pos_adj;
					mat_end_pos[region_cnt-1] = end_pos_adj;
					mat_length[region_cnt-1] = region_length + 1;
					mat_region[region_cnt-1] = region_cnt;
					mat_num_CpGs[region_cnt-1] = num_CpGs;
					int_mat_state[region_cnt-1] = tmp_direction;

				}

			}

		}
		else if((direction[bin_cnt] == 0) && (flag == 1)){
			flag = 0;
			end_pos_adj = obs_end[bin_cnt-1];
			end_bin = bin_cnt - 1;
			DMR_length = end_pos_adj - start_pos_adj;
			num_CpGs = 0;
			for(i=0; i<total_lines; i++){
				if(matrix_pos[i] >= start_pos_adj && matrix_pos[i] <= end_pos_adj){
					num_CpGs = num_CpGs + 1;
				}
			}

			if(DMR_length < min_length || num_CpGs < min_CpGs){
				for(i=start_bin; i<=end_bin; i++){
					direction[i] = 0;
				}
			}
			else{
				region_cnt = region_cnt + 1;
				if(tmp_direction == 1){
					strcpy(region_state, "hyper");
				}
				if(tmp_direction == 2){
					strcpy(region_state, "hypo");
				}
				region_length = end_pos_adj - start_pos_adj;

				mat_region = (int*) realloc (mat_region, region_cnt * sizeof(int));
				mat_start_pos = (int*) realloc (mat_start_pos, region_cnt * sizeof(int));
				mat_end_pos = (int*) realloc (mat_end_pos, region_cnt * sizeof(int));
				int_mat_state = (int*) realloc (int_mat_state, region_cnt * sizeof(int));
				mat_length = (int*) realloc (mat_length, region_cnt * sizeof(int));
				mat_num_CpGs = (int*) realloc (mat_num_CpGs, region_cnt * sizeof(int));

				mat_start_pos[region_cnt-1] = start_pos_adj;
				mat_end_pos[region_cnt-1] = end_pos_adj;
				mat_length[region_cnt-1] = region_length + 1;
				mat_region[region_cnt-1] = region_cnt;
				mat_num_CpGs[region_cnt-1] = num_CpGs;
				int_mat_state[region_cnt-1] = tmp_direction;

			}
		}

	}

	Rprintf("Post Adjustment completed...\n");

	//Display output
	Rprintf("Output DMRs...\n");

	int display_dmrs = 10;
	if(region_cnt < display_dmrs){
		display_dmrs = region_cnt;
	}

	if(region_cnt == 0){
		Rprintf("There are %d DMRs in total. Please try the Post Adjustment algorithm with different input values.\n", region_cnt);
	}
	else{
		Rprintf("There are %d DMRs in total. The first %d are displayed.\n", region_cnt, display_dmrs);
		Rprintf("Region: Start, End, Length, Direction, CpGs\n");
		for(i=0;i<display_dmrs;i++){
			Rprintf("%d: %d, %d, %d, %d, %d\n", mat_region[i], mat_start_pos[i], mat_end_pos[i], mat_length[i], int_mat_state[i], mat_num_CpGs[i]);
			if(i == 9){
				break;
			}
		}

	}

	Rprintf("Saving output...\n");

	SEXP S_postadj = PROTECT(allocMatrix(INTSXP, full_bins, 4));
	for(i=0;i<full_bins;i++){
		INTEGER(S_postadj)[idx2c(i,0,full_bins)] = direction[i];
		INTEGER(S_postadj)[idx2c(i,1,full_bins)] = obs_dist[i];
		INTEGER(S_postadj)[idx2c(i,2,full_bins)] = obs_start[i];
		INTEGER(S_postadj)[idx2c(i,3,full_bins)] = obs_end[i];
	}

	SEXP S_dmrs = PROTECT(allocMatrix(INTSXP, region_cnt, 6));
	for(i=0;i<region_cnt;i++){
		INTEGER(S_dmrs)[idx2c(i,0,region_cnt)] = mat_region[i];
		INTEGER(S_dmrs)[idx2c(i,1,region_cnt)] = mat_start_pos[i];
		INTEGER(S_dmrs)[idx2c(i,2,region_cnt)] = mat_end_pos[i];
		INTEGER(S_dmrs)[idx2c(i,3,region_cnt)] = mat_length[i];
		INTEGER(S_dmrs)[idx2c(i,4,region_cnt)] = int_mat_state[i];
		INTEGER(S_dmrs)[idx2c(i,5,region_cnt)] = mat_num_CpGs[i];
	}

	const char *names[] = {"obsPostAdj", "DMRs", ""};
	SEXP S_res = PROTECT(mkNamed(VECSXP, names));
	SET_VECTOR_ELT(S_res, 0, S_postadj);
	SET_VECTOR_ELT(S_res, 1, S_dmrs);

	UNPROTECT(3);

	// Free memory
	free(direction);
	free(obs_dist);
	free(obs_start);
	free(obs_end);
	free(matrix_pos);

	free(mat_region);
	free(mat_start_pos);
	free(mat_end_pos);
	free(int_mat_state);
	free(mat_length);
	free(mat_num_CpGs);

	Rprintf("****** Program ended. ******\n");

	return (S_res);
}
