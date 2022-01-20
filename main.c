#include <stdio.h>
#include "stdlib.h"
#include "Math.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <time.h>

#define POS_ROW 201
#define POS_COL 1

#define VEL_ROW 201
#define VEL_COL 1

#define POSHAT_ROW 201
#define POSHAT_COL 1

#define VELHAT_ROW 201
#define VELHAT_COL 1

#define POSHATINF_ROW 201
#define POSHATINF_COL 1

#define VELHATINF_ROW 201
#define VELHATINF_COL 1

#define HINFGAINS_ROW 2
#define HINFGAINS_COL 201

#define KALMANGAINS_ROW 2
#define KALMANGAINS_COL 201


double gama = 0.01;
double dt = 0.1;
double duration = 60.0;
double Steadystate = 0.0;

#pragma region Inverse of Matrix

gsl_matrix * invert_a_matrix(gsl_matrix *matrix, size_t size)
{
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}

#pragma endregion

void print_mat_contents(gsl_matrix *matrix, size_t size)
{
    size_t i, j;
    double element;

    for (i = 0; i < size; ++i) {
        for (j = 0; j < size; ++j) {
            element = gsl_matrix_get(matrix, i, j);
            printf("%f ", element);
        }
        printf("\n");
    }
}

double rand_generator(){
	double random_value;

    srand ( time ( NULL));

    random_value = (double)rand()/RAND_MAX*2.0-1.0;//float in range -1 to 1

	return random_value;
}

#pragma region Main Program

void main()
{
	gsl_matrix *pos = gsl_matrix_alloc(POS_ROW, POS_COL);
	gsl_matrix *vel = gsl_matrix_alloc(VEL_ROW, VEL_COL);
	gsl_matrix *velhat = gsl_matrix_alloc(VELHAT_ROW, VELHAT_COL);
	gsl_matrix *poshat = gsl_matrix_alloc(POSHAT_ROW, POSHAT_COL);
	gsl_matrix *poshatinf = gsl_matrix_alloc(POSHAT_ROW, POSHAT_COL);
	gsl_matrix *velhatinf = gsl_matrix_alloc(VELHAT_ROW, VELHAT_COL);
	gsl_matrix *hinfgains = gsl_matrix_alloc(HINFGAINS_ROW, HINFGAINS_COL);
	gsl_matrix *kalmangains = gsl_matrix_alloc(KALMANGAINS_ROW, KALMANGAINS_COL);

	HinfinityInC(gama, dt, duration, Steadystate, pos, vel, poshat, velhat,
				poshatinf, velhatinf, hinfgains, kalmangains);
}
#pragma endregion

void HinfinityInC(double gama,double dt,double duration,double Steadystate,
					gsl_matrix *pos,
					gsl_matrix *vel,
					gsl_matrix *poshat,
					gsl_matrix *velhat,
					gsl_matrix *poshatinf,
					gsl_matrix *velhatinf,
					gsl_matrix *hinfgains,
					gsl_matrix *kalmangains)
{
	#pragma region Initialization
	/*For loop counter declarations */
	int count;

	//double K[2][1] = {0,0};
	//double L[2][2] = {0,0,0,0};
	double measnoise = 2.0;  // nominal noise velocity measurement noise
	double accelnoise = 0.2; // nominal acceleration noise

	double a_d[] = { 	1.0, dt,
                 	0.0, 1.0};

	double b_d[] = { 	pow(dt,2.0)/2.0,
						dt};

	double c_d[] = {  	0.0,
						1.0};

	double x_d[] = {  	0.0,
						0.0};

	double y_d[] = {0.0};

	gsl_matrix_view a = gsl_matrix_view_array(a_d, 2, 2); 	// transition matrix
	gsl_matrix_view b = gsl_matrix_view_array(b_d, 2, 1);		// Input matrix
	gsl_matrix_view c = gsl_matrix_view_array(c_d, 1, 2);		// measurement matrix
	gsl_matrix_view x = gsl_matrix_view_array(x_d, 2, 1);		// initial state vector
	gsl_matrix_view y = gsl_matrix_view_array(y_d, 1, 1);		// initial state vector

	/*y = c * x */
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		1.0, &c.matrix, &x.matrix,
		0.0, &y.matrix);							//Initial measurement

	double identity_d[] = {  	1.0, 0.0,
								0.0, 1.0};
	gsl_matrix_view identity = gsl_matrix_view_array(identity_d, 2, 2);

	/*
		Initialize Kalman filter variables
	*/
	double xhat_d[]  = {0.0, 0.0};
	gsl_matrix_view xhat = gsl_matrix_view_array(xhat_d, 2, 1);

	double Sz_d[]  = {pow(measnoise,2)};
	gsl_matrix_view Sz = gsl_matrix_view_array(Sz_d, 1, 1);

	double Sw_d[] = { ((pow(accelnoise,2))*(pow(dt,4))/4.0),  ((pow(accelnoise,2))*(pow(dt,3))/2.0),
						((pow(accelnoise,2))*(pow(dt,3))/2.0), (pow(accelnoise,2))*(pow(dt,2))};
	gsl_matrix_view Sw = gsl_matrix_view_array(Sw_d, 2, 2);

	double P_d[] = { 	0.0, 0.0,
						0.0, 0.0};
	gsl_matrix_view P = gsl_matrix_view_array(P_d, 2, 2);

	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		1.0, &Sw.matrix, &identity.matrix,
		0.0, &P.matrix);

	double xhatinf_d[] = {  	0.0,
								0.0};
	gsl_matrix_view xhatinf = gsl_matrix_view_array(xhatinf_d, 2, 1);

	double Pinf_d[] = {  	0.01, 0.0,
							0.0, 0.01};
	gsl_matrix_view Pinf = gsl_matrix_view_array(Pinf_d, 2, 2);

	double W_d[] = {  	0.0003/1000.0, 0.0050/1000.0,
						0.0050/1000.0, 0.1000/1000.0};
	gsl_matrix_view W = gsl_matrix_view_array(W_d, 2, 2);

	double V_d[] = {0.01};
	gsl_matrix_view V = gsl_matrix_view_array(V_d, 1, 1);

	double Q_d[] = {  	0.01, 0.0,
						0.0, 0.01};
	gsl_matrix_view Q = gsl_matrix_view_array(Q_d, 2, 2);

	/* Initialize the Pos value */
	for (int i = 0; i < POS_ROW; i++)
		for (int j = 0; j < POS_COL; j++)
			gsl_matrix_set (pos, i, j, x_d[0]);

	/* Initialize the Vel value */
	for (int i = 0; i < VEL_ROW; i++)
		for (int j = 0; j < VEL_COL; j++)
			gsl_matrix_set (vel, i, j, x_d[1]);

	/* Initialize the Poshat value */
	for (int i = 0; i < POSHAT_ROW; i++)
		for (int j = 0; j < POSHAT_COL; j++)
			gsl_matrix_set (poshat, i, j, xhat_d[0]);

	/* Initialize the Velhat value */
	for (int i = 0; i < VELHAT_ROW; i++)
		for (int j = 0; j < VELHAT_COL; j++)
			gsl_matrix_set (velhat, i, j, xhat_d[1]);

	/* Initialize the Poshatinf value */
	for (int i = 0; i < POSHATINF_ROW; i++)
		for (int j = 0; j < POSHATINF_COL; j++)
			gsl_matrix_set (poshatinf, i, j, xhatinf_d[0]);

	/* Initialize the Velhatinf value */
	for (int i = 0; i < VELHAT_ROW; i++)
		for (int j = 0; j < VELHAT_COL; j++)
			gsl_matrix_set (velhatinf, i, j, xhatinf_d[1]);

	/* Initialize the hinfgains value */
	for (int i = 0; i < HINFGAINS_ROW; i++)
		for (int j = 0; j < HINFGAINS_COL; j++)
			gsl_matrix_set (hinfgains, i, j, 0.0);

	/* Initialize the Kalmangains value */
	for (int i = 0; i < KALMANGAINS_ROW; i++)
		for (int j = 0; j < KALMANGAINS_COL; j++)
			gsl_matrix_set (kalmangains, i, j, 0.0);

	int dt_t = (int)(dt*10.0);
	int dur_t = (int)(duration*10.0);
	#pragma endregion

	for(count=0; count<dur_t; count+=dt){
		// use a constant commanded acceleration of 1 foot/sec^2
		double u_d[] = {1.0};
		double k_d[] = {0.0,0.0};
		double l_d[] = {0.0,0.0,0.0,0.0};

		gsl_matrix_view k = gsl_matrix_view_array(k_d, 2, 1);		// Use steady-state H-infinity gains
		gsl_matrix_view l = gsl_matrix_view_array(l_d, 2, 2);
		gsl_matrix_view u = gsl_matrix_view_array(u_d, 1, 1);

		// figure out the H-Infinity gains
		if(Steadystate == 1.0){
			// Use Steady State H-Infinity Gain
			gsl_matrix_set(&k.matrix, 0,0,0.11);
			gsl_matrix_set(&k.matrix, 1,0,0.09);
		}
		else{
			size_t size = 1;
			#pragma region " L = inv(eye(2) - g * Q * Pinf + c' * inv(V) * c * Pinf)"
			#pragma region Part 3 of equation

			double c_transpose_d[] = {  	0.0,
											0.0};
			gsl_matrix_view c_transpose = gsl_matrix_view_array(c_transpose_d, 2, 1);
			gsl_matrix_transpose_memcpy(&c_transpose.matrix, &c.matrix);

			gsl_matrix *inverse_of_V = invert_a_matrix(&V.matrix, size);

			double c_times_Pinf_d[] = { 0.00,
										0.00};

			gsl_matrix_view c_times_Pinf = gsl_matrix_view_array(c_times_Pinf_d, 1, 2);

			double invV_c_time_Pinf_d[] = { 0.00,
										0.00};

			gsl_matrix_view invV_c_time_Pinf = gsl_matrix_view_array(invV_c_time_Pinf_d, 1, 2);

			double till_c_transpose_d[] = { 	0.0, 0.0,
                 								0.0, 0.0};

			gsl_matrix_view till_c_transpose = gsl_matrix_view_array(till_c_transpose_d, 2, 2);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c.matrix, &Pinf.matrix,
							0.0, &c_times_Pinf.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, inverse_of_V, &c_times_Pinf.matrix,
							0.0, &invV_c_time_Pinf.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c_transpose.matrix, &invV_c_time_Pinf.matrix,
							0.0, &till_c_transpose.matrix);
			#pragma endregion
			#pragma region Part 2 of equation

			double part_2_eq_d[] = {	0.0, 0.0,
										0.0, 0.0};

			gsl_matrix_view part_2_eq = gsl_matrix_view_array(part_2_eq_d, 2, 2);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							gama, &Q.matrix, &Pinf.matrix,
							0.0, &part_2_eq.matrix);

			#pragma endregion
			#pragma region Part 1 of equation
			double part_1_eq_d[] = {	1.0, 0.0,
										0.0, 1.0};

			gsl_matrix_view part_1_eq = gsl_matrix_view_array(part_1_eq_d, 2, 2);
			#pragma endregion

			double temp_sum_d[] = 	{	0.0, 0.0,
										0.0, 0.0};

			gsl_matrix_view temp_sum = gsl_matrix_view_array(temp_sum_d, 2, 2);
			gsl_matrix_add(&temp_sum.matrix, &part_1_eq.matrix);
			gsl_matrix_sub(&temp_sum.matrix, &part_2_eq.matrix);
			gsl_matrix_add(&temp_sum.matrix, &till_c_transpose.matrix);

			gsl_matrix_add(&l.matrix, invert_a_matrix(&temp_sum.matrix, 2));

			#pragma endregion

			#pragma region "K = a * Pinf * L * c' * inv(V);"

			double c_time_invV_d[] = 		{ 	0.0,
                 								0.0};

			gsl_matrix_view c_time_invV = gsl_matrix_view_array(c_time_invV_d, 2, 1);

			double L_times_c_time_invV_d[] = { 	0.0,
                 								0.0};

			gsl_matrix_view L_times_c_time_invV = gsl_matrix_view_array(L_times_c_time_invV_d, 2, 1);

			double Pinf_L_times_c_time_invV_d[] = { 	0.0,
                 										0.0};

			gsl_matrix_view Pinf_L_times_c_time_invV = gsl_matrix_view_array(Pinf_L_times_c_time_invV_d, 2, 1);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c_transpose.matrix, inverse_of_V,
							0.0, &c_time_invV.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &l.matrix, &c_time_invV.matrix,
							0.0, &L_times_c_time_invV.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &Pinf.matrix, &L_times_c_time_invV.matrix,
							0.0, &Pinf_L_times_c_time_invV.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &a.matrix, &Pinf_L_times_c_time_invV.matrix,
							0.0, &k.matrix);

			#pragma endregion

			#pragma region "Pinf = a * Pinf * L * a' + W;"

			double a_transpose_d[] = {  	0.0, 0.0,
											0.0, 0.0};
			gsl_matrix_view a_transpose = gsl_matrix_view_array(a_transpose_d, 2, 2);

			gsl_matrix_transpose_memcpy(&a_transpose.matrix, &a.matrix);

			double L_times_a_trans_d[] = {  	0.0, 0.0,
												0.0, 0.0};
			gsl_matrix_view L_times_a_trans = gsl_matrix_view_array(L_times_a_trans_d, 2, 2);

			double Pinf_L_times_a_trans_d[] = {  	0.0, 0.0,
													0.0, 0.0};
			gsl_matrix_view Pinf_L_times_a_trans = gsl_matrix_view_array(Pinf_L_times_a_trans_d, 2, 2);

			double a_Pinf_L_times_a_trans_d[] = {  	0.0, 0.0,
													0.0, 0.0};
			gsl_matrix_view a_Pinf_L_times_a_trans = gsl_matrix_view_array(a_Pinf_L_times_a_trans_d, 2, 2);

			double Pinf_temp_d[] = {  	0.0, 0.0,
										0.0, 0.0};
			gsl_matrix_view Pinf_temp = gsl_matrix_view_array(Pinf_temp_d, 2, 2);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &l.matrix, &a_transpose.matrix,
							0.0, &L_times_a_trans.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &Pinf.matrix, &L_times_a_trans.matrix,
							0.0, &Pinf_L_times_a_trans.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &a.matrix, &Pinf_L_times_a_trans.matrix,
							0.0, &a_Pinf_L_times_a_trans.matrix);

			gsl_matrix_add(&Pinf_temp.matrix, &a_Pinf_L_times_a_trans.matrix);
			gsl_matrix_add(&Pinf_temp.matrix, &W.matrix);
			
			gsl_matrix_memcpy(&Pinf.matrix, &Pinf_temp.matrix);

			#pragma endregion

			// Force Pinf to be symmetric

			#pragma region "Pinf = (Pinf + Pinf') / 2;"
			double Pinf_transpose_d[] = {  	0.0, 0.0,
											0.0, 0.0};
			gsl_matrix_view Pinf_transpose = gsl_matrix_view_array(Pinf_transpose_d, 2, 2);

			gsl_matrix_transpose_memcpy(&Pinf_transpose.matrix, &Pinf.matrix);

			double Pinf_sum_d[] = {  	0.0, 0.0,
										0.0, 0.0};
			gsl_matrix_view Pinf_sum = gsl_matrix_view_array(Pinf_sum_d, 2, 2);


			gsl_matrix_add(&Pinf_sum.matrix, &Pinf.matrix);
			gsl_matrix_add(&Pinf_sum.matrix, &Pinf_transpose.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							0.5, &Pinf_sum.matrix, &identity.matrix,
							0.0, &Pinf.matrix);
			#pragma endregion

			// Make sure the eigenvalues of Pinf are less than 1 in magnitude

			#pragma region "lambda = eig(Pinf);  if (abs(lambda(1)) >= 1) | (abs(lambda(2)) >= 1)  disp('gamma is too large');return; end"
			gsl_vector *eval = gsl_vector_alloc (2);
  			gsl_matrix *evec = gsl_matrix_alloc (2, 2);

			gsl_eigen_symmv_workspace * workspace_eigen = gsl_eigen_symmv_alloc (2);

			gsl_eigen_symmv (&Pinf.matrix, eval, evec, workspace_eigen);
			gsl_eigen_symmv_free (workspace_eigen);
			gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_ASC);

			double eval_0_i = gsl_vector_get (eval, 0);
			double eval_1_i = gsl_vector_get (eval, 1);
			
			if(fabs(eval_0_i) >= 1.0 || fabs(eval_1_i) >= 1.0){
				printf("gamma value is too large");
				return;
			}
			#pragma endregion
		}

		

		#pragma region "xhatinf = a * xhatinf + b * u + K * (y - c * xhatinf);"

		double y_minus_c_times_xhatinf_d[] = {0.0};
		gsl_matrix_view y_minus_c_times_xhatinf = gsl_matrix_view_array(y_minus_c_times_xhatinf_d, 1, 1);

		double c_times_xhatinf_d[] = {0.0};
		gsl_matrix_view c_times_xhatinf = gsl_matrix_view_array(c_times_xhatinf_d, 1, 1);

		double K_y_minus_c_times_xhatinf_d[] = {0.0, 0.0};
		gsl_matrix_view K_y_minus_c_times_xhatinf = gsl_matrix_view_array(K_y_minus_c_times_xhatinf_d, 2, 1);

		double b_times_u_d[] = {0.0, 0.0};
		gsl_matrix_view b_times_u = gsl_matrix_view_array(b_times_u_d, 2, 1);

		double a_times_xhatinf_d[] = {0.0, 0.0};
		gsl_matrix_view a_times_xhatinf = gsl_matrix_view_array(a_times_xhatinf_d, 2, 1);

		double xhatinfsum_d[] = {0.0, 0.0};
		gsl_matrix_view xhatinfsum = gsl_matrix_view_array(xhatinfsum_d, 2, 1);


		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c.matrix, &xhatinf.matrix,
							0.0, &c_times_xhatinf.matrix);

		gsl_matrix_add(&y_minus_c_times_xhatinf.matrix, &y.matrix);
		gsl_matrix_sub(&y_minus_c_times_xhatinf.matrix, &c_times_xhatinf.matrix);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &k.matrix, &y_minus_c_times_xhatinf.matrix,
							0.0, &K_y_minus_c_times_xhatinf.matrix);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &b.matrix, &u.matrix,
							0.0, &b_times_u.matrix);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &a.matrix, &xhatinf.matrix,
							0.0, &a_times_xhatinf.matrix);

		gsl_matrix_add(&xhatinfsum.matrix, &a_times_xhatinf.matrix);
		gsl_matrix_add(&xhatinfsum.matrix, &b_times_u.matrix);
		gsl_matrix_add(&xhatinfsum.matrix, &K_y_minus_c_times_xhatinf.matrix);

		for(int i=0; i<2;i++)
			for(int j=0;j<1;j++)
				gsl_matrix_set(&xhatinf.matrix, i,j, 0.0);

		gsl_matrix_add(&xhatinf.matrix, &xhatinfsum.matrix);
		#pragma endregion

		
		gsl_matrix_set(hinfgains, 0, count, k_d[0]);
		gsl_matrix_set(hinfgains, 1, count, k_d[1]);
		
		gsl_matrix_set(poshatinf, count, 0, xhatinf_d[0]);

		gsl_matrix_set(velhatinf, count, 0, xhatinf_d[1]);

		//printf("Hinf(0) : %f \n", gsl_matrix_get(hinfgains, 0, count));
		//printf("Hinf(1) : %f \n", gsl_matrix_get(hinfgains, 1, count));

		//printf("poshatinf(0) : %f \n", gsl_matrix_get(poshatinf, count, 0));

		

		
		// Simulate the linear system and noisy measurement
		// Note that rand functino in Matlab's gaussian (normal) random number
		// generator; rand is matlab's uniform random number generator
		#pragma region " ProcessNoise = 2 * accelnoise * b .* [randn; randn];"
		double rand_val1 = rand_generator();
		double rand_val2 = rand_generator();
		double random_number_d[] = {rand_val1, rand_val2};
		gsl_matrix_view random_number = gsl_matrix_view_array(random_number_d, 2, 1);

		double ProcessNoise_d[] = {0.0, 0.0};
		gsl_matrix_view ProcessNoise = gsl_matrix_view_array(ProcessNoise_d, 2, 1);

		gsl_matrix_set(&ProcessNoise.matrix, 0, 0, 2 * accelnoise* b_d[0] * random_number_d[0]);
		gsl_matrix_set(&ProcessNoise.matrix, 1, 0, 2 * accelnoise* b_d[1] * random_number_d[1]);
		#pragma endregion

		#pragma region x = a * x + b * u + ProcessNoise;

		double a_times_x_d[] = {0.0, 0.0};
		gsl_matrix_view a_times_x = gsl_matrix_view_array(a_times_x_d, 2, 1);

		double b_times_u_2_d[] = {0.0, 0.0};
		gsl_matrix_view b_times_u_2 = gsl_matrix_view_array(b_times_u_2_d, 2, 1);

		double x_updated_d[] = {0.0, 0.0};
		gsl_matrix_view x_updated = gsl_matrix_view_array(x_updated_d, 2, 1);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &a.matrix, &x.matrix,
							0.0, &a_times_x.matrix);
		
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &b.matrix, &u.matrix,
							0.0, &b_times_u_2.matrix);

		gsl_matrix_add(&x_updated.matrix, &a_times_x.matrix);
		gsl_matrix_add(&x_updated.matrix, &b_times_u_2.matrix);
		gsl_matrix_add(&x_updated.matrix, &ProcessNoise.matrix);

		for(int i=0; i<2;i++)
			for(int j=0;j<1;j++)
				gsl_matrix_set(&x.matrix, i,j, 0.0);

		gsl_matrix_add(&x.matrix, &x_updated.matrix);

		#pragma endregion

		double Measnoise_d[] = {measnoise * (rand_generator() - 0.5)};
		gsl_matrix_view Measnoise = gsl_matrix_view_array(Measnoise_d, 1, 1);

		double c_times_x_d[] = {0.0};
		gsl_matrix_view c_times_x = gsl_matrix_view_array(c_times_x_d, 1, 1);

		double y_updated_d[] = {0.0};
		gsl_matrix_view y_updated = gsl_matrix_view_array(y_updated_d, 1, 1);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c.matrix, &x.matrix,
							0.0, &c_times_x.matrix);

		gsl_matrix_add(&y_updated.matrix, &c_times_x.matrix);
		gsl_matrix_add(&y_updated.matrix, &Measnoise.matrix);

		// Compute the Kalman Filter estimate
		// Exterpolate the most recent state estimate to the present time

		double b_times_u_3_d[] = {0.0, 0.0};
		gsl_matrix_view b_times_u_3 = gsl_matrix_view_array(b_times_u_3_d, 2, 1);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &b.matrix, &u.matrix,
							0.0, &b_times_u_3.matrix);

		double a_times_xhat_d[] = {0.0, 0.0};
		gsl_matrix_view a_times_xhat = gsl_matrix_view_array(a_times_xhat_d, 2, 1);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &a.matrix, &xhat.matrix,
							0.0, &a_times_xhat.matrix);

		double xhat_sum_d[] = {0.0, 0.0};
		gsl_matrix_view xhat_sum = gsl_matrix_view_array(xhat_sum_d, 2, 1);

		gsl_matrix_add(&xhat_sum.matrix, &a_times_xhat.matrix);
		gsl_matrix_add(&xhat_sum.matrix, &b_times_u_3.matrix);

		gsl_matrix_set(poshat, 0,count, xhat_d[0]);
		gsl_matrix_set(velhat, 0,count, xhat_d[1]);

		// Form the Innovation vector.
		double c_times_xhat_d[] = {0.0};
		gsl_matrix_view c_times_xhat = gsl_matrix_view_array(c_times_xhat_d, 1, 1);

		double Inn_d[] = {0.0};
		gsl_matrix_view Inn = gsl_matrix_view_array(Inn_d, 1, 1);

		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c.matrix, &xhat.matrix,
							0.0, &c_times_xhat.matrix);

		
		gsl_matrix_add(&Inn.matrix, &y.matrix);
		gsl_matrix_sub(&Inn.matrix, &c_times_xhat.matrix);

		if(Steadystate == 1.0){
			gsl_matrix_set(&k.matrix, 0,0,0.1);
			gsl_matrix_set(&k.matrix, 1,0,0.01);
		}
		else{
			// compute the covariance of the Innovation
			#pragma region "s = c * P * c' + Sz;"
			
			double P_times_c_tranpose_d[] = {0.0, 0.0};
			gsl_matrix_view P_times_c_tranpose = gsl_matrix_view_array(P_times_c_tranpose_d, 2, 1);

			double c_transpose_d[] = {  	0.0,
											0.0};
			gsl_matrix_view c_transpose = gsl_matrix_view_array(c_transpose_d, 2, 1);
			gsl_matrix_transpose_memcpy(&c_transpose.matrix, &c.matrix);

			double c_P_times_c_tranpose_d[] = {0.0};
			gsl_matrix_view c_P_times_c_tranpose = gsl_matrix_view_array(c_P_times_c_tranpose_d, 1, 1);

			double S_d[] = {0.0};
			gsl_matrix_view S = gsl_matrix_view_array(S_d, 1, 1);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &P.matrix, &c_transpose.matrix,
							0.0, &P_times_c_tranpose.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c.matrix, &P_times_c_tranpose.matrix,
							0.0, &c_P_times_c_tranpose.matrix);

			gsl_matrix_add(&S.matrix, &c_P_times_c_tranpose.matrix);
			gsl_matrix_add(&S.matrix, &Sz.matrix);

			#pragma endregion

			#pragma region "K = a * P * c' * inv(s);"
			
			gsl_matrix *inverse_of_S = invert_a_matrix(&S.matrix, 1);

			double P_times_c_trans_invS_d[] = {0.0, 0.0};
			gsl_matrix_view P_times_c_trans_invS = gsl_matrix_view_array(P_times_c_trans_invS_d, 2, 1);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &P_times_c_tranpose.matrix, inverse_of_S,
							0.0, &P_times_c_trans_invS.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &a.matrix, &P_times_c_trans_invS.matrix,
							0.0, &k.matrix);


			#pragma endregion
		
			#pragma region "P = a * P * a' - a * P * c' * inv(s) * c * P * a' + Sw;"
			double a_transpose_d[] = {  	0.0, 0.0,
											0.0, 0.0};
			gsl_matrix_view a_transpose = gsl_matrix_view_array(a_transpose_d, 2, 2);
			gsl_matrix_transpose_memcpy(&a_transpose.matrix, &a.matrix);

			double P_times_a_tranpose_d[] = {  	0.0, 0.0,
											 	0.0, 0.0};
			gsl_matrix_view P_times_a_tranpose = gsl_matrix_view_array(P_times_a_tranpose_d, 2, 2);

			double c_P_times_a_tranpose_d[] = {  	0.0, 
													0.0};
			gsl_matrix_view c_P_times_a_tranpose = gsl_matrix_view_array(c_P_times_a_tranpose_d, 1, 2);

			double invS_c_P_times_a_tranpose_d[] = {  	0.0, 
													    0.0};
			gsl_matrix_view invS_c_P_times_a_tranpose = gsl_matrix_view_array(invS_c_P_times_a_tranpose_d, 1, 2);

			double ct_invS_c_P_times_a_tranpose_d[] = {  	0.0, 0.0,
													    0.0, 0.0};
			gsl_matrix_view ct_invS_c_P_times_a_tranpose = gsl_matrix_view_array(ct_invS_c_P_times_a_tranpose_d, 2, 2);

			double P_ct_invS_c_P_times_a_tranpose_d[] = {  	0.0, 0.0,
													    0.0, 0.0};
			gsl_matrix_view P_ct_invS_c_P_times_a_tranpose = gsl_matrix_view_array(P_ct_invS_c_P_times_a_tranpose_d, 2, 2);

			double a_P_ct_invS_c_P_times_a_tranpose_d[] = {  	0.0, 0.0,
													    0.0, 0.0};
			gsl_matrix_view a_P_ct_invS_c_P_times_a_tranpose = gsl_matrix_view_array(a_P_ct_invS_c_P_times_a_tranpose_d, 2, 2);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &P.matrix, &a_transpose.matrix,
							0.0, &P_times_a_tranpose.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c.matrix, &P_times_a_tranpose.matrix,
							0.0, &c_P_times_a_tranpose.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, inverse_of_S, &c_P_times_a_tranpose.matrix,
							0.0, &invS_c_P_times_a_tranpose.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &c_transpose.matrix, &invS_c_P_times_a_tranpose.matrix,
							0.0, &ct_invS_c_P_times_a_tranpose.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &P.matrix, &ct_invS_c_P_times_a_tranpose.matrix,
							0.0, &P_ct_invS_c_P_times_a_tranpose.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &a.matrix, &P_ct_invS_c_P_times_a_tranpose.matrix,
							0.0, &a_P_ct_invS_c_P_times_a_tranpose.matrix);

			double P_time_at_d[] = {  	0.0, 0.0,
										0.0, 0.0};
			gsl_matrix_view P_time_at = gsl_matrix_view_array(P_time_at_d, 2, 2);

			double a_P_time_at_d[] = {  	0.0, 0.0,
										0.0, 0.0};
			gsl_matrix_view a_P_time_at = gsl_matrix_view_array(a_P_time_at_d, 2, 2);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &P.matrix, &a_transpose.matrix,
							0.0, &P_time_at.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							1.0, &a.matrix, &P_time_at.matrix,
							0.0, &a_P_time_at.matrix);

			
			double P_sum_d[] = {  	0.0, 0.0,
										0.0, 0.0};
			gsl_matrix_view P_sum = gsl_matrix_view_array(P_sum_d, 2, 2);

			gsl_matrix_add(&P_sum.matrix, &a_P_time_at.matrix);
			gsl_matrix_sub(&P_sum.matrix, &a_P_ct_invS_c_P_times_a_tranpose.matrix);
			gsl_matrix_add(&P_sum.matrix, &Sw.matrix);

			gsl_matrix_memcpy(&P.matrix, &P_sum.matrix);
			#pragma endregion
			
			#pragma region "P = (P + P') / 2;"
			double P_transpose_d[] = {  	0.0, 0.0,
											0.0, 0.0};
			gsl_matrix_view P_transpose = gsl_matrix_view_array(P_transpose_d, 2, 2);

			gsl_matrix_transpose_memcpy(&P_transpose.matrix, &P.matrix);

			double P_sum_2_d[] = {  	0.0, 0.0,
										0.0, 0.0};
			gsl_matrix_view P_sum_2 = gsl_matrix_view_array(P_sum_2_d, 2, 2);


			gsl_matrix_add(&P_sum_2.matrix, &P.matrix);
			gsl_matrix_add(&P_sum_2.matrix, &P_transpose.matrix);

			gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							0.5, &P_sum_2.matrix, &identity.matrix,
							0.0, &P.matrix);
			#pragma endregion

		}

		// Update the Kalman filter state estimate
		double K_times_Inn_d[] = {  	0.0,
										0.0};
		gsl_matrix_view K_times_Inn = gsl_matrix_view_array(K_times_Inn_d, 2, 1);
		gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
							0.5, &k.matrix, &Inn.matrix,
							0.0, &K_times_Inn.matrix);

		gsl_matrix_add(&xhat.matrix, &K_times_Inn.matrix);

		// Save some parameters for plotting later

		gsl_matrix_set(kalmangains, 0,count, k_d[0]);
		gsl_matrix_set(kalmangains, 1,count, k_d[1]);
		
		gsl_matrix_set(pos, count, 0, x_d[0]);

		gsl_matrix_set(vel, count, 0, x_d[1]);

		//printf("kalman gain [0] : %f\n", gsl_matrix_get(kalmangains, 0, count));
		//printf("kalman gain [1] : %f\n", gsl_matrix_get(kalmangains, 1, count));
	}
}


