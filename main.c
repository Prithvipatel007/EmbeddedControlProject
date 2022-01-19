#include <stdio.h>
#include "stdlib.h"
#include "Math.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

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
	int i,j, count, t,u;

	//double K[2][1] = {0,0};
	//double L[2][2] = {0,0,0,0};

	double k_d[] = {0,0};
	double l_d[] = {0,0,0,0};

	gsl_matrix_view k = gsl_matrix_view_array(k_d, 2, 1);		// Use steady-state H-infinity gains
	gsl_matrix_view l = gsl_matrix_view_array(l_d, 2, 2);

	double measnoise = 2.0;  // nominal noise velocity measurement noise
	double accelnoise = 0.2; // nominal acceleration noise

	double a_d[] = { 	1.0, dt,
                 	0.0, 1.0};

	double b_d[] = { 	pow(dt,2.0)/2.0,
						dt};

	double c_d[] = {  	0.0,
						1.0};

	double x_d[] = {  	0.1,
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
		0.0, &y);							//Initial measurement

	/*
		Initialize Kalman filter variables
	*/
	double xhat_d[]  = {0,0};
	gsl_matrix_view xhat = gsl_matrix_view_array(xhat_d, 2, 2);

	double Sw_d[] = { (pow(accelnoise,2))*(pow(dt,4)/4),  (pow(accelnoise,2))*(pow(dt,3)/2), 
						(pow(accelnoise,2))*(pow(dt,3)/2), (pow(accelnoise,2))*(pow(dt,2))};
	gsl_matrix_view Sw = gsl_matrix_view_array(Sw_d, 2, 2);

	double P_d[] = { (pow(accelnoise,2))*(pow(dt,4)/4),  (pow(accelnoise,2))*(pow(dt,3)/2), 
						(pow(accelnoise,2))*(pow(dt,3)/2), (pow(accelnoise,2))*(pow(dt,2))};
	gsl_matrix_view P = gsl_matrix_view_array(P_d, 2, 2);

	double xhatinf_d[] = {  	0.0,
								0.0};
	gsl_matrix_view xhatinf = gsl_matrix_view_array(xhatinf_d, 2, 1);


	double Pinf_d[] = {  	0.01, 0,
							0, 0.01};
	gsl_matrix_view Pinf = gsl_matrix_view_array(Pinf_d, 2, 2);

	double W_d[] = {  	0.0003/1000, 0.0050/1000,
						0.0050/1000, 0.1000/1000};
	gsl_matrix_view W = gsl_matrix_view_array(W_d, 2, 2);

	double V_d[] = {0.01};
	gsl_matrix_view V = gsl_matrix_view_array(V_d, 1, 1);

	double Q_d[] = {  	0.01, 0.0,
						0.0, 0.01};
	gsl_matrix_view Q = gsl_matrix_view_array(Q_d, 2, 2);

	/* Initialize the Pos value */
	for (i = 0; i < POS_ROW; i++)
		for (j = 0; j < POS_COL; j++)
			gsl_matrix_set (pos, i, j, x_d[0]);

	/* Initialize the Vel value */
	for (i = 0; i < VEL_ROW; i++)
		for (j = 0; j < VEL_COL; j++)
			gsl_matrix_set (vel, i, j, x_d[1]);

	/* Initialize the Poshat value */
	for (i = 0; i < POSHAT_ROW; i++)
		for (j = 0; j < POSHAT_COL; j++)
			gsl_matrix_set (poshat, i, j, xhat_d[0]);

	/* Initialize the Velhat value */
	for (i = 0; i < VELHAT_ROW; i++)
		for (j = 0; j < VELHAT_COL; j++)
			gsl_matrix_set (velhat, i, j, xhat_d[1]);

	/* Initialize the Poshatinf value */
	for (i = 0; i < POSHATINF_ROW; i++)
		for (j = 0; j < POSHATINF_COL; j++)
			gsl_matrix_set (poshatinf, i, j, xhatinf_d[0]);

	/* Initialize the Velhatinf value */
	for (i = 0; i < VELHAT_ROW; i++)
		for (j = 0; j < VELHAT_COL; j++)
			gsl_matrix_set (velhatinf, i, j, xhatinf_d[1]);
	
	double hinfgains_d[] = {  	0.0,
								0.0};
	
	/* Initialize the hinfgains value */
	for (i = 0; i < HINFGAINS_ROW; i++)
		for (j = 0; j < HINFGAINS_COL; j++)
			gsl_matrix_set (hinfgains, i, j, 0.0);

	/* Initialize the Kalmangains value */
	for (i = 0; i < KALMANGAINS_ROW; i++)
		for (j = 0; j < KALMANGAINS_COL; j++)
			gsl_matrix_set (hinfgains, i, j, 0.0);

	int dt_t = (int)(dt*10.0);
	int dur_t = (int)(duration*10.0);
	#pragma endregion

	for(i=0; i<dur_t; i+=dt){
		// use a constant commanded acceleration of 1 foot/sec^2
		double u = 1.0;

		// figure out the H-Infinity gains
		if(Steadystate == 1.0){
			// Use Steady State H-Infinity Gain
			gsl_matrix_set(&k.matrix, 0,0,0.11);
			gsl_matrix_set(&k.matrix, 1,0,0.09);
		}
		else{
			//L = inv(eye(2) - g * Q * Pinf + c' * inv(V) * c_time_Pinf);
			size_t size = 1;
			#pragma region  Finding the L equation
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

			gsl_matrix_add(&l.matrix, &part_1_eq.matrix);
			gsl_matrix_sub(&l.matrix, &part_2_eq.matrix);
			gsl_matrix_add(&l.matrix, &till_c_transpose.matrix);
			#pragma endregion


		}
	}

	/*double test_d[] = {1.0,2.0,3.0,4.0};
	gsl_matrix_view test = gsl_matrix_view_array(test_d, 2, 2);
	gsl_matrix *transpose = gsl_matrix_alloc(2,2);
	gsl_matrix_transpose_memcpy(transpose, &test.matrix);
	print_mat_contents(transpose, size);*/

	/*double test_d[] = {1.0,2.0,3.0,4.0};
			gsl_matrix_view test = gsl_matrix_view_array(test_d, 2, 2);
			gsl_matrix *inverse = invert_a_matrix(&test.matrix, size);
			print_mat_contents(inverse, size);*/


/*
	for (i = 0; i < HINFGAINS_ROW; i++)
		for (j = 0; j < HINFGAINS_COL; j++)
			printf ("hinfgains(%d,%d) = %f\n", i, j, gsl_matrix_get (hinfgains, i, j));
	/*printf ("[ %f, %f\n", xhatinf[0], xhatinf[1]);
  	printf ("  %f, %f ]\n", xhatinf[2], xhatinf[3]);*/

	/*
		Initialize H-infinity filter variables
	*/
/*
	double xhatinf_d[]  = {0,0};
	gsl_matrix_view xhatinf = gsl_matrix_view_array(xhatinf_d, 2, 1);

	double Pinf_d[]  = {0.01, 0, 0, 0.01};
	gsl_matrix_view Pinf = gsl_matrix_view_array(Pinf_d, 2, 1);

	double W_d[]  = {0.0003/1000, 0.0050/1000, 0.0050/1000, 0.10001/1000};
	gsl_matrix_view Pinf = gsl_matrix_view_array(Pinf_d, 2, 1);*/




/*

	//Initialize Kalman filter variables
	double xhat[2][1]  = {0,0};		//initial Kalman filter state estimate
	int Sz = 4;					// Sz = meanoise^2, measurement error covariance

	//Sw = accelnoise^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]
	double Sw[2][2] = {0, 0.00002, 0.00002, 0.0004};	// Process noise cov

	/*Initial Kalman filter estimation covariance */
/*	double P[POS_ROW][POS_COL];
	for(i=0; i< POS_ROW ; i++)					//Copy Sw array values
	{
		for(j=0; j< POS_COL ; j++)
		{
			P[i][j] = Sw[i][j];
		}
	}

	/*Initialize Hinfinity filter variables */
/*	double xhatinf[2][1] = {0,0};			//Initial Hinfinity filter state estimate

	//Pinf = 0.01*eye(2);
	double Pinf[2][2] = {0.01, 0 ,  0 , 0.01};
	double w_mat[] = {(3/10000000), (1/200000), (1/200000), (1/10000) };
	double v = 0.01;
	double q_mat[2][2] = {0.01,0,0.01};

	for(count=0; count < 599 ; ++count)
	{
		t = count/ 10.0f ;

		//Use a constant commanded acceleration of 1 foot/sec^2.
		u =1;

		 //Figure out the H-infinity estimate.
		if(Steadystate == 1)
		{
			K[0][0] = 0.11;
			K[1][0] = 0.09;
		}
		else
		{

		}



	}*/


}
