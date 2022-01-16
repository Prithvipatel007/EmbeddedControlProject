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

	printf ("[ %f, %f\n", P_d[0], P_d[1]);
  	printf ("  %f, %f ]\n", P_d[2], P_d[3]);

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
