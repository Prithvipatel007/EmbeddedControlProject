#include <stdio.h>
#include "stdlib.h"
#include "Math.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include "invertor.c"

#define POS_ROW 2
#define POS_COL 2

#define VEL_ROW 2
#define VEL_COL 2

#define POSHAT_ROW 2
#define POSHAT_COL 2

#define VELHAT_ROW 2
#define VELHAT_COL 2

#define POSHATINF_ROW 2
#define POSHATINF_COL 2

#define VELHATINF_ROW 2
#define VELHATINF_COL 2

#define HINFGAINS_ROW 2
#define HINFGAINS_COL 1

#define KALMANGAINS_ROW 2
#define KALMANGAINS_COL 1



void HinfinityInC(float, float, int, int, double[POS_ROW][POS_COL], double[VEL_ROW][VEL_ROW],
					double[POSHAT_ROW][POSHAT_COL],double[VELHAT_ROW][VELHAT_COL],
					double[POSHATINF_ROW][POSHATINF_COL],double[VELHATINF_ROW][VELHATINF_COL],
					double[HINFGAINS_ROW][HINFGAINS_COL], double[KALMANGAINS_ROW][KALMANGAINS_COL]);

float gama = 0.01;
float dt = 0.1;
int duration = 60;
int Steadystate = 0;

void main()
{
	int i,j;

	double pos[POS_ROW][POS_COL] = {0,0,0,0};
	double vel[VEL_ROW][VEL_COL] = {0,0,0,0};
	double poshat[POSHAT_ROW][POSHAT_COL] = {0,0,0,0};
	double velhat[VELHAT_ROW][VELHAT_COL] = {0,0,0,0};
	double poshatinf[POSHATINF_ROW][POSHATINF_COL] = {0,0,0,0};
	double velhatinf[VELHATINF_ROW][VELHATINF_COL]= {0,0,0,0};
	double hinfgains[HINFGAINS_ROW][HINFGAINS_COL] = {0,0};
	double kalmangains[KALMANGAINS_ROW][KALMANGAINS_COL] = {0,0};

	HinfinityInC(gama, dt, duration, Steadystate, pos, vel,poshat, velhat,
				poshatinf, velhatinf, hinfgains, kalmangains);

}

void HinfinityInC(float gama,float dt,int duration,int Steadystate, double pos[POS_ROW][POS_COL],
					double vel[VEL_ROW][VEL_COL],
					double poshat[POSHAT_ROW][POSHAT_COL],
					double velhat[VELHAT_ROW][VELHAT_COL],
					double poshatinf[POSHATINF_ROW][POSHATINF_COL],
					double velhatinf[VELHATINF_ROW][VELHATINF_COL],
					double hinfgains[HINFGAINS_ROW][HINFGAINS_COL],
					double kalmangains[KALMANGAINS_ROW][KALMANGAINS_COL])
{
	/*For loop counter declarations */
	int i,j, count, t,u;
	double K[2][1] = {0,0};
	double L[2][2] = {0,0,0,0};


	int measnoise = 2;  // nominal noise velocity measurement noise
	float accelnoise = 0.2; // nominal acceleration noise
	double a_mat[2][2] = {1,0.1, 0,1}; // transition matrix
	double b_mat[2][1] = {0.005,0.1}; // Input matrix
	double c_mat[] = {0,1}; 			// measurement matrix
	double x_mat[2][1] = {0,0};			// Initial state vector

	/*y = c * x */
	double y_mat = 0;					//Initial measurement


	//Initialize Kalman filter variables
	double xhat[2][1]  = {0,0};		//Initial Kalman filter estimat
	int Sz = 4;					// Sz = meanoise^2, measurement error covariance

	//Sw = accelnoise^2 * [dt^4/4 dt^3/2; dt^3/2 dt^2]
	double Sw[2][2] = {0, 0.00002, 0.00002, 0.0004};	// Process noise cov

	/*Initial Kalman filter estimation covariance */
	double P[POS_ROW][POS_COL];
	for(i=0; i< POS_ROW ; i++)					//Copy Sw array values
	{
		for(j=0; j< POS_COL ; j++)
		{
			P[i][j] = Sw[i][j];
		}
	}

	/*Initialize Hinfinity filter variables */
	double xhatinf[2][1] = {0,0};			//Initial Hinfinity filter state estimate

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
			L = invert_a_matrix(K);
		}



	}


}
