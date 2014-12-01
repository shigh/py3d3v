
#pragma once
#define _USE_MATH_DEFINES
#include <omp.h>
#include <math.h>

void move_par(const int N, const double dt,
			  double *zp, const double *vz,
			  double *yp, const double *vy,
			  double *xp, const double *vx);

void accel_par(const int N, const double dt, const double *qm,
			   const double *Ez, double *vz,
			   const double *Ey, double *vy,
			   const double *Ex, double *vx);

void scale_array(const int N, double *x, double sx);

void scale_array_3_copy(const int N,
						const double *z, const double sz, double *zc,
						const double *y, const double sy, double *yc,
						const double *x, const double sx, double *xc);

void calc_E_short_range_par_gaussian(int N,
									 double* Ezp, const double* zp, double Lz,
									 double* Eyp, const double* yp, double Ly,
									 double* Exp, const double* xp, double Lx,
									 const double* q, double rmax, double beta);

void calc_E_short_range_par_gaussian(int N1, int N2,
									 double* Ezp, const double* zp, double Lz,
									 double* Eyp, const double* yp, double Ly,
									 double* Exp, const double* xp, double Lx,
									 const double* zp2, const double* yp2,
									 const double* xp2, const double* q2,
									 double rmax, double beta);

void calc_E_short_range_par_s2(int N,
							   double* Ezp, const double* zp, double Lz,
							   double* Eyp, const double* yp, double Ly,
							   double* Exp, const double* xp, double Lx,
							   const double* q, double rmax, double beta);

struct Gaussian
{

private:

	// Beta and Beta^2
	double beta, beta2;
	// Constants in E
	double c, d, r2;

public:

	Gaussian(double beta_);

	double E(double r);

};

/*! S2 from Hockney and Eastwood
 */
struct S2
{

private:

	// Beta and Beta^2
	double beta;
	// Constants in E
	double c, d, fpii, b2, r2;

public:

	S2(double beta_);

	double E(double r);

};
