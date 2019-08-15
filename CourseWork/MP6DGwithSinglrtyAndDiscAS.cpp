#include "pch.h"
#include "MP6DGwithSinglrtyAndDiscAS.h"


MP6DGwithSinglrtyAndDiscAS::MP6DGwithSinglrtyAndDiscAS()
{
}

MP6DGwithSinglrtyAndDiscAS::MP6DGwithSinglrtyAndDiscAS(const vector<double>& inData, const vector<double>& domainParams)
{
}


MP6DGwithSinglrtyAndDiscAS::~MP6DGwithSinglrtyAndDiscAS()
{
}

double MP6DGwithSinglrtyAndDiscAS::realSol(const double & x, const double & y)
{
	return 0.0;
}

double MP6DGwithSinglrtyAndDiscAS::BC1(const double & x, const double & y)
{
	return 0.0;
}

double MP6DGwithSinglrtyAndDiscAS::BC2(const double & x, const double & y)
{
	return 0.0;
}

double MP6DGwithSinglrtyAndDiscAS::BC3(const double & x, const double & y)
{
	return 0.0;
}

double MP6DGwithSinglrtyAndDiscAS::source(const double & x, const double & y)
{
	// vector<double> Vxy = V(x, y);
	// lymbda=params [0], eps = params[1], gamma = params[2];
	return 0.;/* -params[0] * (MP2_d2AnSol_dx2(x, y, params[1]) + MP2_d2AnSol_dy2(x, y, params[1]))
		+ params[2] * MP2_AnSol(x, y, params[1]);*/
}

double MP6DGwithSinglrtyAndDiscAS::lymbda(const double & x, const double & y)
{
	if ((y >= 0.) & (y < 1.)) {
		return pow(lyambdaParam ,2);
	}
}

double MP6DGwithSinglrtyAndDiscAS::gamma(const double & x, const double & y)
{
	if ((y >= 0.) & (y < 0.5)) {
		return 2.;
	}
	if ((y >= 0.5) & (y <= 1.0)) {
		return 1.;
	}

}

vector<double> MP6DGwithSinglrtyAndDiscAS::velocity(const double & x, const double & y)
{
	return { 0.,0. };
}
double MP6DGwithSinglrtyAndDiscAS::MP2_U1(const double &x) {
	if ((x >= 0.) & (x <= 1.)) {
		return sin(4. * PI*x) + 2.;
	}
	// else cout error return NULL;
}
double MP6DGwithSinglrtyAndDiscAS::MP2_U2(const double &x, const double &eps) {
	if ((x >= 0.) & (x <= 1.)) {
		return 1. + exp(-1. / eps) - exp((-x) / eps) - exp((x - 1.) / eps);
	}
}

double MP6DGwithSinglrtyAndDiscAS::MP2_dU1(const double &x) {
	if ((x >= 0.) & (x <= 1.)) {
		return 4.*PI*cos(4.*PI*x);
	}
}
double MP6DGwithSinglrtyAndDiscAS::MP2_dU2(const double &x, const double &eps) {
	if ((x >= 0.) & (x <= 1.)) {
		return (1. / eps)*(exp((-x) / eps) - exp((x - 1.) / eps));
	}

}

double MP6DGwithSinglrtyAndDiscAS::MP2_d2U1(const double &x) {
	if ((x >= 0.) & (x <= 1.)) {
		return -16. * PI*PI*sin(4. * PI*x);
	}
}
double MP6DGwithSinglrtyAndDiscAS::MP2_d2U2(const double &x, const double &eps) {
	if ((x >= 0.) & (x <= 1.)) {
		return (1. / (eps*eps))*(-exp((-x) / eps) - exp((x - 1.) / eps));
	}
}
// Ux - x-компонента функции
double MP6DGwithSinglrtyAndDiscAS::MP2_Ux(const double &x, const double &eps) {
	if ((x >= 0.) & (x <= 1.)) {
		return MP2_U1(x)*MP2_U2(x, eps);
	}
}
double MP6DGwithSinglrtyAndDiscAS::MP2_Uy(const double &y, const double &eps) {
	if ((y >= 0.) & (y <= 1.)) {
		return MP2_U2(y, eps);
	}
}

// dUx -полный диф-ал
double MP6DGwithSinglrtyAndDiscAS::MP2_dU_dx(const double &x, const double &y, const double &eps) {
	return (MP2_dU1(x)*MP2_U2(x, eps) +
		MP2_U1(x)*MP2_dU2(x, eps))
		*MP2_Uy(y, eps);
}
double MP6DGwithSinglrtyAndDiscAS::MP2_dU_dy(const double &x, const double &y, const double &eps) {
	return (MP2_U1(x)*MP2_U2(x, eps))
		*MP2_dU2(y, eps);
}

double MP6DGwithSinglrtyAndDiscAS::MP2_d2U_dx2(const double &x, const double &y, const double &eps) {
	return (MP2_d2U1(x)*MP2_U2(x, eps) +
		2.*MP2_dU1(x)*MP2_dU2(x, eps) +
		MP2_U1(x)*MP2_d2U2(x, eps))
		*MP2_U2(y, eps);
}
double MP6DGwithSinglrtyAndDiscAS::MP2_d2U_dy2(const double &x, const double &y, const double &eps) {
	return (MP2_U1(x)*MP2_U2(x, eps))*MP2_d2U2(y, eps);
}

double MP6DGwithSinglrtyAndDiscAS::MP2_velocity(const double &x, const double &y, const double &eps) {
	if ((y >= 0.) & (y < 0.5)) {
		return 3. - exp((2.*y - 1.) / (2.*eps));
	}
	if ((y >= 0.5) & (y <= 1.0)) {
		return 1. + exp((1. - 2.*y) / (2.*eps));
	}
}
double MP6DGwithSinglrtyAndDiscAS::MP2_dV_dx(const double &x, const double &y, const double &eps) { return 0.; }
double MP6DGwithSinglrtyAndDiscAS::MP2_dV_dy(const double &x, const double &y, const double &eps) {
	if ((y >= 0.) & (y < 0.5)) {
		return (-1. / eps) * exp((2.*y - 1.) / (2.*eps));
	}
	if ((y >= 0.5) & (y <= 1.0)) {
		return (-1. / eps) * exp((1. - 2.*y) / (2.*eps));
	}
}

double MP6DGwithSinglrtyAndDiscAS::MP2_d2V_dx2(const double &x, const double &y, const double &eps) { return 0.; }
double MP6DGwithSinglrtyAndDiscAS::MP2_d2V_dy2(const double &x, const double &y, const double &eps) {
	if ((y >= 0.) & (y < 0.5)) {
		return (-1. / (eps*eps)) * exp((2.*y - 1.) / (2.*eps));
	}
	if ((y >= 0.5) & (y <= 1.0)) {
		return (1. / (eps*eps)) * exp((1. - 2.*y) / (2.*eps));
	}
}


double MP6DGwithSinglrtyAndDiscAS::MP2_AnSol(const double &x, const double &y, const double &eps) {
	return 0.25*MP2_velocity(x, y, eps) * MP2_Ux(x, eps) * MP2_Uy(y, eps);
}

double MP6DGwithSinglrtyAndDiscAS::MP2_dAnSol_dx(const double &x, const double &y, const double &eps) {
	return 0.25*(
		MP2_dV_dx(x, y, eps) * (MP2_Ux(x, eps) * MP2_Uy(y, eps))
		+ MP2_velocity(x, y, eps) * MP2_dU_dx(x, y, eps));
}
double MP6DGwithSinglrtyAndDiscAS::MP2_dAnSol_dy(const double &x, const double &y, const double &eps) {
	return 0.25*(
		MP2_dV_dy(x, y, eps) * (MP2_Ux(x, eps) * MP2_Uy(y, eps))
		+ MP2_velocity(x, y, eps) * MP2_dU_dy(x, y, eps));
}

double MP6DGwithSinglrtyAndDiscAS::MP2_d2AnSol_dx2(const double &x, const double &y, const double &eps) {
	return 0.25*(
		MP2_d2V_dx2(x, y, eps) * (MP2_Ux(x, eps) * MP2_Uy(y, eps))
		+ 2.* MP2_dV_dx(x, y, eps) * MP2_dU_dx(x, y, eps)
		+ MP2_velocity(x, y, eps) * MP2_d2U_dx2(x, y, eps));
}
double MP6DGwithSinglrtyAndDiscAS::MP2_d2AnSol_dy2(const double &x, const double &y, const double &eps) {
	return 0.25*(
		MP2_d2V_dy2(x, y, eps) * (MP2_Ux(x, eps) * MP2_Uy(y, eps))
		+ 2.* MP2_dV_dy(x, y, eps) * MP2_dU_dy(x, y, eps)
		+ MP2_velocity(x, y, eps) * MP2_d2U_dy2(x, y, eps));
}

