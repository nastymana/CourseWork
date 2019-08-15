#include "pch.h"
#include "basisLagrangeTrg.h"

LagrangeTrg::LagrangeTrg()
{
	cout << "\nLagrangeTrg::Constructor()\n";
	fi = { Lagr_trg_fixy_1, Lagr_trg_fixy_2, Lagr_trg_fixy_3 };
	grad = { Lagr_trg_grad_fixy_1, Lagr_trg_grad_fixy_2, Lagr_trg_grad_fixy_3 };
	carrier = { 0., 0., 1., 0., 0., 1 };
	GaussPoints = {
				{0.659027622374092,  0.231933368553031},
				{0.659027622374092,  0.109039009072877 },
				{0.231933368553031,  0.659027622374092},
				{0.231933368553031,  0.109039009072877},
				{0.109039009072877,  0.659027622374092},
				{0.109039009072877,  0.231933368553031 } };
	GaussWeights = {
				0.16666666666666666667,
				0.16666666666666666667,
				0.16666666666666666667,
				0.16666666666666666667,
				0.16666666666666666667,
				0.16666666666666666667 };
}


LagrangeTrg::~LagrangeTrg()
{
	cout << "\nLagrangeTrg::DEstructor()\n";
}

// Lagrange basis for triangles для мастер -элемента прямоуг треуг с координатами [00, 01, 10]
double Lagr_trg_fixy_1(const double&x, const double&y) {
	return 1 - x - y;
}
double Lagr_trg_fixy_2(const double&x, const double&y) {
	return x;
}
double Lagr_trg_fixy_3(const double&x, const double &y) {
	return y;
}

double Lagr_trg_dx_fixy_1(const double&x, const double&y) {
	return -1.;
}
double Lagr_trg_dx_fixy_2(const double&x, const double&y) {
	return 1.;
}
double Lagr_trg_dx_fixy_3(const double&x, const double&y) {
	return 0.;
}

double Lagr_trg_dy_fixy_1(const double&x, const double&y) {
	return -1.;
}
double Lagr_trg_dy_fixy_2(const double&x, const double&y) {
	return 0;
}
double Lagr_trg_dy_fixy_3(const double&x, const double&y) {
	return 1.;
}

void Lagr_trg_grad_fixy_1(const double&x, const double&y, double &resultx, double &resulty) {
	resultx = Lagr_trg_dx_fixy_1(x, y);
	resulty = Lagr_trg_dy_fixy_1(x, y);
}
void Lagr_trg_grad_fixy_2(const double&x, const double&y, double &resultx, double &resulty) {
	resultx = Lagr_trg_dx_fixy_2(x, y);
	resulty = Lagr_trg_dy_fixy_2(x, y);
}
void Lagr_trg_grad_fixy_3(const double&x, const double&y, double &resultx, double &resulty) {
	resultx = Lagr_trg_dx_fixy_3(x, y);
	resulty = Lagr_trg_dy_fixy_3(x, y);
}

vector<double> Lagr_trg_grad_fixy_1(const double&x, const double&y) {
	vector<double> result(2, 0.);
	result[0] = Lagr_trg_dx_fixy_1(x, y);
	result[1] = Lagr_trg_dy_fixy_1(x, y);
	return result;
}
vector<double> Lagr_trg_grad_fixy_2(const double&x, const double&y) {
	vector<double> result(2, 0.);
	result[0] = Lagr_trg_dx_fixy_2(x, y);
	result[1] = Lagr_trg_dy_fixy_2(x, y);
	return result;
}
vector<double> Lagr_trg_grad_fixy_3(const double&x, const double &y) {
	vector<double> result(2, 0.);
	result[0] = Lagr_trg_dx_fixy_3(x, y);
	result[1] = Lagr_trg_dy_fixy_3(x, y);
	return result;
}
