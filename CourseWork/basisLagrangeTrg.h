#pragma once
#include "pch.h"
#include "basisFunction.h"
class LagrangeTrg :
	public basisFunction
{ // мноооогое переделать
public:
	LagrangeTrg();
	virtual~LagrangeTrg();
private:
	vector<scalar_function> fi;
	vector<vector_function> grad;
	vector<double> carrier;
	vector<vector<double>> GaussPoints; // protected - открыть доступ для класса SLAE
	vector<double> GaussWeights;
private:

};


// Lagrange basis for triangles для мастер -элемента прямоуг треуг с координатами [00, 01, 10]
double Lagr_trg_fixy_1(const double&x, const double&y);
double Lagr_trg_fixy_2(const double&x, const double&y);
double Lagr_trg_fixy_3(const double&x, const double &y);

double Lagr_trg_dx_fixy_1(const double&x, const double&y);
double Lagr_trg_dx_fixy_2(const double&x, const double&y);
double Lagr_trg_dx_fixy_3(const double&x, const double&y);

double Lagr_trg_dy_fixy_1(const double&x, const double&y);
double Lagr_trg_dy_fixy_2(const double&x, const double&y);
double Lagr_trg_dy_fixy_3(const double&x, const double&y);

void Lagr_trg_grad_fixy_1(const double&x, const double&y, double &resultx, double &resulty);
void Lagr_trg_grad_fixy_2(const double&x, const double&y, double &resultx, double &resulty);
void Lagr_trg_grad_fixy_3(const double&x, const double&y, double &resultx, double &resulty);

vector<double> Lagr_trg_grad_fixy_1(const double&x, const double&y);
vector<double> Lagr_trg_grad_fixy_2(const double&x, const double&y);
vector<double> Lagr_trg_grad_fixy_3(const double&x, const double &y);

