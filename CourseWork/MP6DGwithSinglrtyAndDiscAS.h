#pragma once
#include "AbstractModel.h"

class MP6DGwithSinglrtyAndDiscAS :
		public AbstractModel{
public:
	MP6DGwithSinglrtyAndDiscAS();
	MP6DGwithSinglrtyAndDiscAS(const vector<double> &inData, const vector<double>& domainParams);
	virtual ~MP6DGwithSinglrtyAndDiscAS();

	double realSol(const double &x, const double &y);
	double BC1(const double &x, const double &y);
		double BC2(const double &x, const double &y);
		double BC3(const double &x, const double &y);
		double source(const double &x, const double &y);
		double lymbda(const double &x, const double &y);
		double gamma(const double &x, const double &y);
		virtual double initialCond(const double &x, const double &y) = 0;
		vector<double> velocity(const double &x, const double &y);
	private:

		double MP2_U1(const double &x);
		double MP2_U2(const double &x, const double &eps);
		double MP2_dU1(const double &x);
		double MP2_dU2(const double &x, const double &eps);
		double MP2_d2U1(const double &x);
		double MP2_d2U2(const double &x, const double &eps);
		// Ux - x-компонента функции
		double MP2_Ux(const double &x, const double &eps);
		double MP2_Uy(const double &y, const double &eps);


		// dUx -полный диф-ал
		double MP2_dU_dx(const double &x, const double &y, const double &eps);
		double MP2_dU_dy(const double &x, const double &y, const double &eps);
		double MP2_d2U_dx2(const double &x, const double &y, const double &eps);
		double MP2_d2U_dy2(const double &x, const double &y, const double &eps);

		double MP2_velocity(const double &x, const double &y, const double &eps);

		double MP2_dV_dx(const double &x, const double &y, const double &eps);
		double MP2_dV_dy(const double &x, const double &y, const double &eps);

		double MP2_d2V_dx2(const double &x, const double &y, const double &eps);
		double MP2_d2V_dy2(const double &x, const double &y, const double &eps);

		double MP2_AnSol(const double &x, const double &y, const double &eps);

		double MP2_dAnSol_dx(const double &x, const double &y, const double &eps);

		double MP2_dAnSol_dy(const double &x, const double &y, const double &eps);

		double MP2_d2AnSol_dx2(const double &x, const double &y, const double &eps);

		double MP2_d2AnSol_dy2(const double &x, const double &y, const double &eps);

};