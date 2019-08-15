#pragma once
#include "basisFunction.h"
namespace Lagrange {
	class basisLagrangeQuad :public basisFunction
	{
	public:
		basisLagrangeQuad(const int &rankItem);
		virtual ~basisLagrangeQuad();
		void coutBasisParams();
		double epsToX(const double &x1, const double &x2, const double &eps);
		double XtoEps(const double & x1, const double & x2, const double & x);

	};

	double fiBiline1(const double &ksi, const double & eta);
	double fiBiline2(const double &ksi, const double & eta);

	double fiBiline3(const double &ksi, const double & eta);
	double fiBiline4(const double &ksi, const double & eta);

	vector<double> gradFiBiline1(const double &ksi, const double & eta);
	vector<double> gradFiBiline2(const double &ksi, const double & eta);

	vector<double> gradFiBiline3(const double &ksi, const double & eta);
	vector<double> gradFiBiline4(const double &ksi, const double & eta);

	////////////////////////////////////////////////////////////////////////
	double fiBiquad1(const double &ksi, const double & eta);
	double fiBiquad2(const double &ksi, const double & eta);
	double fiBiquad3(const double &ksi, const double & eta);
	double fiBiquad4(const double &ksi, const double & eta);
	double fiBiquad5(const double &ksi, const double & eta);
	double fiBiquad6(const double &ksi, const double & eta);
	double fiBiquad7(const double &ksi, const double & eta);
	double fiBiquad8(const double &ksi, const double & eta);
	double fiBiquad9(const double &ksi, const double & eta);

	vector<double> gradFiBiquad1(const double &ksi, const double & eta);
	vector<double> gradFiBiquad2(const double &ksi, const double & eta);
	vector<double> gradFiBiquad3(const double &ksi, const double & eta);
	vector<double> gradFiBiquad4(const double &ksi, const double & eta);
	vector<double> gradFiBiquad5(const double &ksi, const double & eta);
	vector<double> gradFiBiquad6(const double &ksi, const double & eta);
	vector<double> gradFiBiquad7(const double &ksi, const double & eta);
	vector<double> gradFiBiquad8(const double &ksi, const double & eta);
	vector<double> gradFiBiquad9(const double &ksi, const double & eta);

	vector<double> laplacianFiBiquad1(const double & ksi, const double & eta);
	vector<double> laplacianFiBiquad2(const double & ksi, const double & eta);
	vector<double> laplacianFiBiquad3(const double & ksi, const double & eta);
	vector<double> laplacianFiBiquad4(const double & ksi, const double & eta);
	vector<double> laplacianFiBiquad5(const double & ksi, const double & eta);
	vector<double> laplacianFiBiquad6(const double & ksi, const double & eta);
	vector<double> laplacianFiBiquad7(const double & ksi, const double & eta);
	vector<double> laplacianFiBiquad8(const double & ksi, const double & eta);
	vector<double> laplacianFiBiquad9(const double & ksi, const double & eta);

	/////////////////////////////////////////////////////////////////////////////////////////

		// вспомогательные функции:
	double PsiLine1(const double &ksi);
	double PsiLine2(const double &ksi);

	double dPsiLine1(const double &ksi);
	double dPsiLine2(const double &ksi);

	////////////////////////////////
	double PsiQuad1(const double &ksi);
	double PsiQuad2(const double &ksi);
	double PsiQuad3(const double &ksi);

	double dPsiQuad1(const double &ksi);
	double dPsiQuad2(const double &ksi);
	double dPsiQuad3(const double &ksi);

	double d2PsiQuad1(const double &ksi);
	double d2PsiQuad2(const double &ksi);
	double d2PsiQuad3(const double &ksi);

	///////////////////////////////
	double PsiCube1(const double &ksi);
	double PsiCube2(const double &ksi);
	double PsiCube3(const double &ksi);
	double PsiCube4(const double &ksi);

	double dPsiCube1(const double &ksi);
	double dPsiCube2(const double &ksi);
	double dPsiCube3(const double &ksi);
	double dPsiCube4(const double &ksi);
}