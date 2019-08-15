#include "pch.h"
#include "basisLagrangeQuad.h"
using namespace Lagrange;
// возможны проблемы с именами. тк в другом базисе я (функции строители) их буду называть также
basisLagrangeQuad::basisLagrangeQuad(const int &rankItem)
{
	cout << "\basisLagrangeQuad::Constructor()\n";
	basisRank = rankItem;

	/*
	vector<vector_function> laplacian; 
	vector<double> carrierX, carrierY;
	vector<double> GaussPoints;
	vector<double> GaussWeights;*/

	if (rankItem == 1) {
		
		fi = { fiBiline1,	fiBiline2, fiBiline3,  fiBiline4 };
		grad = { gradFiBiline1, gradFiBiline2, gradFiBiline3, gradFiBiline4 };
		laplacian = {};

		basisSize = 4;
		GaussRank = 2;

		GaussWeights = { 1., 1. };
		GaussPoints = { 0.211324855, 0.7886751345 };
		
		carrierX = { 0., 1. };
		carrierY = { 0., 1. };

	}
	else if (rankItem == 2) {
		fi = { &fiBiquad1, &fiBiquad2, &fiBiquad3, &fiBiquad3, &fiBiquad4,
					&fiBiquad5, &fiBiquad6, &fiBiquad7, &fiBiquad8, &fiBiquad9 };

		grad = { &gradFiBiquad1, &gradFiBiquad2, &gradFiBiquad3, &gradFiBiquad4,
				&gradFiBiquad5, &gradFiBiquad6, &gradFiBiquad7, &gradFiBiquad8, &gradFiBiquad9 };
		laplacian = { laplacianFiBiquad1, laplacianFiBiquad2, laplacianFiBiquad3, laplacianFiBiquad4,laplacianFiBiquad5,
					 laplacianFiBiquad6, laplacianFiBiquad7, laplacianFiBiquad8, laplacianFiBiquad9 };
		basisSize = 9;
		GaussRank = 3;
		GaussPoints = { 0.1127016655,	 0.5,				0.8872983345};
		GaussWeights = { 0.5555555555,	 0.8888888888889,	0.55555555};
		carrierX = { 0., 0.5, 1. };
		carrierY = { 0., 0.5, 1. };
	}
	else if (rankItem == 3) {
		fi = {};
		grad = {};
		basisSize = 16;
		GaussRank = 5;
		GaussPoints = { 0.049101, 0.23076535, 0.5, 0.76923465, 0.9530899 };
		GaussWeights = { 0.4786287,0.2369269,0.5688888,0.2369269,0.4786287 };
		carrierX = { 0., 0.25, 0.75, 1. };
		carrierY = { 0., 0.25, 0.75, 1. };
	}

	

}

basisLagrangeQuad::~basisLagrangeQuad()
{
	cout << "\basisLagrangeQuad::Destructor()\n";
}

void Lagrange::basisLagrangeQuad::coutBasisParams()
{
	cout << "Lagrange basis for quad element, rank = " << basisRank << endl;
}

double Lagrange::basisLagrangeQuad::epsToX(const double & x1, const double & x2, const double & eps){
		return  (eps * (x2 - x1)+x1);
}
double Lagrange::basisLagrangeQuad::XtoEps(const double & x1, const double & x2, const double & x) {
	return  (x - x1) / (x2 - x1);

}


///////////////////////////////////////////////////////////////////////////
namespace Lagrange {
	double fiBiline1(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiline1()\n";
		return PsiLine1(ksi)*PsiLine1(eta);
	}

	double fiBiline2(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiline2()\n";
		return PsiLine2(ksi)*PsiLine1(eta);
	}

	double fiBiline3(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiline3()\n";
		return PsiLine1(ksi)*PsiLine2(eta);
	}

	double fiBiline4(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiline4()\n";
		return PsiLine2(ksi)*PsiLine2(eta);
	}

	vector<double>  gradFiBiline1(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiline1()\n";
		return { dPsiLine1(ksi)*PsiLine1(eta), PsiLine1(ksi)*dPsiLine1(eta) };
	}

	vector<double>  gradFiBiline2(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiline2()\n";
		return { dPsiLine2(ksi)*PsiLine1(eta), PsiLine2(ksi)*dPsiLine1(eta) };
	}

	vector<double>  gradFiBiline3(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiline3()\n";
		return { dPsiLine1(ksi)*PsiLine2(eta), PsiLine1(ksi)*dPsiLine2(eta) };
	}

	vector<double>  gradFiBiline4(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiline4()\n";
		return { dPsiLine2(ksi)*PsiLine2(eta), PsiLine2(ksi)*dPsiLine2(eta) };
	}


	double  fiBiquad1(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad1()\n";
		return PsiQuad1(ksi)*PsiQuad1(eta);
	}

	double  fiBiquad2(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad2()\n";
		return PsiQuad2(ksi)*PsiQuad1(eta);
	}

	double  fiBiquad3(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad3()\n";
		return PsiQuad3(ksi)*PsiQuad1(eta);
	}

	double  fiBiquad4(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad4()\n";
		return PsiQuad1(ksi)*PsiQuad2(eta);
	}

	double  fiBiquad5(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad5()\n";
		return PsiQuad2(ksi)*PsiQuad2(eta);
	}

	double  fiBiquad6(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad6()\n";
		return PsiQuad3(ksi)*PsiQuad2(eta);
	}

	double  fiBiquad7(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad7()\n";
		return PsiQuad1(ksi)*PsiQuad3(eta);
	}

	double  fiBiquad8(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad8()\n";
		return PsiQuad2(ksi)*PsiQuad3(eta);
	}

	double  fiBiquad9(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::fiBiquad9()\n";
		return PsiQuad3(ksi)*PsiQuad3(eta);
	}

	vector<double>  gradFiBiquad1(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad1()\n";
		return { dPsiQuad1(ksi)*PsiQuad1(eta),  PsiQuad1(ksi)*dPsiQuad1(eta) };
	}

	vector<double> gradFiBiquad2(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad2()\n";
		return { dPsiQuad2(ksi)*PsiQuad1(eta),  PsiQuad2(ksi)*dPsiQuad1(eta) };
	}

	vector<double>  gradFiBiquad3(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad3()\n";
		return { dPsiQuad3(ksi)*PsiQuad1(eta),  PsiQuad3(ksi)*dPsiQuad1(eta) };
	}

	vector<double> gradFiBiquad4(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad4()\n";
		return { dPsiQuad1(ksi)*PsiQuad2(eta),  PsiQuad1(ksi)*dPsiQuad2(eta) };
	}

	vector<double>  gradFiBiquad5(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad5()\n";
		return { dPsiQuad2(ksi)*PsiQuad2(eta),  PsiQuad2(ksi)*dPsiQuad2(eta) };
	}

	vector<double> gradFiBiquad6(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad6()\n";
		return { dPsiQuad3(ksi)*PsiQuad2(eta),  PsiQuad3(ksi)*dPsiQuad2(eta) };
	}

	vector<double>  gradFiBiquad7(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad7()\n";
		return { dPsiQuad1(ksi)*PsiQuad3(eta),  PsiQuad1(ksi)*dPsiQuad3(eta) };
	}

	vector<double> gradFiBiquad8(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad8()\n";
		return { dPsiQuad2(ksi)*PsiQuad3(eta),  PsiQuad2(ksi)*dPsiQuad3(eta) };
	}

	vector<double>  gradFiBiquad9(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::gradFiBiquad9()\n";
		return { dPsiQuad3(ksi)*PsiQuad3(eta),  PsiQuad3(ksi)*dPsiQuad3(eta) };
	}

	vector<double> laplacianFiBiquad1(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::laplacianFiBiquad1()\n";
		return  { d2PsiQuad1(ksi)*PsiQuad1(eta) , PsiQuad1(ksi)*d2PsiQuad1(eta) };
	}

	vector<double> laplacianFiBiquad2(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::laplacianFiBiquad2()\n";
		return  { d2PsiQuad2(ksi)*PsiQuad1(eta), PsiQuad2(ksi)*d2PsiQuad1(eta) };
	}

	vector<double> laplacianFiBiquad3(const double & ksi, const double & eta) {
//		cout << "\basisLagrangeQuad::laplacianFiBiquad3()\n";
		return { d2PsiQuad3(ksi)*PsiQuad1(eta) , PsiQuad3(ksi)*d2PsiQuad1(eta) };
	}

	vector<double> laplacianFiBiquad4(const double & ksi, const double & eta) {
	//	cout << "\basisLagrangeQuad::laplacianFiBiquad4()\n";
		return { d2PsiQuad1(ksi)*PsiQuad2(eta) , PsiQuad1(ksi)*d2PsiQuad2(eta) };
	}

	vector<double> laplacianFiBiquad5(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::laplacianFiBiquad5()\n";
		return { d2PsiQuad2(ksi)*PsiQuad2(eta) , PsiQuad2(ksi)*d2PsiQuad2(eta) };
	}

	vector<double> laplacianFiBiquad6(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::laplacianFiBiquad6()\n";
		return { d2PsiQuad3(ksi)*PsiQuad2(eta) , PsiQuad3(ksi)*d2PsiQuad2(eta) };
	}

	vector<double> laplacianFiBiquad7(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::laplacianFiBiquad7()\n";
		return { d2PsiQuad1(ksi)*PsiQuad3(eta) , PsiQuad1(ksi)*d2PsiQuad3(eta) };
	}

	vector<double> laplacianFiBiquad8(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::laplacianFiBiquad8()\n";
		return { d2PsiQuad2(ksi)*PsiQuad3(eta) , PsiQuad2(ksi)*d2PsiQuad3(eta) };
	}

	vector<double> laplacianFiBiquad9(const double & ksi, const double & eta) {
		// cout << "\basisLagrangeQuad::laplacianFiBiquad9()\n";
		return { d2PsiQuad3(ksi)*PsiQuad3(eta) , PsiQuad3(ksi)*d2PsiQuad3(eta) };
	}

	double PsiLine1(const double & ksi) {	return  (1. - ksi);	}

	double PsiLine2(const double & ksi) {	return ksi;	}

	double  dPsiLine1(const double & ksi) {	return -1.;	}

	double dPsiLine2(const double & ksi) {		return 1.;	}

	double  PsiQuad1(const double & ksi) {		return 2.*(ksi - 0.5)*(ksi - 1.);	}

	double PsiQuad2(const double & ksi) {		return -4.*ksi*(ksi - 1.);	}

	double  PsiQuad3(const double & ksi) {		return 2.*ksi*(ksi - 0.5);	}

	double dPsiQuad1(const double & ksi) {		return 4.*ksi - 3.;	}

	double  dPsiQuad2(const double & ksi) {		return -8.*ksi + 4.;	}

	double dPsiQuad3(const double & ksi) {		return 4.*ksi - 1.;	}
	double d2PsiQuad1(const double & ksi) {		return 4.;	}

	double  d2PsiQuad2(const double & ksi) {		return -8.;	}

	double d2PsiQuad3(const double & ksi) {		return 4.;	}


	double  PsiCube1(const double & ksi)	{		return (-4.5*(ksi - (1. / 3.))*(ksi - (2. / 3.))*(ksi - 1.0));	}

	double  PsiCube2(const double & ksi)	{		return  (13.5*ksi * (ksi - (2. / 3.))*(ksi - 1.0));	}

	double  PsiCube3(const double & ksi)	{		return(-13.5*ksi * (ksi - (1. / 3.))*(ksi - 1));	}

	double  PsiCube4(const double & ksi)	{		return (4.5*ksi * (ksi - (1. / 3.)) * (ksi - (2. / 3.)));	}

	double  dPsiCube1(const double & ksi)	{		return (-13.5*ksi * ksi + 18 * ksi - 5.5);	}

	double  dPsiCube2(const double & ksi)	{		return (40.5*ksi* ksi - 45 * ksi + 9.);	}

	double  dPsiCube3(const double & ksi)	{		return (-40.5 * ksi* ksi + 36 * ksi - 4.5);	}

	double  dPsiCube4(const double & ksi)	{		return (13.5 * ksi* ksi - 9. *ksi + 1.);	}

}