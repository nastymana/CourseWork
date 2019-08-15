#include "pch.h"
#include "MP1biquadAS.h"


MP1biquadAS::MP1biquadAS()
{
	cout << "\n CG_MP1::default Constructor()\n";
	lyambdaParam = 1.;
	gammaParam = 0.;
	velocityParam1 = 1.;
	velocityParam2 = 0.;
	sourceParam = 0.;
	string modelName = "MP1biquadAS";
	addDescription(modelName);
}


MP1biquadAS::MP1biquadAS(const vector<double> &inData, const vector<double>& domainParams)
{
	cout << "\n CG_MP1::Constructor()\n";
	lyambdaParam = domainParams[0];
	gammaParam = domainParams[1];
	velocityParam1 = domainParams[2];
	velocityParam2 = domainParams[3];
	sourceParam = domainParams[4];
	string modelName = "MP1biquadAS";
	
	addDescription(modelName);
}

MP1biquadAS::~MP1biquadAS()
{
	cout << "\n CG_MP1::Destructor()\n";
}

double MP1biquadAS::realSol(const double & x, const double & y)
{
	return (x * (1.0 - x)*y * (1.0 - y));
}

double MP1biquadAS::BC1(const double & x, const double & y){
	if(((x==0.)||(x==1.)) || ((y==0.)||(y==1.))) return 0.0;
}

double MP1biquadAS::BC2(const double & x, const double & y){
	return 0.0;
}

double MP1biquadAS::BC3(const double & x, const double & y){
	return 0.0;
}

double MP1biquadAS::source(const double & x, const double & y){
	double k, m1, m2, g;
	if (((x >= 0.) && (x <= 1.)) && ((y >= 0.) && (y <= 1.))) {
		vector<double> v = velocity(x, y);
	//	coutVector(v, "v", 'H');
	//	cout << lymbda(x, y) << ",  "<< -lymbda(x, y)<<endl;
		k = -lymbda(x, y)*(-2.*y * (1. - y) - 2. * x * (1. - x));
		m1 = v[0] * ((1. - 2. * x)*y * (1 - y));
		m2 = v[1] * ((1. - 2.*y)*x * (1. - x));
		g = gamma(x, y)*realSol(x, y);

			
	}
	//cout << "k=" << k << ", m1=" << m1<<", m2="<<m2 << ", g = " << g << endl;
	return k + m1+m2 + g;
}

double MP1biquadAS::lymbda(const double & x, const double & y){
	
	return lyambdaParam;
}

double MP1biquadAS::gamma(const double & x, const double & y){
	return gammaParam;
}

vector<double> MP1biquadAS::velocity(const double & x, const double & y){
	return { velocityParam1, velocityParam2 };
}

double MP1biquadAS::initialCond(const double & x, const double & y)
{
	return 0.0;
}
