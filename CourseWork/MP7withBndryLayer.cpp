#include "pch.h"
#include "MP7withBndryLayer.h"


MP7withBndryLayer::MP7withBndryLayer()
{
}


MP7withBndryLayer::MP7withBndryLayer(const vector<double>& inData, const vector<double>& domainParams)
{
	cout << "\n MP7withBndryLayer::Constructor()\n";
	lyambdaParam = domainParams[0];
	gammaParam = domainParams[1];
	velocityParam1 = domainParams[2];
	velocityParam2 = domainParams[3];
	sourceParam = lyambdaParam;// domainParams[4]; // ë - lyambda
	string modelName = "MP7withBndryLayer and \n u(x,y) = cos(PI*(x + y))*(1. - exp((x - 1.) / ë))*(1. - exp((y - 1.) / ë));";

	addDescription(modelName);
}

MP7withBndryLayer::~MP7withBndryLayer()
{
}

double MP7withBndryLayer::realSol(const double & x, const double & y)
{
	if (((x >= -1.) && (x <= 1.)) && ((y >= -1.) && (y <= 1.)))
	return cos(PI*(x + y))*(1. - exp((x - 1.) / sourceParam))*(1. - exp((y - 1.) / sourceParam));
}

double MP7withBndryLayer::BC1(const double & x, const double & y)
{
	if (((x >= -1.) && (x <= 1.)) && ((y >= -1.) && (y <= 1.)))
	return cos(PI*(x+y))*(1. - exp((x-1.)/sourceParam))*(1. - exp((y - 1.) / sourceParam));
}

double MP7withBndryLayer::BC2(const double & x, const double & y)
{
	return 0.0;
}

double MP7withBndryLayer::BC3(const double & x, const double & y)
{
	return 0.0;
}

double MP7withBndryLayer::source(const double & x, const double & y)
{
	double valueGradx = 0., valueGrady, valueLaplasX, valueLaplasY, valueU;
	if (((x >= -1.) && (x <= 1.)) && ((y >= -1.) && (y <= 1.))) {
		
		valueLaplasX = (-PI * PI)* cos(PI*(x + y))*(1. - exp((x - 1.) / sourceParam))*(1. - exp((y - 1.) / sourceParam))
			+ (-2.*PI)* sin(PI*(x + y))*((-1. / sourceParam)* exp((x - 1.) / sourceParam))*(1. - exp((y - 1.) / sourceParam))
			+ cos(PI*(x + y))*((-1. / pow(sourceParam, 2))* exp((x - 1.) / sourceParam))*(1. - exp((y - 1.) / sourceParam));
		
		valueLaplasY = (-PI * PI)* cos(PI*(x + y))*(1. - exp((x - 1.) / sourceParam))*(1. - exp((y - 1.) / sourceParam))
			+ (-2.*PI)* sin(PI*(x + y))*((-1. / sourceParam)* exp((y - 1.) / sourceParam))*(1. - exp((x - 1.) / sourceParam))
			+ cos(PI*(x + y))*((-1. / pow(sourceParam, 2))* exp((y - 1.) / sourceParam))*(1. - exp((x - 1.) / sourceParam));
		
		valueGradx = (-PI)* sin(PI*(x + y))*(1. - exp((x - 1.) / sourceParam))*(1. - exp((y - 1.) / sourceParam))
			+ cos(PI*(x + y))*((-1. / sourceParam)* exp((x - 1.) / sourceParam))*(1. - exp((y - 1.) / sourceParam));
		
		valueGrady = (-PI)* sin(PI*(x + y))*(1. - exp((x - 1.) / sourceParam))*(1. - exp((y - 1.) / sourceParam))
			+ cos(PI*(x + y))*((-1. / sourceParam)* exp((y - 1.) / sourceParam))*(1. - exp((x - 1.) / sourceParam));
		valueU = realSol(x, y);


	}
	vector<double> v = velocity(x, y);
			return (-lymbda(x,y)*(valueLaplasX+valueLaplasY)+v[0]*valueGradx+v[1]*valueGrady +gamma(x,y)*valueU);
}

double MP7withBndryLayer::lymbda(const double & x, const double & y)
{
	return lyambdaParam;
}

double MP7withBndryLayer::gamma(const double & x, const double & y)
{
	return gammaParam;
}

vector<double> MP7withBndryLayer::velocity(const double & x, const double & y)
{
	return { velocityParam1, velocityParam2 };
}

double MP7withBndryLayer::initialCond(const double & x, const double & y)
{
	return 0.0;
}
