#include "pch.h"
#include "RP1_filterInPorousMedia.h"


RP1_filterInPorousMedia::RP1_filterInPorousMedia()
{
	lyambdaParam = 1.;
	gammaParam = 0.;
	velocityParam1 = 1.;
	velocityParam2 = 0.;
	sourceParam = 0.;
	string modelName = "RP1_filterInPorousMedia";
	addDescription(modelName);

}

RP1_filterInPorousMedia::RP1_filterInPorousMedia(const vector<double>& inData, const vector<double>& domainParams)
{
	cout << "\n CG_MP1::default Constructor()\n";
	lyambdaParam = domainParams[0];
	gammaParam = domainParams[1];
	velocityParam1 = domainParams[2];
	velocityParam2 = domainParams[3];
	sourceParam = domainParams[4];
	string modelName = "RP1_filterInPorousMedia";
	addDescription(modelName);

}


RP1_filterInPorousMedia::~RP1_filterInPorousMedia()
{
}

double RP1_filterInPorousMedia::realSol(const double & x, const double & y)
{
	return 0.0;
}

double RP1_filterInPorousMedia::initialCond(const double & x, const double & y)
{ 
	// равномерное распределение вещества в начальный момент времени
	if((y >= 0. && y <= 0.01) && (x >= 0. && x <= 0.1)) return 0.1;
}

double RP1_filterInPorousMedia::BC1(const double & x, const double & y)
{
	double BCvalue = 0.;
	if((y==0)&&(x>=0.&&x<=0.1)) BCvalue = 1.;
	else if ((y == 0.01) && (x >= 0.&&x <= 0.1)) BCvalue = 0.3;
	return BCvalue;
}

double RP1_filterInPorousMedia::BC2(const double & x, const double & y)
{
	double BCvalue = 0.;
	// условие непротекания на левой и правой границах
	if ((x == 0) || (x == 0.1) && (y >= 0.&&y <= 0.01)) BCvalue = 0; 

	return BCvalue;
}

double RP1_filterInPorousMedia::BC3(const double & x, const double & y)
{
	return 0.0;
}

double RP1_filterInPorousMedia::source(const double & x, const double & y)
{
	double value = 0.;
	if (x == 0.1)value = 1.;
	return value;
}

double RP1_filterInPorousMedia::lymbda(const double & x, const double & y)
{
	return lyambdaParam;
}

double RP1_filterInPorousMedia::gamma(const double & x, const double & y)
{
	return 0.0;
}

vector<double> RP1_filterInPorousMedia::velocity(const double & x, const double & y)
{
	return {velocityParam1, 0.};
}
