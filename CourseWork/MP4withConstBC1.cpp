#include "pch.h"
#include "MP4withConstBC1.h"


MP4withConstBC1::MP4withConstBC1()
{
}

MP4withConstBC1::MP4withConstBC1(const vector<double>& inData, const vector<double>& domainParams)
{
	cout << "\n MP4withConstBC1::Constructor()\n";
	lyambdaParam = domainParams[0];
	gammaParam = domainParams[1];
	velocityParam1 = domainParams[2];
	velocityParam2 = domainParams[3];
	sourceParam = domainParams[4];
	string modelName = "MP4withConstBC1";

	addDescription(modelName);
}



MP4withConstBC1::~MP4withConstBC1()
{
}

double MP4withConstBC1::realSol(const double & x, const double & y)
{
	return 0.0;
}

double MP4withConstBC1::BC1(const double & x, const double & y)
{
	double BCvalue = 0.;
	if ((y == 0.) && (x >= 0.&&x <= 0.2)) BCvalue = 1.;
	else BCvalue = 0.;
	return BCvalue;
}

double MP4withConstBC1::BC2(const double & x, const double & y)
{
	return 0.0;
}

double MP4withConstBC1::BC3(const double & x, const double & y)
{
	return 0.0;
}

double MP4withConstBC1::source(const double & x, const double & y)
{
	return 0.0;
}

double MP4withConstBC1::lymbda(const double & x, const double & y)
{
	return lyambdaParam;
}

double MP4withConstBC1::gamma(const double & x, const double & y)
{
	return 0.0;
}

vector<double> MP4withConstBC1::velocity(const double & x, const double & y)
{
	return { velocityParam1, velocityParam2 };
}

double MP4withConstBC1::initialCond(const double & x, const double & y)
{
	return 0.0;
}
