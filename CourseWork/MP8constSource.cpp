#include "pch.h"
#include "MP8constSource.h"


MP8constSource::MP8constSource()
{
}


MP8constSource::~MP8constSource()
{
}

MP8constSource::MP8constSource(const vector<double>& inData, const vector<double>& domainParams)
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

double MP8constSource::realSol(const double & x, const double & y)
{
	return 1.0;
}

double MP8constSource::BC1(const double & x, const double & y)
{
	return 1.0;
}

double MP8constSource::BC2(const double & x, const double & y)
{
	return 0.0;
}

double MP8constSource::BC3(const double & x, const double & y)
{
	return 0.0;
}

double MP8constSource::source(const double & x, const double & y)
{
	return 0.0;
}

double MP8constSource::lymbda(const double & x, const double & y)
{
	return lyambdaParam;
}

double MP8constSource::gamma(const double & x, const double & y)
{
	return gammaParam;
}

vector<double> MP8constSource::velocity(const double & x, const double & y)
{
	return {velocityParam1, velocityParam2};
}

double MP8constSource::initialCond(const double & x, const double & y)
{
	return 0.0;
}
