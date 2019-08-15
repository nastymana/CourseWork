#pragma once
#include "AbstractModel.h"
class MP1biquadAS :
	public AbstractModel
{
public:
	MP1biquadAS();
	MP1biquadAS(const vector<double> &inData, const vector<double>& domainParams);
	virtual ~MP1biquadAS();

	virtual double realSol(const double &x, const double &y);
	double BC1(const double &x, const double &y);
	double BC2(const double &x, const double &y);
	double BC3(const double &x, const double &y);
	double source(const double &x, const double &y);
	double lymbda(const double &x, const double &y);
	double gamma(const double &x, const double &y);
	vector<double> velocity(const double &x, const double &y);
	virtual double initialCond(const double &x, const double &y) ;

};

