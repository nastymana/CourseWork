#pragma once
#include "AbstractModel.h"
class MP2withBndryLayers :
	public AbstractModel
{
public:
	MP2withBndryLayers();
	MP2withBndryLayers(const vector<double> &inData, const vector<double>& domainParams);
	~MP2withBndryLayers();
	double realSol(const double &x, const double &y);
	double BC1(const double &x, const double &y);
	double BC2(const double &x, const double &y);
	double BC3(const double &x, const double &y);
	double source(const double &x, const double &y);
	double lymbda(const double &x, const double &y);
	double gamma(const double &x, const double &y);
	vector<double> velocity(const double &x, const double &y);
	virtual double initialCond(const double &x, const double &y);
};

