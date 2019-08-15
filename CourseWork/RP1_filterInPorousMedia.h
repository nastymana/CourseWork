#pragma once
#include "AbstractModel.h"
class RP1_filterInPorousMedia :
	public AbstractModel
{
	/*
	нестационарная задача фильтрация в пористой среде
	область 0.1х0.01 [0., 0.1]x[0., 0.01]
	D = 10^-9 - 10^-7

	*/
public:
	RP1_filterInPorousMedia();
	RP1_filterInPorousMedia(const vector<double> &inData, const vector<double>& domainParams);
	virtual ~RP1_filterInPorousMedia();

	virtual double realSol(const double &x, const double &y);
	virtual double BC1(const double &x, const double &y);
	virtual double BC2(const double &x, const double &y);
	virtual double BC3(const double &x, const double &y);
	virtual double initialCond(const double &x, const double &y);
	virtual double source(const double &x, const double &y);
	virtual double lymbda(const double &x, const double &y);
	virtual double gamma(const double &x, const double &y);
	virtual vector<double> velocity(const double &x, const double &y);


	
};

