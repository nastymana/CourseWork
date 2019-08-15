#pragma once
#include "AbstractModel.h"
class MP4withConstBC1 :
	public AbstractModel{
	// область [0, 1]x[0, 1]
	// Lyambda = 1, gamma= 0., a=(1.,2) 
	// BC1: u=1 (x e(0.,1), y=0) u = 1 на остальной части границы
public:
	MP4withConstBC1();
	MP4withConstBC1(const vector<double> &inData, const vector<double>& domainParams);
	virtual ~MP4withConstBC1();	
	
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

