#pragma once
#include "AbstractModel.h"
class MP7withBndryLayer :
	public AbstractModel
{/*
 для тестирования задачи конвекции-диффузии
 задача с пограничными слоями порядка O(eps`) 
 на правой и верхней границах области [-1,1]x[-1,1]

 ` eps - sourceParam
 */
public:

	MP7withBndryLayer();
	MP7withBndryLayer(const vector<double> &inData, const vector<double>& domainParams);
	virtual ~MP7withBndryLayer();
	
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

