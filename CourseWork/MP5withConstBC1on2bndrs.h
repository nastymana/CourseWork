#pragma once
#include "AbstractModel.h"
class MP5withConstBC1on2bndrs :
	public AbstractModel
{
public:
	MP5withConstBC1on2bndrs();
	~MP5withConstBC1on2bndrs();
	virtual double initialCond(const double &x, const double &y) = 0;
};

