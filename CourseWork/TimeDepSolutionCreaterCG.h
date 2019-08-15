#pragma once
#include "SolutionCreater.h"
class TimeDepSolutionCreaterCG :
	public SolutionCreaterCG
{
public:
	TimeDepSolutionCreaterCG();
	virtual ~TimeDepSolutionCreaterCG();

	SolutionDiscreteAnalog* createDiscreteAnalog(const string & typeItem, 
		const string & basisItem, const int & rankItem,
		AbstractModel * newModel, Mesh * newMesh, const vector<double> &timeParams);
};

