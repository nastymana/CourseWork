#pragma once
#include "SolutionCreater.h"
class SolutionCreaterCG :
	public SolutionCreater
{
public:
	SolutionCreaterCG();
	virtual ~SolutionCreaterCG();
	virtual SolutionDiscreteAnalog* createDiscreteAnalog(const string &typeItem, const string &basisItem, const int &rankItem,  AbstractModel *newModel, Mesh* newMesh);
};

