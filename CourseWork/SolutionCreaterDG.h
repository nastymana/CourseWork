#pragma once
#include "SolutionCreater.h"
class SolutionCreaterDG :
	public SolutionCreater
{
public:
	SolutionCreaterDG();
	~SolutionCreaterDG();
	virtual SolutionDiscreteAnalog* createDiscreteAnalog(const string &typeItem, const string &basisItem, const int &rankItem, AbstractModel *newModel, Mesh* newMesh);
};

