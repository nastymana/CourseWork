#pragma once
class SolutionCreater
{ // аналог - PizzaStore
public:
	SolutionCreater();
	virtual ~SolutionCreater();
	// pure function
	virtual SolutionDiscreteAnalog* createDiscreteAnalog(const string &typeItem, const string &basisItem, const int &rankItem, AbstractModel *newModel, Mesh* newMesh) = 0;
	
};

