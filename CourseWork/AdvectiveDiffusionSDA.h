#pragma once
#include "SolutionDiscreteAnalog.h"
class AdvecDiffusCG_SDA :
	public SolutionDiscreteAnalogCG
{
public:
	AdvecDiffusCG_SDA();
	//*AdvectiveDiffusionSDA(const string &item);
	AdvecDiffusCG_SDA(AbstractModel *newModel, Mesh *newMesh,
		SDAFactory *newDAfitches, const vector<double>& timeParams);
	AdvecDiffusCG_SDA(AbstractModel *newModel, Mesh* newMesh, SDAFactory *newDAfitches);

	virtual ~AdvecDiffusCG_SDA();
		
	void calcLocalMatrices(const currentParamsValue &curParams, const quadrElem &curElem,
		vector<vector<double>> &resultA, vector<double> &result);
	
};

