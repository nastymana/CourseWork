#pragma once
class AdvecDiffusDG_SDA :
	public SolutionDiscreteAnalogDG
{
public:
	AdvecDiffusDG_SDA();
	AdvecDiffusDG_SDA(AbstractModel *newModel, Mesh *newMesh, SDAFactory*  newDAfitches);
	~AdvecDiffusDG_SDA();

	void calcLocalMatrices(const currentParamsValue & curparams, const quadrElem & curElem,
		vector<vector<double>>& resulta, vector<double>& resultB);

};

