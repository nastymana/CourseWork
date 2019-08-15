#include "pch.h"
#include "AdvecDiffusDG_SDA.h"


AdvecDiffusDG_SDA::AdvecDiffusDG_SDA()
{
}

AdvecDiffusDG_SDA::AdvecDiffusDG_SDA(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches)
{
	cout << "/AdvecDiffusDG_SDA:  Constructor/n";
	sdaFactory = newDAfitches;
	mesh = newMesh;
	model = newModel;	
}

AdvecDiffusDG_SDA::~AdvecDiffusDG_SDA()
{
}

void AdvecDiffusDG_SDA::calcLocalMatrices(const currentParamsValue & curParams,
	const quadrElem & curElem, vector<vector<double>>& resultA, vector<double>& resultB)
{
	cout << "\nAdvectiveDiffusionSDA::calcLocalMatrices()\n";
	vector<vector<double>> localStiffMatrix, localAdvecMatrix, LocalMassMatrix;
	vector<double> localbVector, localICvector;

	vector<vector<double>> FI = basis->getFiValuesInGaussPoints();
	vector<vector<vector<double>>> gradFI = basis->getGradFiValuesInGaussPoints();

	SolutionDiscreteAnalog::calcLocalAdvectiveMatrix(curParams.velocity, curElem, FI, gradFI, localAdvecMatrix);
	SolutionDiscreteAnalog::calcLocalStiffnessMatrix(curParams.lymbda, curElem, gradFI, localStiffMatrix);
	SolutionDiscreteAnalog::calcLocalRightVector(curElem, FI, localbVector);


	resultA = localStiffMatrix + localAdvecMatrix;
	resultB = localbVector;

}
