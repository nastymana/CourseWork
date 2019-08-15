#include "pch.h"
#include "AdvectiveDiffusionSDA.h"
#include <map>


AdvecDiffusCG_SDA::AdvecDiffusCG_SDA()
{
}

AdvecDiffusCG_SDA::AdvecDiffusCG_SDA(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches, const vector<double>& timeParams)
{
	T0 = timeParams[0];
	Tn = timeParams[1];
	NtimeSteps = timeParams[2];
	timeStep = (Tn - T0) / NtimeSteps;
	sdaFactory = newDAfitches;
	mesh = newMesh;
	model = newModel;
}

AdvecDiffusCG_SDA::AdvecDiffusCG_SDA(AbstractModel * newModel, Mesh * newMesh, SDAFactory * newDAfitches){
	cout << "/AdvecDiffusÑG_SDA:  Constructor/n";
	sdaFactory = newDAfitches;
	mesh = newMesh;
	model = newModel;

}

AdvecDiffusCG_SDA::~AdvecDiffusCG_SDA()
{
	cout << "\nAdvectiveDiffusionSDA::~Destructor()\n";
}

void AdvecDiffusCG_SDA::calcLocalMatrices( const currentParamsValue &curParams, const quadrElem &curElem,
		vector<vector<double>> &resultA,  vector<double> &resultb) {
	cout << "\nAdvectiveDiffusionSDA::calcLocalMatrices()\n";
	vector<vector<double>> localStiffMatrix, localAdvecMatrix, LocalMassMatrix;
	vector<double> localbVector, localICvector;

	vector<vector<double>> FI = basis->getFiValuesInGaussPoints();
	vector<vector<vector<double>>> gradFI = basis->getGradFiValuesInGaussPoints();
	
	SolutionDiscreteAnalog::calcLocalAdvectiveMatrix(curParams.velocity, curElem, FI, gradFI, localAdvecMatrix);
	SolutionDiscreteAnalog::calcLocalStiffnessMatrix(curParams.lymbda, curElem, gradFI, localStiffMatrix);
	SolutionDiscreteAnalog::calcLocalRightVector(curElem, FI, localbVector);
	
	resultA = localStiffMatrix + localAdvecMatrix; 
	resultb = localbVector;

}



