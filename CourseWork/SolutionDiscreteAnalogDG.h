#pragma once
#include "SolutionDiscreteAnalog.h"
class SolutionDiscreteAnalogDG :
	public SolutionDiscreteAnalog
{
public:
	SolutionDiscreteAnalogDG();
	SolutionDiscreteAnalogDG(AbstractModel *newModel, Mesh *newMesh,
		SDAFactory *newDAfitches, const vector<double>& timeParams);
	SolutionDiscreteAnalogDG(AbstractModel *newModel, Mesh *newMesh, SDAFactory*  newDAfitches);
	
	~SolutionDiscreteAnalogDG();
	
	void defineSolutionProcess(const string &basisItem, const int& rankBasis);

	// сборка SLAE
	virtual void assembleSlae(const int &nTimeStep);
	virtual void  BoundCondDerichlet(const int &nBndry, const int &nTimeStep);
	
	void calcLocalBcDerichletMatrices(const int nBndry, const vector<vector<double>> &Fi, const vector<vector<vector<double>>> &gradfi, quadrElem &elem,
		currentParamsValue &params, vector<vector<double>> &localAd, vector<double> &localBd);
	//virtual void calcLocalMatrices(const currentParamsValue &curParams, const quadrElem &curElem,
	//	vector<vector<double>> &resultA, vector<double> &result);// какое уравнение - конве-диффузии, конв-дифф-реак или др

	void caclLocalNumerFluxMatrices(int nBnd, const quadrElem &elem,
		const currentParamsValue& params, vector<vector<double>> &resultEin, vector<vector<double>> &resultEex);
	
	void addNumFluxToSLAE();
	
	double IPstabParamMU(const int nBnd, const quadrElem &elem, const currentParamsValue& params);
	
};

