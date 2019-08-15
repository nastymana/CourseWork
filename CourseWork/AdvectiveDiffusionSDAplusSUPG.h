#pragma once
#include "SolutionDiscreteAnalog.h"
class AdvecDifSDAplusSUPG :
	public SolutionDiscreteAnalogCG
{
	// возможно и правда не нужно делать два этапа наследования с промежуточным делением на CG | DG
public:
	AdvecDifSDAplusSUPG(AbstractModel *newModel, Mesh* newMesh, SDAFactory *newDAfitches, const vector<double> &timeParams);
	AdvecDifSDAplusSUPG(AbstractModel *newModel, Mesh* newMesh, SDAFactory *newDAfitches);
	~AdvecDifSDAplusSUPG();

	void calcLocalMatrices( const currentParamsValue &Params, const quadrElem &curElem,
		vector<vector<double>> &resultA, vector<double> &result);
private:
	void calcLocalSUPGMatrix1stOrderBasis(const int& nTime, const currentParamsValue &Params, const quadrElem &Elem,
		const vector<vector<vector<double>>> &gradFi,
		vector<vector<double>> &resultSUPGmatrix, vector<double> &resultSUPGvector);
	void calcLocalSUPGMatrix2ndOrderBasis(const int& nTime, const currentParamsValue &Params, const quadrElem &curElem,
		const vector<vector<vector<double>>> &gradFi, const vector<vector<vector<double>>> &laplacFi,
		vector<vector<double>> &resultSUPGmatrix, vector<double> &resultSUPGvector);
	
	void calcRightVectorNewTimeStep(const vector<double>& prevStepU);
	double calcSUPGtau(quadrElem Elem, currentParamsValue Params, double Ctau);


};

