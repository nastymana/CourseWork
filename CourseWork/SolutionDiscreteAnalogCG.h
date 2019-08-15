#pragma once
#include "SolutionDiscreteAnalog.h"
class SolutionDiscreteAnalogCG :
	public SolutionDiscreteAnalog
{
public:
	SolutionDiscreteAnalogCG();
	SolutionDiscreteAnalogCG(AbstractModel *newModel, Mesh *newMesh,
		SDAFactory *newDAfitches, const vector<double>& timeParams);
	SolutionDiscreteAnalogCG(AbstractModel *newModel, Mesh *newMesh, SDAFactory *newDAfitches);

	~SolutionDiscreteAnalogCG();

	// здесь будем определять какой базис использовать, пока только для квадратов.
	// 	virtual SLAE* createSLAE()=0;
	virtual void defineSolutionProcess(const string &basisItem, const int& rankBasis);

	// сборка SLAE. если nTimeStep==0 - то считаем начальное условие, 
	// >0 - пересчитываем только правую часть , <0 - не считаем, задача стацион. 
	virtual void assembleSlae(const int &nTimeStep);
	
	virtual void calcLocalMatrices( const currentParamsValue &curParams, const quadrElem &curElem,
		vector<vector<double>> &resultA, vector<double> &result)=0;// какое уравнение - конве-диффузии, конв-дифф-реак или др
	
	void BoundCondDerichlet(const int &nBoundary, const int&nTimeStep);

private:
};

