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

	// ����� ����� ���������� ����� ����� ������������, ���� ������ ��� ���������.
	// 	virtual SLAE* createSLAE()=0;
	virtual void defineSolutionProcess(const string &basisItem, const int& rankBasis);

	// ������ SLAE. ���� nTimeStep==0 - �� ������� ��������� �������, 
	// >0 - ������������� ������ ������ ����� , <0 - �� �������, ������ �������. 
	virtual void assembleSlae(const int &nTimeStep);
	
	virtual void calcLocalMatrices( const currentParamsValue &curParams, const quadrElem &curElem,
		vector<vector<double>> &resultA, vector<double> &result)=0;// ����� ��������� - �����-��������, ����-����-���� ��� ��
	
	void BoundCondDerichlet(const int &nBoundary, const int&nTimeStep);

private:
};

