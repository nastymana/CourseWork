#pragma once
#include "SolutionDiscreteAnalogFactory.h"
class SDAFactoryCG :
	public SDAFactory
{
public:
	SDAFactoryCG();
	~SDAFactoryCG();
	
	SLAE_CSR createSlaePortret(const Mesh& newMesh);
	basisFunction* createBasis(const string &basisItem, const int &basisRankItem);

private:
	void createComplexBasis(string basisRankItem); // ��� ������������� ������� � DG �� ������ ������
	void createLagrangeBasis(string basisRankItem);
	void createLegandreBasis(string basisRankItem);

};

