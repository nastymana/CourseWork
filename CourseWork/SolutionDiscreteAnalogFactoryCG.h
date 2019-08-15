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
	void createComplexBasis(string basisRankItem); // для иерархических базисов и DG из разных базисо
	void createLagrangeBasis(string basisRankItem);
	void createLegandreBasis(string basisRankItem);

};

