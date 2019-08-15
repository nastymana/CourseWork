#pragma once
#include "SolutionDiscreteAnalogFactory.h"
class SDAFactoryDG:
	public SDAFactory
{
public:
	SDAFactoryDG();
	~SDAFactoryDG();
	
	virtual SLAE_CSR createSlaePortret(const Mesh& newMesh);
	virtual	basisFunction* createBasis(const string &basisItem, const int &basisRankItem);
};

