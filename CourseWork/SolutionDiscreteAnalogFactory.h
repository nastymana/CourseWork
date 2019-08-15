#pragma once
class SDAFactory
{
public:
	SDAFactory();
	virtual ~SDAFactory();
	virtual SLAE_CSR createSlaePortret(const Mesh& newMesh) = 0;
	
	virtual	basisFunction* createBasis(const string &basisItem, const int &basisRankItem) = 0;

private:

	
	
};

