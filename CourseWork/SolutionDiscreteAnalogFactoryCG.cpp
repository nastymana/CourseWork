#include "pch.h"
#include "SolutionDiscreteAnalogFactoryCG.h"


SDAFactoryCG::SDAFactoryCG()
{
	cout << "\nSolutionDiscreteAnalogFactoryCG::Constructor()\n";
}


SDAFactoryCG::~SDAFactoryCG()
{
	cout << "\nSolutionDiscreteAnalogFactoryCG::~Destructor()\n";
}

SLAE_CSR SDAFactoryCG::createSlaePortret(const Mesh& mesh)
{
	cout << "\nSolutionDiscreteAnalogFactoryCG::createSlaePortret()\n";
//	void PortretSLAU(const vector<vector<int>> &num_basis, vector<int> &jg1d, vector<int>  &ig1d, int nnode)
	int nelem = mesh.elemBF.size(),
		nBFinSLAE = mesh.elemBF[nelem - 1].back() + 1;
	//	cout << " = " << nelem << endl;
	//	coutVector(num_basis, "num_basis");
	SLAE_CSR slae;
	vector<set<int>> num_basis_row, jg2d(nBFinSLAE);
	for (int i = 0; i < nelem; i++){
		set<int> help;
		// в этом массиве хранится информация о гранях
		for (int j = 0; j < int(mesh.elemBF[i].size()); j++) {
			help.insert(mesh.elemBF[i][j]);
		}
		//coutSet(help, "help", 'H');
		num_basis_row.push_back(help);
	}

	for (int i = 0; i < nelem; i++){
		int k = 0;
		set<int> help;
		//coutSet(num_basis_row[i], "num_basis", 'H');
		for (set<int>::iterator it = num_basis_row[i].begin(); it != num_basis_row[i].end(); it++){
			help.insert(*it);
		}
		for (set<int>::iterator it = num_basis_row[i].begin(); it != num_basis_row[i].end(); it++){
			jg2d[*it].insert(help.begin(), help.end());
			//		cout << *it;
			//	coutSet(jg2d[*it], "jg2d", 'H');
		}
	}
	///	coutVectorsets(jg2d, " Portret");
	for (int i = 0; i < int(jg2d.size()); i++){
		set<int>::iterator it = jg2d[i].begin();
		slae.iptr.push_back(int(slae.jptr.size()));
		for (it; it != jg2d[i].end(); it++){
			slae.jptr.push_back(*it);
		}
	}
	slae.iptr.push_back(int(slae.jptr.size()));
	slae.A = vector<double>(slae.jptr.size());
	slae.b = vector<double>(slae.iptr.size() - 1);// Проверить!!!
	return slae;
}



basisFunction * SDAFactoryCG::createBasis(const string &basisItem, const int &basisRankItem)
{
	basisFunction* basis;
	cout << "\nSolutionDiscreteAnalogFactoryCG::createBasis()\n";
	if (basisItem==("lagrange")) {
		basis = new Lagrange::basisLagrangeQuad(basisRankItem);
	}
	/*else if (basisItem._Equal("legandre")) {
		
	}
	else if (basisItem._Equal("lagrange")) {
		return nullptr;
	}*/
	else basis = new Lagrange::basisLagrangeQuad(basisRankItem);
	return basis;
}

void SDAFactoryCG::createComplexBasis(string basisRankItem)
{
	cout << "\nSolutionDiscreteAnalogFactoryCG::createComplexBasis()\n";
}

void SDAFactoryCG::createLagrangeBasis(string basisRankItem)
{
	cout << "\nSolutionDiscreteAnalogFactoryCG::createLagrangeBasis()\n";
}

void SDAFactoryCG::createLegandreBasis(string basisRankItem)
{
	cout << "\nSolutionDiscreteAnalogFactoryCG::createLegandreBasis()\n";
}
