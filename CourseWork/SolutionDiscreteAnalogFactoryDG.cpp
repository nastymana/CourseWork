#include "pch.h"
#include "SolutionDiscreteAnalogFactoryDG.h"


SDAFactoryDG::SDAFactoryDG()
{
}


SDAFactoryDG::~SDAFactoryDG()
{
}
SLAE_CSR SDAFactoryDG::createSlaePortret(const Mesh &mesh)
{
	SLAE_CSR newSlae;
// void portretSLAE_forDG(const vector<int> &Elements, const vector<vector<int>> &neighbours, vector<int>& iptr, vector<int> &jptr) 
	int Nelem = mesh.elemBF.size() ;

	int counterElmsBF = 0, // кол-во БФ на эл-те + суммарное кол-во на всех соседях
		NBF_all = mesh.elemBF[Nelem-1].back(); // суммарное кол-во БФ в области - 
	cout << NBF_all << endl;
	newSlae.iptr = vector<int>(NBF_all + 2, 0);
	vector<vector<int>> bf_counter;
	vector<vector<int>> jptr_2D(NBF_all+1);

	for (int i = 0; i < Nelem; i++) {
		vector<int> help;
		// cout << mesh.elemBF[i][0] << endl;
		//counterElmsBF = mesh.elemBF[i + 1][0] - mesh.elemBF[i][0];
		for (int p = mesh.elemBF[i][0]; p <= mesh.elemBF[i].back(); p++) {
			help.push_back(p); // номера БФ самого эл-та
		}
		// coutVector(help, "help just elm", 'H');
		//coutVector(mesh.neighbours[i], "neighbours[i]", 'H');
		for (int k = 1; k < mesh.nghbrs[i].size(); k++) {
			//cout << k << ":: " << endl;

			if (mesh.nghbrs[i][k] != -1) {
			//	counterElmsBF += mesh.elemBF[mesh.nghbrs[i][k+1] + 1][0] - mesh.elemBF[mesh.nghbrs[i][k+1]][0];
				//for (int p = mesh.elemBF[mesh.nghbrs[i][k+1]][0]; p < mesh.elemBF[mesh.nghbrs[i][k+1] + 1][0]; p++) {
				for (int p = 0; p< mesh.elemBF[mesh.nghbrs[i][k]].size();p++){
				//	cout << p;
					help.push_back(p+ mesh.elemBF[mesh.nghbrs[i][k]][0]);
				}
			}
		}
	//	coutVector(help, "help +new neighbour elm", 'H');
		sort(help.begin(), help.end());
		//coutVector(help, "help +new neighbour elm", 'H');
		//cout << "help.size()=" << help.size() << endl;
		for (int j = mesh.elemBF[i][0] + 1; j <= mesh.elemBF[i].back()+1; j++) {
			
			newSlae.iptr[j] = newSlae.iptr[j - 1] + help.size();
			jptr_2D[j - 1] = help;
		}
	}
	//coutVector(newSlae.iptr, "iptr", 'H');
//	coutVector(jptr_2D, "jptr2D");
	
	reshape_vector(jptr_2D, newSlae.jptr);
	
	/*for (int i = 0; i < int(jptr_2D.size()); i++) {
		
		newSlae.iptr.push_back(int(newSlae.jptr.size()));
		for (int j = 0; j<jptr_2D[i].size(); j++) {
			newSlae.jptr.push_back(jptr_2D[i][j]);
		}
	}*/
	//coutVector(newSlae.jptr, "newSlae.jptr", 'H');
//	newSlae.iptr.push_back(int(newSlae.jptr.size()));
	newSlae.A = vector<double>(newSlae.jptr.size());
	newSlae.b = vector<double>(newSlae.iptr.size() - 1);// Проверить!!!

	
	return newSlae;
}

basisFunction * SDAFactoryDG::createBasis(const string &basisItem, const int &basisRankItem)
{
	basisFunction* basis;
	cout << "\nSolutionDiscreteAnalogFactoryCG::createBasis()\n";
	if (basisItem == ("lagrange")) {
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
