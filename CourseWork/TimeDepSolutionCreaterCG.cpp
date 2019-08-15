#include "pch.h"
#include "TimeDepSolutionCreaterCG.h"


TimeDepSolutionCreaterCG::TimeDepSolutionCreaterCG()
{
}


TimeDepSolutionCreaterCG::~TimeDepSolutionCreaterCG()
{
}

SolutionDiscreteAnalog * TimeDepSolutionCreaterCG::createDiscreteAnalog(const string & typeItem, const string & basisItem, const int & rankItem,
	AbstractModel * newModel, Mesh * newMesh, const vector<double> &timeParams)
{
	cout << "\nSolutionCreaterCG::createDiscreteAnalog()\n";
	SolutionDiscreteAnalog* SDA;
	SDAFactory* SDAFactory = new SDAFactoryCG();

	if (typeItem == "AdvDif") {
		cout << "\n There is Advective Diffusion SDA option \n";
		SDA = new AdvecDiffusCG_SDA(newModel, newMesh, SDAFactory, timeParams);
	}
	else if (typeItem == "AdvDifSUPG") {
		cout << "\n There is Advective Diffusion plus SUPG SDA option \n";
		SDA = new AdvecDifSDAplusSUPG(newModel, newMesh, SDAFactory, timeParams);
	}
	else {
		cout << "\n There is not coinsidense with options \n";
		SDA = new AdvecDiffusCG_SDA(newModel, newMesh, SDAFactory);
	}

	SDA->defineSolutionProcess(basisItem, rankItem);
	SDA->assembleSlae(0);
	
	
	return SDA;
}
