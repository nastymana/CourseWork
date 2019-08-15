#include "pch.h"
#include "SolutionCreaterCG.h"


SolutionCreaterCG::SolutionCreaterCG()
{
	cout << "\nSolutionCreaterCG::Constructor()\n";
}


SolutionCreaterCG::~SolutionCreaterCG()
{
	cout << "\nSolutionCreaterCG::Destructor()\n";
}

SolutionDiscreteAnalog* SolutionCreaterCG::createDiscreteAnalog(const string &typeItem, const string &basisItem, const int &rankItem,
	AbstractModel *newModel, Mesh* newMesh)
{
	cout << "\nSolutionCreaterCG::createDiscreteAnalog()\n";

	SolutionDiscreteAnalog* SDA;
	SDAFactory* SDAFactory = new SDAFactoryCG();

	if (typeItem =="AdvDif") {
		cout << "\n There is Advective Diffusion SDA option \n";
		SDA = new AdvecDiffusCG_SDA(newModel, newMesh, SDAFactory);
	}
	else if (typeItem =="AdvDifSUPG") {
		cout << "\n There is Advective Diffusion plus SUPG SDA option \n";
		 SDA = new AdvecDifSDAplusSUPG(newModel, newMesh, SDAFactory);
	}
	else {
		cout << "\n There is not coinsidense with options \n";
		SDA = new AdvecDiffusCG_SDA(newModel, newMesh, SDAFactory);
	}

	SDA->defineSolutionProcess(basisItem, rankItem);
	SDA->assembleSlae(-1);
	return SDA;
}
