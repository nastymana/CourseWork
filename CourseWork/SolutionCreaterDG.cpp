#include "pch.h"
#include "SolutionCreaterDG.h"



SolutionCreaterDG::SolutionCreaterDG()
{
}

SolutionCreaterDG::~SolutionCreaterDG()
{
}

SolutionDiscreteAnalog* SolutionCreaterDG::createDiscreteAnalog(const string &typeItem, const string &basisItem, const int &rankItem,
	AbstractModel *newModel, Mesh* newMesh)
{	cout << "\nSolutionCreaterDG::createDiscreteAnalog()\n";
	SolutionDiscreteAnalog* SDA;
	SDAFactory* SDAfactory = new SDAFactoryDG();
	//if (typeItem=="AdvDif") {
		 SDA = new AdvecDiffusDG_SDA(newModel, newMesh,SDAfactory);
	//}
	//else SDA = new SolutionDiscreteAnalogDG();

	SDA->defineSolutionProcess(basisItem, rankItem);
	SDA->assembleSlae(-1);
	

	return SDA;
}
