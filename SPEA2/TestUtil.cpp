#include "TestUtil.h"



TestUtil::TestUtil()
{
}


TestUtil::~TestUtil()
{
}

void TestUtil::run()
{
	PopulationReader* p = new PopulationReader();	
	vector<Individual*> individuals = p->loadFromArquive("individuals.txt");
	
	/*if (individuals.size() > 0) {
		HyperVolumeCalculator h;
		vector<Individual*>nonDominatedSolutions = nonDominatedSol(individuals);		
		cout << h.calculateForTwoObjectiveWithExtremalReferences(nonDominatedSolutions, individuals) << endl;
		dump(nonDominatedSolutions);
		clearVector(nonDominatedSolutions);
		clearVector(individuals);
	}*/


	if (individuals.size() > 0) {
		HyperVolumeCalculator h;
		vector<Individual*>nonDominatedSolutions = nonDominatedSol(individuals);		
		cout << h.calculateSpread(nonDominatedSolutions) << endl;		
		/*dump(nonDominatedSolutions);*/
		clearVector(nonDominatedSolutions);
		clearVector(individuals);
	}
	delete p;
}

vector<Individual*> TestUtil::nonDominatedSol(vector<Individual*> individuals) {

	vector<Individual*> nonDominated;
	unsigned j = 0;

	for (unsigned i = 0; i < individuals.size(); i++) {

		bool oneDominated = false;

		for (j = 0; j < individuals.size(); j++) {
			if (i != j) {
				if (isDominated(individuals[i], individuals[j])) {
					oneDominated = true;
					break;
				}
			}
		}

		if (!oneDominated) {
			Individual* id = new Individual(*(individuals.at(i)));
			nonDominated.push_back(id);//Verify
		}

	}

	return nonDominated;
}

bool TestUtil::isDominated(Individual * one, Individual * two)
{
	bool anyDominate = false;
	ZDT objective;

	for (unsigned int i = 0; i < one->getAptidao().size(); i++) {

		if (objective.isMinimization(i)) {

			if (one->getAptidao()[i] < two->getAptidao()[i])
				return false;
			else if (one->getAptidao()[i] > two->getAptidao()[i])
				anyDominate = true;
		}
		else {

			if (one->getAptidao()[i] > two->getAptidao()[i])
				return false;
			else if (one->getAptidao()[i] < two->getAptidao()[i])
				anyDominate = true;
		}
	}

	return anyDominate;
}

void TestUtil::clearVector(vector<Individual*>& v)
{
	for (size_t i = 0; i < v.size(); i++) {
		delete v[i];
	}
	v.clear();
}

void TestUtil::dump(vector<Individual*> individuals) {

	/*cout << "\n Population: " << endl;*/
	for (unsigned int i = 0; i < individuals.size(); i++) {

		cout << individuals.at(i)->getAptidao()[0] << ";" << individuals.at(i)->getAptidao()[1] << endl;
	}
}