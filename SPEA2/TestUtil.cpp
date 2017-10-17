#include "TestUtil.h"



TestUtil::TestUtil()
{
}


TestUtil::~TestUtil()
{
}

void TestUtil::run()
{
	vector<Individual*> individuals = new PopulationReader()->loadFromArquive("logs/individuals.txt");
	vector<Individual*>nonDominatedSolutions = nonDominatedSol(individuals);
	HyperVolumeCalculator h;
	cout << h.calculateForTwoObjective(nonDominatedSolutions) << endl;
	clearVector(nonDominatedSolutions);
	clearVector(individuals);
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