#include "HyperVolumeCalculator.h"



HyperVolumeCalculator::HyperVolumeCalculator()
{
}


HyperVolumeCalculator::~HyperVolumeCalculator()
{
}


/*Sort by first function*/
bool compareByFirstObjective(Individual* i, const Individual* j) {

	double f1 = i->getAptidao()[0];
	double f2 = j->getAptidao()[0];

	return (f1 > f2);
}

bool compareBySecondObjective(Individual* i, const Individual* j) {

	double f1 = i->getAptidao()[1];
	double f2 = j->getAptidao()[1];

	return (f1 > f2);
}

double HyperVolumeCalculator::calculateForTwoObjective(vector<Individual*> individuals)
{
	double partialHeight = 0.0;
	this->hypervolume = 0;
	sort(individuals.begin(), individuals.end(), compareByFirstObjective);
	vector<double> extremalPoints = identifyExtremalPoints(individuals);
	this->referenceOne = extremalPoints.at(0);
	this->referenceTwo = extremalPoints.at(1);
	
	for (size_t i = 0; i < individuals.size(); i++) {

		double height = fabs(referenceOne - individuals[i]->getAptidao()[0]);
		double width = fabs(referenceTwo - individuals[i]->getAptidao()[1]);
		height -= partialHeight;
		partialHeight += height;
		double area = height * width;
		hypervolume += area;
	}

	return hypervolume;
}

vector<double> HyperVolumeCalculator::differenceBetweenExtremalPointsAndReference(vector<Individual*> individuals)
{
	size_t firstPosition = 0;
	size_t lastPosition = individuals.size() - 1;
	vector<double> differences;

	sort(individuals.begin(), individuals.end(), compareByFirstObjective);

	differences.push_back(abs((referenceOne - individuals.at(firstPosition)->getAptidao()[0]) +
		(referenceTwo - individuals.at(firstPosition)->getAptidao()[1])));

	differences.push_back(abs((referenceOne - individuals.at(lastPosition)->getAptidao()[0]) +
		(referenceTwo - individuals.at(lastPosition)->getAptidao()[1])));

	return differences;
}

vector<double> HyperVolumeCalculator::identifyExtremalPoints(vector<Individual*> individuals) {

	vector<double> extremalPoints;
	double referenceOne = individuals.at(0)->getAptidao()[0];
	double referenceTwo = individuals.at(0)->getAptidao()[1];

	for (size_t j = 1; j < individuals.size(); j++) {

		if (individuals.at(j)->getAptidao()[0] > referenceOne)
			referenceOne = individuals.at(j)->getAptidao()[0];
		if (individuals.at(j)->getAptidao()[1] > referenceTwo)
			referenceTwo = individuals.at(j)->getAptidao()[1];
	}

	extremalPoints.push_back(referenceOne);
	extremalPoints.push_back(referenceTwo);

	return extremalPoints;
}
