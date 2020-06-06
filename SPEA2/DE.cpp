//#include "stdafx.h"
#include "DE.h"


DE::DE()
{
	this->objective = new ZDT();
}


DE::~DE()
{
	delete this->objective;
}

void DE::run(Population* population)
{
	for (int generations = 0; generations < numberOfIterations; generations++) {
		vector<Individual*> donators;
		this->mutar(donators, population);
		this->cruzar(donators, population);
		this->clearVector(donators);
	}
}

void DE::setLimiteInferior(vector<double> limiteInferior)
{
	this->limiteInferior = limiteInferior;
}

void DE::setLimiteSuperior(vector<double> limiteSuperior)
{
	this->limiteSuperior = limiteSuperior;
}

void DE::setCrossoverRate(double crossoverRate)
{
	this->crossoverRate = crossoverRate;
}

vector<double> DE::getLimiteInferior()
{
	return this->limiteInferior;
}

vector<double> DE::getLimiteSuperior()
{
	return this->limiteSuperior;
}

double DE::getCrossoverRate()
{
	return this->crossoverRate;
}

void DE::setObjective(ZDT * obj)
{
	this->objective = obj;
}

void DE::mutar(vector<Individual*>& donators, Population* population)
{
	size_t size = population->getIndividuals().size();
	size_t solutionSize = population->getIndividuals()[0]->getGenes().size();
	vector<Individual*> lvDonators;

	for (unsigned int i = 0; i < size; i++) {

		int posOne = i;
		Individual* one = population->getIndividuals()[posOne];
		Individual* donator = new Individual();
		Individual* target = new Individual(*one);

		donator->setAptidao(one->getAptidao());
		donator->setGenes(one->getGenes());
		donator->setTarget(target);

		int posTwo = 0;
		Individual* first;

		do {
			posTwo = getRandomIndividual(population);
		} while (posOne == posTwo);

		first = population->getIndividuals()[posTwo];
		int posThree = 0;
		Individual* second;

		do {
			posThree = getRandomIndividual(population);
		} while (posOne == posThree || posTwo == posThree);

		second = population->getIndividuals()[posThree];
		int posFour = 0;
		Individual* third;

		do {
			posFour = getRandomIndividual(population);
		} while ((posOne == posFour) || (posTwo == posFour) || (posThree == posFour));

		third = population->getIndividuals()[posFour];
		double f = 1.0;

		for (int j = 0; j < solutionSize; j++) {

			double newGene = first->getGenes()[j] + (f * (second->getGenes()[j] - third->getGenes()[j]));

			if (newGene > limiteSuperior[j])
				newGene = limiteSuperior[j];
			else if (newGene < limiteInferior[j])
				newGene = limiteInferior[j];

			donator->getGenes()[j] = newGene;
		}

		lvDonators.push_back(donator);
	}

	donators.swap(lvDonators);
}

void DE::cruzar(vector<Individual*> donators, Population* population)
{
	size_t size = population->getIndividuals().size();
	size_t solutionSize = population->getIndividuals()[0]->getGenes().size();

	for (int i = 0; i < size; i++) {

		Individual* donator = donators[i];
		Individual* target = donator->getTarget();
		Individual* candidate = new Individual();
		candidate->setQtdGenes(donator->getQtdGenes());

		for (int j = 0; j < solutionSize; j++) {

			double tentative = ((double)rand() / (double)(RAND_MAX)) * 1;

			if (tentative > crossoverRate)
				candidate->getGenes().push_back(target->getGenes()[j]);
			else
				candidate->getGenes().push_back(donator->getGenes()[j]);
		}

		candidate->setAptidao(this->evaluateIndividual(candidate));

		if (isDominated(target, candidate)) {
			int pos = i;
			population->getIndividuals()[pos]->setAptidao(candidate->getAptidao());
			population->getIndividuals()[pos]->setGenes(candidate->getGenes());
		}
		delete candidate;
		delete target;
	}
}

bool DE::isDominated(Individual* one, Individual*  two)
{
	bool anyDominate = false;

	for (unsigned int i = 0; i < one->getAptidao().size(); i++) {

		if (this->objective->isMinimization(i)) {

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

int DE::getRandomIndividual(Population* population)
{
	return 0 + (rand() % abs((int)population->getIndividuals().size()));
}

bool DE::nonDominated(Individual* ind, Population* pop) {

	unsigned j = 0;
	for (j = 0; j < pop->getIndividuals().size(); j++) {

		if (this->isDominated(ind, pop->getIndividuals()[j])) {
			return false;
		}
	}

	return true;
}

int DE::getWeakPos(vector<Individual*>& individuals, Individual* ind) {

	for (unsigned int i = 0; i < individuals.size(); i++) {
		if (isDominated(individuals[i], ind)) {
			return i;
		}
	}

	return  0 + (rand() % individuals.size());
}


vector<double> DE::evaluateIndividual(Individual * individual)
{
	ZDT* objective = new ZDT();
	vector<double> value = objective->evaluateIndividual(individual);
	delete objective;

	return value;
}

void DE::clearVector(vector<Individual*>& v)
{
	for (size_t i = 0; i < v.size(); i++)
		delete v[i];
	v.clear();
}
