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

void DE::run(Population& population)
{

	/*for (int generations = 0; generations < 60; generations++) {*/
	vector<Individual*> donators;
	this->mutar(donators, population);
	this->cruzar(donators, population);
	this->clearVector(donators);
	/*}*/
}

void DE::setLimiteInferior(double limiteInferior)
{
	this->limiteInferior = limiteInferior;
}

void DE::setLimiteSuperior(double limiteSuperior)
{
	this->limiteSuperior = limiteSuperior;
}

void DE::setCrossoverRate(double crossoverRate)
{
	this->crossoverRate = crossoverRate;
}

double DE::getLimiteInferior()
{
	return this->limiteInferior;
}

double DE::getLimiteSuperior()
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

void DE::mutar(vector<Individual*>& donators, Population & population)
{
	size_t size = population.getIndividuals().size();
	size_t solutionSize = population.getIndividuals()[0]->getGenes().size();
	vector<Individual*> lvDonators;
	double inferiorGenes[] = { 5, 5, 5, 5, 5 , 5 };
	double superiorGenes[] = { 50, 60, 100 , 120, 100, 60 };

	for (unsigned int i = 0; i < size; i++) {

		int posOne = i;
		Individual* one = population.getIndividuals()[posOne];
		Individual* donator = new Individual();
		donator->setAptidao(one->getAptidao());
		donator->setGenes(one->getGenes());
		donator->setTarget(one);

		int posTwo = 0;
		Individual* first;

		do {
			posTwo = getRandomIndividual(population);
		} while (posOne == posTwo);

		first = population.getIndividuals()[posTwo];
		int posThree = 0;
		Individual* second;

		do {
			posThree = getRandomIndividual(population);
		} while (posOne == posThree || posTwo == posThree);

		second = population.getIndividuals()[posThree];
		int posFour = 0;
		Individual* third;

		do {
			posFour = getRandomIndividual(population);
		} while ((posOne == posFour) || (posTwo == posFour) || (posThree == posFour));

		third = population.getIndividuals()[posFour];
		double f = 1.0;

		for (int j = 0; j < solutionSize; j++) {

			double newGene = first->getGenes()[j] + (f * (second->getGenes()[j] - third->getGenes()[j]));
			
			if (newGene > superiorGenes[j])
				newGene = superiorGenes[j];
			else if (newGene < inferiorGenes[j])
				newGene = inferiorGenes[j];

			donator->getGenes()[j] = newGene;
		}

		lvDonators.push_back(donator);
	}

	donators.swap(lvDonators);
}

void DE::cruzar(vector<Individual*> donators, Population& population)
{
	size_t size = population.getIndividuals().size();
	size_t solutionSize = population.getIndividuals()[0]->getGenes().size();

	for (int i = 0; i < size; i++) {

		Individual* donator = donators[i];
		Individual* target = donator->getTarget();
		Individual* experimental = new Individual();
		experimental->setQtdGenes(donator->getQtdGenes());

		for (int j = 0; j < solutionSize; j++) {

			double tentative = ((double)rand() / (double)(RAND_MAX)) * 1;

			if (tentative > crossoverRate)
				experimental->getGenes().push_back(target->getGenes()[j]);
			else
				experimental->getGenes().push_back(donator->getGenes()[j]);
		}

		experimental->setAptidao(this->evaluateIndividual(experimental));

		if (isDominate(population.getIndividuals()[i], experimental)) {
			population.getIndividuals()[i]->setAptidao(experimental->getAptidao());
			population.getIndividuals()[i]->setGenes(experimental->getGenes());
		}
		delete experimental;
	}
}

bool DE::isDominate(Individual* one, Individual*  two)
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

int DE::getRandomIndividual(Population & population)
{
	return 0 + (rand() % abs((int)population.getIndividuals().size()));
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
}
