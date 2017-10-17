#include "PSO.h"



PSO::PSO()
{
	this->objective = new ZDT();
}


PSO::~PSO()
{
	delete this->objective;
}

void PSO::run(Population& population)
{

	double inferiorGenes[] = { 10, 10, 35, 35, 130 , 125 };
	double superiorGenes[] = { 125, 150, 225 , 210, 325, 315 };
	this->initPopulation(population);
	int swarmSize = population.getIndividuals().size();
	int times = (swarmSize);

	for (int i = 0; i < times; i++) {

		Individual* leader = this->getLeader(population);
		vector<double> p = leader->getGenes();
		size_t solutionSize = leader->getGenes().size();
		int position = 0 + (rand() % swarmSize);
		Individual* individual = population.getIndividuals()[position];
		double r1 = ((double)rand() / (double)(RAND_MAX)) * 1;
		double r2 = ((double)rand() / (double)(RAND_MAX)) * 1;

		w = 0.1 + (0.8 / i);
		Individual* trial = new Individual(*individual);
		for (int j = 0; j < solutionSize; j++) {

			double historicalGene = individual->getHistoricalGenes()[j];
			double v = individual->getVelocity()[j];
			double gene = individual->getGenes()[j];

			double equationA = (w * v) + ((c1 * r1) * (historicalGene - gene)) + ((c2 * r2) * (p[j] - gene));
			trial->getVelocity()[j] = (equationA);
			double equationB = (gene + equationA);

			if (equationB > superiorGenes[j] || equationB < inferiorGenes[j]) {

				/*equationB = manterValorNoLimite(equationB);*/
				equationB = manterValorNoLimiteEed(equationB, inferiorGenes[j], superiorGenes[j]);
				trial->getVelocity()[j] *= (-0.5);
			}

			trial->getGenes()[j] = equationB;
		}

		trial->setAptidao(this->evaluateIndividual(individual));

		if (this->isDominated(individual, trial)) {
			individual->setAptidao(trial->getAptidao());
			individual->setHistoricalAptidao(trial->getHistoricalAptidao());
			individual->setGenes(trial->getGenes());
			individual->setHistoricalGenes(trial->getHistoricalGenes());
			individual->setVelocity(trial->getVelocity());
		}
		delete trial;
	}

}

vector<double> PSO::evaluateIndividual(Individual* individual) {

	ZDT* objective = new ZDT();
	vector<double> aptidao = objective->evaluateIndividual(individual);
	individual->setHistoricalGenes(individual->getGenes());
	individual->setHistoricalAptidao(aptidao);
	delete objective;

	return aptidao;
}

void PSO::initPopulation(Population & population)
{
	for (int i = 0; i < population.getIndividuals().size(); i++) {

		Individual* individual = population.getIndividuals()[i];
		individual->setHistoricalGenes(individual->getGenes());
		individual->setVelocity(individual->getVelocity());
		individual->setHistoricalAptidao(individual->getAptidao());
	}
}

double PSO::manterValorNoLimite(double value)
{
	if (value > this->limiteSuperior)
		return limiteSuperior;
	else if (value < this->limiteInferior)
		return limiteInferior;

	return value;
}

double PSO::manterValorNoLimiteEed(double value, double limiteInferior, double limiteSuperior)
{
	if (value > limiteSuperior)
		return limiteSuperior;
	else if (value < limiteInferior)
		return limiteInferior;

	return value;
}

Individual* & PSO::getLeader(Population & population)
{
	int pos = 0 + (rand() % population.getIndividuals().size());
	return (population.getIndividuals()[pos]);
}

Individual* & PSO::getMinFitnessObject(Population & population)
{
	int pos = 0;
	vector<double> lower = population.getIndividuals()[pos]->getAptidao();

	for (size_t j = 1; j < population.getIndividuals().size(); j++) {

		if (population.getIndividuals()[j]->getAptidao() < lower) {
			pos = j;
			lower = population.getIndividuals()[j]->getAptidao();
		}
	}

	return (population.getIndividuals()[pos]);
}

Individual* & PSO::getLeaderMultiObjective(Population & population)
{
	int size = population.getIndividuals().size();
	int choosed = 0, firstWinner = 0, secondWinner = 0, thirdWinner = 0, fourthWinner = 0, finalist_one, finalist_two;
	int pos1 = 0 + (rand() % abs(size));
	int pos2 = 0 + (rand() % abs(size));
	int pos3 = 0 + (rand() % abs(size));
	int pos4 = 0 + (rand() % abs(size));
	int pos5 = 0 + (rand() % abs(size));
	int pos6 = 0 + (rand() % abs(size));
	int pos7 = 0 + (rand() % abs(size));
	int pos8 = 0 + (rand() % abs(size));

	if (isDominated(population.getIndividuals()[pos1], population.getIndividuals()[pos2])) {
		firstWinner = pos2;
	}
	else {
		firstWinner = pos1;
	}

	if (isDominated(population.getIndividuals()[pos3], population.getIndividuals()[pos4])) {
		secondWinner = pos4;
	}
	else {
		secondWinner = pos3;
	}

	if (isDominated(population.getIndividuals()[firstWinner], population.getIndividuals()[secondWinner])) {
		finalist_one = secondWinner;
	}
	else {
		finalist_one = firstWinner;
	}

	//second group

	if (isDominated(population.getIndividuals()[pos5], population.getIndividuals()[pos6])) {
		thirdWinner = pos6;
	}
	else {
		thirdWinner = pos5;
	}

	if (isDominated(population.getIndividuals()[pos7], population.getIndividuals()[pos8])) {
		fourthWinner = pos8;
	}
	else {
		fourthWinner = pos7;
	}

	if (isDominated(population.getIndividuals()[thirdWinner], population.getIndividuals()[fourthWinner])) {
		finalist_two = fourthWinner;
	}
	else {
		finalist_two = thirdWinner;
	}

	if (isDominated(population.getIndividuals()[finalist_two], population.getIndividuals()[finalist_one])) {
		choosed = finalist_one;
	}
	else {
		choosed = finalist_two;
	}

	return population.getIndividuals()[choosed];
}

bool PSO::isDominated(Individual * one, Individual * two)
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


double PSO::getLimiteSuperior()
{
	return this->limiteSuperior;
}

double PSO::getLimiteInferior()
{
	return this->limiteInferior;
}

void PSO::setLimiteSuperior(double limiteSuperior)
{
	this->limiteSuperior = limiteSuperior;
}

void PSO::setLimiteInferior(double limiteInferior)
{
	this->limiteInferior = limiteInferior;
}

