#include "GA.h"


GA::GA()
{
	//ctor
	this->objective = new ZDT();
}

GA::~GA()
{
	//dtor
	delete this->objective;
}

void GA::run(Population* population) {

	for (int i = 0; i < numberOfIterations; i++) {
		this->crossover(population);
		this->mutating(population);		
	}
}

vector<double> GA::getLimiteInferior()
{
	return this->limiteInferior;
}

vector<double> GA::getlimiteSuperior()
{
	return this->limiteSuperior;
}

double GA::getTheta()
{
	return this->theta;
}

void GA::setLimiteInferior(vector<double> limiteInferior)
{
	this->limiteInferior = limiteInferior;
}

void GA::setLimiteSuperior(vector<double> limiteSuperior)
{
	this->limiteSuperior = limiteSuperior;
}

void GA::setTheta(double theta)
{
	this->theta = theta;
}

void GA::setFunction(ZDT * function)
{
	this->objective = function;
}

vector<double> GA::getGenes(int qtdGenes, double limiteInferior, double limiteSuperior) {

	vector<double> genes;
	int i = 0;

	for (i = 0; i < qtdGenes; i++) {
		double gene = (limiteInferior + (((double)rand() / (double)(RAND_MAX)) * limiteSuperior));
		genes.push_back(gene);
	}
	return genes;
}

vector<Individual*> GA::initPopulation(Population* population, double limiteInferior, double limiteSuperior, int qtdGenes) {

	int i = 0;
	int numIndividuals = population->getSize();
	srand((unsigned)time(NULL));

	for (i = 0; i < numIndividuals; i++) {
		Individual* id = new Individual();
		id->setQtdGenes(qtdGenes);
		id->setGenes(getGenes(id->getQtdGenes(), limiteInferior, limiteSuperior));
		id->setAptidao(evaluateIndividual(id));
		population->getIndividuals().push_back(id);
	}

	return population->getIndividuals();
}

vector<Individual*> GA::mutating(Population* population) {

	int numIndividuals = population->getSize();
	int i = 0;

	for (i = 0; i < numIndividuals; i++) {

		double tentative = ((double)rand() / (double)(RAND_MAX)) * 1;
		Individual* id = population->getIndividuals()[i];
		double gene = 0;

		if (tentative < this->getMutationRate()) {

			int genePosition = 0 + (rand() % abs(id->getQtdGenes()));
			gene = (limiteInferior[genePosition] + (((double)rand() / (double)(RAND_MAX)) * fabs(limiteInferior[genePosition] - limiteSuperior[genePosition])));
			id->getGenes()[genePosition] = gene;
			id->setAptidao(evaluateIndividual(id));
		}
	}

	return population->getIndividuals();
}

vector<Individual*> GA::crossover(Population* population) {

	size_t numIndividuals = population->getIndividuals().size();
	int i = 0;
	int times = (numIndividuals);
	double lambda = 1;

	double s = 0;
	double a = 0.5;
	/*for (i = 0; i < numIndividuals; i++) {

		s += population.getIndividuals().at(i)->getAptidao()[0] + population.getIndividuals().at(i)->getAptidao()[1];
	}
*/
	for (i = 0; i < times; i++) {

		double tentative = ((double)rand() / (double)(RAND_MAX)) * 1;

		if (tentative < this->getCrossOverRate()) {

			Individual* son = new Individual();
			Individual* son2 = new Individual();
			Individual* son3 = new Individual();
			/*int father = this->roulette(s, population);
			int mother = this->roulette(s, population);
			int moreDivesity = this->roulette(s, population);*/
			int father = this->tournament(population);
			int mother = this->tournament(population);
			int uncle = this->tournament(population);

			Individual* cromossomeFather = population->getIndividuals()[father];
			Individual* cromossomeMother = population->getIndividuals()[mother];
			Individual* cromossomeUncle = population->getIndividuals()[uncle];
			size_t qtdGenes = cromossomeFather->getGenes().size();
			
			son->setQtdGenes(qtdGenes);
			son2->setQtdGenes(qtdGenes);
			son3->setQtdGenes(qtdGenes);

			double alpha = random(-a, (1.0 + a));

			for (int j = 0; j < qtdGenes; j++) {

				double gene1 = (alpha * cromossomeFather->getGenes()[j]) + ((1.0 - alpha) * cromossomeMother->getGenes()[j]);
				double gene2 = (alpha * cromossomeMother->getGenes()[j]) + ((1.0 - alpha) *  cromossomeFather->getGenes()[j]);
				double gene3 = (alpha * cromossomeUncle->getGenes()[j]) + ((1.0 - alpha) * cromossomeMother->getGenes()[j]);

				son->getGenes().push_back(normalize(gene1, limiteInferior[j], limiteSuperior[j]));
				son2->getGenes().push_back(normalize(gene2, limiteInferior[j], limiteSuperior[j]));
				son3->getGenes().push_back(normalize(gene3, limiteInferior[j], limiteSuperior[j]));
			}

			son->setAptidao(evaluateIndividual(son));
			son2->setAptidao(evaluateIndividual(son2));
			son3->setAptidao(evaluateIndividual(son3));

			if (isDominated(cromossomeFather, son)) {
				population->getIndividuals()[father]->setAptidao(son->getAptidao());
				population->getIndividuals()[father]->setGenes(son->getGenes());
			}
			else if (isDominated(cromossomeMother, son)) {
				population->getIndividuals()[mother]->setAptidao(son->getAptidao());
				population->getIndividuals()[mother]->setGenes(son->getGenes());
			}
			else if (isDominated(cromossomeUncle, son)) {
				population->getIndividuals()[uncle]->setAptidao(son->getAptidao());
				population->getIndividuals()[uncle]->setGenes(son->getGenes());
			}

			if (isDominated(cromossomeMother, son2)) {
				population->getIndividuals()[mother]->setAptidao(son2->getAptidao());
				population->getIndividuals()[mother]->setGenes(son2->getGenes());
			}
			else if (isDominated(cromossomeFather, son2)) {
				population->getIndividuals()[father]->setAptidao(son2->getAptidao());
				population->getIndividuals()[father]->setGenes(son2->getGenes());
			}
			else if (isDominated(cromossomeUncle, son2)) {
				population->getIndividuals()[uncle]->setAptidao(son2->getAptidao());
				population->getIndividuals()[uncle]->setGenes(son2->getGenes());
			}

			if (isDominated(cromossomeFather, son3)) {
				population->getIndividuals()[father]->setAptidao(son3->getAptidao());
				population->getIndividuals()[father]->setGenes(son3->getGenes());
			}
			else if (isDominated(cromossomeMother, son3)) {
				population->getIndividuals()[mother]->setAptidao(son3->getAptidao());
				population->getIndividuals()[mother]->setGenes(son3->getGenes());
			}
			else if (isDominated(cromossomeUncle, son3)) {
				population->getIndividuals()[uncle]->setAptidao(son3->getAptidao());
				population->getIndividuals()[uncle]->setGenes(son3->getGenes());
			}

			delete son;
			delete son2;
			delete son3;
		}
	}

	return population->getIndividuals();
}

vector<Individual*> GA::minimizerCrossover(Population* population) {

	size_t numIndividuals = population->getIndividuals().size();
	int i = 0;
	double gamma = 1;
	double lambda = 1;

	double s = 0;
	double a = 0.1;

	for (i = 0; i < numIndividuals; i++) {

		s += population->getIndividuals().at(i)->getAptidao()[0] + population->getIndividuals().at(i)->getAptidao()[1];
	}

	for (i = 0; i < numIndividuals; i++) {

		double tentative = ((double)rand() / (double)(RAND_MAX)) * 1;

		if (tentative < this->getCrossOverRate()) {

			Individual* son = new Individual();
			Individual* son2 = new Individual();
			Individual* son3 = new Individual();
			int father = this->roulette(s, population);
			int mother = this->roulette(s, population);
			int moreDivesity = this->roulette(s, population);

			Individual* cromossomeFather = population->getIndividuals()[father];
			Individual* cromossomeMother = population->getIndividuals()[mother];
			Individual* cromossomeUncle = population->getIndividuals()[moreDivesity];

			size_t qtdGenes = cromossomeFather->getGenes().size();
			son->setQtdGenes(qtdGenes);
			son2->setQtdGenes(qtdGenes);
			son3->setQtdGenes(qtdGenes);

			double alpha = random(-a, ((double) 1.0 + a));

			for (int j = 0; j < qtdGenes; j++) {

				double gene1 = (alpha * cromossomeFather->getGenes()[j]) + ((1.0 - alpha) * cromossomeMother->getGenes()[j]);
				double gene2 = (alpha * cromossomeMother->getGenes()[j]) + ((1.0 - alpha) *  cromossomeFather->getGenes()[j]);
				double gene3 = (alpha * cromossomeFather->getGenes()[j]) + ((1.0 - alpha) * cromossomeUncle->getGenes()[j]);

				son->getGenes().push_back(normalize(gene1, limiteInferior[j], limiteSuperior[j]));
				son2->getGenes().push_back(normalize(gene2, limiteInferior[j], limiteSuperior[j]));
				son3->getGenes().push_back(normalize(gene3, limiteInferior[j], limiteSuperior[j]));

			}

			son->setAptidao(evaluateIndividual(son));
			son2->setAptidao(evaluateIndividual(son2));
			son3->setAptidao(evaluateIndividual(son3));

			if (son->getAptidao() < population->getIndividuals()[father]->getAptidao()) {

				population->getIndividuals()[father]->setAptidao(son->getAptidao());
				population->getIndividuals()[father]->setGenes(son->getGenes());
			}
			else if (son->getAptidao() < population->getIndividuals()[mother]->getAptidao()) {

				population->getIndividuals()[mother]->setAptidao(son->getAptidao());
				population->getIndividuals()[mother]->setGenes(son->getGenes());
			}

			if (son2->getAptidao() < population->getIndividuals()[mother]->getAptidao()) {

				population->getIndividuals()[mother]->setAptidao(son2->getAptidao());
				population->getIndividuals()[mother]->setGenes(son2->getGenes());
			}
			else if (son2->getAptidao() < population->getIndividuals()[father]->getAptidao()) {

				population->getIndividuals()[father]->setAptidao(son2->getAptidao());
				population->getIndividuals()[father]->setGenes(son2->getGenes());
			}

			if (son3->getAptidao() < population->getIndividuals()[father]->getAptidao()) {

				population->getIndividuals()[father]->setAptidao(son3->getAptidao());
				population->getIndividuals()[father]->setGenes(son3->getGenes());
			}
			else if (son3->getAptidao() < population->getIndividuals()[mother]->getAptidao()) {

				population->getIndividuals()[mother]->setAptidao(son3->getAptidao());
				population->getIndividuals()[mother]->setGenes(son3->getGenes());
			}

			delete son;
			delete son2;
			delete son3;
		}
	}

	return population->getIndividuals();
}

int GA::roulette(double s, Population* population) {

	size_t numIndividuals = population->getIndividuals().size();
	double soma = 0;
	double r = ((double)rand() / (double)(RAND_MAX)) * s;

	for (int i = 0; i < numIndividuals; i++) {

		soma += population->getIndividuals().at(i)->getAptidao()[0] + population->getIndividuals().at(i)->getAptidao()[1];
		if (soma > r)
			return i;
	}
	return (0 + (rand() % abs(population->getIndividuals().at(0)->getQtdGenes())));
}

int GA::tournament(Population* population) {

	int size = population->getIndividuals().size();
	int choosed = 0, firstWinner = 0, secondWinner = 0, thirdWinner = 0, fourthWinner = 0, finalist_one, finalist_two;
	int pos1 = 0 + (rand() % abs(size));
	int pos2 = 0 + (rand() % abs(size));
	int pos3 = 0 + (rand() % abs(size));
	int pos4 = 0 + (rand() % abs(size));
	int pos5 = 0 + (rand() % abs(size));
	int pos6 = 0 + (rand() % abs(size));
	int pos7 = 0 + (rand() % abs(size));
	int pos8 = 0 + (rand() % abs(size));

	if (isDominated(population->getIndividuals()[pos1], population->getIndividuals()[pos2])) {
		firstWinner = pos2;
	}
	else {
		firstWinner = pos1;
	}

	if (isDominated(population->getIndividuals()[pos3], population->getIndividuals()[pos4])) {
		secondWinner = pos4;
	}
	else {
		secondWinner = pos3;
	}

	if (isDominated(population->getIndividuals()[firstWinner], population->getIndividuals()[secondWinner])) {
		finalist_one = secondWinner;
	}
	else {
		finalist_one = firstWinner;
	}

	//second group

	if (isDominated(population->getIndividuals()[pos5], population->getIndividuals()[pos6])) {
		thirdWinner = pos6;
	}
	else {
		thirdWinner = pos5;
	}

	if (isDominated(population->getIndividuals()[pos7], population->getIndividuals()[pos8])) {
		fourthWinner = pos8;
	}
	else {
		fourthWinner = pos7;
	}

	if (isDominated(population->getIndividuals()[thirdWinner], population->getIndividuals()[fourthWinner])) {
		finalist_two = fourthWinner;
	}
	else {
		finalist_two = thirdWinner;
	}

	if (isDominated(population->getIndividuals()[finalist_two], population->getIndividuals()[finalist_one])) {
		choosed = finalist_one;
	}
	else {
		choosed = finalist_two;
	}
		
	return choosed;
}

vector<double> GA::evaluateIndividual(Individual* individual) {

	ZDT* objective = new ZDT();
	vector<double> value = objective->evaluateIndividual(individual);
	delete objective;

	return value;
}

double GA::min(double x1, double x2) {

	return x1 < x2 ? x1 : x2;
}

double GA::max(double x1, double x2) {

	return x1 > x2 ? x1 : x2;
}

double GA::random(double min, double max) {

	return (min + (((double)rand() / (double)(RAND_MAX)) * max));
}

double GA::normalize(double value, double lowerBound, double upperBound) {

	if (value < lowerBound)
		return lowerBound;
	else if (value > upperBound)
		return upperBound;

	return value;
}

double GA::scalar(double x, vector<double> v) {

	double sum = 0;

	for (int i = 0; i < v.size(); i++) {
		sum += (x * v[i]);
	}

	return sum;
}

bool GA::isDominated(Individual * one, Individual * two)
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

void GA::nonDominated(vector<Individual*> inds, Population* pop) {

	unsigned j = 0;
	vector<bool>flag(inds.size());

	for (j = 0; j < pop->getIndividuals().size(); j++) {

		for (unsigned int i = 0; i < inds.size(); i++) {
			if (this->isDominated(inds[i], pop->getIndividuals()[j])) {
				flag[i] = true;
				continue;
			}
		}
	}

	for (unsigned int i = 0; i < inds.size(); i++) {

		if (!inds.at(i)) {
			int pos = getWeakPos(pop->getIndividuals(), inds[i]);
			pop->getIndividuals()[pos]->setAptidao(inds[i]->getAptidao());
			pop->getIndividuals()[pos]->setGenes(inds[i]->getGenes());
		}
	}

}

int GA::getWeakPos(vector<Individual*>& individuals, Individual* ind) {

	for (unsigned int i = 0; i < individuals.size(); i++) {
		if (isDominated(individuals[i], ind)) {
			return i;
		}
	}

	return  0 + (rand() % individuals.size());
}

void GA::setMutationRate(double mutationRate) {
	this->mutationRate = mutationRate;
}

double GA::getMutationRate() {
	return this->mutationRate;
}

void GA::setCrossOverRate(double crossoverRate) {
	this->crossoverRate = crossoverRate;
}

double GA::getCrossOverRate() {
	return this->crossoverRate;
}


