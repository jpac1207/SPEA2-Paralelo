#include "Population.h"

Population::Population()
{

}

Population::~Population()
{
	for (size_t i = 0; i < this->individuals.size(); i++) {
		delete individuals[i];
	}

	individuals.clear();
}


bool compare(Individual* i, const Individual* j) {

	double f1 = i->getFitness();
	double f2 = j->getFitness();

	return (f1 < f2);
}


bool compareMax(Individual * i, Individual* j) {

	return (i->getFitness() > j->getFitness());
}

Individual* Population::getBestIndividual(bool minimization) {
	//verify

	if (minimization) {

		sort(this->getIndividuals().begin(), this->getIndividuals().end(), compare);

	}
	else {
		sort(this->getIndividuals().begin(), this->getIndividuals().end(), compareMax);
	}

	return this->getIndividuals().at(0);
}

const vector<Individual*> Population::getIndividualsCopy() {
	return (individuals);
}

vector<Individual*>& Population::getIndividuals() {
	return this->individuals;
}

void Population::setIndividuals(vector<Individual*> individuals) {
	this->individuals = individuals;
}

int Population::getSize() {
	return this->size;
}

void Population::setSize(int size) {
	this->size = size;
}

vector<double> Population::evaluateIndividual(Individual* individual) {

	ZDT* zdt1 = new ZDT();
	vector<double> value = zdt1->evaluateIndividual(individual);
	delete zdt1;

	return value;
}

vector<double> Population::getGenes(int qtdGenes, double limiteInferior, double limiteSuperior) {

	vector<double> genes;
	int i = 0;
	/*genes.push_back((0 + (((double)rand() / (double)(RAND_MAX)) * (1))));*/
	for (i = 0; i < qtdGenes; i++) {
		double gene = (limiteInferior + (((double)rand() / (double)(RAND_MAX)) * fabs(limiteInferior - limiteSuperior)));
		genes.push_back(gene);
	}
	return genes;
}

vector<double> Population::getGenesEed(int qtdGenes) {
	double inferiorGenes[] = { 5, 5, 5, 5, 5 , 5 };
	double superiorGenes[] = { 50, 60, 100 , 120, 100, 60 };
	vector<double> genes;
	int i = 0;

	for (i = 0; i < qtdGenes; i++) {

		double gene = (inferiorGenes[i] + (((double)rand() / (double)(RAND_MAX)) * (superiorGenes[i] - inferiorGenes[i])));
		genes.push_back(gene);
	}
	return genes;
}

vector<double> Population::getVelocity(int qtdGenes)
{
	vector<double> velocity;
	for (int i = 0; i < qtdGenes; i++)velocity.push_back(0.5);
	return velocity;
}

void Population::initPopulation(double limiteInferior, double limiteSuperior, int qtdGenes) {

	int i = 0;
	int numIndividuals = this->getSize();
	srand((unsigned)time(NULL));

	for (i = 0; i < numIndividuals; i++) {
		Individual* id = new Individual();
		id->setQtdGenes(qtdGenes);
		id->setGenes(getGenes(id->getQtdGenes(), limiteInferior, limiteSuperior));
		id->setVelocity(getVelocity(id->getQtdGenes()));
		id->setAptidao(evaluateIndividual(id));
		this->getIndividuals().push_back(id);
	}

}

/*Environmental Economic dispatch*/
void Population::initPopulationEed(int qtdGenes) {

	int i = 0;
	int numIndividuals = this->getSize();
	srand((unsigned)time(NULL));

	for (i = 0; i < numIndividuals; i++) {
		Individual* id = new Individual();
		id->setQtdGenes(qtdGenes);
		id->setGenes(getGenesEed(id->getQtdGenes()));
		id->setVelocity(getVelocity(id->getQtdGenes()));
		id->setAptidao(evaluateIndividual(id));
		this->getIndividuals().push_back(id);
	}

}

void Population::dump() {

	cout << "Population: " << "\n";
	for (int i = 0; i < this->individuals.size(); i++) {

		cout << this->getIndividuals().at(i)->getAptidao()[0] << ";" << this->getIndividuals().at(i)->getAptidao()[1] << "\n";
	}
}

