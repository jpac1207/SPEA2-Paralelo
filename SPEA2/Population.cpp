#include "Population.h"

Population::Population()
{

}

Population::Population(const Population& population) {
	this->size = population.size;

	for (size_t i = 0; i < this->size; i++)
		this->individuals.push_back(new Individual(*population.individuals.at(i)));
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

	ZDT* zdt = new ZDT();
	vector<double> value = zdt->evaluateIndividual(individual);
	delete zdt;

	return value;
}

vector<double> Population::getGenes(int qtdGenes, vector<double> limiteInferior, vector<double> limiteSuperior) {

	vector<double> genes;
	int i = 0;

	for (i = 0; i < qtdGenes; i++) {
		double gene = (limiteInferior[i] + (((double)rand() / (double)(RAND_MAX)) * fabs(limiteInferior[i] - limiteSuperior[i])));
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

void Population::initPopulation(vector<double> limiteInferior, vector<double> limiteSuperior, int qtdGenes) {

	/*int i = 0;
	srand((unsigned)time(NULL));

	for (i = 0; i < this->size; i++) {
		Individual* id = new Individual();
		id->setQtdGenes(qtdGenes);
		id->setGenes(getGenes(qtdGenes, limiteInferior, limiteSuperior));
		id->setVelocity(getVelocity(qtdGenes));
		id->setAptidao(evaluateIndividual(id));
		this->individuals.push_back(id);
	}*/

	int i = 0;
	srand((unsigned)time(NULL));
	vector<double> genes;
	genes.push_back(113.997453);
	genes.push_back(113.626347);
	genes.push_back(97.399937);
	genes.push_back(179.733101);
	genes.push_back(90.494299);
	genes.push_back(105.400153);
	genes.push_back(259.599877);
	genes.push_back(299.9);
	genes.push_back(284.601078);
	genes.push_back(204.799816);
	genes.push_back(94.1);
	genes.push_back(94.2);
	genes.push_back(214.759791);
	genes.push_back(394.279373);
	genes.push_back(394.279398);
	genes.push_back(394.279381);
	genes.push_back(489.279397);
	genes.push_back(489.27939);
	genes.push_back(511.279371);
	genes.push_back(511.279371);
	genes.push_back(523.27939);
	genes.push_back(523.279437);
	genes.push_back(523.279474);
	genes.push_back(523.279398);
	genes.push_back(523.279375);
	genes.push_back(523.27937);
	genes.push_back(10.1);
	genes.push_back(10.1);
	genes.push_back(10.1);
	genes.push_back(88.297938);
	genes.push_back(189.9);
	genes.push_back(189.9);
	genes.push_back(189.9);
	genes.push_back(164.88839);
	genes.push_back(164.812509);
	genes.push_back(199.9);
	genes.push_back(91.371556);
	genes.push_back(93.306261);
	genes.push_back(109.9);
	genes.push_back(511.279371);
	
	Individual* custom = new Individual();
	custom->setQtdGenes(qtdGenes);
	custom->setGenes(genes);
	custom->setVelocity(getVelocity(qtdGenes));
	custom->setAptidao(evaluateIndividual(custom));
	this->individuals.push_back(custom);

	for (i = 0; i < this->size - 1; i++) {
		Individual* id = new Individual();
		id->setQtdGenes(qtdGenes);
		id->setGenes(getGenes(qtdGenes, limiteInferior, limiteSuperior));
		id->setVelocity(getVelocity(qtdGenes));
		id->setAptidao(evaluateIndividual(id));
		this->individuals.push_back(id);
	}
}

void Population::dump() {

	cout << "Population: " << "\n";
	for (int i = 0; i < this->individuals.size(); i++) {

		cout << this->getIndividuals().at(i)->getAptidao()[0] << ";" << this->getIndividuals().at(i)->getAptidao()[1] << "\n";
	}
}

