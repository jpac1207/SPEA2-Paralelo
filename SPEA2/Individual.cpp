#include "Individual.h"

Individual::Individual()
{
	//ctor
}

Individual::Individual(int qtdGenes)
{
	for (int i = 0; i < qtdGenes; i++) this->velocity.push_back(1);
}

Individual::Individual(const Individual & ind)
{
	for (int i = 0; i < ind.genes.size(); i++) {
		this->velocity.push_back(1);
		this->genes.push_back(ind.genes[i]);
	}
	for (int i = 0; i < ind.historicalGenes.size(); i++) {
		this->historicalGenes.push_back(ind.historicalGenes[i]);
	}
	for (int i = 0; i < ind.aptidao.size(); i++) {
		this->aptidao.push_back(ind.aptidao[i]);
	}
	for (int i = 0; i < ind.historicalAptidao.size(); i++) {
		this->historicalAptidao.push_back(ind.historicalAptidao[i]);
	}

	this->fitness = ind.fitness;
	this->d = ind.d;
	this->raw = ind.raw;
	this->s = ind.s;
	this->qtdGenes = ind.qtdGenes;
}

Individual::~Individual()
{

}

vector<double>& Individual::getGenes() {
	return this->genes;
}

void Individual::setGenes(vector<double> genes) {
	this->genes = genes;
}

double Individual::getFitness() const {
	return this->fitness;
}

void Individual::setFitness(double fitness) {
	this->fitness = fitness;
}

vector<double> Individual::getAptidao() const {
	return this->aptidao;
}

void Individual::setAptidao(vector<double> aptidao) {
	this->aptidao = aptidao;
}

void Individual::setHistoricalGenes(vector<double> historicalGenes)
{
	this->historicalGenes = historicalGenes;
}

void Individual::setHistoricalAptidao(vector<double> historicalAptidao)
{
	this->historicalAptidao = historicalAptidao;
}

void Individual::setVelocity(vector<double> velocity)
{
	this->velocity = velocity;
}

int Individual::getQtdGenes() {
	return this->qtdGenes;
}

void Individual::setQtdGenes(int qtd) {
	this->qtdGenes = qtd;
}

double Individual::getRaw() {
	return this->raw;
}
double Individual::getS() {
	return this->s;
}
double Individual::getD() {
	return this->d;
}
Individual* Individual::getTarget()
{
	return this->target;
}
vector<double> Individual::getHistoricalGenes()
{
	return this->historicalGenes;
}
vector<double> Individual::getHistoricalAptidao()
{
	return this->historicalAptidao;
}
vector<double>& Individual::getVelocity()
{
	return this->velocity;
}
void Individual::setRaw(double raw) {
	this->raw = raw;
}
void Individual::setS(double s) {
	this->s = s;
}
void Individual::setD(double d) {
	this->d = d;
}

void Individual::setTarget(Individual* id) {
	this->target = id;
}
