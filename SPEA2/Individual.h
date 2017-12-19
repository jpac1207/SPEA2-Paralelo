#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <iostream>
#include<vector>

using namespace std;

class Individual
{
public:
	Individual();
	Individual(int qtdGenes);
	Individual(const Individual& ind);
	virtual ~Individual();
	vector<double>& getGenes();
	void setGenes(vector<double>);
	void setFitness(double fitness);
	void setRaw(double raw);
	void setS(double s);
	void setD(double d);
	void setTarget(Individual* id);
	void setAptidao(vector<double> aptidao);
	void setHistoricalGenes(vector<double> historicalGenes);
	void setHistoricalAptidao(vector<double> historicalAptidao);
	void setVelocity(vector<double> velocity);
	int getQtdGenes();
	void setQtdGenes(int qtd);
	double getRaw();
	double getS();
	double getD();
	double getFitness() const;
	vector<double> getAptidao() const;
	Individual* getTarget();
	vector<double> getHistoricalGenes();
	vector<double> getHistoricalAptidao();
	vector<double>& getVelocity();

protected:

private:
	vector<double> genes;
	vector<double> aptidao;
	double fitness;
	int qtdGenes;
	double raw;
	double s;
	double d;
	Individual* target;
	vector<double> historicalGenes;
	vector<double> historicalAptidao;
	vector<double> velocity;

};

#endif // INDIVIDUAL_H
