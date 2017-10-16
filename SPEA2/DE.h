#pragma once

#include<cstdlib>
#include<ctime>
#include<math.h>
#include<vector>
#include"ObjectiveFunction.h"
#include"ZDT.h"
#include"Population.h"
#include"Individual.h"
#include"HyperVolumeCalculator.h"

using namespace std;

class DE
{
public:
	DE();
	~DE();
	void run(Population& population);
	void setLimiteInferior(double limiteInferior);
	void setLimiteSuperior(double limiteSuperior);
	void setCrossoverRate(double crossoverRate);
	double getLimiteInferior();
	double getLimiteSuperior();
	double getCrossoverRate();
	void setObjective(ZDT* obj);
private:
	double limiteInferior;
	double limiteSuperior;
	double crossoverRate;
	ZDT* objective;
	void mutar(vector<Individual*>& donators, Population& population);
	void cruzar(vector<Individual*> donators, Population& population);
	bool isDominate(Individual* one, Individual* two);
	int getRandomIndividual(Population& population);
	vector<double> evaluateIndividual(Individual* individual);
	void clearVector(vector<Individual*> &v);
};

