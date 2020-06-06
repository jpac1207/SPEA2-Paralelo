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
	void run(Population* population);
	void setLimiteInferior(vector<double> limiteInferior);
	void setLimiteSuperior(vector<double> limiteSuperior);
	void setCrossoverRate(double crossoverRate);
	vector<double> getLimiteInferior();
	vector<double> getLimiteSuperior();
	double getCrossoverRate();
	void setObjective(ZDT* obj);
private:
	const int numberOfIterations = 30;
	vector<double> limiteInferior;
	vector<double> limiteSuperior;
	double crossoverRate;
	ZDT* objective;
	void mutar(vector<Individual*>& donators, Population* population);
	void cruzar(vector<Individual*> donators, Population* population);
	bool isDominated(Individual* one, Individual* two);
	int getRandomIndividual(Population* population);
	bool nonDominated(Individual * ind, Population * pop);
	int getWeakPos(vector<Individual*>& individuals, Individual * ind);
	vector<double> evaluateIndividual(Individual* individual);
	void clearVector(vector<Individual*> &v);
};

