#pragma once
#include<cstdlib>
#include<ctime>
#include<math.h>
#include<vector>
#include<stdlib.h>

#include"ObjectiveFunction.h"
#include"ZDT.h"
#include"Population.h"
#include"Individual.h"

using namespace std;

class PSO
{
public:
	PSO();
	~PSO();
	void run(Population& population);	
	double getLimiteSuperior();
	double getLimiteInferior();
	void setLimiteSuperior(double limiteSuperior);
	void setLimiteInferior(double limiteInferior);
	
private:
	double w = 0.9;
	const double c1 = 2;
	const double c2 = 2;
	double limiteSuperior;
	double limiteInferior;
	ZDT* objective;
	void initPopulation(Population& population);
	double manterValorNoLimite(double value);
	double manterValorNoLimiteEed(double value, double limiteInferior, double limiteSuperior);
	Individual* & getLeader(Population& population);
	Individual *& getMinFitnessObject(Population & population);
	Individual *& getLeaderMultiObjective(Population & population);
	bool isDominated(Individual * one, Individual * two);
	vector<double> evaluateIndividual(Individual * individual);
};

