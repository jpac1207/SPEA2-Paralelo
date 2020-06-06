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
	void run(Population* population);	
	vector<double> getLimiteSuperior();
	vector<double> getLimiteInferior();
	void setLimiteSuperior(vector<double> limiteSuperior);
	void setLimiteInferior(vector<double> limiteInferior);
	
private:
	double w = 0.9;
	const double c1 = 0.5;
	const double c2 = 0.5;
	const int numberOfIterations = 30;
	vector<double> limiteSuperior;
	vector<double> limiteInferior;
	ZDT* objective;
	void initPopulation(Population* population);
	double keepValuesInBounds(double value, double limiteInferior, double limiteSuperior);
	Individual *& getLeader(Population* population);
	Individual *& getMinFitnessObject(Population* population);
	Individual *& getLeaderMultiObjective(Population* population);
	bool isDominated(Individual * one, Individual* two);
	bool nonDominated(Individual * ind, Population * pop);
	vector<int> identifyExtremalPositions(vector<Individual*> individuals);
	int getWeakPos(vector<Individual*>& individuals, Individual *ind);
	vector<double> evaluateIndividual(Individual* individual);
};

