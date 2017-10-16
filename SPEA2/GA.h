#ifndef GA_H
#define GA_H

#include<cstdlib>
#include<ctime>
#include<math.h>
#include<stdexcept>
#include "stdafx.h"
#include <iostream>
#include"ObjectiveFunction.h"
#include"ZDT.h"
#include"Population.h"
#include"Individual.h"
#include"HyperVolumeCalculator.h"


using namespace std;

class GA
{
    public:
        GA();
        virtual ~GA();
        vector<Individual*> initPopulation(Population* population, double limiteInferior, double limiteSuperior, int qtdGenes);
        double getMutationRate();
        double getCrossOverRate();
        void setMutationRate(double mutationRate);
        void setCrossOverRate(double crossoverRate);
        vector<double> evaluateIndividual(Individual* individual);		
        vector<Individual*> mutating(Population& population);
		vector<Individual*> mutatingEed(Population & population);
        vector<Individual*> crossover(Population& population);
		vector<Individual*> minimizerCrossover(Population & population);
        void run(Population& population);
		double getLimiteInferior();
		double getlimiteSuperior();
		double getTheta();
		void setLimiteInferior(double limiteInferior);
		void setLimiteSuperior(double limiteSuperior);
		void setTheta(double theta);
		void setFunction(ZDT* function);
		

    protected:

    private:
        double mutationRate;
        double crossoverRate;	
		double limiteInferior;
		double limiteSuperior;
		double theta;
        const int M = 10;
        vector<double> getGenes(int qtdGenes, double limiteInferior, double limiteSuperior);
        int roulette(double s, Population& p);
		int tournament(Population & population);
		double min(double x1, double x2);
		double max(double x1, double x2);
		double random(double min, double max);
		double normalize(double value);
		double normalize(double value, double lowerBound, double upperBound);
		double scalar(double x, vector<double> v);
		bool isDominated(Individual* one, Individual* two);
		ZDT* objective;
};

#endif // GA_H
