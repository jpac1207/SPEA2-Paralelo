#ifndef POPULATION_H
#define POPULATION_H


#include<algorithm>
#include<vector>
#include<cstdlib>
#include<ctime>
#include<iostream>
#include<stdio.h>

#include"Individual.h"
#include"ZDT.h"

using namespace std;

class Population
{
    public:
        Population();		
        virtual ~Population();
        vector<Individual*>& getIndividuals();
        const vector<Individual*> getIndividualsCopy();
        void setIndividuals(vector<Individual*> individuals);
        int getSize();
        void setSize(int size);
        Individual* getBestIndividual(bool minimization);
        void initPopulation(double limiteInferior, double limiteSuperior, int qtdGenes);
		void initPopulationEed(int qtdGenes);
		void dump();
    protected:

    private:
        vector<Individual*> individuals;
        int size;
        vector<double> getGenes(int qtdGenes, double limiteInferior, double limiteSuperior);
		vector<double> getGenesEed(int qtdGenes);		
		vector<double> getVelocity(int qtdGenes);
        vector<double> evaluateIndividual(Individual* individual);

};

#endif // POPULATION_H
