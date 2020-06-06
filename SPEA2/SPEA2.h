#ifndef SPEA2_H
#define SPEA2_H

#define ARRAY_SIZE(array) (sizeof((array[0]))/sizeof((*array[0])))

#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<stdio.h>
#include<vector>
#include<algorithm>
#include<iostream>
//#include<cstdlib>
//#include<ctime>
#include"mpi.h"
#include <time.h>

#include"ObjectiveFunction.h"
#include"ZDT.h"
#include"Individual.h"
#include"Population.h"
#include"GA.h"
#include"DE.h"
#include"PSO.h"
#include"HyperVolumeCalculator.h"
#include"PopulationReader.h"
#include"MatTransform.h"
#include"FuzzyController.h"

using namespace std;

class SPEA2
{	
    public:
        SPEA2();
		SPEA2(int argc, char* argv[]);
		SPEA2(int rank);
        virtual ~SPEA2();
        void run();
		void runSingle();
		

    protected:
		static const bool debug = false;
		
    private:
        static const int numberOfIterations = 300;
        static const int archiveSize = 50;
        static const int populationSize = 1500;
        static const int qtdGenes = 40;		
		static const int qtdThreads = 4;
		static const int numberOfSlaves = 3;
		static const int qtdMetaHeuristcs = 3;		
		vector<double> limiteInferior;
		vector<double> limiteSuperior;
        const int k = (int) sqrt( populationSize );
		double lastHipervolume;
		int argc;
		char** argv;		
		int rank;
        ObjectiveFunction* objective;
        vector<Individual*> archive;
		vector<Individual*> RT;
		vector<Individual*> index;
		vector<double> probabilities;
        Population* population;
        vector<Individual*> nonDominatedSol(vector<Individual*> populacao);
        bool isDominated(Individual* one, Individual* two);
		vector<Individual*> fillRT(vector<Individual*> pop, vector<Individual*> archive);
        void compute_s(vector<Individual*>& RT);
		void compute_raw(vector<Individual*>& RT);
		void compute_d(vector<Individual*>& RT);
        void distances(vector< vector<double> >& matrix, vector<Individual*> individuos);
        void orderMatrix(vector< vector<double> >& matriz);
        void getNonDominatedSolutions(vector<Individual*>& index, vector<Individual*> RT);
        void fillArchive(vector<Individual*>& RT, vector<Individual*>& archive, vector<Individual*>& index);
        void clustering(vector<Individual*>& archive, vector<Individual*>& individuos);
		void genericalClustering(vector<Individual*>& individuos, int nToRemove);
        void getLineFromAMatrix(vector<double>& lineMatrix, vector< vector<double> >& matriz, int line);
		void getLinePieceFromAMatrix(vector<double>& lineMatrix, vector< vector<double> >& matriz, int line, int initColumn);
        int getMinValuePosition(vector<double> vetor);
		int getMaxValuePosition(vector<double> vetor);
        void removeLineAndColumn(vector< vector<double> >& matriz, int numLinha, int numColuna);
        void copyVector(vector<double>& line, vector<double> original, int posInicial, int posFinal);
        bool contains(vector<Individual*>& vetor, Individual* individual);
		void survive(vector<Individual*>& individuals, int n);
        void dump(vector<Individual*> individuals);
		void dumpPopulation();
		void fullDump(vector<Individual*> individuals);
		void dumpIndividual(Individual* ind);		
		int getMetaHeuristc();	
		string getMetaHeuristcaByNumber(int number);
		void competitiveAbordage(vector<double> hipervolumes, vector< vector< Individual*> >& populationsOfEachMetaheuristic);
		void fuzzyAbordage(vector<double> hipervolumes, vector< vector< Individual*> >& populationsOfEachMetaheuristic);
		void cooperativeAbordage(vector< vector<Individual*> >& populationsOfEachMetaheuristic);		
		int comparePopulations(vector<vector<Individual*>>& populationsOfEachMetaheuristic);
		int randomPosition(Population* & pop);
		int randomPosition(int size);
		void modifyPopulation(int rank, Population* population);
		vector< vector<Individual*> > modifyPopulationSingle(Population* population);
		vector<double> calcHiperVolume(vector< vector<Individual*> >& individualsOffAllMetaheuristics);
		bool disperse(vector<Individual*> v, Individual * id);
		void sendMessage(int numOfProcesses);
		double receiveMessage();
		void waitMessageForMe(int tag);
		void sendMatrix(int numOfProcesses, double* matrix, size_t numLines, size_t numCollums);
		double* waitMatrixForMe(int tag, size_t numLines, size_t numCollums);
		double* waitMatrixWithHypervolumeForMe(int tag, size_t numLines, size_t numCollums);
		void clearVector(vector<Individual*>& v);
		void copyVectorIndividuals(vector<Individual*>& one, vector<Individual*>& two);
		void initBounds(double lowerBound, double upperBound);
		void initBounds();
};

#endif // SPEA2_H
