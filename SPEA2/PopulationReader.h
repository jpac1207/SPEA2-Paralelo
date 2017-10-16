#pragma once

#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include <sstream>
#include"ZDT.h"

#include"Individual.h"

using namespace std;

class PopulationReader
{
public:
	PopulationReader();
	~PopulationReader();
	vector<Individual*> loadFromArquive(string path);
	vector<Individual*> loadPopulation(int rank);
	vector<Individual*> getMasterPopulation(int rank);
	vector<string> split(string text, char delimiter);
	void savePopulation(vector<Individual*> pop, int rank);
	void savePopulationForAll(vector<Individual*> pop, int numberOfSlaves);
	vector<double> evaluateIndividual(Individual* individual);
};

