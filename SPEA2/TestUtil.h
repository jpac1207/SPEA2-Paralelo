#pragma once
#include "Individual.h"
#include"ZDT.h"
#include"PopulationReader.h"
#include"HyperVolumeCalculator.h"
#include<vector>

using namespace std;
class TestUtil
{
public:
	TestUtil();
	~TestUtil();
	static void run();
	static vector<Individual*> nonDominatedSol(vector<Individual*> individuals);
	static bool isDominated(Individual * one, Individual * two);
	static void clearVector(vector<Individual*>& v);
	static void dump(vector<Individual*> individuals);
};

