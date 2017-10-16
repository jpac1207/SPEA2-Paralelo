#pragma once
#include "Individual.h"
#include "HyperVolumeCalculator.h"

#define HIGH 1
#define MEDIUM 0.5
#define LOW 0.1

using namespace std;
class FuzzyController
{
public:
	FuzzyController();
	~FuzzyController();
	vector<Individual*> inference(vector<double> lvHipervolumes, vector< vector< Individual*> > populations);
	int getOutputVariables();
private:
	int outputVariable;
	const double MaximumZdt1HiperVolume = 120.657446;
	const double MaximumZdt3HiperVolume = 128.773637;
	int Fuzzyfication(vector<double> hipervolumDifferences, vector< vector<double> > differencesPoints);
	double Aggregation(double inputOne, double inputTwo, double inputThree);
	double getLinguisticVariableHipervolume(double value);
	double getLinguisticVariableFirstPoint(double value);
	double getLinguisticVariableLastPoint(double value);
	double OR(double varOne, double varTwo);
	double AND(double varOne, double varTwo);
	void clearVector(vector<Individual*>& v);
};

