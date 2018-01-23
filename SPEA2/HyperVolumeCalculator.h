#pragma once

#include<vector>
#include<algorithm>
#include<iostream>
#include<math.h>

#include"Individual.h"
#include"Population.h"

using namespace std;

class HyperVolumeCalculator
{
public:
	HyperVolumeCalculator();
	~HyperVolumeCalculator();
	double calculateForTwoObjective(vector<Individual*> individuals);
	double calculateSpread(vector<Individual*> individuals);
	vector<double> differenceBetweenExtremalPointsAndReference(vector<Individual*> individuals);
	vector<double> identifyExtremalPoints(vector<Individual*> individuals);
private:
	double hypervolume = 0.0;
	double referenceOne = 11.0;
	double referenceTwo = 11.0;	
	double euclideanDistance(vector<double> one, vector<double> two);
};

