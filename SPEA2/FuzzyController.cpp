#include "FuzzyController.h"



FuzzyController::FuzzyController()
{
}


FuzzyController::~FuzzyController()
{
}

//Returns the number of choiced population
vector<Individual*> FuzzyController::inference(vector<double> lvHipervolumes, vector<vector<Individual*>> populations)
{
	vector<double> differencesHipervolumes;
	vector< vector<double> > differencesPoints;
	HyperVolumeCalculator h;

	for (size_t i = 0; i < populations.size(); i++) {

		differencesPoints.push_back(h.differenceBetweenExtremalPointsAndReference(populations.at(i)));
		differencesHipervolumes.push_back(MaximumZdt3HiperVolume - lvHipervolumes.at(i));//get difference between real and found hipervolume		
	}

	this->outputVariable = this->Fuzzyfication(differencesHipervolumes, differencesPoints);

	// delete pointers not useds
	for (size_t i = 0; i < populations.size(); i++) {

		if (i != outputVariable) {

			clearVector(populations[i]);
		}
	}

	return populations.at(this->outputVariable);
}

int FuzzyController::getOutputVariables()
{
	return this->outputVariable;
}

int FuzzyController::Fuzzyfication(vector<double> hipervolumDifferences, vector<vector<double>> differencesPoints)
{
	vector<double> values;
	//get only decimal parts of numbers
	for (size_t i = 0; i < hipervolumDifferences.size(); i++) {

		double x = hipervolumDifferences.at(i);
		double xLine = this->getLinguisticVariableHipervolume(x); // get a value from hipervolume difference
		double y = differencesPoints.at(i).at(0); // Tem que ser pequeno
		double yLine = this->getLinguisticVariableFirstPoint(y);
		double z = differencesPoints.at(i).at(1);
		double zLine = this->getLinguisticVariableLastPoint(z);

		/*cout << "Heuristic(" << i << ") --> " << " x -> "<<xLine<<" y -> " << yLine<<" z -> "<<zLine<< endl;*/

		values.push_back(this->Aggregation(xLine, yLine, zLine));
	}

	double data = this->OR(values.at(0), this->OR(values.at(1), values.at(2)));

	for (size_t i = 0; i < values.size(); i++) {

		if (values.at(i) == data) {

			return i;
		}
	}
}
/*Centroide method with a 10 range*/
double FuzzyController::Aggregation(double inputOne, double inputTwo, double inputThree)
{
	int initPoint = 0;
	int finalPoint = 10;
	double dividendo = 0;
	double divisor = 0;
	double variable = inputOne;
	int times = 0;
	bool ok = false;

	for (double i = initPoint; i <= finalPoint; i++) {

		dividendo += i;
		times++;

		if (i == finalPoint)
			ok = !ok; //in last interation

		if ((times == 3) || ok) {
			dividendo *= variable;
			divisor += (variable * 3);

			if (variable == inputOne)
				variable = inputTwo;
			else
				variable = inputThree;

			times = 0;
		}

	}

	return (dividendo / divisor);
}

double FuzzyController::getLinguisticVariableHipervolume(double value) {

	if (value <= 8.5)
		return HIGH;
	else if (value < 15.0)
		return MEDIUM;
	else
		return LOW;
}

double FuzzyController::getLinguisticVariableFirstPoint(double value) {

	if (value >= 22.0)
		return HIGH;
	else if (value > 21.00 && value < 22.00)
		return MEDIUM;
	else
		return LOW;
}

double FuzzyController::getLinguisticVariableLastPoint(double value) {

	if (value >= 21.065)
		return HIGH;
	else if (value > 21.06 && value < 21.065)
		return MEDIUM;
	else
		return LOW;
}

double FuzzyController::OR(double varOne, double varTwo) {

	if (varOne > varTwo) return varOne;
	else return varTwo;
}

double FuzzyController::AND(double varOne, double varTwo) {

	if (varOne < varTwo) return varOne;
	else return varTwo;
}

void FuzzyController::clearVector(vector<Individual*>& v)
{
	for (size_t i = 0; i < v.size(); i++) {
		delete v[i];
	}
	v.clear();
}
