#ifndef ZDT_H
#define ZDT_H
#define PI (3.141592653589793)

#include<iostream>
#include<stdio.h>
#include<math.h>
#include"ObjectiveFunction.h"

using namespace std;

class ZDT : public ObjectiveFunction
{
public:
	ZDT();
	virtual ~ZDT();
	vector<double> evaluateIndividual(Individual* i);
	bool isMinimization(int function);
	int numberOfFunctions();

protected:

private:
};

#endif // ZDT_H
