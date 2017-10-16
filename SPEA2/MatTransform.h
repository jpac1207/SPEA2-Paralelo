#pragma once
#include"Individual.h"
#include"ZDT.h"

#include<iostream>
#include<vector>
#include<stdlib.h>

using namespace std;

class MatTransform
{
public:
	MatTransform();
	~MatTransform();
	double* getIndividualsInMatrix(vector<Individual*> individuals);
	double* getIndividualsInMatrixWithHipervolume(vector<Individual*> individuals, double hipervolume);
	vector<Individual*> getIndividualsInVector(double*  individuals, size_t numLines, size_t numCollums);	
	vector<Individual*> getIndividualsInVectorWithHipervolume(double*  individuals, size_t numLines, size_t numCollums);
	vector<double> evaluateIndividual(Individual* individual);
	double* allocMatrix(size_t numLines, size_t numCollums);
	double * allocMatrixWithHypervolume(size_t numLines, size_t numCollums);
};

