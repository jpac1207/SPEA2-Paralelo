#include "MatTransform.h"



MatTransform::MatTransform()
{
}


MatTransform::~MatTransform()
{
}

/*Return an array with all individual gens*/
double * MatTransform::getIndividualsInMatrix(vector<Individual*> individuals)
{
	size_t numLines = individuals.size();
	size_t col = individuals.at(0)->getGenes().size();	
	double* matrix = (double *) malloc( (numLines * col) * sizeof(double));
	int cont = 0;
	
	for (size_t i = 0; i < numLines; i++) {		

		for (size_t j = 0; j < col; j++) {
			matrix[cont] = individuals.at(i)->getGenes()[j]; //put genes from individual "i" in line "i" of matrix
			cont++;
		}					
	}	

	return matrix;
}

double * MatTransform::getIndividualsInMatrixWithHipervolume(vector<Individual*> individuals, double hipervolume)
{
	size_t numLines = individuals.size();
	size_t col = individuals.at(0)->getGenes().size();
	double* matrix = (double *) malloc( ( (numLines * col) + 1 ) * sizeof(double));
	int cont = 1;

	matrix[0] = hipervolume;

	for (size_t i = 0; i < numLines; i++) {

		for (size_t j = 0; j < col; j++) {
			matrix[cont] = individuals.at(i)->getGenes()[j]; //put genes from individual "i" in line "i" of matrix
			cont++;
		}
	}

	return matrix;
}

vector<Individual*> MatTransform::getIndividualsInVector(double* individuals, size_t numLines, size_t numCollums)
{
	vector<Individual*> individualsVector (numLines);
	int count = 0;

	for (size_t i = 0; i < numLines; i++) {

		vector<double> genes;

		for (size_t j = 0; j < numCollums; j++) {
			genes.push_back(individuals[count]);
			count++;
		}

		Individual* id = new Individual(genes.size());
		id->setGenes(genes);
		id->setQtdGenes(genes.size());
		id->setAptidao(this->evaluateIndividual(id));
		individualsVector[i] = (id);
	}
	return individualsVector;
}

vector<Individual*> MatTransform::getIndividualsInVectorWithHipervolume(double * individuals, size_t numLines, size_t numCollums)
{
	vector<Individual*> individualsVector(numLines);
	int count = 1;

	for (size_t i = 0; i < numLines; i++) {

		vector<double> genes;

		for (size_t j = 0; j < numCollums; j++) {
			genes.push_back(individuals[count]);
			count++;
		}

		Individual* id = new Individual(genes.size());
		id->setGenes(genes);
		id->setQtdGenes(genes.size());
		id->setAptidao(this->evaluateIndividual(id));
		individualsVector[i] = (id);
	}
	return individualsVector;
}

vector<double> MatTransform::evaluateIndividual(Individual* individual) {

	ZDT* zdt1 = new ZDT();
	vector<double> value =  zdt1->evaluateIndividual(individual);
	delete zdt1;

	return value;
}

double * MatTransform::allocMatrix(size_t numLines, size_t numCollums)
{
	double* matrix = (double *) malloc((numLines * numCollums) * sizeof(double));

	return matrix;
}

double * MatTransform::allocMatrixWithHypervolume(size_t numLines, size_t numCollums)
{
	double* matrix = (double *) malloc( ( (numLines * numCollums) + 1 ) * sizeof(double));

	return matrix;
}
