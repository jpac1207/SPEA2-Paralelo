#include <iostream>
#include<stdlib.h>
#include<stdio.h>
//#include<cstdlib>
#include<time.h>
#include"mpi.h"

#include"Individual.h"
#include"Population.h"
#include"PopulationReader.h"
#include"HyperVolumeCalculator.h"
#include"SPEA2.h"
#include <sstream>

using namespace std;

int main(int argc, char* argv[])
{
	/*
		Implementado por Jo�o Pedro Augusto Costa
		Instituto Federal do Maranh�o
		Campus S�o Lu�s, Monte Castelo
		Departamento Acad�mico de Informatica
	*/

	srand((unsigned int)time(NULL));
	/* to identify the current process	*/
	int rank = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (int i = 0; i < 2; i++) {

		SPEA2* spea2 = new SPEA2(rank);
		spea2->run();
		delete spea2;
	}	
	
	MPI_Finalize();
	if(rank == 0)
		std::getchar();
	
	return 0;
}
