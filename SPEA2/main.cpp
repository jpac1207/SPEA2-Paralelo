#include <iostream>
#include<stdlib.h>
#include<stdio.h>
//#include<cstdlib>
#include<time.h>
#include <sstream>

#include"mpi.h"

#include"Individual.h"
#include"Population.h"
#include"PopulationReader.h"
#include"HyperVolumeCalculator.h"
#include"SPEA2.h"
#include"TestUtil.h"


using namespace std;

int main(int argc, char* argv[])
{
	/*
		Implementado por Jo�o Pedro Augusto Costa
		Instituto Federal do Maranh�o
		Campus S�o Lu�s, Monte Castelo
		Departamento Acad�mico de Informatica
	*/

	/*TestUtil t;
	t.run();
	getchar();*/


	srand((unsigned int)time(NULL));

	int rank = 0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for (int i = 0; i < 31; i++) {

		SPEA2* spea2 = new SPEA2(rank);
		spea2->run();
		delete spea2;
	}

	MPI_Finalize();
	if (rank == 0) {
		cout << "FIM" << endl;
		std::getchar();
	}

	return 0;
}
