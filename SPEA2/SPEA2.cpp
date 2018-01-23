#include "SPEA2.h"

SPEA2::SPEA2()
{
	//ctor
}

SPEA2::SPEA2(int argc, char* argv[])
{
	this->argc = argc;
	this->argv = argv;
}

SPEA2::SPEA2(int rank)
{
	this->rank = rank;
}

SPEA2::~SPEA2()
{
	for (size_t j = 0; j < this->archive.size(); j++) delete archive.at(j);
	archive.clear();
	delete objective;
	delete population;
}

bool compare(double d1, double d2) {

	return (d1 < d2);
}

bool compareIndividual(Individual* i, Individual* j) {

	double f1 = i->getFitness();
	double f2 = j->getFitness();

	return (f1 < f2);
}

bool compareDensity(Individual* i, Individual* j) {

	double f1 = i->getD();
	double f2 = j->getD();

	return (f1 < f2);
}

void SPEA2::run() {

	this->lastHipervolume = 0.0;
	this->objective = new ZDT();
	this->population = new Population();
	this->population->setSize(this->populationSize);
	this->initBounds();
	this->population->initPopulation(this->limiteInferior, this->limiteSuperior, this->qtdGenes);
	this->archive = this->nonDominatedSol(this->population->getIndividuals());
	HyperVolumeCalculator* hypervolume = new HyperVolumeCalculator();
	MatTransform* matTransform = new MatTransform();
	this->lastHipervolume = hypervolume->calculateForTwoObjective(this->archive);

	for (int i = 0; i < this->numberOfIterations; i++) {

		if (this->rank == 0) {

			this->clearVector(this->RT);
			RT = this->fillRT(population->getIndividuals(), archive);
			this->compute_s(RT);
			this->compute_raw(RT);
			this->compute_d(RT);
			this->clearVector(archive);
			this->clearVector(index);

			this->getNonDominatedSolutions(index, RT);

			if (index.size() == this->archiveSize) {
				this->copyVectorIndividuals(index, archive);
			}
			else if (index.size() < this->archiveSize) {
				this->copyVectorIndividuals(index, archive);
				this->fillArchive(RT, this->archive, index);
			}
			else {
				this->clustering(this->archive, index);
			}

			if (i == (this->numberOfIterations - 1)) break;

			double* matrixIndividuals = matTransform->getIndividualsInMatrix(this->population->getIndividuals());	// get a matrix with the individuals
			this->sendMatrix(3, matrixIndividuals, (size_t) this->populationSize, (size_t) this->qtdGenes); // send matrix with individuals to all slaves
			free(matrixIndividuals);

			vector<double> lvHipervolumes; //to store the hipervolume of each metaheuristic
			vector< vector<Individual*> > lvIndividuals;// to store individuals of each metaheuristic 

			int lvCount = 0;
			while (lvCount < this->numberOfSlaves) {
				double* matrixOfEachHeuristic = this->waitMatrixWithHypervolumeForMe(rank, this->populationSize, this->qtdGenes);
				lvHipervolumes.push_back(matrixOfEachHeuristic[0]); // get the hipervolume founded			
				lvIndividuals.push_back(matTransform->getIndividualsInVectorWithHipervolume(matrixOfEachHeuristic, this->populationSize, this->qtdGenes));
				free(matrixOfEachHeuristic);
				lvCount++;
			}

			/*this->fuzzyAbordage(lvHipervolumes, lvIndividuals);*/
			this->competitiveAbordage(lvHipervolumes, lvIndividuals);
			/*this->cooperativeAbordage(lvIndividuals);*/

			/*clear copys in RT and index*/
			this->clearVector(RT);
			this->clearVector(index);
		}
		else {

			if (i == (this->numberOfIterations - 1)) break;

			double* matrixIndividual = this->waitMatrixForMe(rank, (size_t) this->populationSize, (size_t) this->qtdGenes); // wait matrix of master
			this->clearVector(this->population->getIndividuals());
			this->population->setIndividuals(matTransform->getIndividualsInVector(matrixIndividual, this->populationSize, this->qtdGenes)); //get population saved for master					
			this->modifyPopulation(rank, population); //Modify the population, save and notify master		
			free(matrixIndividual);
		}
	}

	if (rank == 0) {
		/*double hp = hypervolume->calculateForTwoObjective(this->archive);
		cout << hp << endl;*/
		dump(archive);	
		/*fullDump(archive);*/
	}

	this->clearVector(RT);//clear copys in RT	
	this->clearVector(index);
	delete hypervolume;
	delete matTransform;
}

bool SPEA2::isDominated(Individual* one, Individual* two) {

	bool anyDominate = false;

	for (size_t i = 0; i < one->getAptidao().size(); i++) {

		if (this->objective->isMinimization(i)) {

			if (one->getAptidao()[i] < two->getAptidao()[i])
				return false;
			else if (one->getAptidao()[i] > two->getAptidao()[i])
				anyDominate = true;
		}
		else {

			if (one->getAptidao()[i] > two->getAptidao()[i])
				return false;
			else if (one->getAptidao()[i] < two->getAptidao()[i])
				anyDominate = true;
		}
	}

	return anyDominate;
}

vector<Individual*> SPEA2::fillRT(vector<Individual*> pop, vector<Individual*> archive)
{
	vector<Individual*> RT;

	for (size_t i = 0; i < pop.size(); i++) {
		Individual* ind = new Individual(*pop.at(i));
		RT.push_back(ind);
	}
	for (size_t i = 0; i < archive.size(); i++) {
		Individual* ind = new Individual(*archive.at(i));
		RT.push_back(ind);
	}

	return RT;
}

void SPEA2::compute_s(vector<Individual*>& RT) {

	for (unsigned int i = 0; i < RT.size(); i++) {

		for (unsigned int j = i + 1; j < RT.size(); j++) {

			if (this->isDominated(RT.at(j), RT.at(i))) {//Se j for dominado por i				
				RT[i]->setS((RT.at(i)->getS() + 1));
			}
			else if (this->isDominated(RT.at(i), RT.at(j))) {//Se i for dominado por j				
				RT[j]->setS((RT.at(j)->getS() + 1));
			}
		}
	}
}

void SPEA2::compute_raw(vector<Individual*>& RT) {

	for (unsigned int i = 0; i < RT.size(); i++) {

		for (unsigned int j = i + 1; j < RT.size(); j++) {

			if (this->isDominated(RT.at(i), RT.at(j)))// se i for dominado por j
				RT[i]->setRaw(RT.at(i)->getRaw() + RT.at(j)->getS());
			else if (this->isDominated(RT.at(j), RT.at(i)))// se j for dominado por i
				RT[j]->setRaw(RT.at(j)->getRaw() + RT.at(i)->getS());
		}
	}
}

void SPEA2::compute_d(vector<Individual*>& RT) {

	vector< vector<double> > distanceMatrix;
	size_t sizeOfPopulation = RT.size();
	this->distances(distanceMatrix, RT);
	orderMatrix(distanceMatrix);

	for (unsigned int i = 0; i < sizeOfPopulation; i++) {
		RT[i]->setD(1 / (distanceMatrix[i][this->k] + 2));
		RT[i]->setFitness(RT[i]->getRaw() + RT[i]->getD());
	}
}

void SPEA2::fillArchive(vector<Individual*>& RT, vector<Individual*>& archive, vector<Individual*>& index) {

	int i = 0;

	sort(RT.begin(), RT.end(), compareIndividual);

	/* As soluções com fitness menor que um já foram
	adicionadas no archive pelo index*/
	while (archive.size() < this->archiveSize) {

		/*	if (!this->contains(archive, RT[i])) {*/
		Individual* ind = new Individual(*RT.at(i));
		archive.push_back(ind);
		/*	}*/
		i++;
	}
}

vector<Individual*> SPEA2::nonDominatedSol(vector<Individual*> individuals) {

	vector<Individual*> nonDominated;
	unsigned j = 0;

	for (unsigned i = 0; i < individuals.size(); i++) {

		bool isDominatedByLeastOne = false;

		for (j = 0; j < individuals.size(); j++) {
			if (i != j) {
				if (this->isDominated(individuals[i], individuals[j])) {
					isDominatedByLeastOne = true;
					break;
				}
			}
		}

		if (!isDominatedByLeastOne) {
			Individual* id = new Individual(*individuals.at(i));
			nonDominated.push_back(id);
		}

	}

	return nonDominated;
}

void SPEA2::getNonDominatedSolutions(vector<Individual*>& index, vector<Individual*> RT) {

	vector<Individual*> newVector;

	for (unsigned int j = 0; j < RT.size(); j++) {

		Individual* ind = new Individual(*RT.at(j));

		if (ind->getFitness() < 1)
			newVector.push_back(ind);
		else
			delete ind;
	}

	index.swap(newVector);
}

void SPEA2::clustering(vector<Individual*>& archive, vector<Individual*>& individuos) {

	size_t qtdToRemove = (individuos.size() - this->archiveSize);
	vector< vector<double> > matrix;
	this->distances(matrix, individuos);

	for (int i = 0; i < qtdToRemove; i++) {

		vector<int> indexesToRemove;
		vector<double> value;

		for (unsigned int j = 0; j < individuos.size(); j++) {

			if (j != (individuos.size() - 1)) {
				vector<double> lineMatrix;
				this->getLinePieceFromAMatrix(lineMatrix, matrix, j, (j + 1));
				indexesToRemove.push_back(this->getMinValuePosition(lineMatrix) + (j + 1));
				value.push_back(matrix[j][(indexesToRemove[j])]);
			}
		}
		int lineToRemove = getMinValuePosition(value);
		int colToRemove = indexesToRemove.at(lineToRemove);
		this->removeLineAndColumn(matrix, lineToRemove, colToRemove);
		delete individuos[lineToRemove]; // desalocatte memory
		individuos.erase(individuos.begin() + (lineToRemove));
		indexesToRemove.clear();
		value.clear();
	}

	archive.swap(individuos);
}

void SPEA2::genericalClustering(vector<Individual*>& individuos, int nToRemove) {

	size_t qtdToRemove = nToRemove;
	int finalSize = ((int)individuos.size() - nToRemove);

	while (individuos.size() > (finalSize)) {
		int pos_one = 0 + (rand() % individuos.size());
		int pos_two = 0 + (rand() % individuos.size());

		if (this->isDominated(individuos[pos_one], individuos[pos_two])) {
			individuos.erase(individuos.begin() + pos_one);
		}
		else {
			individuos.erase(individuos.begin() + pos_two);
		}
	}
}

void SPEA2::removeLineAndColumn(vector< vector<double> >& matriz, int numLinha, int numColuna) {

	//quadratic matrix
	if (matriz.size() == matriz.at(numLinha).size()) {
		for (unsigned int i = 0; i < matriz.size(); i++) {
			matriz[numLinha][i] = DBL_MAX;
			matriz[i][numColuna] = DBL_MAX;
		}
	}
	else {
		for (unsigned int i = 0; i < matriz.at(numLinha).size(); i++) {
			matriz[numLinha][i] = DBL_MAX;
		}
		for (unsigned int j = 0; j < matriz.size(); j++) {
			matriz[j][numColuna] = DBL_MAX;
		}
	}
}

int SPEA2::getMinValuePosition(vector<double> vetor) {

	double minValue = vetor.at(0);
	int pos = 0;

	for (unsigned int i = 1; i < vetor.size(); i++) {

		if (vetor.at(i) < minValue) {
			minValue = vetor.at(i);
			pos = i;
		}
	}
	return pos;
}

int SPEA2::getMaxValuePosition(vector<double> vetor)
{
	double maxValue = vetor.at(0);
	int pos = 0;

	for (unsigned int i = 1; i < vetor.size(); i++) {

		if (vetor.at(i) > maxValue) {
			maxValue = vetor.at(i);
			pos = i;
		}
	}
	return pos;
}

void SPEA2::getLineFromAMatrix(vector<double>& lineMatrix, vector< vector<double> >& matriz, int line) {

	vector<double> vetor(matriz.at(line).size());

	for (unsigned int i = 0; i < matriz.at(line).size(); i++) {
		vetor.push_back(matriz.at(line).at(i));
	}

	lineMatrix.swap(vetor);
}

void SPEA2::getLinePieceFromAMatrix(vector<double>& lineMatrix, vector< vector<double> >& matriz, int line, int initColumn) {

	vector<double> vetor(matriz.at(line).size());

	for (unsigned int i = initColumn; i < matriz.at(line).size(); i++) {
		vetor.push_back(matriz.at(line).at(i));
	}

	lineMatrix.swap(vetor);
}

void SPEA2::distances(vector< vector<double> >& matrix, vector<Individual*> individuos) {

	size_t individualsSize = individuos.size();
	vector< vector<double> > distanceMatrix;

	/*for (int i = 0; i < individualsSize; i++) {

		vector<double> row(individualsSize);
		distanceMatrix.push_back(row);
		distanceMatrix[i][i] = DBL_MAX;
	}*/

	/****   line 0  ****/
	vector<double> row(individualsSize);
	distanceMatrix.push_back(row);
	distanceMatrix[0][0] = DBL_MAX;

	for (int i = 0; i < individualsSize; i++) {

		Individual* one = individuos.at(i);

		for (unsigned int j = i + 1; j < individualsSize; j++) {

			if (i == 0) {
				vector<double> row(individualsSize);
				distanceMatrix.push_back(row);
				distanceMatrix[j][j] = DBL_MAX;
			}

			Individual* two = individuos.at(j);

			double somatorioDistanciaPontos = 0;

			for (unsigned int k = 0; k < one->getAptidao().size(); k++) {

				somatorioDistanciaPontos += pow(fabs(one->getAptidao()[k] - two->getAptidao()[k]), 2);
			}

			double differenceSquad = sqrt(somatorioDistanciaPontos);
			distanceMatrix[i][j] = differenceSquad;
			distanceMatrix[j][i] = differenceSquad;
		}
	}

	matrix.swap(distanceMatrix);
}

void SPEA2::orderMatrix(vector< vector<double> >& matriz) {

	int count = 0;
	vector<double> arrayAux;
	size_t length = matriz.size();

	for (int i = 0; i < length; i++) {
		vector<double> arrayaux;
		for (int j = 0; j < length; j++) {
			arrayaux.push_back(matriz[i][j]);
		}
		std::sort(arrayaux.begin(), arrayaux.end(), compare);
		for (int j = 0; j < length; j++) {
			matriz[i][j] = arrayaux.at(j);
		}
		arrayaux.clear();
	}

	/*for (unsigned int i = 0; i < length; i++) {
		for (unsigned int j = 0; j < length; j++) {
			arrayAux.push_back(matriz[i][j]);
		}
	}

	std::sort(arrayAux.begin(), arrayAux.end(), compare);

	for (unsigned int i = 0; i < length; i++) {
		for (unsigned int j = 0; j < length; j++) {
			matriz[i][j] = arrayAux.at(count);
			count++;
		}
	}*/
}

void SPEA2::copyVector(vector<double>& line, vector<double> original, int posInicial, int posFinal) {

	vector<double> copy(posFinal - posInicial);

	for (int i = posInicial; i < posFinal; i++) {

		copy.push_back(original[i]);
	}

	line.swap(copy);
}

bool SPEA2::contains(vector<Individual*>& vetor, Individual* individual) {

	for (unsigned int i = 0; i < vetor.size(); i++) {
		bool flag = true;

		for (unsigned int j = 0; j < individual->getAptidao().size(); j++) {
			if (vetor[i]->getAptidao()[j] != individual->getAptidao()[j]) {
				flag = false;
				break;
			}
		}

		if (flag)
			return true;
	}

	return false;
}

void SPEA2::survive(vector<Individual*>& individuals, int n) {

	sort(individuals.begin(), individuals.end(), compareDensity);

	while (individuals.size() > n) {

		individuals.pop_back();
	}
}

void SPEA2::dump(vector<Individual*> individuals) {

	/*cout << "\n Population: " << endl;*/
	for (unsigned int i = 0; i < individuals.size(); i++) {

		cout << individuals.at(i)->getAptidao()[0] << ";" << individuals.at(i)->getAptidao()[1] << endl;
	}
}

void SPEA2::dumpPopulation() {

	cout << "\n Population: " << endl;
	for (unsigned int i = 0; i < this->population->getIndividuals().size(); i++) {

		cout << this->population->getIndividuals()[i]->getAptidao()[0] << ";" << this->population->getIndividuals()[i]->getAptidao()[1] << endl;
	}

}

void SPEA2::fullDump(vector<Individual*> individuals) {

	for (size_t i = 0; i < individuals.size(); i++) {

		this->dumpIndividual(individuals[i]);
	}
}

void SPEA2::dumpIndividual(Individual* ind) {
	/*cout << "Individual: " << "";*/
	for (unsigned int i = 0; i < ind->getGenes().size(); i++) {

		cout << ind->getGenes()[i] << ";";
	}
	cout << "" << endl;
}

int SPEA2::getMetaHeuristc()
{
	if ((probabilities[0] == probabilities[1]) && (probabilities[1] == probabilities[2])) {
		int number = (rand() % abs(qtdMetaHeuristcs));
		return number;
	}
	else {
		double number = (0 + (((double)rand() / (double)(RAND_MAX)) * 100));

		if (number < probabilities[0])
			return 0;
		else if (number >= probabilities[0] && number <= (probabilities[0] + probabilities[1])) {
			return 1;
		}
		else {
			return 2;
		}
	}
}

string SPEA2::getMetaHeuristcaByNumber(int number)
{
	switch (number)
	{
	case 0: {
		return "GA";
	}
	case 1: {
		return "DE";
	}
	case 2: {
		return "PSO";
	}
	default:
		break;
	}
}

void SPEA2::competitiveAbordage(vector<double> hipervolumes, vector< vector< Individual*> >& populationsOfEachMetaheuristic)
{
	int pos = this->getMaxValuePosition(hipervolumes);
	double bestHipervolumeFounded = hipervolumes.at(pos);
	/*cout << bestHipervolumeFounded << " | " << this->lastHipervolume << endl;*/
	if (bestHipervolumeFounded > this->lastHipervolume) {
		vector<Individual*> newPop = populationsOfEachMetaheuristic.at(pos);
		this->clearVector(this->population->getIndividuals());
		this->population->setIndividuals(newPop);
		this->lastHipervolume = bestHipervolumeFounded;
	}
	else {
		pos = -1; 	//to delete all populations obtained by heuristics if their are not used
	}

	//desallocating space of not choosen heuristics
	for (size_t i = 0; i < populationsOfEachMetaheuristic.size(); i++) {

		if (i != pos) {
			this->clearVector(populationsOfEachMetaheuristic[i]);
		}
	}
	if (pos == -1)
		populationsOfEachMetaheuristic.clear();
}

void SPEA2::fuzzyAbordage(vector<double> hipervolumes, vector<vector<Individual*>>& populationsOfEachMetaheuristic)
{
	FuzzyController* controller = new FuzzyController();
	vector<Individual*> newPop = controller->inference(hipervolumes, populationsOfEachMetaheuristic);
	this->clearVector(this->population->getIndividuals());
	this->population->setIndividuals(newPop);
	delete controller;
}

void SPEA2::cooperativeAbordage(vector< vector<Individual*> >& populationsOfEachMetaheuristic)
{
	HyperVolumeCalculator* calc = new HyperVolumeCalculator();
	vector<Individual*> allIndividuals;

	for (size_t i = 0; i < populationsOfEachMetaheuristic.size(); i++) {
		this->copyVectorIndividuals(populationsOfEachMetaheuristic.at(i), allIndividuals);
	}

	this->genericalClustering(allIndividuals, (allIndividuals.size() - this->populationSize));	
	vector<Individual*> nonDominated = this->nonDominatedSol(allIndividuals);
	double value = calc->calculateForTwoObjective(nonDominated);
	
	if (value > lastHipervolume) {
		this->clearVector(this->population->getIndividuals());
		this->population->setIndividuals(allIndividuals);
		this->lastHipervolume = value;
	}

	/*desallocating memory of populationsOfEachMetaheuristic,
	all of used individuals	of this vectors are passed just like copy*/
	for (size_t i = 0; i < populationsOfEachMetaheuristic.size(); i++) {
		this->clearVector(populationsOfEachMetaheuristic[i]);
	}
	populationsOfEachMetaheuristic.clear();
	this->clearVector(nonDominated);
	delete calc;
	
}

int SPEA2::comparePopulations(vector< vector< Individual*> >& populationsOfEachMetaheuristic) {

	vector<double> nonDominatedForEach;

	for (size_t i = 0; i < populationsOfEachMetaheuristic.size(); i++) {
		vector<Individual*> lvAux = this->nonDominatedSol(populationsOfEachMetaheuristic.at(i));
		nonDominatedForEach.push_back(lvAux.size());
		clearVector(lvAux);
	}

	return this->getMaxValuePosition(nonDominatedForEach);
}

int SPEA2::randomPosition(Population *& pop)
{
	return 0 + (rand() % abs((int)pop->getIndividuals().size()));
}

int SPEA2::randomPosition(int size)
{
	return 0 + (rand() % abs(size));;
}

void SPEA2::modifyPopulation(int rank, Population* population) {

	HyperVolumeCalculator* calc = new HyperVolumeCalculator();
	PopulationReader* p = new PopulationReader();
	MatTransform* mat = new MatTransform();

	if (rank == 1) {

		Population* toBeModifiedByGA = population;
		GA* ga = new GA();
		ga->setCrossOverRate(0.8);
		ga->setMutationRate(0.2);
		ga->setTheta(0.01);
		ga->setLimiteInferior(this->limiteInferior);
		ga->setLimiteSuperior(this->limiteSuperior);
		ga->run(toBeModifiedByGA);
		vector<Individual*> nonDominated = this->nonDominatedSol(toBeModifiedByGA->getIndividuals());
		double hypervolumeFound = calc->calculateForTwoObjective(nonDominated);
		this->clearVector(nonDominated);
		if (debug) cout << "Finished GA " << hypervolumeFound << endl;
		double* matrix = mat->getIndividualsInMatrixWithHipervolume(toBeModifiedByGA->getIndividuals(), hypervolumeFound); //first position is hypervolume		
		MPI_Send(matrix, ((this->populationSize * this->qtdGenes) + 1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		this->clearVector(toBeModifiedByGA->getIndividuals());
		free(matrix);
		delete ga;
	}

	if (rank == 2) {

		Population* toBeModifiedByDE = population;
		DE* de = new DE();
		de->setCrossoverRate(0.5);
		de->setLimiteSuperior(this->limiteSuperior);
		de->setLimiteInferior(this->limiteInferior);
		de->run(toBeModifiedByDE);
		vector<Individual*> nonDominated = this->nonDominatedSol(toBeModifiedByDE->getIndividuals());
		double hypervolumeFound = calc->calculateForTwoObjective(nonDominated);
		this->clearVector(nonDominated);
		if (debug) cout << "Finished DE " << hypervolumeFound << endl;
		double* matrix = mat->getIndividualsInMatrixWithHipervolume(toBeModifiedByDE->getIndividuals(), hypervolumeFound); //first position is hypervolume		
		MPI_Send(matrix, ((this->populationSize * this->qtdGenes) + 1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		this->clearVector(toBeModifiedByDE->getIndividuals());
		free(matrix);
		delete de;
	}

	if (rank == 3) {

		Population* toBeModifiedByPSO = population;
		PSO* pso = new PSO();
		pso->setLimiteInferior(this->limiteInferior);
		pso->setLimiteSuperior(this->limiteSuperior);
		pso->run(toBeModifiedByPSO);
		vector<Individual*> nonDominated = this->nonDominatedSol(toBeModifiedByPSO->getIndividuals());
		double hypervolumeFound = calc->calculateForTwoObjective(nonDominated);
		this->clearVector(nonDominated);
		if (debug) cout << "Finished PSO " << hypervolumeFound << endl;
		double* matrix = mat->getIndividualsInMatrixWithHipervolume(toBeModifiedByPSO->getIndividuals(), hypervolumeFound); //first position is hypervolume		
		MPI_Send(matrix, ((this->populationSize * this->qtdGenes) + 1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		this->clearVector(toBeModifiedByPSO->getIndividuals());
		free(matrix);
		delete pso;
	}

	delete mat;
	delete p;
	delete calc;
}

bool SPEA2::disperse(vector<Individual*> v, Individual* id) {

	for (size_t i = 0; i < v.size(); i++) {

		if (v.at(i)->getAptidao()[0] == id->getAptidao()[0] && v.at(i)->getAptidao()[1] == id->getAptidao()[1]) {
			return false;
		}
	}

	return true;
}

/*--------------------------------------- Communication methods -------------------------------------- */

/*Send a simple message with a integer 1 for n processes*/
void SPEA2::sendMessage(int numOfProcesses)
{
	double message = 1;

	for (int count = 1; count <= numOfProcesses; count++) {
		if (debug)
			cout << "Send ok to ---> " << count << endl;
		MPI_Send(&message, 1, MPI_DOUBLE, count, count, MPI_COMM_WORLD);
	}
}

/*Receive a double as a messagefrom any source*/
double SPEA2::receiveMessage()
{
	MPI_Request request;
	MPI_Status status;
	double message;
	int flag = -1;

	while (1) {
		if (flag != 0) {
			MPI_Irecv(&message, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &request);
			flag = 0;
		}
		MPI_Test(&request, &flag, &status);

		if (flag != 0) {
			if (debug)
				cout << "Message: " << message << " received from " << status.MPI_SOURCE << endl;
			return message;
		}
	}
}

void SPEA2::waitMessageForMe(int tag)
{
	MPI_Request request;
	MPI_Status status;
	double message;
	int flag = -1;

	while (1) {
		if (flag != 0) {
			MPI_Irecv(&message, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &request);
			flag = 0;
		}
		MPI_Test(&request, &flag, &status);

		if (flag != 0) {
			if (debug)
				cout << "Message to rank: " << tag << " received from rank " << status.MPI_SOURCE << endl;
			break;
		}
	}
}

void SPEA2::sendMatrix(int numOfProcesses, double * matrix, size_t numLines, size_t numCollums)
{

	for (int count = 1; count <= numOfProcesses; count++) {
		if (debug)
			cout << "Send matrix to ---> " << count << endl;
		MPI_Send(matrix, ((int)(numLines * numCollums)), MPI_DOUBLE, count, count, MPI_COMM_WORLD);
	}
}

double * SPEA2::waitMatrixForMe(int tag, size_t numLines, size_t numCollums)
{
	MatTransform* matT = new MatTransform();
	MPI_Request request;
	MPI_Status status;
	double* matrix = matT->allocMatrix(numLines, numCollums);
	int flag = -1;
	delete matT;

	while (1) {
		if (flag != 0) {
			MPI_Irecv(matrix, (numLines * numCollums), MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &request);
			flag = 0;
		}
		MPI_Test(&request, &flag, &status);

		if (flag != 0) {
			if (debug) {
				cout << "Matrix to rank: " << tag << " received from rank " << status.MPI_SOURCE << endl;
			}
			break;
		}
	}

	return matrix;
}

double * SPEA2::waitMatrixWithHypervolumeForMe(int tag, size_t numLines, size_t numCollums)
{
	MatTransform* matT = new MatTransform();
	MPI_Request request;
	MPI_Status status;
	double* matrix = matT->allocMatrixWithHypervolume(numLines, numCollums);
	int flag = -1;
	delete matT;

	while (1) {
		if (flag != 0) {
			MPI_Irecv(matrix, ((numLines * numCollums) + 1), MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &request);
			flag = 0;
		}
		MPI_Test(&request, &flag, &status);

		if (flag != 0) {
			if (debug)
				cout << "Matrix to rank: " << tag << " received from rank " << status.MPI_SOURCE << endl;

			break;
		}
	}

	return matrix;
}

void SPEA2::clearVector(vector<Individual*>& v)
{
	for (size_t i = 0; i < v.size(); i++) {
		delete v[i];
	}
	v.clear();
}

void SPEA2::copyVectorIndividuals(vector<Individual*>& one, vector<Individual*>& two)
{
	for (size_t i = 0; i < one.size(); i++) {

		Individual * id = new Individual(*one.at(i));
		two.push_back(id);
	}
}

void SPEA2::initBounds(double lowerBound, double upperBound)
{
	for (int i = 0; i < qtdGenes; i++) {
		limiteInferior.push_back(lowerBound);
		limiteSuperior.push_back(upperBound);
	}
}

void SPEA2::initBounds()
{
	limiteInferior = { 0, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5};
	limiteSuperior = { 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
}
