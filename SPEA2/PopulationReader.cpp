#include "PopulationReader.h"



PopulationReader::PopulationReader()
{
}


PopulationReader::~PopulationReader()
{
}
/*Pega as aptidões em txt*/
vector<Individual*> PopulationReader::loadFromArquive(string path)
{
	//"../logs/individuals.txt"
	string line;
	ifstream file(path);
	vector<Individual*> individuals;

	if (file.is_open()) {		
		while (!file.eof()) {

			getline(file, line);
			/*cout << line << endl;*/
			vector<string> data = this->split(line, ',');
			vector<double> aptidoes;
			aptidoes.push_back(std::stod(data.at(0)));
			aptidoes.push_back(std::stod(data.at(1)));
			
			Individual* id = new Individual();
			id->setAptidao(aptidoes);
			individuals.push_back(id);
		}

		file.close();
	}
	else {
		cout << "Erro ao tentar abrir o arquivo!" << endl;
	}

	return individuals;
}

vector<Individual*> PopulationReader::loadPopulation(int rank)
{
	string line;
	ifstream file("logs/population_rank_" + to_string(rank) + ".bin");
	vector<Individual*> individuals;

	if (file.is_open()) {

		while (!file.eof()) {

			getline(file, line);			
			vector<string> data = this->split(line, ' ');
			vector<double> genes;
			for (int j = 0; j < data.size(); j++ )
				genes.push_back(std::stod(data.at(j)));			

			Individual* id = new Individual();
			id->setGenes(genes);
			id->setQtdGenes(genes.size());
			id->setAptidao(this->evaluateIndividual(id));
			individuals.push_back(id);			
		}		
		file.close();
	}
	else {
		cout << "Erro ao tentar abrir o arquivo!" << endl;
	}

	return individuals;
}

vector<Individual*> PopulationReader::getMasterPopulation(int rank)
{
	string line;
	ifstream file("logs/population_master_rank_" + to_string(rank) + ".bin");
	vector<Individual*> individuals;

	if (file.is_open()) {

		while (!file.eof()) {

			getline(file, line);
			vector<string> data = this->split(line, ' ');
			vector<double> genes;
			for (int j = 0; j < data.size(); j++)
				genes.push_back(std::stod(data.at(j)));

			Individual* id = new Individual(genes.size());
			id->setGenes(genes);
			id->setQtdGenes(genes.size());
			id->setAptidao(this->evaluateIndividual(id));
			individuals.push_back(id);
		}
		file.close();
	}
	else {
		cout << "Erro ao tentar abrir o arquivo!" << endl;
	}

	return individuals;
}

/*
	http://code.runnable.com/VHb0hWMZp-ws1gAr/splitting-a-string-into-a-vector-for-c%2B%2B
*/
vector<string> PopulationReader::split(string text, char delimiter)
{
	vector<string> internal;
	stringstream ss(text); // Turn the string into a stream.
	string tok;

	while (getline(ss, tok, delimiter)) {
		internal.push_back(tok);
	}

	return internal;
}

void PopulationReader::savePopulation(vector<Individual*> pop, int rank)
{
	ofstream file;
	string fileName = "logs/population_rank_" + to_string(rank) + ".bin";	
	file.open(fileName, std::ios::out | std::ios::binary);

	if (file.is_open()) {
		for (int i = 0; i < pop.size(); i++) {
			vector<double> id = pop.at(i)->getGenes();

			for (int j = 0; j < id.size(); j++)
			{
				file << to_string(id.at(j)) + " ";
			}
			if (i < pop.size() - 1)
				file << std::endl;
		}
	}
	else
	{
		cout << "no!" << endl;
	}

	file.close();
}

void PopulationReader::savePopulationForAll(vector<Individual*> pop, int numberOfSlaves)
{
	for (int i = 1; i <= numberOfSlaves; i++) {

		ofstream file;
		string fileName = "logs/population_master_rank_" + to_string(i) + ".bin";
		file.open(fileName, std::ios::out | std::ios::binary);

		if (file.is_open()) {
			for (int i = 0; i < pop.size(); i++) {
				
				vector<double> id = pop.at(i)->getGenes();//get genes

				for (int j = 0; j < id.size(); j++)
				{
					file << to_string(id.at(j)) + " "; // save each gene
				}

				if (i < pop.size() - 1)
					file << std::endl; //fim do arquivo
			}
		}
		else
		{
			cout << "no!" << endl;
		}

		file.close();
	}
}

vector<double> PopulationReader:: evaluateIndividual(Individual* individual) {

	ZDT* zdt1 = new ZDT();
	return zdt1->evaluateIndividual(individual);
	delete zdt1;
}
