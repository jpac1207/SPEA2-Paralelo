#include "ZDT.h"


ZDT::ZDT()
{
	//ctor
}

ZDT::~ZDT()
{
	//dtor
}

vector<double> ZDT::evaluateIndividual(Individual* i) {

	/* ------------------------ zdt1 ---------------------- */
	/*double f1 = i->getGenes()[0];

	size_t m = i->getGenes().size();

	double somatoria = 0;

	for (int count = 1; count < m; count++)
		somatoria += i->getGenes()[count];

	double g = 1 + (9 / (m - 1) * (somatoria));
	double h = (1 - sqrt((f1 / g)));

	double f2 = g * h;*/

	/* ------------------------ zdt2 ---------------------- */
	/*double f1 = i->getGenes()[0];
	int n = i->getGenes().size();
	double somatoria = 0;

	for (int j = 1; j < i->getGenes().size(); j++) {
		somatoria += i->getGenes()[j];
	}
	double g = 1 + ((9 / (n - 1)) * (somatoria));
	double f2 = g * (1 - pow((f1 / g), 2));*/


	/* ------------------------ zdt3 ---------------------- */
	/*double f1 = i->getGenes()[0];
	double somatoria = 0;
	int n = i->getGenes().size();

	for (int j = 1; j < n; j++) {

		somatoria += i->getGenes()[j];
	}
	double g = 1 + ((9 / (n - 1)) * (somatoria));
	double f2 = g * (1 - sqrt((f1 / g)) - (f1 / g) *  sin((10 * PI * f1)));*/


	/* ------------------------ zdt4 ---------------------- */

	//double f1 = i->getGenes()[0];
	//double somatoria = 0.0;
	//int n = i->getGenes().size();

	//for (int j = 1; j < n; j++) {
	//	double x = i->getGenes()[j];
	//	somatoria += pow(x, 2) - (10.0 * cos(4 * PI * x));
	//}
	//double g = 91 + somatoria;
	//double  h = 1 - std::sqrt(f1 / g);
	//double f2 = g * h;

	/* ------------------------ zdt6 ---------------------- */

	/*double f1 = (1 - (std::exp(-4.0 * i->getGenes()[0]) * std::pow(std::sin(6 * PI * i->getGenes()[0]), 6)));
	double somatoria = 0.0;
	int n = i->getGenes().size();

	for (int j = 1; j < n; j++) {
		somatoria += i->getGenes()[j];
	}

	double g = 1 + 9 * std::pow(somatoria / (n - 1), 0.25);

	double f2 =  g * (1 - (pow((f1 / g), 2)));*/

	/* ------------------------ Schaffer 1 ---------------------- */
	/*double f1 = (pow(i->getGenes()[0], 2));
	double f2  = (pow(i->getGenes()[0] - 2, 2));*/


	/* ------------------------ Schaffer 2 ---------------------- */
	/*double f1 = 0.0;

	if (i->getGenes()[0] <= 1.0)
		f1 = i->getGenes()[0] * (-1.0);
	else if (i->getGenes()[0] > 1.0 && i->getGenes()[0] <= 3.0)
		f1 = i->getGenes()[0] - 2.0;
	else if (i->getGenes()[0] > 3.0 && i->getGenes()[0] <= 4.0)
		f1 = 4.0 - i->getGenes()[0];
	else if (i->getGenes()[0] > 4.0)
		f1 = i->getGenes()[0] -4.0;

	double f2 = pow(i->getGenes()[0] - 5.0, 2);*/

	/* ------------------------ Fonseca & Flemming ---------------------- */
	/*double n = i->getGenes().size();
	double somatorio = 0;

	for (int j = 0; j < n; j++) {

		somatorio += pow(i->getGenes()[j] - (1 / sqrt(n)), 2);
	}

	somatorio *= (-1);

	double f1 = 1 - exp(somatorio);

	somatorio = 0;
	for (int j = 0; j < n; j++) {

		somatorio += pow(i->getGenes()[j] + (1 / sqrt(n)), 2);
	}
	somatorio *= (-1);
	double f2 = 1 - exp(somatorio);*/

	/* ------------------------ Environmental Economic Dispatch ---------------------- */
	size_t n = i->getGenes().size();
	double a[] = { 0.01, 0.012, 0.004, 0.006, 0.004, 0.01 };
	double b[] = { 2, 1.5, 1.8, 1.0, 1.8, 1.5 };
	double c[] = { 10, 10, 20, 10, 20, 10 };

	double alpha[] = { 0.00419, 0.00419, 0.00683, 0.00683, 0.00461, 0.00461 };
	double beta[] = { 0.32767, 0.32767, -0.545514, -0.54551, -0.51116, -0.51116 };
	double gamma[] = { 13.85932, 13.85932, 40.26690, 40.26690, 42.89553, 42.89553 };

	double f1 = 0;
	double f2 = 0;
	double sumP = 0;

	for (size_t j = 0; j < n; j++) {
		sumP += i->getGenes()[j];
		f1 += (a[j] * pow(i->getGenes()[j], 2)) + (b[j] * i->getGenes()[j]) + c[j];
		f2 += (alpha[j] * pow(i->getGenes()[j], 2)) + (beta[j] * i->getGenes()[j]) + gamma[j];
	}

	//constraint
	if (sumP < 500) {
		f1 += 100000;
		f2 += 100000;
	}	
	/* ------------------------ to be returned (f1 and f2 must be created and filled) ---------------------- */
	vector<double> aptidoes;
	aptidoes.push_back(f1);
	aptidoes.push_back(f2);
	return aptidoes;
}

bool ZDT::isMinimization(int function) {

	return true;
}
int ZDT::numberOfFunctions() {
	return 2;
}
