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

	/*double f1 = i->getGenes()[0];
	double somatoria = 0.0;
	int n = i->getGenes().size();

	for (int j = 1; j < n; j++) {
		double x = i->getGenes()[j];
		somatoria += pow(x, 2) - (10.0 * cos(4 * PI * x));
	}
	double g = 91 + somatoria;
	double  h = 1 - std::sqrt(f1 / g);
	double f2 = g * h;*/

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

	//double a[] = { 0.01, 0.012, 0.004, 0.006, 0.004, 0.01 };
	//double b[] = { 2, 1.5, 1.8, 1.0, 1.8, 1.5 };
	//double c[] = { 10, 10, 20, 10, 20, 10 };

	//double alpha[] = { 0.00419, 0.00419, 0.00683, 0.00683, 0.00461, 0.00461 };
	//double beta[] = { 0.32767, 0.32767, -0.545514, -0.54551, -0.51116, -0.51116 };
	//double gamma[] = { 13.85932, 13.85932, 40.26690, 40.26690, 42.89553, 42.89553 };

	//double f1 = 0;
	//double f2 = 0;
	//double sumP = 0;

	//for (size_t j = 0; j < n; j++) {
	//	sumP += i->getGenes()[j];
	//	f1 += (a[j] * pow(i->getGenes()[j], 2)) + (b[j] * i->getGenes()[j]) + c[j];
	//	f2 += (alpha[j] * pow(i->getGenes()[j], 2)) + (beta[j] * i->getGenes()[j]) + gamma[j];
	//}

	////constraint
	//if (sumP < 700) {
	//f1 += 100000;
	//f2 += 100000;
	//}

	double f1 = 0;
	double f2 = 0;
	double sumP = 0;

	double a[] = { 0.0069 , 0.0069 , 0.02028, 0.00942, 0.0114, 0.01142, 0.00357, 0.00492, 0.00573, 0.00605,
						0.00515, 0.00569, 0.00421, 0.00752, 0.00708, 0.00708, 0.00313, 0.00313, 0.00313, 0.00313,
						0.00298 , 0.00298, 0.00284, 0.00284, 0.00277, 0.00277, 0.52124, 0.52124, 0.52124, 0.0114,
						0.0016, 0.0016, 0.0016, 0.0001, 0.0001, 0.0001, 0.0161, 0.0161, 0.0161, 0.00313 };

	double b[] = { 6.73 , 6.73 , 7.07, 8.18, 5.35, 8.05, 6.99, 6.6, 6.6, 12.9,
		12.9, 12.8, 12.5, 8.84, 9.15, 9.15, 7.97, 7.97, 7.97, 7.97,
		6.63 , 6.63, 6.66, 6.66, 7.1, 7.1, 3.33, 3.33, 3.33, 5.35,
		6.43, 6.43, 6.43, 8.95, 8.62, 8.62, 5.88, 5.88, 5.88, 7.97 };

	double c[] = { 94.705 , 94.705 , 309.54, 369.03, 148.89, 222.33, 278.71, 391.98, 455.76, 722.82,
		635.2, 654.69, 913.4, 1760.4, 1728.3, 1728.3, 647.85, 649.69, 647.83, 647.81,
		785.96 , 785.96, 794.53, 794.53, 801.32, 801.32, 1055.1, 1055.1, 1055.1, 148.89,
		222.92, 222.92, 222.92, 107.87, 116.58, 116.58, 307.45, 307.45, 307.45, 647.83 };


	double alpha[] = { 0.0057, 0.0046, 0.0025, 0.0028, 0.0058, 0.0053, 0.0052, 0.0056, 0.0057, 0.0052,
						0.0033, 0.0059, 0.0047, 0.0047, 0.004, 0.0056, 0.0059, 0.0043, 0.0051, 0.0049,
						0.0024, 0.004, 0.005, 0.0036, 0.0027, 0.0038, 0.0056, 0.006, 0.0025, 0.0024,
						0.0029, 0.0049, 0.0051, 0.0042, 0.005, 0.006, 0.0058, 0.0022, 0.0056, 0.0026 };

	double beta[] = { 0.033, 0.0458, 0.0469, -0.0446, 0.0008, 0.0481, 0.0167, 0.0478, 0.0499, 0.0411,
						-0.0553, 0.0281, 0.01, -0.0319, 0.0498, 0.046, -0.0208, -0.0417, -0.0034, 0.0463,
						0.0092, 0.0387, 0.0479, 0.0462, 0.0497, 0.0356, 0.0054, 0.0088, 0.0472, -0.0435,
						0.0491, -0.0328, 0.0311, -0.0313, 0.0069, -0.0009, 0.03, 0.0423, 0.0327, -0.0408 };

	double gamma[] = { 7.248, 19.834, 18.317, 19.22, 10.18, 14.774, 6.007, 17.934, 14.468, 17.984,
						11.002, 21.727, 16.742, 5.492, 17.754, 19.684, 13.608, 6.374, 17.277, 6.81,
						20.634, 11.574, 9.36, 19.848, 12.101, 18.162, 21.305, 18.734, 19.399, 14.765,
						5.914, 7.28, 7.546, 20.767, 22, 9.143, 7.102, 11.21, 11.206, 6.195 };

	for (size_t j = 0; j < n; j++) {
		sumP += i->getGenes()[j];
		f1 += (a[j] * pow(i->getGenes()[j], 2)) + (b[j] * i->getGenes()[j]) + c[j];
		f2 += (alpha[j] * pow(i->getGenes()[j], 2)) + (beta[j] * i->getGenes()[j]) + gamma[j];
	}

	if (sumP < 10500) {
		f1 += FLT_MAX;
		f2 += FLT_MAX;
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
