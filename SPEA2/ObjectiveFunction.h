#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include"Individual.h"

using namespace std;

class ObjectiveFunction
{
    public:
        ObjectiveFunction();
        virtual ~ObjectiveFunction();
        virtual vector<double> evaluateIndividual(Individual* i) = 0;
        virtual bool isMinimization(int function) = 0;
        virtual int numberOfFunctions() = 0;

    protected:

    private:
};

#endif // OBJECTIVEFUNCTION_H
