#ifndef COMPARATOR_H
#define COMPARATOR_H

#include<algorithm>

class Comparator
{
    public:
        Comparator();
        virtual ~Comparator();
        bool compare(float f1, float f2);
        bool compareMax(float f1, float f2);
    protected:

    private:
};

#endif // COMPARATOR_H
