#include "Comparator.h"

Comparator::Comparator()
{
    //ctor
}

Comparator::~Comparator()
{
    //dtor
}

bool Comparator:: compare(float f1, float f2){

   return (f1 < f2);
}

bool Comparator:: compareMax(float f1, float f2){

   return (f1 > f2);
}

