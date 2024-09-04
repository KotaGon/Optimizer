#include <iostream>
#include "SA.h"

using namespace myOptimizeName;

int main()
{

  SAModel model;

  model.importData();
  model.optimize();

  model.free();

  return 0;
}
