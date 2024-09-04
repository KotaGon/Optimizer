#include "LK.h"

using namespace myOptimizeName;

int main()
{

   LKModel model;
   model.importData();
   model.optimize();
   model.exportData();
   model.free();

   return 0;
}
