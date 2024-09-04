#include "LK.h"

namespace myOptimizeName
{

  void LKDataModel::readProblem()
  {
    int *ID, *X, *Y;
    //const string filename = "./a280.tsp";
    //const string filename = "E3k.0.tsp";
    const string filename = "./problem/test_problem1500.tsp";
    string line;
    vector<string> elems;
    ifstream file(filename);

    cout << "filename = " << filename << endl;

    while(getline(file, line))
    {
      if(line == "")
      { continue; }
      elems = split(line);
      if(elems[0] == "DIMENSION")
      {
       	this->Dimension = intParse(elems[1]); 
	cout << this->Dimension << endl;
	ID = new int[this->Dimension+1];
	X = new int[this->Dimension+1];
	Y = new int[this->Dimension+1];
      }

      if(elems[0] == "NODE_COORD_SECTION")
      {
	int i = 1;
	while(getline(file, line))
	{
	  if(line == "")
	  { continue; }
	  if(line == "EOF")
	  { break; }
	  
	  elems = split(line);

	  ID[i] = intParse(elems[0]);
	  X[i] = intParse(elems[1]); 
	  Y[i] = intParse(elems[2]);
	  ++i;
	}
      }
    }

    createNodes(ID, X, Y);

    file.close();

    delete [] ID;
    delete [] X;
    delete [] Y;
  }
  void LKDataModel::createNodes(int *ID, int *X, int *Y)
  {
    cout << "create nodes..." << this->Dimension <<endl;
    NodeClass *Prev, *N;
    NodeSet = new NodeClass[Dimension+1];

    for(int i = 1; i <= Dimension; ++i, Prev = N)
    {
      N = &NodeSet[i]; 
      N->V = 1;
      if(i == 1)
      { FirstNode = N; }
      else 
      { Link(Prev, N); }
      N->Id = ID[i];
      N->X = (double) X[i];
      N->Y = (double) Y[i];
      //N->Z = (double) N->Z;
    }
    Link(N, FirstNode);
#if 0
    N = FirstNode;
    do
        if (!N->V && N->Id <= Dimension)
            break;
    while ((N = N->Suc) != FirstNode);
#endif
  }
};
