#include "LK.h"

namespace myOptimizeName
{

  void LKDataModel::readProblem()
  {
    int *ID;
    long  *X, *Y;
    //const string filename = "./a280.tsp";
    //const string filename = "E3k.0.tsp";
    const string filename = "./problem/test_problem1000.tsp";
    //const string filename = "./problem/E10k.0";
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
	X = new long[this->Dimension+1];
	Y = new long[this->Dimension+1];
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
	  X[i] = longParse(elems[1]); 
	  Y[i] = longParse(elems[2]);
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
  void LKDataModel::createNodes(int *ID, long *X, long *Y)
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
    }
    Link(N, FirstNode);
  }
};
