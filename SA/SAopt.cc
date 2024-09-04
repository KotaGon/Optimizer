#include "SA.h"
#include <sstream>
#include <iomanip>

namespace myOptimizeName
{

  void SADataModel::readProblem()
  {
    int *ID;
    long *X, *Y;
//    const string filename = "./test_problem5000.tsp";
    const string filename = "./E10k.0";
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
  void SADataModel::createNodes(int *ID, long *X, long *Y)
  {
    cout << "create nodes..." << this->Dimension <<endl;
    NodeClass *Pred, *N;
    NodeSet = new NodeClass[Dimension+1];

    for(int i = 1; i <= Dimension; ++i, Pred = N)
    {
      N = &NodeSet[i]; 
      if(i == 1)
      { FirstNode = N; }
      else 
      { Link(Pred, N); }
      N->Id = ID[i];
      N->X = (long) X[i];
      N->Y = (long) Y[i];
    }
    Link(N, FirstNode);
  }

  void SAModel::importData()
  {
    cout << "importData... " << endl;

    SADataModel *data = (SADataModel*) DataModel; 
    SASolver *solver = (SASolver*) OptModel;

    data->readProblem();

    const int Dim = data->Dimension;
    BestTour = new int[Dim];
    BetterTour = new int[Dim];

    NodeClass *Ni, *Nj, *FirstNode = data->FirstNode;

    data->CostMatrix = new long[Dim * (Dim - 1) / 2];
    Ni = FirstNode->Suc;

    do {
      Ni->C = &data->CostMatrix[(Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
      for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
	Ni->C[Nj->Id] = (long) (sqrt(sq(Ni->X - Nj->X) + sq(Ni->Y - Nj->Y)) + 0.5);
    }
    while ((Ni = Ni->Suc) != FirstNode);

    return;
  }

  void SAModel::optimize()
  {
    cout << "optimize..."  << endl;

    SASolver *Solver = (SASolver*)OptModel;
    Solver->SimulatedAnnealing();

    return;
  }

  void SAModel::exportData()
  {
    cout << "exportData..." << endl; 
  }

  void SAModel::free()
  {
    cout << "free..." << endl;
    DataModel->free();
    OptModel->free();
    delete DataModel;
    delete OptModel;
    delete [] BestTour;
    delete [] BetterTour;
  }

  void SASolver::CreateInitialTour()
  {

     cout << "Create Initial Tour" << endl;

     SADataModel *Data = (SADataModel*) _myModel->getDataModel();
     NodeClass *N = Data->FirstNode = &Data->NodeSet[1 + rand() % Data->Dimension], *Next;

     do
       N->Activate = 0;
     while((N = N->Suc) != Data->FirstNode);

     Data->FirstNode->Activate = 1;

     int c =0;

     while(N->Suc != Data->FirstNode)
     {
       NodeClass *NextN = Next = N->Suc;
       long a = LONG_MAX;

       while(NextN->Suc != Data->FirstNode)
       {
	 //if(NextN->Activate)
	    //continue; 

	 if(Data->C(N, NextN) < a)
	 {
	   a = Data->C(N, NextN);
	   Next = NextN;
	 }
	 NextN = NextN->Suc;
       }

       Follow(Next, N);

       N = Next;
       N->Activate = 1;
     }

     return;

  }

  void SASolver::SimulatedAnnealing()
  {

    SADataModel *Data = (SADataModel*) _myModel->getDataModel();

    NodeClass *FirstNode, *t1, *t2;
    long Gain, Cost=0, InitialCost;
    int nOut = 0;
    bool update = false;

    CreateInitialTour();

    t1 = FirstNode = Data->FirstNode;
    do
    {
      t2 = SUC(t1);
      Cost += Data->C(t1, t2); 
      t1->Activate = 1;
      ActiveNodes.push(t1);
    }
    while((t1 = t2) != FirstNode);

    cout << "InitialCost = " << (InitialCost=Cost) << endl;

    for(int Iterator = 0; Iterator < MaxIterator; ++Iterator)
    {

      while(ActiveNodes.size() > 0)
      {
	t1 = ActiveNodes.front(); ActiveNodes.pop();
	t2 = SUC(t1);

	t1->Activate = 0;

	Best2OptMove(t1, t2, &Gain, &(update=false));

	//if(Gain > 0)
	if(update)
	{
	  Cost -= Gain;

	  if(Iterator % 100 == 0)
	    out_tour("sa", nOut++);

	  if(!t1->Activate)
	  {
	    ActiveNodes.push(t1); 
	    t1->Activate = 1;
	  }

	  break;
	}

      }

      if(ActiveNodes.size() == 0)
      {
	t1 = FirstNode;
	do
	{
	  t2 = SUC(t1);
	  t1->Activate = 1;
	  ActiveNodes.push(t1);
	}
	while((t1 = t2) != FirstNode);
      }
      cout << (Gain > 0 ? "*" : "") << "[Iterate" << Iterator << "] " << "Cost = " << Cost << ", temp = " << temp << endl;

      //temp *= 0.99;
      //if(temp <= 1.0e-6)
	//temp = 1.0e-6;
      
      temp = Temp(Iterator);
      //temp = 1000 + (1 - 1000) * (double) Iterator / (double) (MaxIterator);
      //if(temp <= 0.0) temp = 1.0e-6;
      //temp = InitialCost/1000.0 - pow(InitialCost/1000.0, (double) Iterator / MaxIterator);

    }

    myGetChar();

    return ;

  }

  double SASolver::Temp(const int Iterator)
  {
     const int StartTemp = 10;
     const int EndTemp = 1;

     //double t = StartTemp * pow(0.999, Iterator);
     double t = StartTemp + (EndTemp - StartTemp) * (double) Iterator / (double) (MaxIterator);
     
     if(t <= 1.0e-6)
       t = 1.0e-6;

     return t; 

  }

  void SASolver::Best2OptMove(NodeClass *t1, NodeClass *t2, long *Gain, bool *update)
  {

    SADataModel *Data = (SADataModel*) _myModel->getDataModel();
    NodeClass *FirstNode = Data->FirstNode;

    NodeClass *t3, *t4;

    t4 = SUC(t2); 
    t3 = SUC(t4);

    long G0 = Data->C(t1, t2);

    do
    {
      t3 = SUC(t4);

      if(t1 == t3 || t1 == t4 || t2 == t3 || t3 == t4)
	continue;

      long G1 = G0 - Data->C(t2, t3);
      long G2 = G1 + Data->C(t3, t4);

      *Gain = G2 - Data->C(t1, t4);

      double rnd = (double) rand() / RAND_MAX;
      //cout << rnd << " " << exp(*Gain / temp) << endl;

      if(*Gain > 0 || (*Gain < 0 && rnd < exp( *Gain / temp)))
      //if(*Gain > 0)
      {
	//cout << "Gain " << *Gain << endl;
	//Swap 
	t2->Pred = NULL;
	NodeClass *u1, *u2 = t4;

	while((u1 = u2))
	{
	  u2 = u1->Pred;
	  u1->Pred = u1->Suc;
	  u1->Suc = u2;
	}

	t1->Suc = t4;
	t2->Suc = t3;
	t3->Pred = t2;
	t4->Pred = t1;

	if(!t1->Activate)
	{ 
	  t1->Activate = 1;
	  ActiveNodes.push(t1);
	}
	if(!t2->Activate)
	{ 
	  t2->Activate = 1;
	  ActiveNodes.push(t2);
	}
	if(!t3->Activate)
	{ 
	  t3->Activate = 1;
	  ActiveNodes.push(t3);
	}
	if(!t4->Activate)
	{ 
	  t4->Activate = 1;
	  ActiveNodes.push(t4);
	}

	*update = true;

	return;

      }

    }
    while((t4 = t3) != FirstNode);


    return;  
  }

  void SASolver::out_tour(const string &filename, const int out)
  {
    SADataModel *Data = (SADataModel*) _myModel->getDataModel();
    NodeClass *t = Data->FirstNode;

    std::ostringstream ss;
    ss << std::setw(4) << std::setfill('0') << out;
    std::string tour_number(ss.str());

    ofstream tour("output/"+filename + tour_number + ".gnu");
    tour << "#!/bin/gnuplot" << endl;
    tour << "set terminal png" << endl;
    tour << "set output 'png/tour" << tour_number << ".png'" << endl; 
    tour << "plot '-' w l," << endl;

    int Cost = 0;

    do
    {
      tour << t->X << " " << t->Y << endl;
      tour << SUC(t)->X << " " << SUC(t)->Y << endl;

      tour << endl;
      tour << endl; 
    }
    while((t = SUC(t)) != Data->FirstNode);
    tour << "e" << endl;

    tour.close();
  }
};
