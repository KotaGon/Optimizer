#include "LK.h"

namespace myOptimizeName
{
  void LKModel::importData()
  {
    cout << "importData... " << endl;

    LKDataModel *data = (LKDataModel*) DataModel; 
    LKSolver *solver = (LKSolver*) OptModel;

    solver->Runs = 100;
    solver->Seed = 1;
    solver->MaxTrials = 10;
    solver->MaxSwaps = -1;
    solver->MaxCandidates = 100;
    solver->AscentCandidates = 50;
    solver->Precision = 100;
    solver->Subgradient = 1;
    solver->Excess = 0.0;
    solver->RestrictedSearch = 1;
    solver->InitialStepSize = 1;
    solver->kopt = 5;

    data->readProblem();

    const int Dim = data->Dimension;
    BestTour = new int[Dim];
    BetterTour = new int[Dim];
    //SwapStack = ...

    /*
       assert(Rand = (int *) malloc((Dimension + 1) * sizeof(int)));
       if (Seed == 0) { Seed = 1; }
       srand(Seed);
       for (i = 1; i <= Dimension; i++)
       Rand[i] = rand();
     */

    srand(solver->Seed);
    solver->Swaps = 0;
    //MakeHeap(Dimension);
    if (solver->MaxCandidates > Dim)
      solver->MaxCandidates = Dim;
    if (solver->AscentCandidates > Dim)
      solver->AscentCandidates = Dim;
    if (solver->InitialPeriod == 0) 
    {
      solver->InitialPeriod = Dim / 2;
      if (solver->InitialPeriod < 1000)
	solver->InitialPeriod = 1000;
    }
    if (solver->Excess == 0.0)
      solver->Excess = 1.0 / Dim;
    if (solver->MaxTrials == 0)
      solver->MaxTrials = Dim;
    if (solver->MaxSwaps < 0)
      solver->MaxSwaps = Dim;

    solver->SwapStack = new SwapRecord[solver->MaxSwaps + 30];
    solver->FlipSegs = new SegmentClass[solver->kopt + 1];

    solver->MakeHeap(Dim);
    
    //NodeClass *(LKSolver::solver->BestMove = &(LKSolver::Best2OptMove);
    //NodeClass *(LKSolver::*BestMove)(NodeClass *, NodeClass*, int *, int *) = &LKSolver::Best2OptMove;

    NodeClass *Ni, *Nj, *FirstNode = data->FirstNode;

    data->CostMatrix = new int[Dim * (Dim - 1) / 2];
    Ni = FirstNode->Suc;

    do {
      Ni->C = &data->CostMatrix[(Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
      for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
	//Ni->C[Nj->Id] = Fixed(Ni, Nj) ? 0 : sqrt(sq(Ni->X - Nj->X) + sq(Ni->Y - Nj->Y));
	Ni->C[Nj->Id] = (int) (sqrt(sq(Ni->X - Nj->X) + sq(Ni->Y - Nj->Y)) + 0.5);
    }
    while ((Ni = Ni->Suc) != FirstNode);
    //WeightType = EXPLICIT;
    //c = 0;
    //C = WeightType == EXPLICIT ? C_EXPLICIT : C_FUNCTION;
    //D = WeightType == EXPLICIT ? D_EXPLICIT : D_FUNCTION;

    int i;
    SegmentClass *SPrev, *S;
    data->GroupSize = sqrt(1.0 * Dim);
    solver->Groups = 0;
    for (i = Dim, SPrev = 0; i > 0; i -= data->GroupSize, SPrev = S) 
    {
      S = new SegmentClass; 
      S->Rank = ++(solver->Groups);
      if (!SPrev)
	data->FirstSegment = S;
      else
	Link(SPrev, S);

    }
    Link(S, data->FirstSegment);

    return;
  }

  void LKModel::optimize()
  {
    cout << "optimize..."  << endl;

    int Succeses = 0;

    LKSolver *solver = (LKSolver*)OptModel;
    solver->CreateCandidateSet();

    if(solver->Norm != 0)
    {
      BestCost = DBL_MAX;
      WorstCost = -DBL_MAX; 
      Succeses = 0;
    }
    else 
    { 
      Succeses = 1;
      solver->Runs = 0;
      solver->RecordBetterTour();
      solver->RecordBestTour();
      BestCost = WorstCost = LowerBound;
    }

    if(Succeses)
      return;
   
    solver->MakeKoptFlip();
    solver->RunLoop();

    return;
  }

  void LKModel::exportData()
  {
    cout << "exportData..." << endl; 
  }

  void LKModel::free()
  {
    cout << "free..." << endl;
    DataModel->free();
    OptModel->free();
    delete DataModel;
    delete OptModel;
    delete [] BestTour;
    delete [] BetterTour;
  }


};
