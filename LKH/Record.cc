#include "LK.h"

namespace myOptimizeName
{

  void LKSolver::RecordBetterTour()
  {
    LKModel *myModel = _myModel;
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();

    NodeClass *N = Data->FirstNode;

    for(int i = 1, k = 0; i <= Data->Dimension; ++i, N = N->Suc)
    {
      myModel->BetterTour[++k] = N->Id; 
      N->NextBestSuc = N->BestSuc;
      N->BestSuc = N->Suc;
    }
  }

  void LKSolver::RecordBestTour()
  {
    LKModel *myModel = _myModel;
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();

    for(int i = 1; i <= Data->Dimension; ++i)
    { myModel->BestTour[i] = myModel->BetterTour[i]; }
  }

  void LKSolver::StoreTour()
  {
    NodeClass *u;
    CandidateClass *Nt;

    while (Swaps > 0) {
      Swaps--;
      for (int i = 1; i <= 4; i++) {
	NodeClass *t = i == 1 ? SwapStack[Swaps].t1 :
	  i == 2 ? SwapStack[Swaps].t2 :
	  i == 3 ? SwapStack[Swaps].t3 : SwapStack[Swaps].t4;
	//Activate(t);
	if(!t->V)
	{
	  ActiveNodes.push(t);
	  t->V = 1;
	}
	t->OldPred = t->Pred;
	t->OldSuc = t->Suc;
	t->OldPredExcluded = t->OldSucExcluded = 0;
	t->Cost = LONG_MAX;
	for(size_t j = 0; j < t->CandidateSet.size(); ++j)
	{
	  Nt = &t->CandidateSet[j];
	  u = Nt->To;
	  if (u != t->Pred && u != t->Suc && Nt->Cost < t->Cost)
	    t->Cost = Nt->Cost;
	}
	
	//for (Nt = t->CandidateSet; u = Nt->To; Nt++)
	//if (u != t->Pred && u != t->Suc && Nt->Cost < t->Cost)
	//t->Cost = Nt->Cost;
      }
    }
  }

  void LKSolver::RestoreTour()
  {
    NodeClass *t1, *t2, *t3, *t4;

    /* Loop as long as the stack is not empty */
    while (Swaps > 0) {
      /* Undo topmost 2-opt move */
      Swaps--;
      t1 = SwapStack[Swaps].t1;
      t2 = SwapStack[Swaps].t2;
      t3 = SwapStack[Swaps].t3;
      t4 = SwapStack[Swaps].t4;
      Flip_SL(t3, t2, t1);
      Swaps--;
      /* Make edges (t1,t2) and (t2,t3) excludable again */
      t1->OldPredExcluded = t1->OldSucExcluded = 0;
      t2->OldPredExcluded = t2->OldSucExcluded = 0;
      t3->OldPredExcluded = t3->OldSucExcluded = 0;
      t4->OldPredExcluded = t4->OldSucExcluded = 0;
    }
  }

};
