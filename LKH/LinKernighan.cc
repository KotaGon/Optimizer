#include "LK.h"

namespace myOptimizeName
{
  void LKSolver::RunLoop()
  {
    LKModel *myModel = _myModel;
    double Cost = 0.0;
    double *BestCost = &myModel->BestCost, *WorstCost = &myModel->WorstCost;
    double Last_time = (double) clock()/CLOCKS_PER_SEC, Now;
   
    OutTour = 0;

    for(int Run = 0; Run <= Runs; ++Run)
    {
      Cost = FindTour();

      if(Cost < *BestCost)
      {
        RecordBestTour();
	*BestCost = Cost;
      }
      if(Cost > *WorstCost)
	*WorstCost = Cost;

      Now = (double) clock()/CLOCKS_PER_SEC - Last_time;
      cout << "Run[" << Run << "], Cost = " << Cost << ", BestCost = " << *BestCost << ", Time = " << Now << endl;
    }

    return ;
  }

  double LKSolver::FindTour()
  {
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    double Cost = 0.0, BetterCost = DBL_MAX;
    double Last_time = (double) clock()/CLOCKS_PER_SEC, Now;
    NodeClass *t;

    t = Data->FirstNode;
    do 
     t->OldPred = t->OldSuc = t->NextBestSuc = t->BestSuc = 0; 
    while((t = t->Suc) != Data->FirstNode);

    Hash.clear();

    for(int Trial = 1; Trial <= MaxTrials; ++Trial)
    {
      CreateInitialTour(); 
      Cost = LinKernighan();

      Now = (double) clock() / CLOCKS_PER_SEC - Last_time;

      if(Data->FirstNode->BestSuc)
      {
        t = Data->FirstNode;
	while((t = t->Next = t->BestSuc) != Data->FirstNode);
	Cost = MergeWithTour();
      }

      if(Cost < BetterCost)
      {

        BetterCost = Cost;
	RecordBetterTour();
	AdjustCandidateSets();
	Hash.clear();
	Hash[(long) Cost*Precision] = Hash.size();
	cout << "*" << Trial << ":Cost = " << Cost << ", Time = " << Now << "[sec] " << endl;

      }
      else 
	cout << " " << Trial << ":Cost = " << Cost << ", Time = " << Now << "[sec] " << endl;
    }

    ResetCandidateSets();
    return BetterCost;

  }

  void LKSolver::CreateInitialTour()
  {
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
   
    NodeClass *N, *NextN;
    //int i = 1 + rand() % Data->Dimension;
    //N = Data->FirstNode = &Data->NodeSet[1 + rand() % Data->Dimension];
    //int i = 1 + xor128() % Data->Dimension;
    N = Data->FirstNode = &Data->NodeSet[1 + xor128() % Data->Dimension];
     
    do 
    {
      N->V = 0;
    }
    while((N = N->Suc) != Data->FirstNode);

    Data->FirstNode->V = 1;
    N = Data->FirstNode;

    vector<NodeClass*> nextNodes; 
    nextNodes.reserve(Data->Dimension);

    while(N->Suc != Data->FirstNode)
    {
      NextN = NULL;
      nextNodes.clear();

      for(size_t i = 0; i < N->CandidateSet.size(); ++i)
      {
	CandidateClass *NN = &N->CandidateSet[i];
	NextN = NN->To;

	//if(!NextN->V && NN->Alpha == 0 && (InBestTour(N, NextN) || InNextBestTour(N, NextN)))
	if(!NextN->V && NN->Alpha == 0 && (InBestTour(N, NextN)))
	{
	  nextNodes.push_back(NextN); 
	}
      }

      if(nextNodes.size() == 0)
      {
        for(size_t i = 0; i < N->CandidateSet.size(); ++i)
	{
	  CandidateClass *NN = &N->CandidateSet[i];
	  NextN = NN->To;
	  if(!NextN->V)
	    nextNodes.push_back(NextN);
	} 
      }

      if(nextNodes.size() == 0)
      {
        NextN = N->Suc;	
	//do
	//  nextNodes.push_back(NextN);
	//while((NextN = NextN->Suc) != Data->FirstNode);
	//while ( NextN->Suc != Data->FirstNode)
	  //NextN = NextN->Suc;
      }
      else
      {
        //size_t pos = rand() % nextNodes.size();
	size_t pos = xor128() % nextNodes.size();
	NextN = nextNodes[pos];
      }
      
      Follow(NextN, N);
      N = NextN;
      N->V = 1;
    }

    return;

  }

  double LKSolver::LinKernighan()
  {
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    NodeClass *t1, *t2, *SUCt1;
    SegmentClass *S;

    Reversed = 0;
    S = Data->FirstSegment;
    int i = 0;

    do
    {
      S->Reversed = S->Size = 0;
      S->Rank = ++i;
      S->First = S->Last = NULL; 
    }
    while((S = S->Suc) != Data->FirstSegment);

    Swaps = i = 0;
    double Cost = 0.0;
    std::queue<NodeClass*>().swap(ActiveNodes); //ActiveNodes size is 0
   
    t1 = Data->FirstNode;
    do
    {
      t2 = t1->OldSuc = t1->Next = t1->Suc; 
      t1->OldPred = t1->Pred;
      t1->Rank = ++i;
      Cost += Data->C(t1, t2) - t1->Pi - t2->Pi;

      t1->Cost = LONG_MAX;
      for(size_t j = 0; j < t1->CandidateSet.size(); ++j)
      {
        CandidateClass *Nt1 = &t1->CandidateSet[j];
	t2 = Nt1->To;
	if(t2 != t1->Pred && t2 != t1->Suc && Nt1->Cost < t1->Cost)
	  t1->Cost = Nt1->Cost;
      }

      t1->Parent = S;
      ++(S->Size);

      if(S->Size == 1)
        S->First = t1;
      S->Last = t1;
      if(S->Size == Data->GroupSize)
	S = S->Suc;

      t1->OldPredExcluded = t1->OldSucExcluded = 0;
      t1->Next = NULL;
      ActiveNodes.push(t1);
      t1->V = 1; //Active!
      t1->PreSide = t1->PostSide = NULL;
      t1->SelectCnt = 0;

    }
    while((t1 = t1->Suc) != Data->FirstNode);

    long Gain = 0;

    do
    {
      while(ActiveNodes.size() > 0)
      {
	t1 = ActiveNodes.front();
	ActiveNodes.pop();
	t1->V = 0;

	SUCt1 = SUC(t1);

	for(int X2 = 1; X2 <= 2; ++X2)
	{
	  t2 = X2 == 1 ? PRED(t1) : SUCt1;

	  long G0 = Data->C(t1, t2);
	  while((t2 = BestKOptMove(t1, t2, &G0, &Gain)))
	  {
	    if(Gain > 0)
	    {
	      Cost -= Gain;
	      StoreTour();
	      if(!t1->V)
	      { 
		ActiveNodes.push(t1);
		t1->V = 1;
	      }

	      if(Hash.find((long) Cost) != Hash.end())
		goto End_LinKernighan;
	      goto Next_t1;
	    }
	  }
	  RestoreTour();
	}
Next_t1:;
      }
      if(Hash.find((long) Cost) != Hash.end())
	goto End_LinKernighan;

      Hash[(long) Cost] = Hash.size();

      Gain = Gain23();
      if(Gain > 0)
      {
	cout << "  ==> By Gain23 " << Cost / Precision << " " << (Cost-Gain) / Precision << " " << Gain / Precision << endl;
	Cost -= Gain;
	StoreTour();

	if(Hash.find((long) Cost) != Hash.end())
	  goto End_LinKernighan;

      }

    }
    while(Gain > 0);

End_LinKernighan:;
    //cout << "output => tour" << OutTour << ".gnu" << endl;
    //out_tour("tour", OutTour++);
    
    NormalizeNodeList();
    
    return Cost / Precision;
  }

  void LKSolver::NormalizeNodeList()
  {
    NodeClass *t1, *t2;
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    t1 = Data->FirstNode;
    do 
    {
      t2 = SUC(t1);
      t1->Pred = PRED(t1);
      t1->Suc = t2;
    }
    while ((t1 = t2) != Data->FirstNode);

    return;
  }

};
