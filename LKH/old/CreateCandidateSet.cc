#include "LK.h"

namespace myOptimizeName
{

  static long HeapCount;   

  void LKSolver::MakeHeap(const long Size)
  {
    (Heap = (NodeClass **) malloc((Size + 1) * sizeof(NodeClass *)));
    HeapCount = 0;
  }

  void LKSolver::SiftUp(NodeClass * N)
  {
    long Loc = N->Loc, P = Loc / 2;

    while (P && N->Rank < Heap[P]->Rank) {
      Heap[Loc] = Heap[P];
      Heap[Loc]->Loc = Loc;
      Loc = P;
      P /= 2;
    }
    Heap[Loc] = N;
    Heap[Loc]->Loc = Loc;
  }

  NodeClass *LKSolver::DeleteMin()
  {
    NodeClass *Remove, *Item;
    long Ch, Loc;

    if (!HeapCount)
      return 0;
    Remove = Heap[1];
    Item = Heap[HeapCount];
    HeapCount--;
    Loc = 1;
    Ch = 2 * Loc;

    while (Ch <= HeapCount) {
      if (Ch < HeapCount && Heap[Ch + 1]->Rank < Heap[Ch]->Rank)
	Ch++;
      if (Heap[Ch]->Rank >= Item->Rank)
	break;
      Heap[Loc] = Heap[Ch];
      Heap[Loc]->Loc = Loc;
      Loc = Ch;
      Ch *= 2;
    }
    Heap[Loc] = Item;
    Item->Loc = Loc;
    Remove->Loc = 0;
    return Remove;
  }

  void LKSolver::Insert(NodeClass * N)
  {
    long Ch, P;

    Ch = ++HeapCount;
    P = Ch / 2;
    while (P && N->Rank < Heap[P]->Rank) {
      Heap[Ch] = Heap[P];
      Heap[Ch]->Loc = Ch;
      Ch = P;
      P /= 2;
    }
    Heap[Ch] = N;
    N->Loc = Ch;
  }
  void LKSolver::CreateCandidateSet()
  {
    LKModel *myModel = _myModel;
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();

    double Cost;
    NodeClass *Na, *Nb;

    if(1)
    { 
      Na = Data->FirstNode;
      do
      {
	for(int i = 1; i < Na->Id; ++i)
	{ Na->C[i] *= Precision; } 
      }
      while((Na = Na->Suc) != Data->FirstNode);
    }

    Cost = Ascent();
    myModel->LowerBound = Cost / Precision; 

    cout << "Lower Bound = " << myModel->LowerBound << endl;

    if(myModel->Optimum != -DBL_MAX && myModel->Optimum != 0)
    { cout << "Gap = " << 100 * (myModel->Optimum - myModel->LowerBound) / myModel->Optimum; }
    
    GenerateCandidates(MaxCandidates, fabs(Excess * Cost), 0);

    if (1) 
    {
      Na = Data->FirstNode;
      do {
	Nb = Na;
	while ((Nb = Nb->Suc) != Data->FirstNode) {
	  if (Na->Id > Nb->Id)
	    Na->C[Nb->Id] += Na->Pi + Nb->Pi;
	  else
	    Nb->C[Na->Id] += Na->Pi + Nb->Pi;
	}
      }
      while ((Na = Na->Suc) != Data->FirstNode);
    }

    return;

  }

  double LKSolver::Ascent()
  {
    //cout << "ascent" << endl;

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();

    NodeClass *t;
    double BestW, W, W0;
    int T, Period;
    int InitialPhase;

Start:
    t = Data->FirstNode;
    do 
    	t->BestPi = t->Pi = t->V = 0;
    while((t = t->Suc) != Data->FirstNode);

    W = Minimum1TreeCost();

    //if (!Subgradient || !Norm || W / Precision == Optimum)    
    if(!Subgradient || !Norm)
      return W;

    GenerateCandidates(AscentCandidates, INT_MAX, 1);
    
    t = Data->FirstNode;
    do 
      t->LastV = t->V;
    while((t = t->Suc) != Data->FirstNode);

    BestW = W0 = W;
    InitialPhase = 1;

    for(Period = InitialPeriod, T = InitialStepSize * Precision;
	Period > 0 && T > 0 && Norm != 0; Period /= 2, T /= 2)
    {

      for(int P = 1; T && P <= Period && Norm > 0; ++P)
      {
	t = Data->FirstNode;
	do
	{
	  if(t->V != 0)
	    t->Pi += T * (7 * t->V + 3 * t->LastV) / 10;
	  t->LastV = t->V;
	}	  
	while((t = t->Suc) != Data->FirstNode);

	W = Minimum1TreeCost();

	t = Data->FirstNode;

	if(W > BestW)
	{
          if( W > INT_MAX && W > 2 * W0)
	  { 
	    T = 0;
	    break;
	  }	    
	  BestW = W;
	  cout << "* T = " << T 
	       << ", Period = " << Period 
	       << ", P = " << P 
	       << ", BestW = " << BestW 
	       << ", Norm = " << Norm << endl;
	  t = Data->FirstNode;
	  do 
	    t->BestPi = t->Pi;
	  while((t = t->Suc) != Data->FirstNode);
	  if(InitialPhase)
	    T *= 2;

	  if(P == Period && (Period *= 2) > InitialPeriod)
	    Period = InitialPeriod;
	}
	else if(InitialPhase && P > Period / 2)
	{
	  InitialPhase = 0;
	  P = 0;
	  T = 3 * T / 4; 

	  /*
	  cout << "  T = " << T 
	       << ", Period = " << Period 
	       << ", P = " << P 
	       << ", BestW = " << BestW 
	       << ", Norm = " << Norm << endl;
	*/
	}
	else 
	{
	  /*
	  cout << "  T = " << T 
	       << ", Period = " << Period 
	       << ", P = " << P 
	       << ", BestW = " << BestW 
	       << ", Norm = " << Norm << endl;
	*/
	}
      }
    }

    t = Data->FirstNode;
    do
    {
      t->CandidateSet.clear();
      t->Pi = t->BestPi;
    }
    while((t = t->Suc) != Data->FirstNode);

    W = Minimum1TreeCost();

    if(W < W0)
    {
      do 
	t->Pi = 0;
      while((t = t->Suc) != Data->FirstNode);
      W = Minimum1TreeCost();
      if(AscentCandidates < Data->Dimension)
      {
	if((AscentCandidates *= 2) > Data->Dimension)
	  AscentCandidates = Data->Dimension;
	goto Start;
      }

    }

    return W;
  }

  double LKSolver::Minimum1TreeCost()
  {
    //cout << "1tree" << endl;
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    NodeClass *N = NULL, *N1 = NULL;
    double Sum = 0.0;
    int Max;

    MinimumSpanningTree();
    N = Data->FirstNode;
    do 
    { 
      N->V = -2;
      Sum += N->Pi; 
    }
    while((N = N->Suc) != Data->FirstNode);

    Sum *= -2;
    while((N = N->Suc) != Data->FirstNode)
    {
      N->V++;
      N->Dad->V++;
      Sum += N->Cost;
      N->Next = 0;
    }
    Data->FirstNode->Dad = Data->FirstNode->Suc;
    Data->FirstNode->Cost = Data->FirstNode->Suc->Cost;
    Max = INT_MIN;

    do
    {
      if(N->V == -1)
      {
	Connect(N, Max);	
	if(N->NextCost > Max)
	{ N1 = N; Max = N->NextCost; }

      }
    }
    while((N = N->Suc) != Data->FirstNode);

    N1->Next->V++;
    N1->V++;
    Sum += N1->NextCost;
    Norm = 0;
    do
      Norm += N->V * N->V;
    while((N = N->Suc) != Data->FirstNode);

    if(N1 == Data->FirstNode)
    { N1->Suc->Dad = 0; }
    else 
    {
      Data->FirstNode->Dad = 0;
      Precede(N1, Data->FirstNode);
      Data->FirstNode = N1; 
    }
    
    if(Norm == 0)
    { 
      for(N = Data->FirstNode->Dad; N; N1 = N, N = N->Dad)
      { Follow(N, N1); }
    }
    return Sum; 
  }

  void LKSolver::MinimumSpanningTree()
  {
    //cout << "spanning tree" << endl;
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    
    NodeClass *Blue = NULL, *NextBlue = NULL, *N = NULL;
    multimap<int, NodeClass*> mp;
    int Min, d;

    Blue = N = Data->FirstNode;
    Blue->Dad = NULL;
    Blue->Loc = 0;
    
    //if(Blue->CandidateSet)
    if(Blue->CandidateSet.size() > 0)
    {
      while(( N = N->Suc) != Data->FirstNode)
      {
        N->Dad = Blue;
	N->Cost = N->Rank = INT_MAX;
	//INSERT(N);
	//multimap<int, NodeClass*>::iterator loc_itr  = mp.insert(make_pair(N->Cost, N));
	//locItr loc_itr = mp.insert(make_pair(N->Cost, N));
	//N->LocItr = loc_itr;
	//N->Loc = 1;
	Insert(N);
      }	

      CandidateClass *NBlue;
      //for( NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++)
      for(size_t i = 0; i < Blue->CandidateSet.size(); ++i)
      {
	NBlue = &Blue->CandidateSet[i];
	N = NBlue->To;
	N->Dad = Blue;

	//mp.erase(N->LocItr);
	N->Cost = N->Rank = NBlue->Cost + N->Pi + Blue->Pi;

	//locItr loc_itr = mp.insert(make_pair(N->Cost, N));
	//N->LocItr = loc_itr;
	SiftUp(N);
      }

      //locItr begin = mp.begin(); //end = mp.end();
      //NextBlue = begin->second;
      //NextBlue->Loc = 0;
      //mp.erase(begin);

      while( NextBlue = DeleteMin())
      //while(NextBlue != NULL)
      {
        Follow(NextBlue, Blue);
	Blue = NextBlue;

	//for(NBlue = Blue->CandidateSet; (N = NBlue->To); NBlue++)
	for(size_t i = 0; i < Blue->CandidateSet.size(); ++i)
	{
	  NBlue = &Blue->CandidateSet[i];
	  N = NBlue->To;

	  if(!N->Loc)
	  { continue; }
	  d = NBlue->Cost + N->Pi + Blue->Pi;
	  if( d < N->Cost)
	  {
	    N->Dad = Blue;

	    //mp.erase(N->LocItr);
	    N->Cost = N->Rank = d;
	    SiftUp(N);
	    //locItr loc_itr = mp.insert(make_pair(N->Cost, N));
	    //N->LocItr = loc_itr;
	  }
	}
	/*
	NextBlue = NULL;
	if(mp.size() != 0)
	{
	  begin = mp.begin();
	  NextBlue = begin->second;
	  NextBlue->Loc = 0;
	  mp.erase(begin);
	}
	*/
      }
    }
    else 
    {
      while((N = N->Suc) != Data->FirstNode)
      { N->Cost = INT_MAX; }

      while(( N = Blue->Suc) != Data->FirstNode)
      {
	Min = INT_MAX;
	do
	{
	  d = Data->D(Blue, N);
	  if(d < N->Cost)
	  { N->Cost = d; N->Dad = Blue; }
	  if(N->Cost < Min)
	  { Min = N->Cost; NextBlue = N; }
	}
	while((N = N->Suc) != Data->FirstNode);
	Follow(NextBlue, Blue);
	Blue = NextBlue;
      }
    }
    //cout << "end spanning tree" << endl;
  
    return;
  }

  void LKSolver::GenerateCandidates(const long MaxCandidates, const long MaxAlpha, const int Symmetric)
  {
    //cout << "generate " << endl;
    NodeClass *From, *To;
    CandidateClass *NFrom, *NN, *NTo;
    int a, d, Count;

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();

//#define Mark Next
//#define Beta NextCost
    
    From = Data->FirstNode;
    do 
    {
      From->CandidateSet.clear();
      //if(From->CandidateSet != NULL)
      //{ delete[] From->CandidateSet; }
      //From->CandidateSet = NULL;
      From->Next = NULL;
    } 
    while((From = From->Suc) != Data->FirstNode);

    From = Data->FirstNode;
    do 
    {
      //From->CandidateSet = new CandidateClass[MaxCandidates + 1]; 
      From->CandidateSet.reserve(Data->Dimension*2);
    }
    while((From = From->Suc) != Data->FirstNode);

    do
    {
      //NFrom = From->CandidateSet;
      //NFrom->To = 0;

      if(From != Data->FirstNode)
      {
        From->NextCost = INT_MIN; 
	for(To = From; To->Dad != 0; To = To->Dad)
	{
  	  To->Dad->NextCost = max(To->NextCost, To->Cost);	  
	  To->Dad->Next = From;
	}
      }
      Count = 0;
      To = Data->FirstNode;

      do
      {

	if(To == From)
	{ continue; }  
	d = Data->D(From, To);
	if(From == Data->FirstNode)
	  a = To == From->Dad ? 0 : d - From->NextCost;
	else if(To == Data->FirstNode)
	  a = From == To->Dad ? 0 : d - To->NextCost;
	else 
	{
	  if(To->Next != From)
	    To->NextCost = max(To->Dad->NextCost, To->Cost);
	  a = d - To->NextCost;
	}

	if(To->NextCost == INT_MIN )
	  continue;

	if(InOptimumTour(From, To))
	{ a = 0; }

	if(a <= MaxAlpha)
	{
	  CandidateClass Cand;
	  Cand.To = To;
	  Cand.Cost = d;
	  Cand.Alpha = a;

	  //if(Count < MaxCandidates)
	  //{ Count++;}
	  From->CandidateSet.push_back(Cand);

	  /*
	  NN = NFrom;
	  while(--NN >= From->CandidateSet)
	  {
	  	if(a > NN->Alpha || (a == NN->Alpha && d >= NN->Cost))
		  break;
		*(NN + 1) = *NN;
	  }
	  NN++;
	  NN->To = To;
	  NN->Cost = d;
	  NN->Alpha = a;
	  if(Count < MaxCandidates)
	  { Count++; NFrom++;}
	  NFrom->To = 0;
	  */
	}

      }
      while((To = To->Suc) != Data->FirstNode);
      
    }
    while((From = From->Suc)!= Data->FirstNode);

    From = Data->FirstNode;
    do 
    {
       sort(From->CandidateSet.begin(), From->CandidateSet.end()); 

       if(From->CandidateSet.size() > MaxCandidates)
       {  From->CandidateSet.erase(From->CandidateSet.begin() + MaxCandidates, From->CandidateSet.end()); }
    }
    while((From = From->Suc) != Data->FirstNode);

    if(!Symmetric)
      return;

    To = Data->FirstNode;
    do
    {
      for(size_t i = 0; i < To->CandidateSet.size(); ++i)
      {
        NTo = &To->CandidateSet[i];
	From = NTo->To;
	bool found_flag = false;
	for(size_t j = 0; j < From->CandidateSet.size(); ++j)
	{
	  NN = NFrom = &From->CandidateSet[j];
	  if(NN->To == To)
	  {
	    found_flag = true;
	    break;
	    /*
	    a = NTo->Alpha;
	    d = NTo->Cost;
	    CandidateClass Cand;
	    Cand.To = To;
	    Cand.Cost = d;
	    Cand.Alpha = a;
	    From->CandidateSet.push_back(Cand); 

	    break;
	     */
	  }
	}
	if(!found_flag)
	{
	  a = NTo->Alpha;
	  d = NTo->Cost;
	  CandidateClass Cand;
	  Cand.To = To;
	  Cand.Cost = d;
	  Cand.Alpha = a;
	  From->CandidateSet.push_back(Cand); 
	}
      }
    
    }
    while((To = To->Suc) != Data->FirstNode);

    From = Data->FirstNode;
    do 
       sort(From->CandidateSet.begin(), From->CandidateSet.end()); 
    while((From = From->Suc) != Data->FirstNode);

    return;
  }

  void LKSolver::Connect(NodeClass *N1, int Max)
  {
    int d = 0;
    NodeClass *N;
    CandidateClass *NN1;
    N1->Next = 0;
    N1->NextCost = INT_MAX;

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();

    if(N1->CandidateSet.size() == 0)
    {
      N = Data->FirstNode;
      do
      {
	if(N == N1 || N == N1->Dad || N1 == N->Dad)
	  continue;
	d = Data->D(N1, N);
        //if( (d=Data->D(N1, N)) < N1->NextCost)
	if(d < N1->NextCost)
	{ 
	  N1->NextCost = d;
	  if(d <= Max)
	    //break;
	    return;
	  N1->Next = N;
	}	  
      }
      while((N = N->Suc) != Data->FirstNode);
    }
    else
    {
      for(size_t i = 0; i < N1->CandidateSet.size(); ++i)
      {
        NN1 = &N1->CandidateSet[i];
	N = NN1->To;
	if(N == N1->Dad || N1 == N->Dad)
	  continue;
	d = NN1->Cost + N1->Pi + N->Pi;
	//if( ( d = NN1->Cost + N1->Pi + N->Pi) < N1->NextCost)
	if(d < N1->NextCost)
	{ 
	  N1->NextCost = d;
	  if(d <= Max)
	    return;
	  N1->Next = N;
	}
      }	
    }
    return;
  }

  void LKSolver::ResetCandidateSets()
  {

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    NodeClass *From = Data->FirstNode;

    do 
    {
      sort(From->CandidateSet.begin(), From->CandidateSet.end());

      size_t pos;
      while( From->CandidateSet[pos = From->CandidateSet.size() - 1].Alpha == INT_MAX)
        From->CandidateSet.pop_back(); 
    } 
    while((From = From->Suc) != Data->FirstNode);

    return;

  }


  void LKSolver::AdjustCandidateSets()
  {

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    NodeClass *From = Data->FirstNode;

    do
    {
       NodeClass *Nodes[2] = {From->Pred, From->Suc };
       for(int i = 0; i < 2; ++i)
       {
         NodeClass *To = Nodes[i];

	 for(size_t j = 0; j < From->CandidateSet.size(); ++j)
	 {
	   CandidateClass *NFrom = &From->CandidateSet[j];
	   if(NFrom->To == To)
	     goto next_from;
	 }

	 CandidateClass Cand;
	 Cand.Cost = Data->C(From, To);
	 Cand.To = To;
	 Cand.Alpha = INT_MAX;
	 From->CandidateSet.push_back(Cand);

       }

       next_from:;
    }
    while((From = From->Suc) != Data->FirstNode);

    do
    {
      int count = 0;
      for(size_t i = 0; i < From->CandidateSet.size(); ++i)
      {
        CandidateClass *NFrom = &From->CandidateSet[i];
	NodeClass *To = NFrom->To;
	if(InBestTour(From, To) && InNextBestTour(From, To))
	{
	  iter_swap(From->CandidateSet.begin() + count++, From->CandidateSet.begin() + i);
	}
      }

    }
    while((From = From->Suc) != Data->FirstNode);
  }

};


