#include "LK.h"

namespace myOptimizeName
{
  int LKSolver::Excludable(const NodeClass *ta, const NodeClass *tb)
  {
    if (ta == tb->OldPred)
      return !tb->OldPredExcluded;
    if (ta == tb->OldSuc)
      return !tb->OldSucExcluded;
    return 0;
  }
  void LKSolver::Exclude(NodeClass *ta, NodeClass *tb)
  {
    if (ta == tb->Pred || ta == tb->Suc)
      return;
    if (ta == tb->OldPred)
      tb->OldPredExcluded = 1;
    else if (ta == tb->OldSuc)
      tb->OldSucExcluded = 1;
    if (tb == ta->OldPred)
      ta->OldPredExcluded = 1;
    else if (tb == ta->OldSuc)  
      ta->OldSucExcluded = 1;
  }

  NodeClass *LKSolver::BestKOptMove(NodeClass *t1, NodeClass *t2, int *G0, int *Gain)
  {

    NodeClass *rval = NULL;
    vector<NodeClass*> Nodes(2*kopt+1), BestNodes(2*kopt+1);
    int BestGk = INT_MIN;
    long int BestFlipKey;

    if(SUC(t1) != t2)
      Reversed ^= 1;

    SegmentClass *FirstFlipSegment = &FlipSegs[1];

    Nodes[1] = t1; Nodes[2] = t2;
    BestNodes[1] = t1; BestNodes[2] = t2;
    BestNodes[2*kopt] = NULL;

    rval = BestKOptMoveR(1, Nodes, BestNodes, FirstFlipSegment, 12, &BestFlipKey, *G0, &BestGk, G0, Gain);

    if(rval != NULL)
      return rval; 
    *Gain = 0;
    if(BestNodes[2*kopt])
    {
      //MakeOptMove(BestNodes, BestFlips);
      MakeOptMove(BestNodes, *FlipHash[BestFlipKey].GetFlipStacks());
      for(int k = 1; k <= kopt; ++k)
	Exclude(BestNodes[2*k-1], BestNodes[2*k]);
      *G0 = BestGk;
      rval = BestNodes[2*kopt];
    }

    return rval;  
  }

  NodeClass *LKSolver::BestKOptMoveR(
      int k, 
      vector<NodeClass*> &Nodes, vector<NodeClass*> &BestNodes, 
      SegmentClass *FirstFlipSegment, 
      long int ParentKey,
      long int *BestFlipKey,
      int Gk, int *BestGk, int *G0, int *Gain, 
      bool brige, const int klimit)
  {
    if(k + 1 > kopt || k + 1 > klimit) 
      return NULL;

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    NodeClass *rval = NULL;
    int gain_ = 0;
    const int _2k = 2*k, _2k_1 = 2*k+1, _2k_2 = 2*k + 2;
    NodeClass *Node2k, *Node2k_1, *Node2k_2;
    Node2k = Nodes[_2k];

    FlipLogClass *Pnt = &FlipHash[ParentKey];

    for(size_t i = 0; i < Node2k->CandidateSet.size(); ++i)
    {
      CandidateClass *Nt = &Node2k->CandidateSet[i];
      Node2k_1 = Nodes[_2k_1] = Nt->To;
      const int G2k_1 = Gk - Nt->Cost;

      if(Node2k->Pred == Node2k_1 || Node2k->Suc == Node2k_1 || G2k_1 <= 0)
	continue;

      if( ( k + 1 >= 3 && Node2k == Nodes[2] && Node2k_1 == Nodes[3]) || 
	  ( k + 1 >= 3 && Node2k == Nodes[3] && Node2k_1 == Nodes[2]) || 
	  ( k + 1 >= 4 && Node2k == Nodes[4] && Node2k_1 == Nodes[5]) || 
	  ( k + 1 >= 4 && Node2k == Nodes[5] && Node2k_1 == Nodes[4]) ||
	  ( k + 1 >= 5 && Node2k == Nodes[6] && Node2k_1 == Nodes[7]) || 
	  ( k + 1 >= 5 && Node2k == Nodes[7] && Node2k_1 == Nodes[6]))
	continue;

      if(Node2k_1 == Nodes[1] && Node2k_1 == Nodes[2])
	continue;

      NextFlipClass *NextCnd = NULL;
      if(!InsertSegR3(_2k_1, Nodes, Pnt, NextCnd))
        continue; 
      
      for(int x = 0; x <= 1; ++x)
      {
	int Xk = NextCnd->Xks[x] + 1;
	
	if(NextCnd->Xks[x] < 0)
	  continue;

	Node2k_2 = Nodes[_2k_2] = Xk == 1 ? SUC(Node2k_1) : PRED(Node2k_1);

	if(( k + 1 >= 2 && Node2k_2 == Nodes[1] && Node2k_1 == Nodes[2]) || 
	    ( k + 1 >= 2 && Node2k_2 == Nodes[2] && Node2k_1 == Nodes[1]) || 
	    ( k + 1 >= 3 && Node2k_2 == Nodes[3] && Node2k_1 == Nodes[4]) || 
	    ( k + 1 >= 3 && Node2k_2 == Nodes[4] && Node2k_1 == Nodes[3]) || 
	    ( k + 1 >= 4 && Node2k_2 == Nodes[5] && Node2k_1 == Nodes[6]) || 
	    ( k + 1 >= 4 && Node2k_2 == Nodes[6] && Node2k_1 == Nodes[5]) ||
	    ( k + 1 >= 5 && Node2k_2 == Nodes[8] && Node2k_1 == Nodes[7]) || 
	    ( k + 1 >= 5 && Node2k_2 == Nodes[7] && Node2k_1 == Nodes[8]))
	  continue;

	if (Node2k_2 == Nodes[1])
	  continue;

	const int G2k = G2k_1 + Data->C(Node2k_1, Node2k_2);

	long int LogKey = NextCnd->LogKeys[x];
	bool feasibleFlip = NextCnd->Feasibles[x];

	if(brige && !feasibleFlip && Nodes[_2k_2] != Nodes[1] && !(klimit > 2 && k + 1 == 2) )
	{
	  gain_ = BrigeMove(Nodes, FirstFlipSegment, LogKey, G2k - Data->C(Nodes[1], Nodes[_2k_2]), k+1);
	  if(gain_ > 0)
	  { 
	    *Gain = gain_;
	    return NULL;
	  }
	}

	if(Node2k_2 != Nodes[1] && feasibleFlip)
	  gain_ = G2k - Data->C(Node2k_2, Nodes[1]);

	if(feasibleFlip && gain_ > 0 && Nodes[1] != Node2k_2 && !brige)
	{
	  MakeOptMove(Nodes, *FlipHash[LogKey].GetFlipStacks());
	  *G0 = G2k;
	  *Gain = gain_;

	  return Node2k_2;
	}

	if(k + 2 <= kopt)
	  rval = BestKOptMoveR(k+1, Nodes, BestNodes, FirstFlipSegment, LogKey, BestFlipKey, G2k, BestGk, G0, Gain, brige, klimit);

	if(rval != NULL)
	  return rval;

	if(k + 1 == kopt && feasibleFlip && Nodes[1] != Node2k_2)
	{
	  if(G2k >= *BestGk && Swaps < MaxSwaps && G2k - Precision >= Node2k_2->Cost 
	      && Excludable(Node2k_1, Node2k_2))
	  {
	    for(int ii = 1; ii <= 2*kopt; ++ii)
	      BestNodes[ii] = &Data->NodeSet[Nodes[ii]->Id];
	    *BestFlipKey = LogKey;
	    *BestGk = G2k;
	  }
	}

	Nodes[_2k_2] = NULL;
      }

      Nodes[_2k_1] = NULL;
    }

    return rval;

  }

  bool LKSolver::find_(vector<NodeClass*> &Nodes, const int until, NodeClass* Node)
  {
    for(size_t i = 0; i < until; ++i)
      if(Nodes[i] != NULL && Nodes[i] == Node)
	return true;
    return false; 
  }

  void LKSolver::MakeOptMove(vector<NodeClass*> &Nodes, vector<SwapRecord> &Flips)
  {
    for(size_t f = 0; f < Flips.size(); ++f)
    {
      NodeClass *f1 = Nodes[Flips[f].id1], *f2 = Nodes[Flips[f].id2], *f3 = Nodes[Flips[f].id3];
      Swap1(f1, f2, f3);
    }  
  }

  bool LKSolver::InsertSegR3(int n, vector<NodeClass*> &Nodes, FlipLogClass *Pnt, NextFlipClass *&NextCnd)
  {
    //const int X = tSeg->GetFrom() < tSeg->GetTo() ? 0 : 1;
    //vector<NextFlipClass> *CandidateBtws = Pnt->GetNextBetweens(X);
    vector<NextFlipClass> *CandidateBtws = Pnt->GetNextBetweens();

    for(size_t i = 0; i < CandidateBtws->size(); ++i)
    {
      NextFlipClass *CndBtw = &CandidateBtws->at(i);

      int From = CndBtw->GetFrom(), To = CndBtw->GetTo();
      //int K1 = (From+1) / 2, K2 = (To+1) / 2;
      NodeClass *t1 = Nodes[From], *t2 = Nodes[n], *t3 = Nodes[To];

      if(BETWEEN(t1, t2, t3))
      {
	//cout << From << " " << To << endl;
	//Key = CndBtw->GetKey();
	//Feasible = CndBtw->Feasible();
	NextCnd = CndBtw;
	return true;
      }
    }

    return false;
  }


  /*
  bool LKSolver::InsertSegR2(SegmentClass *First, SegmentClass *tSeg, vector<NodeClass*> &Nodes, FlipLogClass *Pnt, long int &Key, bool &Feasible)
  {
    //FlipLogClass *Pnt = &FlipHash[ParentKey];

    const int X = tSeg->GetFrom() < tSeg->GetTo() ? 0 : 1;
    vector<NextFlipClass> *CandidateBtws = Pnt->GetNextBetweens(X);

    for(size_t i = 0; i < CandidateBtws->size(); ++i)
    {
      NextFlipClass *CndBtw = &CandidateBtws->at(i);

      int From = CndBtw->GetFrom(), To = CndBtw->GetTo();
      int K1 = (From+1) / 2, K2 = (To+1) / 2;

      SegmentClass *Seg1 = &FlipSegs[K1], *Seg2 = &FlipSegs[K2];
      NodeClass *t1 = Nodes[From], *t2 = Nodes[tSeg->GetFrom()], *t3 = Nodes[To];

      if(BETWEEN(t1, t2, t3))
      {

	//SegmentClass *tmp1 = Seg1->GetSuc();

	//Seg1->SetSuc(tSeg);
	//tSeg->SetSuc(tmp1);

	//tmp1->SetPred(tSeg);
	//tSeg->SetPred(Seg1);
	Key = CndBtw->GetKey();
	Feasible = CndBtw->Feasible();

	return true;
      }
    }

    return false;
  }

  bool LKSolver::InsertSeg(SegmentClass *First, SegmentClass *tSeg, vector<NodeClass*> &Nodes, long int ParentKey, bool brige)
  {
    if(ParentKey == 0)
    {
      SegmentClass *Seg = First;
      do
      {
	SegmentClass *tmp1;
	NodeClass *t1 = Nodes[Seg->GetTo()], *t2 = Nodes[tSeg->GetFrom()], *t3 = Nodes[Seg->GetSuc()->GetFrom()];
	if(BETWEEN(t1, t2, t3))
	{
	  tmp1 = Seg->GetSuc();

	  Seg->SetSuc(tSeg);
	  tSeg->SetSuc(tmp1);

	  tmp1->SetPred(tSeg);
	  tSeg->SetPred(Seg);

	  return true;
	}
      }
      while((Seg = Seg->GetSuc()) != First);
    }
    else 
    {
      FlipLogClass *Pnt = !brige ? &FlipHash[ParentKey] : &BrigeHash[ParentKey];

      const int X = tSeg->GetFrom() < tSeg->GetTo() ? 0 : 1;
      vector<NextFlipClass> *CandidateBtws = Pnt->GetNextBetweens(X);

      for(size_t i = 0; i < CandidateBtws->size(); ++i)
      {
	NextFlipClass *CndBtw = &CandidateBtws->at(i);

	int From = CndBtw->GetFrom(), To = CndBtw->GetTo();
	int K1 = (From+1) / 2, K2 = (To+1) / 2;

	SegmentClass *Seg1 = &FlipSegs[K1], *Seg2 = &FlipSegs[K2];
	NodeClass *t1 = Nodes[From], *t2 = Nodes[tSeg->GetFrom()], *t3 = Nodes[To];

	if(BETWEEN(t1, t2, t3))
	{
	  SegmentClass *tmp1 = Seg1->GetSuc();

	  Seg1->SetSuc(tSeg);
	  tSeg->SetSuc(tmp1);

	  tmp1->SetPred(tSeg);
	  tSeg->SetPred(Seg1);

	  return true;
	}
      }
    }

    return false;
  }
*/
  void LKSolver::InsertFlipLog(const long int Logkey, const bool feasible, vector<SwapRecord> &Flips)
  { 
    if(FlipHash.find(Logkey) != FlipHash.end())
    {
      FlipHash[Logkey].SetLogKey(Logkey);
      FlipHash[Logkey].SetFeasible(feasible);
      FlipHash[Logkey].SetFlipStacks(Flips);   
    } 
    else 
      FlipHash[Logkey] = FlipLogClass(Logkey, feasible, Flips); 
  }
  void LKSolver::InsertBrigeFlipLog(const long int Logkey, const bool feasible, vector<SwapRecord> &Flips)
  {
    if(BrigeHash.find(Logkey) != BrigeHash.end()) 
    {
      BrigeHash[Logkey].SetLogKey(Logkey);
      BrigeHash[Logkey].SetFeasible(feasible);
      BrigeHash[Logkey].SetFlipStacks(Flips);
    }
    else 
      BrigeHash[Logkey] = FlipLogClass(Logkey, feasible, Flips); 
  }

  bool LKSolver::SearchFlipLogKey(long int hashkey)
  { return FlipHash.find(hashkey) != FlipHash.end(); }

  long int LKSolver::GetFlipLogKey(SegmentClass *FirstSegment)
  {
    long int key = 0;
    SegmentClass *Seg = FirstSegment;
    bool nonSequence = false;

    do
    {
      key = 10 * key + Seg->GetFrom();
      key = 10 * key + Seg->GetTo();

      if(Seg->NonSequence)
	nonSequence = true;
    }
    while((Seg = Seg->GetSuc()) != FirstSegment);

    if(nonSequence) 
      key *= -1;

    return key;
  }

};
