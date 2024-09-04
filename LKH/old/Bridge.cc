#include "LK.h"

namespace myOptimizeName
{
  int LKSolver::Gain23()
  {

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    int Gain = 0;
    Reversed = 0;

    for(int X2 = 1; X2 <= 2; ++X2)
    {
      Reversed ^= 1;
      NodeClass *s1 = Data->FirstNode, *s2;

      do
      {
	vector<NodeClass*> Nodes(2*kopt+1), BestNodes(2*kopt+1);
	int BestGk = INT_MIN;
	long int BestFlipKey;

	SegmentClass *FirstFlipSegment = &FlipSegs[1];
	FirstFlipSegment->SetPred(FirstFlipSegment); FirstFlipSegment->SetSuc(FirstFlipSegment);
	FirstFlipSegment->SetFrom(1); FirstFlipSegment->SetTo(2);
	FirstFlipSegment->NonSequence = true;

	s2 = SUC(s1);
	int G0 = Data->C(s1, s2);
	Gain = 0;
	Nodes[1] = s1; Nodes[2] = s2;
	//BestKOptMoveR(1, Nodes, BestNodes, FirstFlipSegment, BestFlips, G0, &BestGk, &G0, &Gain, true, 2);

	BestKOptMoveR(1, Nodes, BestNodes, FirstFlipSegment, -12, &BestFlipKey, G0, &BestGk, &G0, &Gain, true, 2);

	if(Gain > 0)
	  return Gain;
      }
      while((s1 = s2) != Data->FirstNode);
    }
  
    return Gain;
  }

  int LKSolver::BrigeMove(vector<NodeClass*> &Nodes, SegmentClass *FirstFlipSegment, long int ParentKey, int G, const int k)
  {
  
    int Gain = 0;
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
   
    //NodeClass *u2 = Nodes[FirstFlipSegment->GetTo()];
    //NodeClass *u3 = Nodes[FirstFlipSegment->GetSuc()->GetFrom()];
    NodeClass *u2 = Nodes[2];
    NodeClass *u3 = Nodes[3];
    long int ParentKey2, Key;
    bool feasibleFlip;

    FlipLogClass *Pnt = &FlipHash[ParentKey];

    for(NodeClass *t1 = u2; t1 != u3; t1 = SUC(t1))
    {
      NodeClass *t2 = SUC(t1);
      Nodes[2*k+1] = t1; Nodes[2*k+2] = t2;
      //if(find_(Nodes, 2*k, t1) && find_(Nodes,2*k, t2))
//	continue;

      if( (Nodes[1] == t1 && Nodes[2] == t2) || 
	  (Nodes[1] == t2 && Nodes[2] == t1) || 
	  (Nodes[3] == t1 && Nodes[4] == t2) || 
	  (Nodes[3] == t2 && Nodes[4] == t1))
	continue;

      int G0 = G + Data->C(t1, t2);

      for(size_t i = 0; i < t2->CandidateSet.size(); ++i)
      {
	CandidateClass *Nt2 = &t2->CandidateSet[i];
	NodeClass *t3 = Nt2->To;

	if(t3 == t2->Pred || t3 == t2->Suc)
	  continue;
	int G1 = G0 - Nt2->Cost;

	for(int X4 = 1; X4 <= 2; ++X4)
	{
	  NodeClass *t4 = X4 == 1 ? SUC(t3) : PRED(t3);

	  if( (Nodes[1] == t3 && Nodes[2] == t4) || 
	      (Nodes[1] == t4 && Nodes[2] == t3) || 
	      (Nodes[3] == t3 && Nodes[4] == t4) || 
	      (Nodes[3] == t4 && Nodes[4] == t3))
	    continue;
	  if(t4 == t2 || t4 == t1)
	    continue;
	  
	  //if(t4 == t2 || t4 == t1 || (find_(Nodes, 2*k, t3) && find_(Nodes, 2*k, t4)))
	    //continue;
	  int G2 = G1 + Data->C(t3, t4);
	  Gain = G2 - Data->C(t4, t1); 


	  if(Gain > 0)
	  {
	    Nodes[2*k+3] = t3; Nodes[2*k+4] = t4;

	    /*
	    SegmentClass *tSeg = &FlipSegs[k+1], *tSeg2 = &FlipSegs[k+2];
	    tSeg->SetFrom(2*k+1); tSeg->SetTo(2*k+2);
	    tSeg2->SetFrom(X4 == 1 ? 2*k+3 : 2*k+4); tSeg2->SetTo(X4 == 1 ? 2*k+4 : 2*k+3);
	    tSeg2->NonSequence = tSeg->NonSequence = false;
	    */
	    int x1, x2;
	    //vector<NextFlipClass> *CandidateBtws;

	    // ============== //

	    NextFlipClass *NextCnd = NULL, *NextCnd2 = NULL;

	    //if(!(InsertSegR3(2*k+1, Nodes, &FlipHash[ParentKey], NextCnd)))
	    if(!(InsertSegR3(2*k+1, Nodes, Pnt, NextCnd)))
	       goto NEXTNODE1;

	    x1 = NextCnd->Xks[0] == 0 ? 0 : 1;
	    ParentKey2 = NextCnd->LogKeys[x1];

	    if(!(InsertSegR3(2*k+3, Nodes, &BrigeHash[ParentKey2], NextCnd2)))
	      goto NEXTNODE2;

	    x2 = NextCnd2->Xks[0] + 1 == X4 ? NextCnd2->Xks[0] : NextCnd2->Xks[1];

	    if(x2 < 0)
	      continue;

	    Key = NextCnd2->LogKeys[x2];
	    feasibleFlip = NextCnd2->Feasibles[x2];

	    /*
	    cout << Key << " " << (feasibleFlip) << endl;
	    myGetChar();

	    continue;
	    */
	    if(feasibleFlip)
	    {
	      MakeOptMove(Nodes, *BrigeHash[Key].GetFlipStacks()); 
	      return Gain; 
	    }

	    continue;
	    // =============== // 

	    /*
	    if(!(link1=InsertSeg(FirstFlipSegment, tSeg, Nodes, ParentKey, false)))
	      goto NEXTNODE1;
	    ParentKey2 = GetFlipLogKey(FirstFlipSegment);
	    if(!(link2=InsertSeg(FirstFlipSegment, tSeg2, Nodes, ParentKey2, true)))
	      goto NEXTNODE1;
	      */

	    Key = GetFlipLogKey(FirstFlipSegment)*2;
	    feasibleFlip = BrigeHash[Key].IsFeasibleLog();

	    if(feasibleFlip)
	    {
	      //MakeOptMove(Nodes, FlipStacks);
	      MakeOptMove(Nodes, *BrigeHash[Key].GetFlipStacks()); 
	      return Gain; 
	    }

NEXTNODE2:;
	  Nodes[2*k+3] = Nodes[2*k+4] = NULL;
	  }
	}

      }

NEXTNODE1:;

	  Nodes[2*k+1] = Nodes[2*k+2] = NULL;
    }

    return 0;
  }
};
