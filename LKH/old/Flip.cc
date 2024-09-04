#include "LK.h"

namespace myOptimizeName
{

#if 0
  using Node = NodeClass;
  using Segment = SegmentClass;

  void LKSolver::Flip_SL(Node * t1, Node * t2, Node * t3)
  {
    Node *t4, *a, *b, *c, *d;
    Segment *P1, *P2, *P3, *P4, *Q1, *Q2;
    Node *s1, *s2;
    int i, Temp;

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();

    if (t3 == t2->Pred || t3 == t2->Suc)
      return;
    //if (Groups == 1) {
    //  Flip(t1, t2, t3);
     // return;
   // }
    t4 = t2 == SUC(t1) ? PRED(t3) : SUC(t3);
    P1 = t1->Parent;
    P2 = t2->Parent;
    P3 = t3->Parent;
    P4 = t4->Parent;

    /* Split segments if needed */
    if (P1 != P3 && P2 != P4) {
      if (P1 == P2) {
	SegmentSplit(t1, t2);
	P1 = t1->Parent;
	P2 = t2->Parent;
      }
      if (P3 == P4 && P1 != P3 && P2 != P4) {
	SegmentSplit(t3, t4);
	P3 = t3->Parent;
	P4 = t4->Parent;
      }
    } else if ((P1 == P3
	  && abs(t3->Rank - t1->Rank) > 0.75 * Data->GroupSize)
	|| (P2 == P4
	  && abs(t4->Rank - t2->Rank) > 0.75 * Data->GroupSize)) {
      if (P1 == P2) {
	SegmentSplit(t1, t2);
	P1 = t1->Parent;
	P2 = t2->Parent;
	P3 = t3->Parent;
	P4 = t4->Parent;
      }
      if (P3 == P4) {
	SegmentSplit(t3, t4);
	P1 = t1->Parent;
	P2 = t2->Parent;
	P3 = t3->Parent;
	P4 = t4->Parent;
      }
    }
    /* Check if it is possible to flip locally within a segment */
    b = 0;
    if (P1 == P3) {
      /* Either the t1 --> t3 path or the t2 --> t4 path lies 
	 within one segment */
      if (t1->Rank < t3->Rank) {
	if (P1 == P2 && P1 == P4 && t2->Rank > t1->Rank) {
	  a = t1;
	  b = t2;
	  c = t3;
	  d = t4;
	} else {
	  a = t2;
	  b = t1;
	  c = t4;
	  d = t3;
	}
      } else {
	if (P1 == P2 && P1 == P4 && t2->Rank < t1->Rank) {
	  a = t3;
	  b = t4;
	  c = t1;
	  d = t2;
	} else {
	  a = t4;
	  b = t3;
	  c = t2;
	  d = t1;
	}
      }
    } else if (P2 == P4) {
      /* The t2 --> t4 path lies within one segment */
      if (t4->Rank < t2->Rank) {
	a = t3;
	b = t4;
	c = t1;
	d = t2;
      } else {
	a = t1;
	b = t2;
	c = t3;
	d = t4;
      }
    }
    if (b) {
      //int Cbc = Data->C(b, c), Cda = Data->C(d, a);
      /* Flip locally (b --> d) within a segment */
      i = d->Rank;
      d->Suc = 0;
      s2 = b;
      while ((s1 = s2)) {
	s2 = s1->Suc;
	s1->Suc = s1->Pred;
	s1->Pred = s2;
	s1->Rank = i--;
	//Temp = s1->SucCost;
	//s1->SucCost = s1->PredCost;
	//s1->PredCost = Temp;
      }
      d->Pred = a;
      b->Suc = c;
     // d->PredCost = Cda;
      //b->SucCost = Cbc;
      if (a->Suc == b) {
	a->Suc = d;
	//a->SucCost = d->PredCost;
      } else {
	a->Pred = d;
	//a->PredCost = d->PredCost;
      }
      if (c->Pred == d) {
	c->Pred = b;
//	c->PredCost = b->SucCost;
      } else {
	c->Suc = b;
//	c->SucCost = b->SucCost;
      }
      if (b->Parent->First == b)
	b->Parent->First = d;
      else if (d->Parent->First == d)
	d->Parent->First = b;
      if (b->Parent->Last == b)
	b->Parent->Last = d;
      else if (d->Parent->Last == d)
	d->Parent->Last = b;
    } else {
      int Ct2t3, Ct4t1;
      /* Reverse a sequence of segments */
      if (P1->Suc != P2) {
	a = t1;
	t1 = t2;
	t2 = a;
	a = t3;
	t3 = t4;
	t4 = a;
	Q1 = P1;
	P1 = P2;
	P2 = Q1;
	Q1 = P3;
	P3 = P4;
	P4 = Q1;
      }
      /* Find the sequence with the smallest number of segments */

      if ((i = P2->Rank - P3->Rank) < 0)
	i += Groups;
      if (2 * i > Groups) {
	a = t3;
	t3 = t2;
	t2 = a;
	a = t1;
	t1 = t4;
	t4 = a;
	Q1 = P3;
	P3 = P2;
	P2 = Q1;
	Q1 = P1;
	P1 = P4;
	P4 = Q1;
      }
      //Ct2t3 = C(t2, t3);
      //Ct4t1 = C(t4, t1);
      /* Reverse the sequence of segments (P3 --> P1). 
	 Mirrors the corresponding code in the Flip function */
      i = P1->Rank;
      P1->Suc = 0;
      Q2 = P3;
      while ((Q1 = Q2)) {
	Q2 = Q1->Suc;
	Q1->Suc = Q1->Pred;
	Q1->Pred = Q2;
	Q1->Rank = i--;
	Q1->Reversed ^= 1;

      }
      P3->Suc = P2;
      P2->Pred = P3;
      P1->Pred = P4;
      P4->Suc = P1;
      if (t3->Suc == t4) {
	t3->Suc = t2;
	//t3->SucCost = Ct2t3;
      } else {
	t3->Pred = t2;
	//t3->PredCost = Ct2t3;
      }
      if (t2->Suc == t1) {
	t2->Suc = t3;
//	t2->SucCost = Ct2t3;
      } else {
	t2->Pred = t3;
//	t2->PredCost = Ct2t3;
      }
      if (t1->Pred == t2) {
	t1->Pred = t4;
//	t1->PredCost = Ct4t1;
      } else {
	t1->Suc = t4;
//	t1->SucCost = Ct4t1;
     }
      if (t4->Pred == t3) {
	t4->Pred = t1;
//	t4->PredCost = Ct4t1;
      } else {
	t4->Suc = t1;
//	t4->SucCost = Ct4t1;
      }
    }
    SwapStack[Swaps].t1 = t1;
    SwapStack[Swaps].t2 = t2;
    SwapStack[Swaps].t3 = t3;
    SwapStack[Swaps].t4 = t4;
    Swaps++;

  }
#endif
#if 1
  void LKSolver::Flip_SL(NodeClass *t1, NodeClass *t2, NodeClass *t3)
  {
    NodeClass *t4 = t2 == SUC(t1) ? PRED(t3) : SUC(t3);
    SegmentClass *P1 = t1->Parent, *P2 = t2->Parent, 
		 *P3 = t3->Parent, *P4 = t4->Parent;
    SegmentClass *Q1 = NULL, *Q2 = NULL, *Q3 = NULL;

    if(t3 == t2->Pred || t3 == t2->Suc)
      return;
    //cout << Reversed <<endl;
    //cout << "ID = " << t1->Id << " " << t2->Id << " " << t3->Id << " " << t4->Id << endl;
    //cout << "RANK = " << t1->Rank << " " << t2->Rank << "  " << t3->Rank << " " << t4->Rank << endl;
    //cout << "PARENT RANK = " << P1->Rank << " " << P2->Rank << " " << P3->Rank << " " << P4->Rank << endl;
    //cout << "ready flip" << endl;

    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();

    /*
       NodeClass *k = Data->FirstNode;
       int n = 0;
       k = Data->FirstNode;
       do
       {
       cout << k->Id << " " << k->Rank << " " << k->Parent->Rank << "/ " << k->Pred->Id << ". " << k->Suc->Id << endl; 
       ++n;
       }
       while((k = SUC(k)) != Data->FirstNode);

       cout << "BEFORE SPLIT SEG" << endl;
       if(n != Data->Dimension)
       { cout << "n = " << n << endl;}
     */
    if(P1 == P2)
    {
      SegmentSplit(t1, t2);
      P1 = t1->Parent;
      P2 = t2->Parent;
      P3 = t3->Parent;
      P4 = t4->Parent;
    }
    if(P3 == P4)
    {
      SegmentSplit(t3, t4);
      P1 = t1->Parent;
      P2 = t2->Parent;
      P3 = t3->Parent;
      P4 = t4->Parent;
    }


    /*
       cout << endl;
       k = Data->FirstNode;
       n = 0;
       do
       {
       cout << k->Id << " " << k->Rank << " " << k->Parent->Rank << "/ " << k->Pred->Id << ". " << k->Suc->Id << endl; 
       ++n;
       }
       while((k = SUC(k)) != Data->FirstNode);
       cout << "AFTER SPLIT SEGMENT " << endl;
       cout << "NEW PARENT ID = " << P1->Rank << " " << P2->Rank << " " << P3->Rank << " " << P4->Rank << endl;
       if(n != Data->Dimension)
       { cout << "n = " << n << endl;}
     */

    NodeClass *a = NULL, *b = NULL, *c = NULL, *d = NULL;

    if(P1 == P3)
    {
      if(t1->Rank < t3->Rank)
      {
	a = t2;
	b = t1;
	c = t4;
	d = t3;	
      }
      else 
      {
	a = t4;
	b = t3;
	c = t2;
	d = t1;
      }
    }
    else if(P2 == P4)
    {
      //cout << "p2 == p4" << endl;
      if(t4->Rank < t2->Rank)
      {
	//cout << "t4 < t2" << endl;
	a = t3;
	b = t4;
	c = t1;
	d = t2;
      }
      else 
      {
	//cout << "t2 > t4" << endl;
	a = t1;
	b = t2;
	c = t3;
	d = t4;
      }
    }

    if(a)
    {

      int i = b->Rank;
      NodeClass *s1, *s2;

      b->Pred = NULL;
      s2 = d;

      while((s1 = s2))
      {
	s2 = s1->Pred;
	s1->Pred = s1->Suc;
	s1->Suc = s2;
	s1->Rank = i++;
      }

      d->Pred = a;
      b->Suc = c;

      if(a->Suc == b)
	a->Suc = d;
      else 
	a->Pred = d;
      if(c->Pred == d)
	c->Pred = b;
      else 
	c->Suc = b;
      if(b->Parent->First == b)
	b->Parent->First = d;
      else if(d->Parent->First == d)
	d->Parent->First = b;
      if(b->Parent->Last == b)
	b->Parent->Last = d;
      else if(d->Parent->Last == d)
	d->Parent->Last = b;

    }
    else
    {

      if(P1->Rank > P2->Rank)
      {
	a = t1;
	t1 = t2;
	t2 = a;
	a = t3;
	t3 = t4;
	t4 = a;

	P1 = t1->Parent;
	P2 = t2->Parent;
	P3 = t3->Parent;
	P4 = t4->Parent;
      }

      if(P2->Rank > P4->Rank)
      {
	a = t3;
	t3 = t1;
	t1 = a;
	a = t2;
	t2 = t4;
	t4 = a;

	P1 = t1->Parent;
	P2 = t2->Parent;
	P3 = t3->Parent;
	P4 = t4->Parent;
      }

      if(P1->Rank > P2->Rank)
      {
	a = t1;
	t1 = t2;
	t2 = a;
	a = t3;
	t3 = t4;
	t4 = a;

	P1 = t1->Parent;
	P2 = t2->Parent;
	P3 = t3->Parent;
	P4 = t4->Parent;
      }

      int i = P2->Rank;
      P2->Pred = NULL; 
      Q2 = P4;

      while((Q1 = Q2))
      {
	Q2 = Q1->Pred;
	Q1->Pred = Q1->Suc;
	Q1->Suc = Q2;
	Q1->Rank = i++;
	Q1->Reversed ^= 1;
      }
      P2->Suc = P3;                                                                                                                    
      P3->Pred = P2;
      P1->Suc = P4;
      P4->Pred = P1;
      if (t3->Suc == t4)
	t3->Suc = t2;
      else
	t3->Pred = t2;
      if (t2->Suc == t1)
	t2->Suc = t3;
      else
	t2->Pred = t3;
      if (t1->Pred == t2)
	t1->Pred = t4;
      else
	t1->Suc = t4;
      if (t4->Pred == t3)
	t4->Pred = t1;
      else
	t4->Suc = t1; 
    }

    /*
       cout << endl;
       k = Data->FirstNode;

       n = 0;
       do
       {
       string label = "";
       if(k == t1)
       label = "[NODE1] ";
       else if(k == t2)
       label = "[NODE2] ";
       else if(k == t3)
       label = "[NODE3] ";
       else if(k == t4)
       label = "[NODE4] ";
       cout << label << k->Id << " " << k->Rank << " " << k->Parent->Rank << "/ " << k->Pred->Id << ". " << k->Suc->Id << endl; 
       ++n;
       }
       while((k = SUC(k)) != Data->FirstNode);

       cout << "AFTER FLIP" << endl;
       if(n != Data->Dimension)
       { cout << "n = " << n << endl; }
     */
    SwapStack[Swaps].t1 = t1;
    SwapStack[Swaps].t2 = t2;
    SwapStack[Swaps].t3 = t3;
    SwapStack[Swaps].t4 = t4;
    ++Swaps;

    return; 
  }
#endif 

  void LKSolver::SegmentSplit(NodeClass *t1, NodeClass *t2)
  {
    SegmentClass *P = t1->Parent, *Q;

    if(t1->Rank > t2->Rank)
    {
      NodeClass *t = t1;
      t1 = t2;
      t2 = t;
    }

    int Count = t1->Rank - P->First->Rank + 1;
    NodeClass *u = NULL;

    if(Count  < P->Size / 2)
    {
      Q = !P->Reversed ? P->Pred : P->Suc;
      NodeClass *t = P->Reversed == Q->Reversed ? Q->Last : Q->First;

      //===== P =====  ===== Q =====
      //Last --- First First ---- Last 

      int i = t->Rank;
      if(t == Q->Last)
      {
	if( t == Q->First && t->Suc != P->First)
	{
	  Q->Reversed ^= 1;
	  u = t->Suc;
	  t->Suc = t->Pred;
	  t->Pred = u;
	}

	for(NodeClass *t = P->First; t != t2; t = t->Suc)
	{
	  t->Parent = Q;
	  t->Rank = ++i;
	}

	Q->Last = t1;
      }
      else 
      {
	for(NodeClass *t = P->First; t != t2; t = u)
	{
	  t->Parent = Q;	  
	  t->Rank = --i;
	  u = t->Suc; 
	  t->Suc = t->Pred;
	  t->Pred = u;
	}
	Q->First = t1;
      }
      P->First = t2;
    }
    else 
    {
      Q = !P->Reversed ? P->Suc : P->Pred;
      NodeClass *t = P->Reversed == Q->Reversed ? Q->First : Q->Last;
      int i = t->Rank;
      //===== P ========  ====== Q =====
      //First ---t- Last  First ---- Last 
      //===== Q =======  ====== P =========
      //First ---- Last  Last--t ----- First 

      if(t == Q->First)
      {
	if(t == Q->Last && t->Pred != P->Last)
	{
	  u = t->Suc;                                       
	  t->Suc = t->Pred;                                 
	  t->Pred = u;
	  Q->Reversed ^= 1;
	}

	for (NodeClass *t = P->Last; t != t1; t = t->Pred) 
	{
	  t->Parent = Q;
	  t->Rank = --i;
	}
	Q->First = t2;
      }
      else 
      {
	//for(NodeClass *t = t2; t != P->Last; t = t->Suc)
	for(NodeClass *t = P->Last; t != t1; t = u)
	{
	  u = t->Pred;
	  t->Pred = t->Suc;
	  t->Suc = u;
	  t->Rank = ++i;
	  t->Parent = Q;
	}
	Q->Last = t2;

      }

      Count = P->Size - Count;
      P->Last = t1;
    }

    P->Size -= Count;
    Q->Size += Count;

    return;
  }

  bool LKSolver::Flipable(const int kOpt, vector<SegmentClass> &Segs, vector<SwapRecord> &FlipStacks)
  {
    SegmentClass *First = &Segs[0];
    bool feasible = true, restart = false;
    int n = 0, b[11] = {0, 2*kOpt,3,2,5,4,7,6,9,8,1};
    b[2*kOpt] = 1;

RESTART:

    FlipStacks.clear();

    for(size_t i = 0; i < kOpt; ++i)
    {  Segs[i].Reversed = 0; Segs[i].Done = false; }
    for(size_t i = 0; i < kOpt; ++i)
      Segs[i % kOpt].Suc = &Segs[(i + 1) % kOpt];
    for(size_t i = kOpt; i > 0; --i)
      Segs[i % kOpt].Pred = &Segs[(i - 1) % kOpt]; 

    //start: sub tour check
    if(!restart)
    {

      SegmentClass *Parent[kOpt*2];
      int Mark[kOpt*2];

      for(size_t i = 0;  i < kOpt; ++i)
      {
	Parent[Segs[i].GetFrom()-1] = &Segs[i];
	Parent[Segs[i].GetTo()-1] = &Segs[i];
	Mark[2*i] = Mark[2*i+1] = 0;
      }

      int p = 1;
      SegmentClass *Par = Parent[p];
      Mark[p-1] = 1;
      ++n;
      do
      {
	if(p == Par->GetFrom())
	{
	  if(!Mark[Par->GetPred()->GetTo()-1])
	  {  Par = Par->GetPred(); p = Par->GetTo(); }
	  else if(!Mark[b[p]-1])
	  { Par = Parent[b[p]-1]; p = b[p]; }
	  else 
	    break;
	}
	else 
	{
	  if(!Mark[Par->GetSuc()->GetFrom()-1])  
	  { Par = Par->GetSuc(); p = Par->GetFrom(); }
	  else if(!Mark[b[p]-1])
	  { Par = Parent[b[p]-1]; p = b[p]; }
	  else 
	    break;
	}

	Mark[p-1] = 1; 
      }
      while(++n);

      feasible = n == 2*kOpt ? true : false;

      if(!feasible)
	return feasible;
    }
    //end : sub tour chek

    SegmentClass *s1 = First;

    do
    {
      if(s1->Done)
	continue;
      if(s1->GetFrom() == b[s1->GetTo()])
      { s1->Done = true; continue; }

      SegmentClass *s2 = !restart ? s1->GetSuc() : s1->GetPred() ;
      do
      {
	bool case1=false, case2=false;
	if(((case1 = (s2->GetFrom() == b[s1->GetFrom()])) || (case2 = (s2->GetTo() == b[s1->GetTo()]))) || restart)
	{
	  restart = false;
	  SwapRecord record; 
	  record.id1 = s1->GetFrom();
	  record.id2 = s1->GetTo();
	  record.id3 = s2->GetTo();
	  FlipStacks.push_back(record);

	  int p = s1->GetTo();
	  s1->SetTo(s2->GetFrom());
	  s2->SetFrom(p);
	  if(case1)
	    s1->Done = true;
	  if(case2)
	    s2->Done = true;

	  if(s2 == s1->GetSuc())
	  { s2 = s1 = First; continue; }

	  SegmentClass *s4;
	  for(SegmentClass *s3 = s1->GetSuc(); s3 != s2; s3 = s4)
	  {
	    s4 = s3->GetSuc();
	    s3->Reversed ^= 1;
	  }

	  SegmentClass *SUC1 = s1->GetSuc(), *PRED2 = s2->GetPred();

	  s1->SetSuc(PRED2);
	  s2->SetPred(SUC1);

	  SUC1 = s1->GetSuc();
	  SUC1->SetPred(s1);
	  PRED2 = s2->GetPred();
	  PRED2->SetSuc(s2);

	  s2 = s1 = First;

	}
      }
      while((s2 = s2->GetSuc()) != First);
    }
    while((s1 = s1->GetSuc()) != First);

    s1 = First;
    feasible = true;
    do
    { 
      feasible *= s1->Done; 
      //cout << s1->GetFrom() << s1->GetTo() << "-"; 
    }
    while((s1 = s1->GetSuc()) != First);

    if(!feasible)
    { restart = true; goto RESTART; }

    return feasible;
  }

  bool LKSolver::Flipable2(const int kOpt, SegmentClass *First, vector<SwapRecord> &FlipStacks)
  {
    bool feasible = true, restart = false;
    int n = 0, b[2*kOpt + 1], SeqCnt = 0, NonSeqCnt = 0; 
    SegmentClass *Seg = First;
    do
    {
      if(Seg->NonSequence) 
	++NonSeqCnt;
      else 
	++SeqCnt;
    }
    while((Seg = Seg->GetSuc()) != First);


    //make flip from-to vertex indexs
    //b[0] = 0; b[1] = 2*kOpt; b[2*kOpt] = 1;
    //for(int i = 2; i < 2*kOpt; i += 2)
    //{ b[i] = i + 1; b[i+1] = i; }

    b[0] = 0; 
    if(NonSeqCnt)
    { 
      b[1] = 2*NonSeqCnt; b[2*NonSeqCnt] = 1;
      for(int i = 2; i < 2*NonSeqCnt; i += 2)
      { b[i] = i + 1; b[i+1] = i; }
    }
    if(SeqCnt)
    {
      b[2*NonSeqCnt + 1] = 2 * kOpt; b[2*kOpt] = 2*NonSeqCnt + 1;
      for(int i = 2 + 2*NonSeqCnt; i < 2*kOpt; i += 2)
      { b[i] = i + 1; b[i+1] = i; }
    }

    FlipStacks.clear();
    do
    {
      Seg->Reversed = 0;
      Seg->Done = false;
      Seg->FOld = Seg->F;
      Seg->TOld = Seg->T;
      Seg->SucOld = Seg->Suc;
      Seg->PredOld = Seg->Pred;
    }
    while((Seg = Seg->GetSuc()) != First);

    //start: sub tour check
    if(!restart)
    {

      SegmentClass *Parent[kOpt*2];
      int Mark[kOpt*2];

      Seg = First; 
      int i = 0;
      do
      {
	Parent[Seg->GetFrom()-1] = Seg;
	Parent[Seg->GetTo()-1] = Seg;
	Mark[2*i] = Mark[2*i+1] = 0;
	++i;
      }
      while((Seg = Seg->GetSuc()) != First);

      int p = 1;
      SegmentClass *Par = Parent[p];
      Mark[p-1] = 1;
      ++n;
      do
      {

	if(p == Par->GetFrom())
	{
	  if(!Mark[Par->GetPred()->GetTo()-1])
	  {  Par = Par->GetPred(); p = Par->GetTo(); }
	  else if(!Mark[b[p]-1])
	  { Par = Parent[b[p]-1]; p = b[p]; }
	  else 
	    break;
	}
	else 
	{
	  if(!Mark[Par->GetSuc()->GetFrom()-1])  
	  { Par = Par->GetSuc(); p = Par->GetFrom(); }
	  else if(!Mark[b[p]-1])
	  { Par = Parent[b[p]-1]; p = b[p]; }
	  else 
	    break;
	}

	Mark[p-1] = 1; 
      }
      while(++n);

      long int LogKey = GetFlipLogKey(First);
      feasible = n == 2*kOpt ? true : false;

      if(!feasible)
	return feasible;
    }
    //end : sub tour chek

    SegmentClass *s1 = First;

RESTART:
    do
    {
      if(s1->Done)
	continue;
      if(s1->GetFrom() == b[s1->GetTo()])
      { s1->Done = true; continue; }

      //SegmentClass *s2 = !restart ? s1->GetSuc() : s1->GetPred() ;
      SegmentClass *s2 = !restart ? s1->GetSuc() : First->GetPred();
      do
      {
	bool case1=false, case2=false;
	if(((case1 = (s2->GetFrom() == b[s1->GetFrom()])) || (case2 = (s2->GetTo() == b[s1->GetTo()]))) || restart)
	{
	  restart = false;
	  SwapRecord record; 
	  record.id1 = s1->GetFrom();
	  record.id2 = s1->GetTo();
	  record.id3 = s2->GetTo();
	  FlipStacks.push_back(record);

	  int p = s1->GetTo();
	  s1->SetTo(s2->GetFrom());
	  s2->SetFrom(p);
	  if(case1)
	    s1->Done = true;
	  if(case2)
	    s2->Done = true;

	  if(s2 == s1->GetSuc())
	  { s2 = s1 = First; continue; }

	  SegmentClass *s4;
	  for(SegmentClass *s3 = s1->GetSuc(); s3 != s2; s3 = s4)
	  {
	    s4 = s3->GetSuc();
	    s3->Reversed ^= 1;
	  }

	  SegmentClass *SUC1 = s1->GetSuc(), *PRED2 = s2->GetPred();

	  s1->SetSuc(PRED2);
	  s2->SetPred(SUC1);

	  SUC1 = s1->GetSuc();
	  SUC1->SetPred(s1);
	  PRED2 = s2->GetPred();
	  PRED2->SetSuc(s2);

	  s2 = s1 = First;

	}
      }
      while((s2 = s2->GetSuc()) != First);
    }
    while((s1 = s1->GetSuc()) != First);

    feasible = true;
    s1 = First;
    do
      feasible *= s1->Done; 
    while((s1 = s1->GetSuc()) != First && feasible);

    if(!feasible)
    { restart = true; goto RESTART; }

    s1 = First;
    do
    {
      s1->Reversed = 0;
      s1->F = s1->FOld;
      s1->T = s1->TOld;
      s1->Suc = s1->SucOld;
      s1->Pred = s1->PredOld;
    }
    while((s1 = s1->GetSuc()) != First);

    return feasible;
  }

  void LKSolver::MakeKoptFlip()
  {

    SegmentClass *First = &FlipSegs[1];
    First->SetPred(First);
    First->SetSuc(First);

    First->SetFrom(1);
    First->SetTo(2);

    First->NonSequence = false;
    MakeKoptFlipR(12, 2, First, false);
    cout << "end sequnce" << endl;
    First->NonSequence = true;
    MakeKoptFlipR(-12, 2, First, true);
    cout << "end nonsequnce" << endl;

#if 1
    for(auto Ite = FlipHash.begin(); Ite != FlipHash.end(); ++Ite)
    {
       const long int BaseKey = Ite->first;
       FlipLogClass *aLog = &Ite->second;

#if 0
       for(int X = 0; X <= 1; ++X)
       {
         vector<NextFlipClass> *aVec = aLog->GetNextBetweens(X); 

	 for(size_t i = 0; i < aVec->size(); ++i)
	 {
	   NextFlipClass *NextFlip = &aVec->at(i);
	   bool found = false;

	   vector<NextFlipClass> *cands = aLog->GetNextBetweens(); 

	   for(size_t j = 0; j < cands->size(); ++j)
	   {
	     NextFlipClass *NextFlip2 = &cands->at(j);

	     if(NextFlip2->GetFrom() == NextFlip->GetFrom() && 
		 NextFlip2->GetTo() == NextFlip->GetTo())
	     {
	       found = true;
	       NextFlip2->Xks[X] = X;
	       NextFlip2->Feasibles[X] = FlipHash[NextFlip->GetKey()].IsFeasibleLog();
	       NextFlip2->LogKeys[X] = NextFlip->GetKey();

	       break;
	     }

	   }

	   if(!found)
	   {
	     NextFlipClass cand = NextFlipClass(NextFlip->GetKey(), NextFlip->Feasible(), X, NextFlip->GetFrom(), NextFlip->GetTo()); 
	     cand.Xks[X] = X;
	     cand.Feasibles[X] = FlipHash[NextFlip->GetKey()].IsFeasibleLog();
	     cand.LogKeys[X] = NextFlip->GetKey();
	     aLog->AddBtwCond(cand);
	   }
	   
	 }
       }
#endif 
#if 1
       vector<NextFlipClass> *cands = aLog->GetNextBetweens(); 
       for(size_t j = 0; j < cands->size(); ++j)
       {
	 NextFlipClass *NextFlip2 = &cands->at(j);

	 int *Xs = NextFlip2->Xks;
	 bool *feas = NextFlip2->Feasibles;

	 if((Xs[0] == -1 && (Xs[1] >= 0 && feas[1])) || 
	    (Xs[0] >= 0 && !feas[0] && Xs[1] >= 0 && feas[1]))
	 {
	   auto tmp1 = NextFlip2->Xks[1];
	   NextFlip2->Xks[1] = NextFlip2->Xks[0];
	   NextFlip2->Xks[0] = tmp1;

	   auto tmp2 = NextFlip2->Feasibles[1];
	   NextFlip2->Feasibles[1] = NextFlip2->Feasibles[0];
	   NextFlip2->Feasibles[0] = tmp2;

	   auto tmp3 = NextFlip2->LogKeys[1];
	   NextFlip2->LogKeys[1] = NextFlip2->LogKeys[0];
	   NextFlip2->LogKeys[0] = tmp3;
	 }

       }
#endif 
    }

    for(auto Ite = FlipHash.begin(); Ite != FlipHash.end(); ++Ite)
    {
      const long int BaseKey = Ite->first;
      //if(BaseKey < 0)
      if(BaseKey > 0)
	continue;
      if(BaseKey < -100000000)
	continue;
      FlipLogClass *aLog = &Ite->second;
      vector<NextFlipClass> *cands = aLog->GetNextBetweens(); 

      cout << "BASEKEY = " << BaseKey << " my Feasible = " << (aLog->IsFeasibleLog()) << ", ";
      //if(BrigeHash.find(BaseKey) != BrigeHash.end())
	//cout << (BrigeHash[BaseKey].IsFeasibleLog()) << endl;
      if(BrigeHash.find(2*BaseKey) != BrigeHash.end())
	cout << (BrigeHash[BaseKey*2].IsFeasibleLog()) << endl;
      else 
	cout << "NOT FOUND" << endl;

      for(size_t j = 0; j < cands->size(); ++j)
       {
	 NextFlipClass *NextFlip2 = &cands->at(j);

	 /*
	 cout << "basekey = " << BaseKey << " " 
	      << "from, to = " << NextFlip2->GetFrom() << " " << NextFlip2->GetTo() << " "
	      << "Xs = " <<  NextFlip2->Xks[0] << " " << NextFlip2->Xks[1]  << " "
	      << "feas = " << NextFlip2->Feasibles[0] << " " << NextFlip2->Feasibles[1] << " "
	      << "kes = " << NextFlip2->LogKeys[0] << " " << NextFlip2->LogKeys[1] << endl;
*/
       }
       //if(cands->size() == 0)
	 //cout << "basekey = " << BaseKey << " " << (aLog->IsFeasibleLog()) << endl;
    }

    myGetChar();
#endif 
    return;
  }

  bool LKSolver::MakeKoptFlipR(const long int ParentKey, const int k, SegmentClass *First, const bool Brige, const int BrigeDepth)
  {
    bool rval = false;
    SegmentClass *Seg = First;
    vector<SwapRecord> FlipStacks;

    do
    {
      SegmentClass *NewSeg = &FlipSegs[k];    

      NewSeg->NonSequence = Brige;
      NextFlipClass cand(0, false, -1, Seg->GetTo(), Seg->GetSuc()->GetFrom());

      for(int X = 1; X <= 2; ++X)
      {
	//About New Seg..
	NewSeg->SetFrom(X == 1 ? 2*k-1 : 2*k);
	NewSeg->SetTo(X == 1 ? 2*k : 2*k-1);
	NewSeg->NonSequence = Brige && BrigeDepth == 0 ? Brige : false;

	SegmentClass *tmp1 = Seg->GetSuc();

	//insert new set between Seg to SUC(Seg)
	Seg->SetSuc(NewSeg);
	NewSeg->SetSuc(tmp1);

	tmp1->SetPred(NewSeg);
	NewSeg->SetPred(Seg);

	FlipStacks.clear();
	long int LogKey = 0;
	bool feasible = false, ValidFlip = false;

	LogKey = GetFlipLogKey(First);
	feasible = Flipable2(k, First, FlipStacks);

	if((Brige && !feasible) || (1 <= BrigeDepth && BrigeDepth <= 2))
	{
	  if(!BrigeDepth)
	    InsertFlipLog(LogKey, feasible, FlipStacks);
	  else 
	    InsertBrigeFlipLog(LogKey*BrigeDepth, feasible, FlipStacks);


	  if(BrigeDepth < 2 && k+1 <= kopt)
	    ValidFlip = MakeKoptFlipR(LogKey, k+1, First, Brige, BrigeDepth+1);

	  if(feasible || ValidFlip)
	  {
	    const int x = X - 1;
	    cand.Xks[x] = x;
	    cand.Feasibles[x] = BrigeHash[LogKey*BrigeDepth].IsFeasibleLog(); //FlipHash[LogKey].IsFeasibleLog();
	    cand.LogKeys[x] = LogKey*BrigeDepth; 
	    //BrigeHash[ParentKey].Add_NextFlip(LogKey*BrigeDepth, BrigeHash[LogKey*BrigeDepth].IsFeasibleLog(), X-1, Seg->GetTo(), tmp1->GetFrom());
	  }

	  if(1 <= BrigeDepth)
	    goto NEXTLOOP; 
	}
	else 
	  InsertFlipLog(LogKey, feasible, FlipStacks);

	if(feasible)
	  ValidFlip = true;

	if(k + 1 <= kopt)
	  if(MakeKoptFlipR(LogKey, k+1, First, Brige))
	    ValidFlip = true;

	if(ValidFlip && ParentKey != 0) 
	{
	  const int x = X - 1;
	  cand.Xks[x] = x;
	  cand.Feasibles[x] = FlipHash[LogKey].IsFeasibleLog();
	  cand.LogKeys[x] = LogKey; 
	  //FlipHash[ParentKey].Add_NextFlip(LogKey, FlipHash[LogKey].IsFeasibleLog(), X-1, Seg->GetTo(), tmp1->GetFrom());
	}

NEXTLOOP:;

	 //remove new seg 
	 NewSeg->GetSuc()->SetPred(NewSeg->GetPred());
	 NewSeg->GetPred()->SetSuc(NewSeg->GetSuc());
	 NewSeg->SetSuc(NULL);
	 NewSeg->SetPred(NULL);

	 if(ValidFlip)
	   rval = true;

      }

      if(cand.Xks[0] >= 0 || cand.Xks[1] >= 0)
      {
	if(!BrigeDepth)
	  FlipHash[ParentKey].AddBtwCond(cand);
	else 
	  BrigeHash[ParentKey].AddBtwCond(cand);
      }

    }
    while((Seg = Seg->GetSuc()) != First);

    return rval;
  }


};
