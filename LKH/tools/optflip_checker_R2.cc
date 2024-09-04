#include <iostream>
#include <assert.h>
#include <vector>

#define OPT (4)

using namespace std;

  
bool feasible = true;
bool restart = false;

class SegmentClass
{
  public:
    SegmentClass(const int _id, const int _s, const int _e)
    { id = _id; s = _s; e = _e; Done = false; reversed = 0;}
    int id;
    int s, e;
    int Done;
    int reversed;

    int GetFrom(){ return (!reversed ? s : e); }
    int GetTo(){ return (!reversed ? e : s); }
    void SetFrom(const int _s) { (!reversed ? s : e) = _s; }
    void SetTo(const int _e) { (!reversed ? e : s) = _e; }

    SegmentClass *Pred, *Suc;

    SegmentClass *GetSuc(){ return !reversed ? Suc : Pred; }
    SegmentClass *GetPred(){ return !reversed ? Pred : Suc; }
    void SetSuc(SegmentClass *s){ (!reversed ? Suc : Pred) = s; }
    void SetPred(SegmentClass *s) { (!reversed ? Pred : Suc) = s; }

};

void display(SegmentClass *FirstSegment)
{
  int loop = 0;
  SegmentClass *s = FirstSegment;
  do
  {
  
    cout << s->GetFrom() << s->GetTo() << "(" << s->id << ")-"; 
    ++loop;
    if(loop >= OPT + 1)
    { cout << "stop" << endl; getchar(); }
  }
  while((s = s->GetSuc()) != FirstSegment);
  cout << endl;
  return;
}

void init_ab(int *a, int *b)
{

  for(int i = 0; i < 2*OPT; i += 2)
  { 
    a[i] = i + 1;
    a[i+1] = i + 2;
  }

  a[0] = 1;
  a[1] = 2;
  a[2] = 5;
  a[3] = 6;
  a[4] = 3;
  a[5] = 4;
  a[6] = 7;
  a[7] = 8;

  /*
  b[0] = 0;
  b[1] = 2*OPT;
  for(int i = 2; i < 2*OPT; i += 2)
  {
   b[i] = i + 1;
   b[i+1] = i;  
  }
  b[2*OPT] = 1;*/

  b[0] = 0;
  b[1] = 4;
  b[2] = 3;
  b[3] = 2;
  b[4] = 1;
  b[5] = 7;
  b[6] = 8;
  b[7] = 5;
  b[8] = 6;
}

void swaps(vector<SegmentClass> &Segs)
{

  for(int i = 0; i < 100; ++i)
  {
    int j = rand() % Segs.size();
    int k = rand() % Segs.size();

    if(j == 0 || k == 0)
      continue;

    SegmentClass tmp = Segs[j];
    Segs[j] = Segs[k];
    Segs[k] = tmp;

    double rnd = (double) rand() / RAND_MAX;

    if(rnd < 0.5)
    {
      int tmp1 = Segs[j].s;
      Segs[j].s = Segs[j].e;
      Segs[j].e = tmp1;
    }
  }
}

bool getFeasibility(vector<SegmentClass> &Segs, int *b)
{

  if(!restart)
  {
    SegmentClass *Parent[OPT*2];
    int n = 0, Mark[OPT*2];

    for(size_t i = 0; i < Segs.size(); ++i)
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

    cout << n << endl;
    feasible = n == 2*OPT ? true : false;
    cout << (feasible == true ? "feasible" : "infeasible") << endl;

    if(!feasible)
      return feasible;
  }
  else 
    return true;
  return true;
}

// main //
int main()
{ 
  vector<SegmentClass> Segs;
  Segs.reserve(100);

  int *a = new int[2*OPT], *b = new int[2*OPT+1];

  //init random
  srand((unsigned int)time(NULL));

  init_ab(a, b);

  Segs.clear();
  for(int i = 0; i < OPT; ++i)
    Segs.push_back(SegmentClass(i, a[2*i], a[2*i+1]));

  SegmentClass *First = &Segs[0];
  SegmentClass *s1 = First;

  if(!restart)
  { 
    //swaps(Segs);
    for(int i = 0; i < OPT; ++i)
    {
      a[2*i] = Segs[i].GetFrom();
      a[2*i + 1] = Segs[i].GetTo(); 
    }
  }

  for(int i = 0; i < OPT; ++i)
    Segs[i % OPT].Suc = &Segs[(i + 1) % OPT];
  for(int i = OPT; i > 0; --i)
    Segs[i % OPT].Pred = &Segs[(i - 1) % OPT];

  display(First);

  //start: sub tour check
  feasible = getFeasibility(Segs, b);
  if(!feasible)
    goto CLEANUP;
  //end : sub tour chek

RESTART:

  do
  {
    if(s1->Done)
      continue;
    if(s1->GetFrom() == b[s1->GetTo()])
    { s1->Done = true; continue; }

    //SegmentClass *s2 = !restart ? s1->GetSuc() : s1->GetPred();
    SegmentClass *s2 = !restart ? s1->GetSuc() : First->GetPred();

    do
    {
      bool case1=false, case2=false;
      if(((case1 = (s2->GetFrom() == b[s1->GetFrom()])) || (case2 = (s2->GetTo() == b[s1->GetTo()]))) || restart)
      {
	restart = false;

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
	  s3->reversed ^= 1;
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
    feasible *= s1->Done; 
  while((s1 = s1->GetSuc()) != First && feasible);

  display(First);
  cout << (feasible == true ? "feasible" : "infeasible => wrong ans?") << endl;

  if(!feasible)
  { restart = true; goto RESTART; }

CLEANUP:

  delete [] b;
  delete [] a;

  return 0;
}

