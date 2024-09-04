#include <iostream>
#include <assert.h>
#include <vector>

#define OPT (4)

using namespace std;

class SegmentClass
{
  public:
    SegmentClass(const int _id, const int _s, const int _e)
    { id = _id; s = _s; e = _e; inf = false; reversed = 0;}
    int id;
    int s, e;
    int inf;
    int reversed;

    int getStart(){ return (!reversed ? s : e); }
    int getEnd(){ return (!reversed ? e : s); }
    void setStart(const int _s) { (!reversed ? s : e) = _s; }
    void setEnd(const int _e) { (!reversed ? e : s) = _e; }

    SegmentClass *Pred, *Suc;

    SegmentClass *SUC(){ return !reversed ? Suc : Pred; }
    SegmentClass *PRED(){ return !reversed ? Pred : Suc; }
    void setSUC(SegmentClass *s){ (!reversed ? Suc : Pred) = s; }
    void setPred(SegmentClass *s) { (!reversed ? Pred : Suc) = s; }

};

void display(SegmentClass *FirstSegment)
{
  int loop = 0;

  SegmentClass *s = FirstSegment;

  do
  {
    cout << s->getStart() << s->getEnd() << "(" << s->id << ")-"; 
  }
  while((s = s->SUC()) != FirstSegment);
  cout << endl;

  return;
}

int main()
{

  vector<SegmentClass> Segs;
  SegmentClass *FirstSegment;

#if 1
  int a[10] = {
    1,2,
    5,6,
    8,7,
    4,3,
    10,9,
  };
#endif 
#if 0
int a[10] = {
    1,2,
    9,10,
    8,7,
    3,4,
    6,5,
  };
#endif 
#if 0
int a[10] = {
    1,2,
    6,5,
    4,3,
    //9,10,
    -1,-1,
    -1,-1,
  };
#endif 
  //p5:12-910-78-56-43
    //p6:12-43-78-65-
  for(int i = 0; i < OPT; ++i)
    Segs.push_back(SegmentClass(i, a[2*i], a[2*i+1]));

  for(int i = 0; i < OPT; ++i)
    Segs[i % OPT].Suc = &Segs[(i + 1) % OPT];

  for(int i = OPT; i > 0; --i)
    Segs[i % OPT].Pred = &Segs[(i - 1) % OPT];

  FirstSegment = &Segs[0];

  //p1:12-910-87-34-65-
  //p2:12-910-34-87-65-
  //p3::12-87-910-34-65-
//p4:12-65-109-78-43-

  int b[11] = {0, 
    2*OPT,3,
    2,5,
    4,7,
    6,9,
    8,1
  }; b[2*OPT] = 1;

  display(FirstSegment);
  
  do
  {
    SegmentClass *s1 = FirstSegment;

    do
    {
      if(s1->inf)
	continue;

      if(s1->getStart() == b[s1->getEnd()])
      { s1->inf = true; continue; }

      SegmentClass *s2 = s1->Suc;
      do
      {
	bool case1, case2;
	if((case1 = (s2->getStart() == b[s1->getStart()])) ||
	   (case2 = (s2->getEnd() == b[s1->getEnd()])))
	{

	  int p = s1->getEnd();
	  s1->setEnd(s2->getStart());
	  s2->setStart(p);
	  if(case1)
	    s1->inf = true;
	  if(case2)
	    s2->inf = true;

	  if(s2 == s1->SUC())
	    goto FLIP;

	  SegmentClass *s4;
	  for(SegmentClass *s3 = s1->SUC(); s3 != s2; s3 = s4)
	  {
	    s4 = s3->SUC();
	    s3->reversed ^= 1;
	  }

	  SegmentClass *SUC1 = s1->SUC(), *PRED2 = s2->PRED();

	  s1->setSUC(PRED2);
	  s2->setPred(SUC1);

	  SUC1 = s1->SUC();
	  SUC1->setPred(s1);
	  PRED2 = s2->PRED();
	  PRED2->setSUC(s2);

	  s1=FirstSegment;
	  s2=FirstSegment;
	  goto FLIP;
	}
      }
      while((s2 = s2->SUC()) != FirstSegment);
    }
    while((s1 = s1->SUC()) != FirstSegment);

    break;

FLIP:;

     display(FirstSegment);

  }
  while(1);

  bool feasible = true; 

  SegmentClass *s1 = FirstSegment;
  do
  {
   cout << s1->getStart() << s1->getEnd() << "-"; 
   feasible *= s1->inf;
  }
  while((s1 = s1->SUC()) != FirstSegment);
  cout << endl;
      
  cout << (feasible ? "feasible" : "infeasible") << endl;

  return 0;
}
