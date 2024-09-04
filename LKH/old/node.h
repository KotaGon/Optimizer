#ifndef _NODE_H
#define _NODE_H

#include "LK.h"

namespace myOptimizeName
{

  /* Macro definitions */
#define Fixed(a,b) ((a)->FixedTo1 == (b) || (a)->FixedTo2 == (b))
#define Follow(b,a)\
  ((a)->Suc != (b) ?\
   Link((b)->Pred,(b)->Suc), Link(b,(a)->Suc), Link(a,b) : 0) 
#define InBestTour(a,b) ((a)->BestSuc == (b) || (b)->BestSuc == (a))
#define InNextBestTour(a,b) ((a)->NextBestSuc == (b) || (b)->NextBestSuc == (a))
#define Link(a,b) ((a)->Suc = (b), (b)->Pred = (a))
#define Near(a,b) ((a)->BestSuc ? InBestTour(a,b) : (a)->Dad == (b) || (b)->Dad == (a))
#define InOptimumTour(a,b) ((a)->OptimumSuc == (b) || (b)->OptimumSuc == (a))
#define IsCommonEdge(a,b) (((a)->MergeSuc[0] == (b) || (b)->MergeSuc[0] == (a)) &&\
    ((a)->MergeSuc[1] == (b) || (b)->MergeSuc[1] == (a)))
#define Precede(a,b)\
  ((b)->Pred != (a) ?\
   Link((a)->Pred,(a)->Suc), Link((b)->Pred,a), Link(a,b) : 0)

#define TWO_LEVEL_TREE
  
#ifdef TWO_LEVEL_TREE
#define PRED(a) (Reversed == (a)->Parent->Reversed ? (a)->Pred : (a)->Suc)
#define SUC(a) (Reversed == (a)->Parent->Reversed ? (a)->Suc : (a)->Pred)
#define BETWEEN(a, b, c) Between_SL(a, b, c)
#define FLIP(a, b, c, d) Flip_SL(a, b, c)
#endif

#define Swap1(a1,a2,a3)\
        FLIP(a1,a2,a3,0)
#define Swap2(a1,a2,a3, b1,b2,b3)\
        (Swap1(a1,a2,a3), Swap1(b1,b2,b3))
#define Swap3(a1,a2,a3, b1,b2,b3, c1,c2,c3)\
        (Swap2(a1,a2,a3, b1,b2,b3), Swap1(c1,c2,c3))
#define Swap4(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3)\
        (Swap3(a1,a2,a3, b1,b2,b3, c1,c2,c3), Swap1(d1,d2,d3))
#define Swap5(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3, e1,e2,e3)\
        (Swap4(a1,a2,a3, b1,b2,b3, c1,c2,c3, d1,d2,d3), Swap1(e1,e2,e3))

  class LKSolver;
  class SegmentClass;
  class CandidateClass;

  class NodeClass
  {
    friend class LKSolver; 

    private:
    public:
      NodeClass() = default;
      int Id;
      int V;
      int LastV;
      int Pi;
      int BestPi;
      int *C;
      int Cost;
      int Rank;
      int NextCost;
      int OldPredExcluded, OldSucExcluded;
      int Loc;
      int SelectCnt;
      double X,Y,Z;
      NodeClass *Pred, *Suc;
      NodeClass *OldPred, *OldSuc;
      NodeClass *Dad;
      NodeClass *Next;
      NodeClass *OptimumSuc;
      NodeClass *NextBestSuc;
      NodeClass *BestSuc;
      NodeClass *PreSide, *PostSide;
      SegmentClass *Parent, *ParentFlip;
      std::vector<CandidateClass> CandidateSet;
      std::multimap<int, NodeClass*>::iterator LocItr;
  };
  using locItr = std::multimap<int, NodeClass*>::iterator;
  
  class SegmentClass
  {
    private:
      
    public:

      SegmentClass()
      {
        Rank = Size = F = T = FOld = TOld = Done = Reversed = 0;
	First = Last = NULL;
	PredOld = SucOld = Pred = Suc = NULL;
      };
      int Reversed;
      NodeClass *First, *Last;
      SegmentClass *Pred, *Suc;
      SegmentClass *PredOld, *SucOld;
      int Rank;
      int Size;
      int F, T;
      int FOld, TOld;
      int Done;
      bool NonSequence;

      inline int GetFrom(){ return (!Reversed ? F : T); }
      inline int GetTo(){ return (!Reversed ? T : F); }
      inline void SetFrom(const int _F){ (!Reversed ? F : T) = _F; }
      inline void SetTo(const int _T){ (!Reversed ? T : F) = _T; }

      inline SegmentClass *GetSuc(){ return !Reversed ? Suc : Pred; }
      inline SegmentClass *GetPred(){ return !Reversed ? Pred : Suc; }
      inline void SetSuc(SegmentClass *s){ (!Reversed ? Suc : Pred) = s; }
      inline void SetPred(SegmentClass *s) { (!Reversed ? Pred : Suc) = s; }

  };
  
  class CandidateClass
  {
    private:
    public:
      NodeClass *To;
      int Cost;
      int Alpha;

      inline bool operator < (const CandidateClass &rhs) const
      { return Alpha != rhs.Alpha ? Alpha < rhs.Alpha : Cost < rhs.Cost; }
  };
};

#endif
