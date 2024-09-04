#ifndef _NODE_H
#define _NODE_H
 
namespace myOptimizeName
{
#define Link(a,b) ((a)->Suc = (b), (b)->Pred = (a))
  
#define Follow(b,a)\
  ((a)->Suc != (b) ?\
   Link((b)->Pred,(b)->Suc), Link(b,(a)->Suc), Link(a,b) : 0) 
#define PRED(a) (!Reversed ? (a)->Pred : (a)->Suc)
#define SUC(a) (!Reversed ? (a)->Suc : (a)->Pred)

  class NodeClass
  {
    public:
      NodeClass() = default;
      int Id;
      int Rank;
      long *C;
      long X, Y;
      int Activate;
      NodeClass *Suc, *Pred;
  };

};

#endif
