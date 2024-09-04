#ifndef _FLIP_H
#define _FLIP_H

#include "LK.h"

namespace myOptimizeName
{
  class SwapRecord;

  class NextFlipClass
  {
    private:
      long int LogKey;
      bool Feas;
      int X;
      int From, To;

    public:
      NextFlipClass() = default;
      NextFlipClass(long int _LogKey, bool _Feas, int _X, int _From, int _To)
	: LogKey(_LogKey), Feas(_Feas), X(_X), From(_From), To(_To) 
      { 
	//Xks[0] = _X; Xks[1] = -1; 
	Xks[0] = Xks[1] = -1;
      }
      ~NextFlipClass(){};
      void clear(){ Xks[0] = Xks[1] = -1; }

      int Xks[2];
      long int LogKeys[2];
      bool Feasibles[2];

      inline long int GetKey(){ return LogKey; }
      inline int GetX() { return X; }
      inline int GetFrom(){ return From; }
      inline int GetTo() { return To; }
      inline bool Feasible(){ return Feas; }


      /*
      inline void SetXks(const int i, const int x)
      { Xks[i] = x; }
      inline int* GetXks(){ return Xks; }
      */
  };

  class FlipLogClass
  {
    private:

      long int Logkey;
      bool Feasible;
      std::vector<SwapRecord> FlipStacks;
      std::vector<NextFlipClass> Children[2];
      std::vector<NextFlipClass> BtwCands;

    public:

      FlipLogClass() = default;
      FlipLogClass(long int _Logkey, const bool _Feasible, std::vector<SwapRecord> &_Flips) 
	: Logkey(_Logkey), Feasible(_Feasible), FlipStacks(_Flips) {}

      inline void SetLogKey(long int value)
      { Logkey = value; }
      inline void SetFeasible(bool value)
      { Feasible = value; }
      inline void SetFlipStacks(std::vector<SwapRecord> &value)
      { FlipStacks = value; }

      inline bool IsFeasibleLog(){ return Feasible; }
      inline std::vector<SwapRecord> *GetFlipStacks(){ return &FlipStacks; }
      inline std::vector<NextFlipClass> *GetNextBetweens(const int X){ return &Children[X]; }
      inline std::vector<NextFlipClass> *GetNextBetweens(){ return &BtwCands; }

      inline void Add_NextFlip(const long int LogKey, const bool Feas, const int X, const int From, const int To)
      { Children[X].push_back(NextFlipClass(LogKey, Feas, X, From, To)); }
      inline void AddBtwCond(NextFlipClass value)
      {
        BtwCands.push_back(value);
      }

  };

};
#endif 
