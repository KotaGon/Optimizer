#ifndef _LK_H
#define _LK_H

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <map>
#include <unordered_map>
#include <cfloat>
#include <time.h>

#include "model.h"
#include "node.h"
#include "flip.h"
#include "utils.h"

using namespace std;

namespace myOptimizeName
{

  class LKDataModel;
  class LKSolver;
  class LKModel;
  class SwapRecord;
  class FlipLogClass;

  class SwapRecord
  {
     private:
     public:
       NodeClass *t1, *t2, *t3, *t4;
       int id1, id2, id3;
  };

  class LKDataModel : public DataModelClass
  {
    private:

      int Dimension;
      int GroupSize;
      int *CostMatrix;
      NodeClass *NodeSet, *FirstNode;
      SegmentClass *FirstSegment;

      LKModel *_myModel;

    public:
      LKDataModel() = default;
      LKDataModel(LKModel *model) :  _myModel(model) {}
      ~LKDataModel(){}

      void readProblem();
      void createNodes(int *ID, int *X, int *Y);

      inline int D(NodeClass *Na, NodeClass *Nb)
      { return (Na->Id < Nb->Id ? Nb->C[Na->Id] : Na->C[Nb->Id]) + Na->Pi + Nb->Pi; }
      inline int C(NodeClass *Na, NodeClass *Nb)
      { return Na->Id < Nb->Id ? Nb->C[Na->Id] : Na->C[Nb->Id]; }

      void free()
      { 
	delete [] NodeSet; 
	SegmentClass *s = FirstSegment, *snext = s->Suc;
	FirstSegment->Pred->Suc = NULL;

	do
	{ 
	  snext = s->Suc;
  	  delete s;
	  s = snext;
	}
	while(snext != NULL);
	//while(snext != FirstSegment);
      }
     
      friend class LKModel;
      friend class LKSolver;

  };

  class LKSolver : public OptModelClass
  {

    private:
      int Runs;
      int Seed;
      int MaxTrials;
      int MaxSwaps;
      int MaxCandidates;
      int AscentCandidates;
      int InitialPeriod;
      int InitialStepSize;
      int Precision;
      int Subgradient;
      double Excess;
      int RestrictedSearch;
      int Swaps;
      int Groups;
      int Norm;
      int Reversed;
      int OutTour;
      long int ParentFlipKey;
      queue<NodeClass*> ActiveNodes;
      unordered_map<int, int> Hash;
      unordered_map<long int, FlipLogClass> FlipHash;
      unordered_map<long int, FlipLogClass> BrigeHash;
      map<pair<NodeClass*, NodeClass*>, int> PreEdges, PostEdges;

      SwapRecord *SwapStack;
      int kopt;
      SegmentClass *FlipSegs; 
      NodeClass **Heap;

      LKModel *_myModel;

      double check_time = 0;

    public:

      LKSolver() = default;
      LKSolver(LKModel*model) : _myModel(model) {}
      ~LKSolver(){}

      void SiftUp(NodeClass * N);
      void MakeHeap(const long Size);
      NodeClass *DeleteMin();
      void Insert(NodeClass * N);

      void CreateCandidateSet();
      double Ascent();
      double Minimum1TreeCost();
      void MinimumSpanningTree();
      void GenerateCandidates(const long MaxCandidates, const long MaxAlpha, const int Symmetric=1);
      void Connect(NodeClass *N1, int Max);
      void RecordBetterTour();
      void RecordBestTour();
      void RunLoop();
      double FindTour();
      void CreateInitialTour();
      double LinKernighan();
      void AdjustCandidateSets();
      void ResetCandidateSets();
      void StoreTour();
      void RestoreTour();
      void NormalizeNodeList();
      void Flip_SL(NodeClass *t1, NodeClass *t2, NodeClass *t3);
      int Excludable(const NodeClass *ta, const NodeClass *tb);
      void Exclude(NodeClass *ta, NodeClass *tb);
      void SegmentSplit(NodeClass *t1, NodeClass *t2);
      bool Between_SL(NodeClass *t1, NodeClass *t2, NodeClass *t3);
      bool Flipable(const int kOpt, vector<SegmentClass> &Segs, vector<SwapRecord> &FlipStacks);
      bool Flipable2(const int kOpt, SegmentClass *First, vector<SwapRecord> &FlipStacks);
      bool find_(vector<NodeClass*> &Nodes, const int until, NodeClass *Node);
      int Gain23();
      int BrigeMove(vector<NodeClass*> &Nodes, SegmentClass *FirstSegment, long int ParentKey, int G, const int k);

      void InsertFlipLog(const long int Logkey, const bool feasible, vector<SwapRecord> &Flips);
      void InsertBrigeFlipLog(const long int Logkey, const bool feasible, vector<SwapRecord> &Flips);
      long int GetFlipLogKey(SegmentClass *FirstSegment);
      bool SearchFlipLogKey(long int hashkey);

      NodeClass* (LKSolver::*BestMove) (NodeClass *t1, NodeClass *t2, int *G0, int *Gain);
      NodeClass* BestKOptMove (NodeClass *t1, NodeClass *t2, int *G0, int *Gain);
      //NodeClass* BestKOptMoveR(int k, vector<NodeClass*> &Nodes, vector<NodeClass*> &BestNodes, SegmentClass *FirstFlipSegment, 
	//  vector<SwapRecord> &BestFlipStacks, int Gk, int *BestGk, int *G0, int *Gain, bool bridge = false, const int klimit=INT_MAX);

      NodeClass* BestKOptMoveR(int k, vector<NodeClass*> &Nodes, vector<NodeClass*> &BestNodes, SegmentClass *FirstFlipSegment, 
	  long int ParentKey, long int *BestFlipKey, int Gk, int *BestGk, int *G0, int *Gain, bool bridge = false, const int klimit=INT_MAX);
      void MakeOptMove(vector<NodeClass*> &Nodes, vector<SwapRecord> &Flips);
      
      //bool InsertSeg(SegmentClass *First, SegmentClass *Seg, vector<NodeClass*> &Nodes, long int ParentKey=0, bool brige = false);
      //bool InsertSegR2(SegmentClass *First, SegmentClass *Seg, vector<NodeClass*> &Nodes, FlipLogClass *Pnt, long int &Key, bool &Feas);

      bool InsertSegR3(int n, vector<NodeClass*> &Nodes, FlipLogClass *Pnt, NextFlipClass *&NextCnd);
      void MakeKoptFlip();
      bool MakeKoptFlipR(const long int ParentKey, const int k, SegmentClass *First, const bool Brige=false, const int BrigeDepth=0);

      void out_tour(const string &filename = "tour", const int out = 0);
      void out_minimum_spanning_tree();
      void out_minimum_1_tree();
      void check_route(vector<NodeClass*> &NowNodes, int k);

      void free()
      { delete [] SwapStack; delete [] FlipSegs; }

      friend class LKModel;
      friend class LKDataModel;

  };

  class LKModel : public modelClass
  {
    private:

      int *BestTour;
      int *BetterTour;
      double LowerBound;
      double Optimum;
      double BestCost, WorstCost;
      //SwapRecord *SwapStack
    public: 
      LKModel()
      {
        DataModel = new LKDataModel(this);
	OptModel = new LKSolver(this);
	Optimum = -DBL_MAX;
      };
      ~LKModel(){}
      void importData();
      void optimize();
      void exportData();
      void free();

      friend class LKDataModel;
      friend class LKSolver;
      
  };

};

#endif
