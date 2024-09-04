#ifndef _SA_H
#define _SA_H

#include "model.h"
#include "node.h"
#include "utils.h"
#include <queue>

namespace myOptimizeName
{

  class SADataModel;
  class SASolver;
  class SAModel;

  class SADataModel : public DataModelClass
  {
    private:

      int Dimension;
      long *CostMatrix;
      NodeClass *NodeSet, *FirstNode;

      SAModel *_myModel;

    public:
      SADataModel() = default;
      SADataModel(SAModel *model) :  _myModel(model) {}
      ~SADataModel(){}

      void readProblem();
      void createNodes(int *ID, long *X, long *Y);

      inline long C(NodeClass *Na, NodeClass *Nb)
      { return Na->Id < Nb->Id ? Nb->C[Na->Id] : Na->C[Nb->Id]; }

      void free()
      { 
	delete [] NodeSet; 
      }
     
      friend class SAModel;
      friend class SASolver;

  };

  class SASolver : public OptModelClass
  {

    private:

      SAModel *_myModel;

    public:

      SASolver() = default;
      SASolver(SAModel*model) : _myModel(model) {}
      ~SASolver(){}

      int MaxIterator = 100000;
      int Reversed;
      double temp;
      queue<NodeClass*> ActiveNodes;

      void free()
      { }

      void CreateInitialTour();
      void SimulatedAnnealing();
      void Best2OptMove(NodeClass *t1, NodeClass *t2, long *Gain, bool *update);
      void out_tour(const string &filename, const int out);
      double Temp(const int Iterator);

      friend class SAModel;
      friend class SADataModel;

  };

  class SAModel : public modelClass
  {
    private:

      int *BestTour;
      int *BetterTour;
      long BestCost;
    public: 
      SAModel()
      {
        DataModel = new SADataModel(this);
	OptModel = new SASolver(this);
      };
      ~SAModel(){}
      void importData();
      void optimize();
      void exportData();
      void free();

      friend class LKDataModel;
      friend class LKSolver;
      
  };


};


#endif
