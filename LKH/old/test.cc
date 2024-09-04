#include "LK.h"

#include <sstream>
#include <iomanip>

namespace myOptimizeName
{

  void LKSolver::out_tour(const string &filename, const int out)
  {
     LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
     NodeClass *t = Data->FirstNode;
     
     std::ostringstream ss;
     ss << std::setw(4) << std::setfill('0') << out;
     std::string tour_number(ss.str());

     ofstream tour("output/"+filename + tour_number + ".gnu");
     tour << "#!/bin/gnuplot" << endl;
     tour << "set terminal png" << endl;
     tour << "set output 'png/tour" << tour_number << ".png'" << endl; 
     tour << "plot '-' w l," << endl;
     
     int Cost = 0;

     do
     {
       tour << t->X << " " << t->Y << endl;
       tour << SUC(t)->X << " " << SUC(t)->Y << endl;

       tour << endl;
       tour << endl; 
       Cost += Data->C(t, SUC(t)) - t->Pi - SUC(t)->Pi;
     }
     while((t = SUC(t)) != Data->FirstNode);
     tour << "e" << endl;

     //cout << " Cost = " << Cost << endl;
    
     tour.close();
  }

  void LKSolver::out_minimum_spanning_tree()
  {
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    ofstream mst("output/mst.gnu");
    NodeClass *t = Data->FirstNode;

    mst << "#!/bin/gnuplot" << endl;
    mst << "plot '-' w l," << endl;

    do
    {
      if(t->Dad == NULL)
      { continue; }
      mst << t->Dad->X << " " << t->Dad->Y << endl;
      mst << t->X << " " << t->Y << endl;
      mst << endl;
      mst << endl;
    }
    while((t = t->Suc) != Data->FirstNode);

    mst << "e" << endl;

    mst.close();


  }

  void LKSolver::check_route(vector<NodeClass*> &NowNodes, int k)
  {

    k *= 2;
    LKDataModel *Data = (LKDataModel*) _myModel->getDataModel();
    NodeClass *Node = Data->FirstNode;

    do
    {
      SegmentClass *Parent = Node->Parent;
      if( k >= 1 && Node == NowNodes[1])
	cout << "[NODE01] " << Node->Id << " " << Node->Rank << " " << Parent->Rank << endl;
      else if( k >= 2 && Node == NowNodes[2])
	cout << "[NODE02]" << Node->Id << " " << Node->Rank << " " << Parent->Rank << endl;
      else if( k >= 3 && Node == NowNodes[3])
	cout << "[NODE03]" << Node->Id << " " << Node->Rank << " " <<  Parent->Rank << endl;
      else if( k >= 4 && Node == NowNodes[4])
	cout << "[NODE04]" << Node->Id << " " <<Node->Rank << " " << Parent->Rank << endl;
      else if( k >= 5 && Node == NowNodes[5])
	cout << "[NODE05]" << Node->Id << " " <<Node->Rank << " " << Parent->Rank << endl;
      else if( k >= 6 && Node == NowNodes[6])
	cout << "[NODE06]" << Node->Id << " " <<Node->Rank << " " << Parent->Rank << endl ;
      else if( k >= 7 && Node == NowNodes[7])
	cout << "[NODE07]" << Node->Id << " " <<Node->Rank << " " << Parent->Rank << endl;
      else if( k >= 8 && Node == NowNodes[8])
	cout << "[NODE08]" << Node->Id << " " <<Node->Rank << " " << Parent->Rank << endl;
      else if( k >= 9 && Node == NowNodes[9])
	cout << "[NODE09]" << Node->Id << " " <<Node->Rank << " " << Parent->Rank << endl;
      else if( k >= 10 && Node == NowNodes[10])
	cout << "[NODE10]" << Node->Id << " " << Node->Rank << " " << Parent->Rank << endl;
      else 
	cout << Node->Id << " " << Node->Rank << " " << Parent->Rank << endl;
    }
    while((Node = SUC(Node)) != Data->FirstNode);

    cout << "in check_route" << k << endl;

    return;
  }
};
