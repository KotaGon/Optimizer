#ifndef _GRAPH_H_
#define _GRAPH_H_
#include "simplex.h"

namespace OptimizeName
{
  class variable;

  class node
  {
    private:
    public:
      double x, y;
      node() = default;
      node(double x, double y) : x(x), y(y) { }
  };

  class graph
  {
    public:
      graph() = default;
      vector<vector<double>> dist;
      vector<node> nodes;
      vector<vector<variable>> vars;
      vector<int> basisIds;
      int model = 0;
  };

  class greedy_tour_solver
  {
    public:
      vector<int> getTour(graph &G)
      {
        int N = G.nodes.size();
        vector<int> tour;
        vector<int> visited(G.nodes.size(), 0);
        int from = 0, i = 0;
        int init_sol = 0.0;

        do
        {
          visited[from] = 1;
          tour.push_back(from);
          i = -1;

          double dist_min = 1e10;
          for(int to = 0; to < N; ++to)
            if(!visited[to] && dist_min > G.dist[from][to])
            { dist_min = G.dist[from][to]; i = to; }
          from = i;
        }
        while(i >= 0);

        for(int i = 0; i < tour.size(); ++i)
          init_sol += G.dist[tour[i%N]][tour[(i+1)%N]];

        cout << "initial tour sol = " << init_sol << endl;

        return tour;
      }
  };

};
#endif
