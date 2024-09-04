#include <iostream>
#include <fstream>
#include "simplex.h"

using namespace OptimizeName;

double sq(const double x) { return x * x; }

int main()
{
  Timer t;
  simplex *lp_solver = new simplex();
  linexpr expr;

  vector<variable> vars;
  int N = 2000;
  vector<double> x, y;
  vector<vector<double>> dist(N, vector<double>(N, 0));
  vector<vector<variable>> gvars(N, vector<variable>(N));
 
  int ii = 0;
  while(ii < N)
  {
    double xp = (double) rand() / (double) RAND_MAX;
    double yp = (double) rand() / (double) RAND_MAX;
    double l = sqrt(sq(xp - 0.5) + sq(yp - 0.5));
    if(!(0.42 < l && l < 0.48)) 
    {
      continue;
    }
    x.push_back(xp);
    y.push_back(yp);
    ++ii;
  }

#if 0
  x.clear();
  y.clear();
  for(int i = 0; i < N; ++i)
    x.push_back((double) rand() / RAND_MAX);
  for(int i = 0; i < N; ++i)
    y.push_back((double) rand() / RAND_MAX);
#endif 

  for(int i = 0; i < N; ++i)
    for(int j = 0; j < N; ++j)
    {
      dist[i][j] = sqrt((x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]));
      if(i == j)
        continue;
      double d = dist[i][j];
      dist[i][j] *= 1000000.0;
      dist[i][j] = (int) dist[i][j];
      dist[i][j] /= 1000000.0;
    }

    cout << "create variable..";

  vector<double> ci;
  map<int, map<int, variable>> varmap;
  for(int i = 0; i < N; ++i)
  {
    for(int j = i + 1; j < N; ++j)
    {
      if(dist[i][j] == 0.0) continue;

      variable var = lp_solver->addVar(continuos, 0.0, 1.0);
      varmap[i][j] = var;
      gvars[i][j] = var;
      vars.push_back(var);
      ci.push_back(dist[i][j]);

      //expr.clear();
      //expr += var;
      //lp_solver->addConstr(expr <= 1);
    }
  }

  cout << "done" << endl;

  cout << "set objective..";
  expr.clear();
  for(int i = 0; i < ci.size(); ++i)
    expr.fast_add(ci[i], vars[i]);
    //expr += (ci[i] * vars[i]);
  lp_solver->setObj(expr, minimization);
  cout << "done" << endl;

  cout << "create constraints..";
  
  for(int i = 0; i < N; ++i)
  {
    linexpr expr;

    for(int j = 0; j <= i; ++j)
    {
      if(dist[j][i] == 0.0) continue;
      expr += varmap[j][i];
    }

    for(int j = i + 1; j < N; ++j)
    {
      if(dist[i][j] == 0.0) continue;
      expr += varmap[i][j];
    }

    lp_solver->addConstr(expr == 2.0);
  }
  cout << "done" << endl;

  ofstream file("./out/dist.dat");
  for(int i = 0; i < N; ++i)
    for(int j = 0; j < N; ++j)
      file << dist[i][j] << endl;
  file.close();

  ofstream pos_file("./out/pos.dat");
  for(int i = 0; i < N; ++i)
    pos_file << x[i] << "\t" << y[i] << endl;
  pos_file.close();

    lp_solver->setGraph(x, y, dist, gvars);
  lp_solver->optimize();

  ofstream sol_file("./out/sol.dat");

  double opt_val = 0.0;
  for(int i = 0; i < N; ++i)
  {
    for(int  j = i + 1; j < N; ++j)
    {
      if(dist[i][j] == 0.0) continue;
      if(varmap[i][j].getSol() > 1.0e-16) 
      {
        sol_file << i << "\t" << j << "\t" << varmap[i][j].getSol() << "\t" << dist[i][j] << endl;
        opt_val += dist[i][j] * varmap[i][j].getSol();
      }
    }
  }
  sol_file.close();

  cout << opt_val << endl;


  return 0;
}

