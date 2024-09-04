#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
using namespace std;

#define BOXSIZE 50000

int main(int argc, char *argv[])
{

  try
  {
    string msg = "";

    if(argc < 3)
    {
      msg +=  "run fault\n";
      msg +=  "arguments: [DIMENSION] [FILENAME]";
      throw invalid_argument(msg);
    }

    srand((unsigned) time(NULL));

    int N = atoi(argv[1]);
    const char *filename = argv[2];
    map<int, int> XMap, YMap;

    ofstream problem(filename);

    problem << "NAME : TEST" << endl;
    problem << "TYPE : TSP" << endl;
    problem << "DIMENSION : " << N << endl;
    problem << "NODE_COORD_SECTION" << endl;

    int i = 1;

    while( i <= N)
    {
      int X = rand() % BOXSIZE;
      int Y = rand() % BOXSIZE;

      if(XMap.find(X) != XMap.end() && YMap.find(Y) != YMap.end())
	continue;

      /*}
      double L2 = (X - BOXSIZE/2) * (X - BOXSIZE/2) + (Y - BOXSIZE/2)*(Y - BOXSIZE/2);
      double r2 = (double) BOXSIZE*BOXSIZE/4.0;
      if(!(0.7*0.7 * r2 <= L2 && L2 <= r2))
	continue;
	*/
      XMap[X] = X;
      YMap[Y] = Y;
    
      //problem << i << " " << rand() % BOXSIZE << " " << rand() % BOXSIZE << endl;
      problem << i << " " << X << " " << Y << endl;
     
      ++i;
    }

    /*
    for(int i = 1; i <= N; ++i)
      problem << i << " " << rand() % BOXSIZE << " " << rand() % BOXSIZE << endl;
    */
    problem << "EOF" << endl;

    problem.close();

  }
  catch (exception &e)
  {
    cerr << e.what() << endl;
    return 1;
  }

  return 0;
}
