#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

#define BOXSIZE 1000

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

    ofstream problem(filename);

    problem << "NAME : TEST" << endl;
    problem << "TYPE : TSP" << endl;
    problem << "DIMENSION : " << N << endl;
    problem << "NODE_COORD_SECTION" << endl;

    for(int i = 1; i <= N; ++i)
      problem << i << " " << rand() % BOXSIZE << " " << rand() % BOXSIZE << endl;
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
