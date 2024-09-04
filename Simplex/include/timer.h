#ifndef _TIMER_H_
#define _TIMER_H_

using namespace std;

class Timer{
  chrono::system_clock::time_point _start, _end;
  long long int _sum = 0, _count = 0;

  public:
  void start(){
    _start = chrono::system_clock::now();
  }

  void stop(){
    _end = chrono::system_clock::now();
  }

  void add(){
    const chrono::system_clock::time_point now = chrono::system_clock::now();
    _sum += static_cast<double>(chrono::duration_cast<chrono::nanoseconds>(now - _start).count());
    _count++;
  }

  long long int sum(){
    return _sum / 1000;
  }

  string average(){
    if(_count == 0){
      return "NaN";
    }
    return to_string(_sum / 1000 / _count);
  }

  void reset(){
    _start = chrono::system_clock::now();
    _sum = 0;
    _count = 0;
  }

  inline int ms() const{
    const chrono::system_clock::time_point now = chrono::system_clock::now();
    return static_cast<double>(chrono::duration_cast<chrono::microseconds>(now - _start).count() / 1000);
  }

  inline int ns() const{
    const chrono::system_clock::time_point now = chrono::system_clock::now();
    return static_cast<double>(chrono::duration_cast<chrono::microseconds>(now - _start).count());
  }
};

#endif 


