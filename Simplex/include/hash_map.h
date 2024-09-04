#ifndef _HASH_MAP_H_
#define _HASH_MAP_H_

#include "simplex.h"

using namespace std;

namespace OptimizeName
{
  class entry;

  class HashMap
  {
    private:
      int free_           = 0;
      entry **Tables      = 0;
      long int SIZE1      = 0;
      long int TABLE_SIZE = 100'000'103;

    public:
      HashMap()
      {
	Tables = new entry*[TABLE_SIZE];
	for(int i = 0; i < TABLE_SIZE; ++i)
	  Tables[i] = 0;
      }
      ~HashMap(){ free(); }

      
      void free();
      void set_size(long int size);
      void insert(entry *ent);
      void erase(entry *ent);
      long int hash_test(long int i, long int j);
      entry* get(long int i, long int j);
  };
};

/* old version => tables has intptr class.. objects is created every insert...
   class HashMap
   {
   class intptr
   {
   public:
   intptr *next = 0;
   long int hash;
   long int i, j, val; 
   };

   private :
   long int TABLE_SIZE = 100000000;
   intptr **Tables, *recycleStore = 0;
   int free_ = 0;
   public:
   HashMap()
   { 
   Tables = new intptr*[TABLE_SIZE]; 
   for(int i = 0; i < TABLE_SIZE; ++i)
   Tables[i] = 0;
   recycleStore = 0;
   }
   ~HashMap()
   { free(); }

   long int hash_test(long int i, long int j)
   {
//auto hash1 = hash<long int>{}((i));
//auto hash2 = hash<long int>{}((j));
auto hash1 = abs(i), hash2 = abs(j);
//重複しないようにハッシュ処理
long int seed = 0;
seed ^= hash1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
seed ^= hash2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
seed = seed % TABLE_SIZE;
return seed;
}

void insert(long int i, long int j, long int val)
{
auto seed = hash_test(i, j);
intptr *p = 0;

if(recycleStore)
{ 
p = recycleStore;
recycleStore = recycleStore->next;
}
else
p = new intptr();
p->i = i;
p->j = j;
p->val = val;
p->hash = seed;

p->next = Tables[p->hash];
Tables[p->hash] = p;
}

long int getval(long int i, long int j)
{
auto seed = hash_test(i, j);
auto *ret = Tables[seed];

while(ret && (ret->i != i || ret->j != j))
ret = ret->next;

return (ret ? ret->val : -1);
}

void erase(long int i, long int j)
{
  auto seed = hash_test(i, j);
  if(!Tables[seed]) return;

  intptr *ret = Tables[seed], *old = 0;
  while(ret && (ret->i != i || ret->j != j))
  { old = ret; ret = ret->next; }
  auto *next = ret->next;

  ret->next = recycleStore;
  recycleStore = ret;
  if(old) old->next = next;
  if(ret == Tables[seed]) Tables[seed] = next;
}

void free()
{
  if(free_++) return ;

  int c1 = 0, c2 = 0;
  intptr *next = 0;
  for(long int i = 0; i < TABLE_SIZE; ++i)
  {
    if(!Tables[i]) continue;
    if(Tables[i]) c1++;
    for(intptr *p = Tables[i]; p != 0; p = next)
    {
      next = p->next;
      delete p;
    }
  }
  for(intptr *p = recycleStore; p != 0; p = next, c2++)
  {
    next = p->next;
    delete p;
  }
  cout << "?" << (double) c1 / TABLE_SIZE << " " << c2 << endl;
  delete [] Tables;
}


};
*/
#endif 
