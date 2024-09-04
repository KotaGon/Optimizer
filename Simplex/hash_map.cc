#include "simplex.h"
#include "hash_map.h"

namespace OptimizeName
{
  void HashMap::set_size(long int size) 
  { SIZE1 = size; }
  long int HashMap::hash_test(long int i, long int j)
  {
    //long int hash1 = abs(i), hash2 = abs(j);
    //auto hash1 = i, hash2 = j;
    //auto hash1 = hash<long int>{}(i);
    //auto hash2 = hash<long int>{}(j);
    long int seed = 0;

    seed = j + SIZE1 * (i + 2);

    //seed ^= hash1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    //seed ^= hash2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    //seed = hash1 ^ hash2;

    seed = seed % TABLE_SIZE;
    return seed;
  }

  void HashMap::free()
  {
    int cnt1 = 0, cnt2 = 0;
    for(int i = 0; i < TABLE_SIZE; ++i)
    {
      if(!Tables[i]) continue;
      ++cnt1;
      for(entry *p = Tables[i]; p != 0; p = p->next_h)
	++cnt2;
    }

    cout << "hash free.. [ratio, tot] = " << (double) cnt1 / TABLE_SIZE << ", " << cnt2 << endl;

    if(free_++) return;
    delete [] Tables;
  }

  void HashMap::insert(entry *ent)
  {
    long int i = ent->row, j = ent->var_id; 
    long int seed = hash_test(i, j);
    ent->hash = seed;
    ent->next_h = Tables[seed];
    ent->pre_h = 0;

    //if(Tables[seed])
      //Tables[seed]->pre_h = ent;
    if(ent->next_h)
      ent->next_h->pre_h = ent;
    Tables[seed] = ent;
  }
  entry *HashMap::get(long int i, long int j)
  {
    long int  seed = hash_test(i, j);
    entry *ret = Tables[seed];
    if(!ret) return ret;

    while(ret && (ret->row != i || ret->var_id != j))
    { ret = ret->next_h; }

    return ret;
  }
  void HashMap::erase(entry *ent)
  {
    if(ent->pre_h)  ent->pre_h->next_h = ent->next_h;
    if(ent->next_h) ent->next_h->pre_h = ent->pre_h;
    if(ent == Tables[ent->hash]) Tables[ent->hash] = ent->next_h;
    ent->hash = -1;
    ent->pre_h = ent->next_h = 0;
  }
};
