#ifndef _SIMPLEX_H_
#define _SIMPLEX_H_
#include <iostream>
#include <ostream>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <cassert>
#include "timer.h"
#include "graph.h"
#include "hash_map.h"
/* --- class and struct --- */

//#define pragma(x) _Pragma(#x)
//#pragma GCC target("sse, sse2, sse3, ssse3, sse4, popcnt, abm, mmx, avx, tune=nagative");
//#pragma GCC optimize("O3")

#define eps 1.0e-6

inline bool equal(double a, double b)
{ return fabs(a - b) < eps; }

constexpr double inf = std::numeric_limits<double>::infinity();

namespace OptimizeName
{
  using namespace std;

  static size_t var_id = 0;
  static size_t const_id = 0;

  class simplex;
  class linexpr;
  class variable;
  class constraint;
  class simplexDictionary;
  class entry;

  enum simplex_phase
  {
    phaseI,
    phaseII,
  };
  enum operation_type
  {
    pivot,
    upper_trans1,
    upper_trans2,
  };

  enum vtype
  {
    binary,
    integer,
    continuos,
  };

  class variable
  {
    friend class simplex;

    private:
    size_t id;
    string name;
    vtype type;
    double lb, ub;
    double sol = 0.0;
    void *var_ptr;
    public:
    variable() = default;
    variable(size_t id, string name, vtype type, double lb, double ub) :
      id(var_id++), name(name), type(type), lb(lb), ub(ub) { }

    size_t getId() const 
    { return id; }

    double getSol()
    { 
      variable *var = static_cast<variable*>(var_ptr);
      return var->sol; 
    }

    bool operator < (const variable& right) const 
    { return id == right.id ? id < right.id : id < right.id; } 
    bool operator == (const variable& right) const
    { return id == right.id; }
    bool operator > (const variable& right) const
    { return id == right.id ? id > right.id : id > right.id; }

  };

  using term = pair<double, variable>;

  class linexpr
  {

    public:
      friend inline linexpr operator * (const double c, const linexpr &x);
      friend inline linexpr operator / ( const double c, const linexpr &x);
      friend inline linexpr operator + (const double c, const linexpr &x);
      friend inline linexpr operator - (const double c, const linexpr &x);
      friend inline linexpr operator + (const linexpr &x, const linexpr &y);
      friend inline linexpr operator - (const linexpr &x, const linexpr &y);
      friend inline linexpr operator * (const double c, variable &x);
      friend inline ostream& operator << (ostream &stream, linexpr &expr);

    private:
      double constant;
      vector<term> terms;
    public:
      linexpr(double val = 0.0) : constant(val){ }

      double getConst() { return constant; }
      vector<term> getTerms(){ return terms; }

      void clear()
      { constant = 0.0; terms.clear();} 

      void addterm(const term &aterm) 
      { terms.push_back(aterm); }

      size_t size(){ return terms.size(); }

      void fast_add(const double c, variable &x)
      { addterm({c, x}); }

      inline linexpr operator += (const double c)
      { constant += c; return *this; }
      inline linexpr operator -= (const double c)
      { constant -= c; return *this; }
      inline linexpr operator *= (const double c)
      {
	constant *= c;
	for(auto &aterm : terms) aterm.first *= c;
	return *this;
      }
      inline linexpr operator /= (const double c)
      { 
	constant /= c;
	for(auto &aterm : terms) aterm.first /= c;
	return *this; 
      }
      inline linexpr operator += (const linexpr &x)
      {
	auto y = x;
	constant += y.constant;
	for(auto &y_term : y.terms) addterm(y_term);
	return *this;
      }
      inline linexpr operator -= (const linexpr &x)
      {
	auto y = x; y *= -1.0;
	constant += y.constant;
	for(auto &y_term : y.terms) addterm(y_term);
	return *this;
      }
      inline linexpr operator += (variable &x)
      {
	addterm({1.0, x});
	return *this;
      }
      inline linexpr operator -= (variable &x)
      {
	addterm({-1.0, x});
	return *this;
      }

  };

  class constraint
  {
    public:
      friend class variable;
      friend class linexpr;
      friend class simplex;
      friend inline constraint operator <= (linexpr x, linexpr y);
      friend inline constraint operator == (linexpr x, linexpr y);
      friend inline ostream& operator << (ostream &stream, constraint &constr);

    private:
      size_t id;
      string name;
      string sense;
      linexpr expr;
      double constant;
    public:
      constraint() = default;
      constraint(string sense, linexpr &expr) : sense(sense), expr(expr) { }
      constraint(linexpr &expr, double val) : expr(expr), constant(val) { }
      constraint(size_t id, string name, string sense, const linexpr &expr, double constant) : 
	id(const_id++), name(name), sense(sense), expr(expr), constant(constant) { }
  };

  enum objsense
  {
    minimization,
    maximize,
  };

  class entry
  {
    public:
      friend class simplex;
      friend class simplexDictionary;
      friend class HashMap;

    private:
      size_t id       = -1;
      double c        = 0.0;
      int row         = -1;
      int var_id      = -1;
      long int hash   = -1;
      entry *pre_col  = 0;
      entry *next_col = 0;
      entry *pre_row  = 0;
      entry *next_row = 0;
      entry *next     = 0;
      entry *next_h   = 0;
      entry *pre_h    = 0;

      void erase_col()
      {
	if(pre_col) pre_col->next_col = next_col;
	if(next_col) next_col->pre_col = pre_col;
	pre_col = next_col = 0;
      }
      void erase_row()
      {
	if(pre_row) pre_row->next_row = next_row;
	if(next_row) next_row->pre_row = pre_row;
	pre_row = next_row = 0;
      }

    public:
      entry() = default;
      entry(size_t id, double c) : id(id), c(c) { }
      entry(size_t id, double c, int row) : id(id), c(c), row(row) { }
      entry(size_t id, double c, int row, int var_id) : id(id), c(c), row(row), var_id(var_id){ }
      void print()
      { cout << "c = " << c << ", xi = " << var_id << endl; }

      inline bool operator <= (const entry& right) const 
      { return c == right.c ? var_id < right.var_id : c < right.c; } 
      bool operator == (const entry& right) const
      { return c == right.c; }
      inline bool operator >= (const entry& right) const
      { return c == right.c ? var_id < right.var_id : c > right.c; } 
  };

  class simplexDictionary
  {
    public:
      friend class simplex;

    private:

      objsense senses[2];

      entry *recycle_enrty = 0;
      entry *first_entry_objI,  *last_entry_objI,  *const_entry_objI;
      entry *first_entry_objII, *last_entry_objII, *const_entry_objII;
      vector<int>    basisIds;
      vector<int>    is_upper_transforms;
      vector<entry>  entries;
      vector<entry*> first_entry_rows;
      vector<entry*> first_entry_cols;
      vector<entry*> const_entry_rows;
      HashMap *HashTable;

      void erase_from_obj(entry *ent, simplex_phase phase)
      {
	if(phase == phaseII && ent == first_entry_objII) first_entry_objII = ent->next_col;
	else if(phase == phaseI && ent == first_entry_objI) first_entry_objI = ent->next_col;

	if(phase == phaseII && ent == last_entry_objII) last_entry_objII = ent->pre_col;
	else if(phase == phaseI && ent == last_entry_objI) last_entry_objI = ent->pre_col;

	if(ent == first_entry_cols[ent->var_id]) first_entry_cols[ent->var_id] = ent->next_row;

	ent->erase_col();
	ent->erase_row();

	ent->next = recycle_enrty;
	recycle_enrty = ent;

	HashTable->erase(ent);

	//erase(ent);
      }

      void erase_from_constr(entry *ent)
      {
	if(ent == first_entry_rows[ent->row]) first_entry_rows[ent->row] = ent->next_col;
	if(ent == first_entry_cols[ent->var_id]) first_entry_cols[ent->var_id] = ent->next_row;

	ent->erase_col();
	ent->erase_row();

	ent->next = recycle_enrty;
	recycle_enrty = ent;

	HashTable->erase(ent);

	//erase(ent);
      }

      /*
	 void erase(entry*ent)
	 {
	 size_t id = ent->id;
	 size_t nsize = entries.size();
	 entry *ent_b = &entries[nsize-1];

	 if(id != nsize - 1)
	 {
	 if(ent_b == first_entry_rows[ent_b->row]) first_entry_rows[ent_b->row] = &entries[id];
	 if(ent_b == first_entry_cols[ent_b->var_id]) first_entry_cols[ent_b->var_id] = &entries[id];
	 if(ent_b == first_entry_objII) first_entry_objII = &entries[id];
	 if(ent_b == first_entry_objI) first_entry_objI = &entries[id];
	 if(ent_b->row >= 0 && ent_b == const_entry_rows[ent_b->row]) const_entry_rows[ent_b->row] = &entries[id];
	 if(ent_b->row < 0 && ent_b == const_entry_objII) const_entry_objII = &entries[id];
	 if(ent_b->row < 0 && ent_b == const_entry_objI) const_entry_objI = &entries[id];
	 }
	 else 
	 {
	 if(ent_b == first_entry_rows[ent_b->row]) first_entry_rows[ent_b->row] = ent_b->next_col;
	 if(ent_b == first_entry_cols[ent_b->var_id]) first_entry_cols[ent_b->var_id] = ent_b->next_col;
	 if(ent_b == first_entry_objII) first_entry_objII = ent_b->next_col;
	 if(ent_b == first_entry_objI) first_entry_objI = ent_b->next_col;
	 }
      //ent_b->erase_col();
      //ent_b->erase_row();

      if(id != nsize - 1)
      {
      entry *ptr = id != nsize - 1 ? &entries[id] : 0;
      if(ent_b->pre_col) ent_b->pre_col->next_col = ptr;//ent_b->next_col;
      if(ent_b->next_col) ent_b->next_col->pre_col = ptr;// ent_b->pre_col;
      if(ent_b->pre_row) ent_b->pre_row->next_row = ptr; //ent_b->next_row;
      if(ent_b->next_row) ent_b->next_row->pre_row = ptr; //ent_b->pre_row;
      }

      //if(nsize -1 != id)
      {
      //HashTable->erase(ent_b->row, ent_b->var_id);
      //HashTable->insert(ent_b->row, ent_b->var_id, id);
      }

      ent_b->id = id;
      //swap(entries.back(), entries[id]);
      //swap(entries[nsize-1], entries[id]);
      entries[id] = entries[nsize-1];
      //entries[id] = move(entries[nsize-1]);

      ent_b = &entries.back();
      ent_b->pre_col = 0;
      ent_b->next_col = 0;
      ent_b->pre_row = 0;
      ent_b->next_row = 0;
      ent_b->row = -1;
      ent_b->var_id = -1;

      entries.resize(nsize - 1);
      }
       */

      void add_list_obj(entry *ent, simplex_phase phase)
      {
	double coeff = ent->c;

	if(phase == phaseII)
	{
	  if((senses[1] == minimization && coeff > 0.0) || (senses[1] == maximize && coeff < 0.0))
	  {
	    ent->pre_col = last_entry_objII;
	    if(last_entry_objII) last_entry_objII->next_col = ent;
	    last_entry_objII = ent;
	  }
	  else 
	  {
	    ent->next_col = first_entry_objII;
	    if(first_entry_objII) first_entry_objII->pre_col = ent;
	    first_entry_objII = ent;
	  }

	  if(!first_entry_objII) first_entry_objII = ent;
	  if(!last_entry_objII)  last_entry_objII  = ent;
	}
	else 
	{
	  if((senses[0] == minimization && coeff < 0.0) || (senses[0] == maximize && coeff > 0.0))
	  {
	    ent->pre_col = last_entry_objI;
	    if(last_entry_objI) last_entry_objI->next_col = ent;
	    last_entry_objI = ent;
	  }
	  else 
	  {
	    ent->next_col = first_entry_objI;
	    if(first_entry_objI) first_entry_objI->pre_col = ent;
	    first_entry_objI = ent;
	  }

	  if(!first_entry_objI) first_entry_objI = ent;
	  if(!last_entry_objI)  last_entry_objI  = ent;
	}
	return;
      };

      void push_to_obj(double coeff, int var_id = -1, simplex_phase phase = phaseII)
      {
	size_t nsize = entries.size();
	if(var_id < 0)
	{
	  if(recycle_enrty)
	  {
	    nsize = recycle_enrty->id;
	    recycle_enrty->c = coeff;
	    recycle_enrty->var_id = var_id;
	    recycle_enrty = recycle_enrty->next;
	  }
	  else 
	  {
	    //entries.push_back({nsize, coeff});
	    entries.emplace_back(nsize, coeff);
	  }
	  if(phase == phaseII)
	    const_entry_objII = &entries[nsize];
	  else 
	    const_entry_objI = &entries[nsize];
	}
	else 
	{
	  if(phase == phaseII)
	  {

	    if(recycle_enrty)
	    {
	      nsize = recycle_enrty->id;
	      recycle_enrty->c = coeff;
	      recycle_enrty->row = -2;
	      recycle_enrty->var_id = var_id;
	      recycle_enrty = recycle_enrty->next;
	    }
	    else 
	      //entries.push_back({nsize, coeff, -2, var_id});
	      entries.emplace_back(nsize, coeff, -2, var_id);
	  
	  }
	  else 
	  {
	    if(recycle_enrty)
	    {
	      nsize = recycle_enrty->id;
	      recycle_enrty->c = coeff;
	      recycle_enrty->row = -1;
	      recycle_enrty->var_id = var_id;
	      recycle_enrty = recycle_enrty->next;
	    }
	    else 
	      //entries.push_back({nsize, coeff, -1, var_id}) ;
	      entries.emplace_back(nsize, coeff, -1, var_id);
	  }

	  add_list_obj(&entries[nsize], phase);

	  HashTable->insert(&entries[nsize]);

	  entries[nsize].next_row = first_entry_cols[var_id];
	  if(first_entry_cols[var_id]) first_entry_cols[var_id]->pre_row = &entries[nsize];
	  first_entry_cols[var_id] = &entries[nsize];
	}
      }

      void push_to_constr(double coeff, int row, int var_id = -1)
      {
	size_t nsize = entries.size();
	if(var_id < 0)
	{
	  if(recycle_enrty)
	  {
	    nsize = recycle_enrty->id;
	    recycle_enrty->row = row;
	    recycle_enrty->c = coeff;
	    recycle_enrty->var_id = -1;
	    recycle_enrty = recycle_enrty->next;
	  }
	  else 
	    entries.emplace_back(nsize, coeff, row);
	    //entries.push_back({ nsize, coeff, row});
	  const_entry_rows[row] = &entries[nsize];
	}
	else 
	{
	  if(recycle_enrty)
	  {
	    nsize = recycle_enrty->id;
	    recycle_enrty->row = row;
	    recycle_enrty->c = coeff;
	    recycle_enrty->var_id = var_id;

	    recycle_enrty = recycle_enrty->next;
	  }
	  else 
	    entries.emplace_back(nsize, coeff, row, var_id);
	    //entries.push_back( { nsize, coeff, row, var_id} );

	  entries[nsize].next_row = first_entry_cols[var_id];
	  if(first_entry_cols[var_id]) first_entry_cols[var_id]->pre_row = &entries[nsize];
	  first_entry_cols[var_id] = &entries[nsize];

	  entries[nsize].next_col = first_entry_rows[row];
	  if(first_entry_rows[row]) first_entry_rows[row]->pre_col = &entries[nsize];
	  first_entry_rows[row] = &entries[nsize];

	  HashTable->insert(&entries[nsize]);
	}
      }

    public:
      simplexDictionary()
      { entries.reserve(INT_MAX); HashTable = new HashMap; }

      ~simplexDictionary() { HashTable->free(); }

      void print(int conv = 0)
      {
	vector<pair<int, double>> cxij;
	cout << endl;
	cout << "objective: " << const_entry_objII->c;
	for(auto *ent = first_entry_objII; ent != 0; ent = ent->next_col)
	  cxij.push_back({ ent->var_id + conv, ent->c} );
	sort(cxij.begin(), cxij.end());
	int j = 0;
	for(auto &pp : cxij)
	{
	  ++j;
	  if( cxij.size() > 10 && j > 5 && j < cxij.size() - 2)
	  { if(j == 6) cout << "..."; continue; }
	  cout << " + " << pp.second << "x" << pp.first;
	}
	cout << "[" << cxij.size() << "]" << endl;
	for(int i = 0; i < first_entry_rows.size(); ++i)
	{
	  if(first_entry_rows.size() > 10 && i > 5 && i < first_entry_rows.size() - 2)
	  { if(i == 6) cout << "..." << endl; continue; }
	  cout << "constr" << i << ": x" << basisIds[i] + conv << " = " << const_entry_rows[i]->c;
	  cxij.clear();
	  for(auto *ent = first_entry_rows[i]; ent != 0; ent = ent->next_col)
	    cxij.push_back({ ent->var_id + conv, ent->c} );

	  sort(cxij.begin(), cxij.end());
	  int j = 0;
	  for(auto &pp : cxij)
	  {
	    ++j;
	    if( cxij.size() > 10 && j > 5 && j < cxij.size() - 2)
	    { if(j == 6) cout << "..."; continue; }
	    cout << " + " << pp.second << "x" << pp.first;
	  }

	  cout << "[" << cxij.size() << "]" << endl;
	}
      }
  };

  class simplex
  {
    private:

      objsense sense, senseI, senseII;
      linexpr objective;
      vector<variable> vars;
      unordered_map<string, size_t> varid_map;
      vector<constraint> constraints;
      unordered_map<string, size_t> constraintid_map;

      size_t iter;
      simplex_phase phase; 
      simplexDictionary dictionary;

      graph G;

      Timer timer;

      void init();
      int getPivotColumn();
      entry* getPivotRow(int col, operation_type &otype);
      void pivotOperation(entry *);
      void enter(int row, entry *, entry *);
      void setOriginalObjective();
      void upperBound_trans(int index, operation_type otype);

    public:
      simplex()
      { vars.reserve(INT_MAX); constraints.reserve(INT_MAX); } 

      void initGraph();
      void setGraph(vector<double> &x, vector<double> &y, vector<vector<double>> &dist, vector<vector<variable>> &gvars)
      { 
        G.model = 1;
        const int N = dist.size();
        for(int i = 0; i < N; ++i) G.nodes.push_back({ x[i], y[i]});
        G.dist = dist;
        G.vars = gvars;
      }

      //add variable
      variable addVar(vtype type = continuos, double lb = -inf, double ub = inf, string name = "")
      {
	size_t nsize = vars.size();
	if(name != "") varid_map[name] = nsize;
	vars.push_back({ var_id, name, type, lb, ub });
	vars[nsize].var_ptr = &vars[nsize];
	return vars[nsize];
      }

      //add constraints
      constraint addConstr(const constraint &constr, string name = "")
      {
	size_t nsize = constraints.size();
	if(name != "") constraintid_map[name] = nsize;
	constraints.push_back({const_id, name, constr.sense, constr.expr, constr.constant});
	if(constr.sense == "=")
	  constraints.push_back({const_id, name, constr.sense, -1.0 * constr.expr, -constr.constant});

	return constraints[nsize];
      }

      //set objective function 
      void setObj(linexpr &expr, objsense obj_sense = minimization)
      { objective = expr; senseII = obj_sense; }

      //getter 
      variable getVar(int id) { return vars[id]; }
      variable getVar(string name) { return vars[varid_map[name]]; }

      void optimize();
      void print();
      void setSol();

  };

  //operators 
  inline linexpr operator * (const double c, const linexpr &x)
  {
    auto ret = x;
    ret *= c;
    return ret;
  }
  inline linexpr operator / ( const double c, const linexpr &x)
  {
    auto ret = x;
    ret /= c;
    return ret;
  }
  inline linexpr operator + (const double c, const linexpr &x)
  {
    auto ret = x;
    ret.constant += c;
    return ret;
  }
  inline linexpr operator - (const double c, const linexpr &x)
  {
    auto ret = x;
    ret.constant -= c;
    return ret;
  }
  inline linexpr operator + (const linexpr &x, const linexpr &y)
  {
    auto ret = x;
    ret += y;
    return ret;
  }
  inline linexpr operator - (const linexpr &x, const linexpr &y)
  {
    auto ret = x;
    ret -= y;
    return ret;
  }
  inline linexpr operator * (const double c, variable &x)
  {
    linexpr ret;
    ret.addterm({c, x});
    return ret;
  }
  inline ostream& operator << (ostream &stream, linexpr &expr)
  {
    auto tempterms = expr.terms;
    sort(tempterms.begin(), tempterms.end(), [&](const term &l, const term &r){
	return l.second.getId() < r.second.getId();
	});
    stream << expr.constant;
    int i = 0;
    for(auto &aTerm : tempterms)
    {
      ++i;
      if(tempterms.size() > 10 && i > 5 && i < tempterms.size() - 2)
      { if(i == 6) cout << "..." << endl; continue; }
      stream << " + " << aTerm.first << "x" << aTerm.second.getId();
    }
    return stream;
  }

  inline constraint operator <= (linexpr x, linexpr y)
  {
    constraint constr = { "<=", x };
    constr.expr -= y;
    constr.constant = -1.0 * constr.expr.getConst();

    return constr;
  }

  inline constraint operator == (linexpr x, linexpr y)
  {
    constraint constr = { "=", x };
    constr.expr -= y;
    constr.constant = -1.0 * constr.expr.getConst();

    return constr;
  }

  inline ostream& operator << (ostream &stream, constraint &constr)
  {
    auto &expr = constr.expr;
    auto tempterms = expr.getTerms();
    sort(tempterms.begin(), tempterms.end(), [&](const term &l, const term &r){
	return l.second.getId() < r.second.getId();
	});
    int i = 0;
    for(auto &aTerm : tempterms)
    {

      ++i;
      if(tempterms.size() > 10 && i > 5 && i < tempterms.size() - 2)
      { if(i == 6) cout << "..." << endl; continue; }

      stream << aTerm.first << "x" << aTerm.second.getId() << " + ";
    }
    stream << " <= " << constr.constant;
    return stream;
  }

};

#endif
