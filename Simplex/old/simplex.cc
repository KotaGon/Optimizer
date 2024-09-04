#include "simplex.h"

using namespace OptimizeName;

Timer t1, t2, t3, t4, t5;
double tt1=0.0, tt2=0.0, tt3=0.0;
double tt4 = 0.0, tt5 = 0.0;
int ppp = 0;

void simplex::init()
{

  cout << "initialize..." ;

  /* all clear */
  dictionary.entries.clear();
  dictionary.first_entry_rows.clear();
  dictionary.first_entry_cols.clear();
  dictionary.basisIds.clear();
  dictionary.const_entry_rows.clear();

  /* alloc size. based resize method */
  long int nvar = vars.size(), nconstr = constraints.size();
  dictionary.first_entry_rows.resize(nconstr, 0);
  dictionary.first_entry_cols.resize(nvar + nconstr + 1, 0);
  dictionary.basisIds.resize(nconstr);
  dictionary.const_entry_rows.resize(nconstr);
  dictionary.is_upper_transforms.resize(nvar, 0);

  dictionary.HashTable->set_size(nvar + nconstr + 2);
  
  /* senses each phase */
  dictionary.senses[0] = senseI;
  dictionary.senses[1] = senseII;

  /* if model is graphical(tsp), then set nodes id to variable id. 
     And get initial solution tour by heuristic method */
  initGraph();

  /* set simplex dict. at same time, set init tour if model is graphical */

  /* object(phaseII)'s dict */
  double constant = objective.getConst();
  dictionary.first_entry_objII = 0;
  dictionary.push_to_obj(constant);

  for(auto &aTerm : objective.getTerms())
  {
    double c = aTerm.first; 
    variable &var = aTerm.second;
    int var_id = var.getId();

    if(!equal(c, 0.0))
    {
      /*
      if(G.model && equal(var.getSol(), 1.0))
      { dictionary.const_entry_objII->c += c; dictionary.push_to_obj(-c, G.basisIds[var_id]); }
      else 
      */
	dictionary.push_to_obj(c, var_id);
    }
  }

  /* object(phaseI)'s dict */ 
  dictionary.first_entry_objI = 0;
  dictionary.push_to_obj(0.0, -1.0, phaseI);
  dictionary.push_to_obj(-1.0, nvar + nconstr, phaseI);

  /* constraint' dict */
  for(int i = 0; i < nconstr; ++i)
  {
    auto &constr = constraints[i];
    double constant = constr.constant;
    dictionary.first_entry_rows[i] = 0;
    dictionary.push_to_constr(constant, i); 

    dictionary.basisIds[i] = nvar + i;

    for(auto &aTerm : constr.expr.getTerms())
    {
      double c = aTerm.first;
      variable &var = aTerm.second;
      int var_id = var.getId();

      /*
      if(G.model && constr.expr.size() == 1 && equal(var.getSol(), 1.0))
      {
	dictionary.basisIds[i] = var_id;
      }
      */

      if(!equal(c, 0.0))
      {
	/*
	if(G.model && equal(var.getSol(), 1.0))
	{
	  if(dictionary.basisIds[i] != var_id)
	  {
	    dictionary.const_entry_rows[i]->c += -c;
	    dictionary.push_to_constr(c, i, G.basisIds[var_id]);
	  }
	  else 
	    dictionary.push_to_constr(-c, i, G.basisIds[var_id]);
	}
	else 
	*/
	  dictionary.push_to_constr(-c, i, var_id);
      }
    }
  }

  /* set initial solution */ 
  for(auto &var : vars)
    if(equal(var.getSol(), 1.0))
      upperBound_trans(var.getId(), upper_trans2);

  /* add dummy variable */
  for(int i = 0; i < nconstr; ++i)
    if(dictionary.const_entry_rows[i]->c < 0.0)
      dictionary.push_to_constr(1.0, i, nvar + nconstr);

  /* test print */
  dictionary.print(0);

  sense = senseI = maximize;

  cout << "done" << endl;
}

void simplex::initGraph()
{
  int nvar = vars.size(), nconstr = constraints.size();

  if(G.model)
  {
    G.basisIds.resize(vars.size()); //変換する基底変数のIdを格納

    for(int i = 0; i < nconstr; ++i)
    {
      auto &constr = constraints[i];
      if(equal(constr.constant, 1.0) && constr.expr.size() == 1)
      {
	int var_id = -1;
	for(auto &aTerm : constr.expr.getTerms())
	{
	  variable &var = aTerm.second;
	  var_id = var.getId();
	  G.basisIds[var_id] = nvar + i;
	  break;
	}

      }
    }

    greedy_tour_solver heuristic;
    vector<int> tour = heuristic.getTour(G);

    assert(tour.size() == G.nodes.size());

    for(int i = 0; i < tour.size(); ++i)
    {
      int from = tour[i], to = tour[(i + 1) % tour.size()];
      variable &var = G.vars[min(from, to)][max(from, to)];
      int id = var.getId();
      vars[id].sol = 1.0;
    }
  }
}

int simplex::getPivotColumn()
{
  /* largest coffiecint rule */
  entry *largest = 0;
  entry *first = phase == phaseII ? dictionary.first_entry_objII : dictionary.first_entry_objI;

  for(entry *ent = first; ent != 0; ent = ent->next_col)
  {
    if(phase == phaseI && iter == 0)
    { largest = ent; break; }

    if(sense == minimization && ent->c > 0.0) break;
    else if(sense == maximize && ent->c < 0.0) break;

    if(largest == 0)
    {
      if( (ent->c <= 0.0 && sense == minimization) || (ent->c >= 0.0 && sense == maximize) )
	largest = ent;
    }
    else 
    {
      if( (*ent <= *largest && sense == minimization) || (*ent >= *largest && sense == maximize) )
	largest = ent;
    }
  }

  //if(largest) cout  << "col : coeff = " << largest->c << " id = " << largest->var_id << endl; 

  return (largest ? largest->var_id : -1);
}

entry *simplex::getPivotRow(int col, operation_type &otype)
{
  double theta_min = inf;
  entry *leaving = 0;

  for(auto *ent = dictionary.first_entry_cols[col]; ent != 0; ent = ent->next_row)
  {
    //if(ent->c >= 0.0 && !(phase == phaseI && iter == 0)) 
      //continue;
    int row = ent->row;
    if(row < 0) continue;

    double c = !(phase == phaseI && iter == 0) ? -ent->c : ent->c;
    //double theta = dictionary.const_entry_rows[row]->c / c;
    double theta = inf;

    operation_type type;
    if(c > 0.0)
    {
      type = pivot;
      theta = dictionary.const_entry_rows[row]->c / c;       
    }
    else if(c < 0.0 && dictionary.basisIds[row] < vars.size())
    {
      int id = dictionary.basisIds[row];
      if(vars[id].ub < inf)
      {
	theta = (vars[id].ub - dictionary.const_entry_rows[row]->c) / -c; 
	type = upper_trans1;
      }
    }

    if(theta < theta_min)
    { 
      theta_min = theta;
      leaving = ent;
      otype = type;
    }
  }
  
  //if(leaving) cout << "row = " << leaving->row << " aij = " << leaving->c << endl;
  
  if(leaving && leaving->var_id < vars.size())
  {
    double ub = vars[leaving->var_id].ub;
    if(theta_min > ub)
      otype = upper_trans2;
  }

  return leaving;
}

void simplex::enter(int row, entry *entering, entry *constant)
{
  int to_row = entering->row;
  /*
  entry *next_col;
  bool flag = false;
  int ncol1 = 0, ncol2 = 0;
  int ncol3 = 0, ncol4 = 0;
  for(auto *ent = first; ent != 0; ent = ent->next_col)
  {
    if(ent->var_id == dictionary.basisIds[row])
    {
      a = ent->c;
      flag = true;

      if(kind < 0)  dictionary.erase_from_obj(ent, phaseI);
      else if(kind == 0) dictionary.erase_from_obj(ent, phaseII);
      else dictionary.erase_from_constr(ent);
      ent->c = 0;
      break;
    }
  }

  if(!flag) 
  {
    cout << "not found" << endl;
    return ;
  }
  */
  double a = entering->c;
  constant->c += a * dictionary.const_entry_rows[row]->c;

#if 0
  t4.start();
  int ncol5 = 0;
  for(auto *ent = first; ent != 0; ent = next_col)
  {
    //int test_val = dictionary.HashTable->getval(ent->row, ent->var_id);
    //long int test_val = dictionary.HashTable->get(ent->row

    assert(!equal(ent->c, 0.0));
    ncol1++;
    next_col = ent->next_col;

    if(equal(coeffs[ent->var_id], 0.0))
    ++ncol5;

    if(ent->var_id != dictionary.basisIds[row])
      ent->c += a * coeffs[ent->var_id];
    else 
      ent->c = 0;
    
    coeffs[ent->var_id] = 0.0;
    
    if(equal(ent->c, 0.0))
    {
      int id = ent->id;
      //if(next_col == &dictionary.entries.back())
	//next_col = &dictionary.entries[id];

      ++ncol3;
      if(kind < 0)  dictionary.erase_from_obj(ent, phaseI);
      else if(kind == 0) dictionary.erase_from_obj(ent, phaseII);
      else dictionary.erase_from_constr(ent);
    }
  }

  tt4 += t4.ns() / 1.0e+6;

  t5.start(); 
  //t5.start();
  for(auto *ent = dictionary.first_entry_rows[row]; ent != 0; ent = ent->next_col)
  {
    ncol2++;
    double *c = &coeffs[ent->var_id];
    if(!equal(*c, 0.0))
    {
      ncol4++;
      if(kind < 0) dictionary.push_to_obj(a * *c, ent->var_id, phaseI);
      else if(kind == 0) dictionary.push_to_obj(a * *c, ent->var_id, phaseII);
      else dictionary.push_to_constr(a * *c, to_row, ent->var_id);
    }
    *c = ent->c;
  }
#endif 
  
#if 1
  t5.start(); 
  for(auto *ent = dictionary.first_entry_rows[row]; ent != 0; ent = ent->next_col)
  {
    const double c = ent->c;
    entry *ent2 = dictionary.HashTable->get(to_row, ent->var_id);

    double c_new = a * c;
    if(ent2 != NULL) c_new += ent2->c;

    if(ent2 == NULL && !equal(c_new, 0.0))
    {
      if(to_row == -1) dictionary.push_to_obj(a * c, ent->var_id, phaseI);
      else if(to_row == -2) dictionary.push_to_obj(a * c, ent->var_id, phaseII);
      else dictionary.push_to_constr(a * c, to_row, ent->var_id);
    }
    else if(equal(c_new, 0.0))
    {
      if(to_row == -1)  dictionary.erase_from_obj(ent2, phaseI);
      else if(to_row == -2) dictionary.erase_from_obj(ent2, phaseII);
      else dictionary.erase_from_constr(ent2);
    }
    else 
    {
      double mul = ent2->c * c_new;

      ent2->c = c_new;
      if(to_row < 0)
      {
	if(mul < 0.0)
	{
	  if(to_row == -2 && ent2 == dictionary.first_entry_objII) dictionary.first_entry_objII = ent2->next_col;
	  else if(to_row == -1 && ent2 == dictionary.first_entry_objI) dictionary.first_entry_objI = ent2->next_col;

	  if(to_row == -2 && ent2 == dictionary.last_entry_objII) dictionary.last_entry_objII = ent2->pre_col;
	  else if(to_row == -1 && ent2 == dictionary.last_entry_objI) dictionary.last_entry_objI = ent2->pre_col;

	  //ent2->c = c_new;
	  ent2->erase_col();
	  dictionary.add_list_obj(ent2, (to_row == -1 ? phaseI : phaseII));
	}
      }
    }
  }
#endif 
  tt5 += t5.ns() / 1.0e+6;
  //cout << ncol1 << "(" << ncol3 << ") "  << ncol2 << "(" << ncol4 << ")" << " " << ncol5 << endl;
}
 
void simplex::pivotOperation(entry *entering)
{
  int row = entering->row, col = entering->var_id;
  int leaving_id = dictionary.basisIds[row];
  double c = entering->c;
  entry *next_ent = 0;

  /* trans pivot row */
  dictionary.basisIds[row] = entering->var_id;
  dictionary.erase_from_constr(entering);
  dictionary.push_to_constr(-1.0 , row, leaving_id);
  dictionary.const_entry_rows[row]->c /= -c;
  for(auto *ent = dictionary.first_entry_rows[row]; ent != 0; ent = ent->next_col)
    ent->c /= -c;

  /* trans all row, leaving var is basisIds[row] */
  entry *next_row = 0;
  for(auto *ent = dictionary.first_entry_cols[col]; ent != 0; ent = next_row)
  {
    next_row = ent->next_row;
    if(ent->row == row) continue;

    if(ent->row == -2) 
    {
      enter(row, ent, dictionary.const_entry_objII);
      dictionary.erase_from_obj(ent, phaseII);
    }
    else if(ent->row == -1)
    {
      enter(row, ent, dictionary.const_entry_objI);
      dictionary.erase_from_obj(ent, phaseI);
    }
    else 
    {
      enter(row, ent, dictionary.const_entry_rows[ent->row]);
      dictionary.erase_from_constr(ent);
    }
    ent->c = 0;
  }

  int basis = dictionary.basisIds[row];
  if(basis < vars.size() && dictionary.is_upper_transforms[basis])
  {
    dictionary.is_upper_transforms[basis] = 0;
    dictionary.const_entry_rows[row]->c -= vars[basis].ub;
    dictionary.const_entry_rows[row]->c *= -1;

    cout << "upper" << " " << ppp++ << " " << dictionary.const_entry_rows[row]->c << endl;
    for(auto *ent = dictionary.first_entry_rows[row]; ent != 0; ent = ent->next_col)
      ent->c *= -1;
  }
}

void simplex::upperBound_trans(int index, operation_type otype)
{
  int var_id = otype == upper_trans1 ? dictionary.basisIds[index] : index;
  double ub = vars[var_id].ub;
  dictionary.is_upper_transforms[var_id] ^= 1;

  if(otype == upper_trans1)
  {
    int row = index;
    dictionary.const_entry_rows[row]->c -= ub;
    dictionary.const_entry_rows[row]->c *= -1;
    for(auto *ent = dictionary.first_entry_rows[row]; ent != 0; ent = ent->next_col)
      ent->c *= -1;
  }
  else 
  {
    entry *next_row = 0;
    //cout << "UPPER" << " " << dictionary.is_upper_transforms[var_id] << " " << var_id << " "<< ppp++ << endl;

    for(auto *ent = dictionary.first_entry_cols[var_id]; ent != 0; ent = next_row)
    {
      next_row = ent->next_row;
      int row = ent->row;
      double a = ent->c;
      double mul = ent->c * (-a);
      ent->c = -a;

      if(mul < 0.0 && row < 0)
      {
	if(row == -2 && ent == dictionary.first_entry_objII) dictionary.first_entry_objII = ent->next_col;
	else if(row == -1 && ent == dictionary.first_entry_objI) dictionary.first_entry_objI = ent->next_col;

	if(row == -2 && ent == dictionary.last_entry_objII) dictionary.last_entry_objII = ent->pre_col;
	else if(row == -1 && ent == dictionary.last_entry_objI) dictionary.last_entry_objI = ent->pre_col;

	ent->erase_col();
	dictionary.add_list_obj(ent, (row == -1 ? phaseI : phaseII));
      }

      if(row == -2) dictionary.const_entry_objII->c += a * ub;
      else if(row == -1) dictionary.const_entry_objI->c += a * ub;
      else dictionary.const_entry_rows[row]->c += a * ub;
      //if(row >= 0 && dictionary.const_entry_rows[row]->c < 0.0)
      //cout << "upper => " << var_id << " " << dictionary.const_entry_rows[row]->c << endl;
    }
  }

}

void simplex::setOriginalObjective()
{
  cout << "ready phaseII .." << endl; 
  /*
  vector<int> indexes(vars.size(), -1);
  
  for(int i = 0; i < dictionary.basisIds.size(); ++i)
  {
    int bId = dictionary.basisIds[i];
    if(bId < vars.size()) 
      indexes[bId] = i;
  }

  dictionary.print();

  entry *next_entry = 0;
  for(auto *ent = dictionary.first_entry_objII; ent != 0; ent = next_entry)
  {	  
    next_entry = ent->next_col;
    int row = indexes[ent->var_id];
    if(row < 0) continue;
    
    for(auto *ent2 = dictionary.first_entry_rows[row]; ent2 != 0; ent2 = ent2->next_col)
      coeffs[ent2->var_id] = ent2->c;
      
    //cout << "start" << endl;
    enter(row, dictionary.first_entry_objII, dictionary.const_entry_objII, 0);
    //cout << "?" << endl;

    cout << "row = " << row << " var_id = " << ent->var_id << endl;
    for(auto *ent2 = dictionary.first_entry_rows[row]; ent2 != 0; ent2 = ent2->next_col)
    {
      cout << "var_id = " << ent2->var_id << endl;
      coeffs[ent2->var_id] = 0;
    }
    //cout << "end " << endl;
  }	 

  cout << "hoge" << endl;
  */
  const int dummy_id = vars.size() + constraints.size();
  
  entry *leaving = 0, *next_col = 0;
  for(int i = 0; i < dictionary.basisIds.size(); ++i)
    if(dictionary.basisIds[i] == dummy_id)
    { leaving = dictionary.first_entry_rows[i]; break; }

  if(leaving) 
    pivotOperation(leaving);

  for(auto *ent = dictionary.first_entry_objII; ent != 0; ent = next_col)
  {
    next_col = ent->next_col;
    if(ent->var_id != dummy_id)
      continue;
    //int id = ent->id;
    //if(next_col == &dictionary.entries.back())
      //next_col = &dictionary.entries[id];
    dictionary.erase_from_obj(ent, phase);
  }

  for(int i = 0; i < dictionary.first_entry_rows.size(); ++i)
  {
    auto *first = dictionary.first_entry_rows[i];

    for(auto *ent = first; ent != 0; ent = next_col)
    {
      next_col = ent->next_col;
      if(ent->var_id != dummy_id)
	continue;
      //int id = ent->id;
      //if(next_col == &dictionary.entries.back())
	//next_col = &dictionary.entries[id];
      dictionary.erase_from_constr(ent);
    }
  }

  dictionary.print(0);

  sense = senseII;
  cout << "initial feasible sol = " << dictionary.const_entry_objII->c << endl;
  cout << "elapsed time = " << timer.ns() / 1.0e6 << "[sec]" << endl;

  //getchar();
}

void simplex::setSol()
{
  int i = 0;
  for(auto &var : vars)
    var.sol = 0.0;

  /*
  for(auto &basis : dictionary.basisIds)
  {
    if(basis < vars.size())
    {
      if(dictionary.is_upper_transforms[var.getId()])
	vars[basis].sol = vars[basis].ub - dictionary.const_entry_rows[i]->c;
      else 
	vars[basis].sol = dictionary.const_entry_rows[i]->c;
    }
    ++i;
  }

  for(auto &var : vars)
    if(dictionary.is_upper_transforms[var.getId()])
      var.sol = var.ub;
  */
  map<int, int> basis_map;
  for(auto &basis : dictionary.basisIds)
    if(basis < vars.size())
      basis_map[basis] = i++;

  for(auto &var : vars)
  {
    auto id = var.getId();
    if(basis_map[id])
    {
      i = basis_map[id];
      if(dictionary.is_upper_transforms[id])
	var.sol = var.ub - dictionary.const_entry_rows[i]->c;
      else 
	var.sol = dictionary.const_entry_rows[i]->c;
    }
    else 
    {
      if(dictionary.is_upper_transforms[id])
	var.sol = var.ub;
    }
  }
}

void simplex::print()
{
  cout << "objective:" << objective << endl;
  for(auto &constr : constraints)
    cout << "constr" << constr.id << ": " << constr << endl;
}

void simplex::optimize()
{
  const size_t niter = INT_MAX;

  init();

  ppp = 0;

  timer.start();
  for(int p = phaseII; p <= phaseII; ++p)
  {
    phase = (simplex_phase) p;
    sense = minimization;
    //if(phase == phaseII) setOriginalObjective();

    for(iter = 0; iter < niter; ++iter)
    {
      t1.start();

      if(iter % 30 == 0)
      {
	cout << "[ITER " << iter << "] " << dictionary.const_entry_objII->c << endl;
      }

      int pivot_col = getPivotColumn(); 
      if(pivot_col < 0)
      { cout << "found optimal sol.. finish simple method" << endl; break;}
      tt1 += t1.ns() / 1.0e6;

      t2.start();
      operation_type otype;
      entry *entering = getPivotRow(pivot_col, otype);
      
      if(!entering) 
      { cout << "no leaving variable." << endl; break; }
      tt2 += t2.ns() / 1.0e6;

      t3.start();

      if(otype != pivot)
	upperBound_trans(otype == upper_trans1 ? entering->row : entering->var_id, otype);
      if(otype != upper_trans2)
	pivotOperation(entering);
      tt3 += t3.ns() / 1.0e+6;
    }

    cout << tt1 << " " << tt2 << " " << tt3 << endl;
    cout << tt4 << " " << tt5 << endl;
    
    if(p == phaseI)
      cout << "finished phaseI.. " << endl;
  }

  dictionary.print(1);

  cout << "optimal solution = " << dictionary.const_entry_objII->c << endl;
  cout << "elapsed time = " << timer.ns() / 1.0e6 << "[sec]" << endl;
  
  setSol();

  dictionary.HashTable->free();

  return;
}


