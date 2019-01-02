#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <string>
// [[Rcpp::export]]
Rcpp::NumericVector get_vit_index(Rcpp::CharacterVector t_id_key, Rcpp::CharacterVector control_treatment_t_ids, Rcpp::NumericVector control_treatment_set_nums)
{
  std::unordered_map<std::string, int> indexMap;
  for (int i = 0; i < t_id_key.size(); ++i)
  {
    std::string key;
    key = t_id_key[i];
    indexMap[key] = i + 1;
  }
  Rcpp::NumericVector intdex(control_treatment_t_ids.size());
  for (int i = 0; i < control_treatment_t_ids.size(); ++i)
  {
    std::string key;
    key = control_treatment_t_ids[i];
    int idx = indexMap[key];
    int addnum = 0;
    if(control_treatment_set_nums[i])
    {
      addnum = control_treatment_set_nums[i] * t_id_key.size();
    }
    intdex[i] = idx + addnum;
  }
  return intdex;
}

// [[Rcpp::export]]
Rcpp::NumericVector get_dits(Rcpp::CharacterVector t_id_key, Rcpp::CharacterVector nonempty_t_ids)
{
  Rcpp::NumericVector dit_vector(t_id_key.size());
  for(int i = 0; i < t_id_key.size(); ++i)
  {
    for (int j = 0; j < nonempty_t_ids.size(); ++j)
    {
    	if (t_id_key[i] == nonempty_t_ids[j])
    	{
    		dit_vector[i] = 1;
    	}
    }
  }
  return dit_vector;
}

// [[Rcpp::export]]
Rcpp::List prep_lead_years(Rcpp::NumericVector ts, Rcpp::NumericVector lead_window)
{
	Rcpp::List numericvecs(ts.size());
	for (int i = 0; i < ts.size(); ++i)
	{
		Rcpp::NumericVector temp(lead_window.size());
		for (int j = 0; j < lead_window.size(); ++j)
		{
			temp[j] = ts[i] + lead_window[j];
		}
		numericvecs[i] = temp;
	}
	return(numericvecs);	
}


// [[Rcpp::export()]]
Rcpp::NumericVector sumwits(int nrow, Rcpp::NumericVector vit_vect)
{
	Rcpp::NumericVector WitVector(nrow);

	for (int i = 0; i < nrow; ++i)
	{
		double sumWit = 0;
		for (int j = i; j < vit_vect.size(); j+= nrow)
		{
			sumWit = sumWit + vit_vect[j];
		}
		WitVector[i] = sumWit;
	}
	return WitVector;
}


// [[Rcpp::export]]
Rcpp::List re_norm_index(Rcpp::NumericMatrix compmat, Rcpp::NumericVector compmat_row_units, Rcpp::NumericVector compmat_cols, int lead, Rcpp::List sets, Rcpp::NumericVector control_start_years)
{
	std::unordered_map<int, int> rowmap;
	for (int i = 0; i < compmat_row_units.size(); ++i)
	{
		rowmap[compmat_row_units[i]] = i;
		//Rcpp::Rcout << "compmat_row_units[i]: " <<compmat_row_units[i] << std::endl;
		//Rcpp::Rcout << "map1 i: " << i << std::endl;
	}
	std::unordered_map<int, int> colmap;
	for (int i = 0; i < compmat_cols.size(); ++i)
	{
		colmap[compmat_cols[i]] = i + 1;
		//Rcpp::Rcout << "compmat_cols[i]: " <<compmat_cols[i] << std::endl;
		//Rcpp::Rcout << "map2 i: " << i << std::endl;
	}
	
	Rcpp::List set_index_list(sets.size());

	for (int i = 0; i < sets.size(); ++i) //should be same size as control_start_years 
	{
		Rcpp::NumericVector controls = sets[i];
		int key;
		key = control_start_years[i];
		int st_year_col = colmap[key];
		Rcpp::LogicalVector control_index(controls.size());
		for (int j = 0; j < controls.size(); ++j)
		{
			int key;
			key = controls[j];
			int idx = rowmap[key];
			//Rcpp::LogicalVector v(lead + 1); //change this to + 2
			Rcpp::LogicalVector v(lead + 2); //change this to + 2
			//for (int k = 0; k <= lead; ++k)
			for (int k = 1; k <= lead + 1; ++k)// init. k to 1, go to lead + 1, and have first check be for t - 1
			{
				//add in separate check for t - 1, rest of these checks should still apply then.
				if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col - 1)))
				{
				  v[0] = false;
				}
				else
				{
				  v[0] = true;
				}
				
				if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col + k)))
				{
					v[k] = false;
				}
				else
				{
					v[k] = true;
				}
			}
			if (Rcpp::is_true(Rcpp::any(v == false)))
			{
				control_index[j] = false;
			}
			else
			{
				control_index[j] = true;
			}
		}
		set_index_list[i] = control_index;	
	}

	return set_index_list;
}



// [[Rcpp::export]]
Rcpp::LogicalVector check_treated_units(Rcpp::NumericMatrix compmat, Rcpp::NumericVector compmat_row_units, Rcpp::NumericVector compmat_cols, int lead, Rcpp::NumericVector treated_ids, Rcpp::NumericVector treated_ts)
{
	
	std::unordered_map<int, int> rowmap;
	for (int i = 0; i < compmat_row_units.size(); ++i)
	{
		rowmap[compmat_row_units[i]] = i;
	}
	std::unordered_map<int, int> colmap;
	for (int i = 0; i < compmat_cols.size(); ++i)
	{
		colmap[compmat_cols[i]] = i + 1;
	}

	Rcpp::LogicalVector set_index(treated_ids.size());
	for (int i = 0; i < treated_ids.size(); ++i)
	{
		int key;
		key = treated_ts[i];
		int st_year_col = colmap[key];
		key = treated_ids[i];
		int idx = rowmap[key];
		// Rcpp::LogicalVector v(lead + 1); // to lead + 2
		Rcpp::LogicalVector v(lead + 2); 
		//for (int k = 0; k <= lead; ++k) // again, thinking init k to 1, iterate to lead + 1, then first check is for t-1
		for (int k = 1; k <= lead + 1; ++k)
		{
			//checking t-1
			if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col - 1))) //assuming that the order holds
			{
			  v[0] = false;
			}
			else
			{
			  v[0] = true;
			}
			
			if( ( (st_year_col + k) > compmat_cols.size() ) || Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col + k)))
			{
				v[k] = false;
			}
			else
			{
				v[k] = true;
			}
		}
		if (Rcpp::is_true(Rcpp::any(v == false)))
		{
			set_index[i] = false;
		}
		else
		{
			set_index[i] = true;
		}
	}
	return set_index;
}

// [[Rcpp::export]]
Rcpp::LogicalVector needs_renormalization(Rcpp::List set_index_list)
{
	
	Rcpp::LogicalVector rewt(set_index_list.size());
	for (int i = 0; i < set_index_list.size(); ++i)
	{
		Rcpp::LogicalVector lv = set_index_list[i];
		if (Rcpp::is_true(Rcpp::any(lv == false)))
		{
			rewt[i] = true;
		}
	}
	return rewt;
}

// [[Rcpp::export]]
Rcpp::List renormalize(Rcpp::List control_index, Rcpp::List sets_to_be_updated)
{
  Rcpp::List updated_sets(sets_to_be_updated.size());
  for(int i = 0; i < sets_to_be_updated.size(); i++)
  {
    Rcpp::NumericVector ctls = sets_to_be_updated[i];
    Rcpp::LogicalVector idx = control_index[i];
    Rcpp::NumericVector oldwts = ctls.attr("weights");
    Rcpp::NumericVector subctl = ctls[idx];
    Rcpp::NumericVector newwts = oldwts[idx];
    double denom = Rcpp::sum(newwts);
    Rcpp::NumericVector results = newwts / denom;  //take advantage of the rcpp sugar?
    subctl.attr("weights") = results;
    updated_sets[i] = subctl;
  }
  return updated_sets;
}

// equality_four <- function(x, y, z){
//   return(sum(x*y)/sum(z))
// }
// [[Rcpp::export]]
Rcpp::NumericVector equality_four_cpp(Rcpp::NumericMatrix Wit_vals, Rcpp::NumericVector y, Rcpp::NumericVector z)
{
  Rcpp::NumericVector results(Wit_vals.ncol());
  for(int i=0; i < Wit_vals.ncol(); i++)
  {
    Rcpp::NumericVector x = Wit_vals(Rcpp::_, i);
    results[i] = Rcpp::sum(x * y) / Rcpp::sum(z);
  }
  return results;
}

