#include <Rcpp.h>
#include <vector>
#include <unordered_map>
#include <string>
/*
 * Function that provides the index of treated/control units so we know where to store the vit values in the large vector which we will eventually use for finding Wits
 */
// [[Rcpp::export]]
Rcpp::NumericVector get_vit_index(Rcpp::CharacterVector t_id_key, 
                                  Rcpp::CharacterVector control_treatment_t_ids, 
                                  Rcpp::NumericVector control_treatment_set_nums)
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
/*
 * This function assigns dits, according to the paper
 */
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


//Rcpp::NumericVector sumwits(int nrow, Rcpp::NumericVector vit_vect)
Rcpp::NumericVector sumwits(int nrow, std::vector<double> &vit_vect)
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

/*
 * Both this function and check_treated_units carry out the dirty work of checking from t-1 to t + lead to see if data is mising. This function looks at control units. 
 * It returns a list of logical vectors, indicating which control units should be dropped/kept in a particular matched set for the following calculations
 */
// [[Rcpp::export]]
Rcpp::List re_norm_index(Rcpp::NumericMatrix compmat, Rcpp::NumericVector compmat_row_units, Rcpp::NumericVector compmat_cols, int lead, 
                         Rcpp::List sets, Rcpp::NumericVector control_start_years)
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
			for (int k = 0; k <= lead; ++k)// init. k to 1, go to lead + 1, and have first check be for t - 1
			{
				//add in separate check for t - 1, rest of these checks should still apply then.
				if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col - 1)))
				{
				  // Rcpp::Rcout << idx << std::endl;
				  // Rcpp::Rcout << st_year_col - 1<< std::endl;
				  // Rcpp::Rcout << compmat(idx, st_year_col - 1)<< std::endl;
				  v[0] = false;
				  
				}
				else
				{
				  // Rcpp::Rcout << idx << std::endl;
				  // Rcpp::Rcout << st_year_col - 1<< std::endl;
				  // Rcpp::Rcout << compmat(idx, st_year_col - 1)<< std::endl;
				  v[0] = true;
				}
				
				if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col + k)))
				{
				  // Rcpp::Rcout << idx << std::endl;
				  // Rcpp::Rcout << st_year_col + k<< std::endl;
				  // Rcpp::Rcout << compmat(idx, st_year_col + k)<< std::endl;
					v[k + 1] = false;
				}
				else
				{
				  // Rcpp::Rcout << idx << std::endl;
				  // Rcpp::Rcout << st_year_col + k<< std::endl;
				  // Rcpp::Rcout << compmat(idx, st_year_col + k)<< std::endl;
					v[k + 1] = true;
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


/*
 * A slightly less complicated version of the re_norm_index function, which looks at treated units for missing data. It returns a logical vector, indicating which treated units can/cannot be used.
 * The entire matched set corresponding with a treated unit with missing data must be dropped entirely from the calculations.
 */

// [[Rcpp::export]]
Rcpp::LogicalVector check_treated_units(Rcpp::NumericMatrix compmat, Rcpp::NumericVector compmat_row_units, 
                                        Rcpp::NumericVector compmat_cols, int lead, Rcpp::NumericVector treated_ids, Rcpp::NumericVector treated_ts)
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
		for (int k = 0; k <= lead; ++k)
		{
			//checking t-1
			if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col - 1))) //assuming that the order holds
			{
			  // Rcpp::Rcout << idx << std::endl;
			  // Rcpp::Rcout << st_year_col - 1<< std::endl;
			  // Rcpp::Rcout << compmat(idx, st_year_col - 1)<< std::endl;
			  v[0] = false;
			}
			else
			{
			  // Rcpp::Rcout << idx << std::endl;
			  // Rcpp::Rcout << st_year_col - 1<< std::endl;
			  // Rcpp::Rcout << compmat(idx, st_year_col - 1)<< std::endl;
			  v[0] = true;
			}
			
			if( ( (st_year_col + k) > compmat_cols.size() ) || Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col + k)))
			{
			  // Rcpp::Rcout << idx << std::endl;
			  // Rcpp::Rcout << st_year_col + k<< std::endl;
			  // Rcpp::Rcout << compmat(idx, st_year_col + k)<< std::endl;
				v[k + 1] = false;
			}
			else
			{
			  // Rcpp::Rcout << idx << std::endl;
			  // Rcpp::Rcout << st_year_col + k<< std::endl;
			  // Rcpp::Rcout << compmat(idx, st_year_col + k)<< std::endl;
				v[k + 1] = true;
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
/*
 * 
 * this function provides us an index of which matched sets need to have their weights updated. It just iterates over the more detailed index found in earlier functions so that we aren't
 * pointlessly iterating over matched sets that don't need to be adjusted.
 */
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

//not deprecated necessarily, but is slower than the R implementation currently
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




/*
 * Function that provides the index of treated/control units so we know where to store the vit values in the large vector which we will eventually use for finding Wits
 */
// [[Rcpp::export]]
std::vector<unsigned int> get_vit_index_unsigned(Rcpp::CharacterVector t_id_key, 
                                                 Rcpp::CharacterVector control_treatment_t_ids, 
                                                 Rcpp::NumericVector control_treatment_set_nums)
{
  std::unordered_map<std::string, int> indexMap;
  for (int i = 0; i < t_id_key.size(); ++i)
  {
    std::string key;
    key = t_id_key[i];
    indexMap[key] = i + 1;
  }
  std::vector<unsigned int> intdex(control_treatment_t_ids.size());
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
Rcpp::NumericVector handle_vits(unsigned int nrow_data, unsigned int mset_size, 
                                unsigned int num_empty, Rcpp::NumericVector weights,
                                Rcpp::CharacterVector tidkey,
                                Rcpp::CharacterVector control_treatment_tids,
                                Rcpp::NumericVector ct_set_nums)
{
  
  std::vector<unsigned int> idxs = get_vit_index_unsigned(tidkey, control_treatment_tids, ct_set_nums);
  unsigned vec_size = (nrow_data * mset_size) - num_empty;
  std::vector<double> vit_vect(vec_size);
  for(int i = 0; i < idxs.size(); i++)
  {
    vit_vect[idxs[i] - 1] = weights[i];
  }
  // vit_vect[idxs] = weights;
  return sumwits(nrow_data, vit_vect);
}




