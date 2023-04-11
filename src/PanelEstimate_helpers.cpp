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
 * 
 * this function provides us an index of which matched sets need to be adjusted because they have missing data
 */
// [[Rcpp::export]]
Rcpp::LogicalVector needs_adjustment(Rcpp::List set_index_list)
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




