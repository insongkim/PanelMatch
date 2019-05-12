//#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <RcppArmadillo.h>
using namespace Rcpp;


/*
This function will return a list of lists. The outer list length is equal to the number of matched sets. Each nested list is of length lag + 1 and each element in this nested list 
 contains indices of control units in the provided data set for a particular year in the lag window
*/ 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List get_yearly_dmats(NumericMatrix expanded_data, NumericVector treated_ids, List ts_to_fetch, CharacterVector row_key, List matched_sets, int lag) 
{
  std::unordered_map<std::string, int> indexMap;
  for(int i = 0; i < row_key.size(); i++)
  {
    std::string key;
    key = row_key[i];
    indexMap[key] = i;
  } 
  
  List df_list(matched_sets.size());
  for(int i = 0; i < matched_sets.size(); i++)
  {
    NumericVector matched_set = matched_sets[i];
    
    List temp_list(lag + 1);
    NumericVector tts = ts_to_fetch[i];
    //NumericVector treated_times_repeated = treated_times[i];
    for(int k = 0; k < tts.size(); k++) // should be equal to lag + 1, same as treated_unit_repeated.size()
    {
      NumericVector index_vector(matched_set.size() + 1); //!!!???
      //int t = treated_unit_repeated[k];
      int t = tts[k];
      std::string tstring = std::to_string(t);
      for(int j = 0; j < matched_set.size(); j++)
      {
        int temp = matched_set[j];
        std::string tempstring = std::to_string(temp);
        std::string lookup = tempstring + "." + tstring;
        index_vector[j] = indexMap[lookup];
      }
      int t_id = treated_ids[i];
      std::string t_id_string = std::to_string(t_id);
      std::string t_lookup = t_id_string + "." + tstring;
      index_vector[index_vector.size() - 1] = indexMap[t_lookup];
      temp_list[k] = index_vector + 1;
    }
    
    df_list[i] = temp_list;
    
  }
  
  return df_list;
}

// [[Rcpp::export]]
Rcpp::LogicalVector check_treated_units_for_treatment_reversion(Rcpp::NumericMatrix compmat, Rcpp::NumericVector compmat_row_units, 
                                                                Rcpp::NumericVector compmat_cols, int lead, Rcpp::NumericVector treated_ids,
                                                                Rcpp::NumericVector treated_ts)
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
      if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col - 1)) || 
         compmat(idx, st_year_col - 1) != 0) //assuming that the order holds
      {
        v[0] = false;
      }
      else
      {
        v[0] = true;
      }
      
      if( ( (st_year_col + k) > compmat_cols.size() ) || Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col + k)) 
            || compmat(idx, st_year_col + k) != 1 )
      {
        v[k + 1] = false;
      }
      else
      {
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

// [[Rcpp::export]]
Rcpp::List check_control_units_for_treatment_restriction(Rcpp::NumericMatrix compmat, Rcpp::NumericVector compmat_row_units, 
                                                         Rcpp::NumericVector compmat_cols, int lead, Rcpp::List sets, 
                                                         Rcpp::NumericVector control_start_years)
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
        if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col - 1)) || compmat(idx, st_year_col - 1) != 0)
        {
          v[0] = false;
        }
        else
        {
          v[0] = true;
        }
        
        if(Rcpp::internal::Rcpp_IsNA(compmat(idx, st_year_col + k)) ||
           compmat(idx, st_year_col + k) != 0)
        {
          v[k + 1] = false;
        }
        else
        {
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

// [[Rcpp::export]]
Rcpp::List multiply_weights_msm(Rcpp::List weights, int number_of_sets)
{
  Rcpp::List final_weights(number_of_sets);
  for(int i = 0; i < number_of_sets; i++)
  {
    Rcpp::NumericVector base_mult = weights[i];
    Rcpp::NumericVector temp2;
    for(int j = i; j < weights.size(); j+=number_of_sets)
    {
      if(j != i)
      {
        Rcpp::NumericVector temp = weights[j];
        temp2 = base_mult * temp;
        base_mult = temp2;
        //base_mult = base_mult * weights[j];
      }

    }
    final_weights[i] = base_mult;
  }
  return final_weights;
}


// [[Rcpp::export]]
Rcpp::List do_exact_matching_refinement(Rcpp::NumericMatrix balanced_data, 
                                        int lag, Rcpp::CharacterVector row_key,
                                        Rcpp::List control_data, Rcpp::CharacterVector treatment_data,
                                        Rcpp::IntegerVector exact_match_variable_column_index)
{
  
  std::unordered_map<std::string, int> indexMap;
  for(int i = 0; i < row_key.size(); i++)
  {
    std::string key;
    key = row_key[i];
    indexMap[key] = i;
  }
  Rcpp::List exact_match_control_lists(exact_match_variable_column_index.size());
  for(int i = 0; i < exact_match_variable_column_index.size(); i++)
  {
    int idx = exact_match_variable_column_index[i];
    Rcpp::List control_idx(treatment_data.size());
    for(int j = 0; j < treatment_data.size(); j ++)
    {
      std::string key;
      key = treatment_data[j];
      int endpoint = indexMap[key];
      Rcpp::NumericVector temptreatmentdata(lag + 1);
      for(int k = 0; k < lag + 1; k++)
      {
        int idxcheck = endpoint - k;
        // Rcpp::Rcout << idxcheck << std::endl;
        // Rcpp::Rcout << idx << std::endl;
        temptreatmentdata[k] = balanced_data(idxcheck, idx);
        //Rcpp::Rcout << balanced_data[idxcheck, 0] << std::endl <<balanced_data[idxcheck, 1] <<std::endl;
      }
      Rcpp::CharacterVector control_strings = control_data[j];
      Rcpp::LogicalVector keep_control_idx(control_strings.size());
      for(int z = 0; z < control_strings.size(); z++)
      {
        std::string ckey;
        ckey = control_strings[z];
        int cendpoint = indexMap[ckey];
        Rcpp::NumericVector tempcontroldata(lag+1);
        for(int a = 0; a < lag + 1; a++)
        {
          int cidxcheck = cendpoint - a;
          // Rcpp::Rcout << balanced_data[cidxcheck, idx] << std::endl;
          tempcontroldata[a] = balanced_data(cidxcheck, idx);
          //Rcpp::Rcout << balanced_data(cidxcheck, 0) << std::endl << balanced_data(cidxcheck, 1) << std::endl << balanced_data(cidxcheck, 2) << std::endl;
        }
        // for(int q = 0; q < tempcontroldata.size(); q++)
        // {
        //   Rcpp::Rcout << tempcontroldata[q];
        // }
        // Rcpp::Rcout << std::endl;
        // for(int q = 0; q < tempcontroldata.size(); q++)
        // {
        //   Rcpp::Rcout << temptreatmentdata[q];
        // }
        //Rcpp::Rcout << std::endl;
        keep_control_idx[z] = Rcpp::is_true(Rcpp::all(tempcontroldata == temptreatmentdata));
      }
      // Rcpp::Rcout << keep_control_idx.size() << std::endl;
      control_idx[j] = keep_control_idx;
    }
    exact_match_control_lists[i] = control_idx;
  }
  return(exact_match_control_lists); //then will need to take the combination of all of these (everything in corresponding indices must be TRUE)
}




