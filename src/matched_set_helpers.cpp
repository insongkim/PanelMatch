// #include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <vector>
using namespace Rcpp ;

/*
 ordered_df: this is basically the same as the data matrix provided by the user, but is ordered by t and unit id
 treated_indices: just a list of the indices (row numbers) of units that have received treatment. should be zero indexed, matching the ordered_df here
 treat_col_idx: column number corresponding to the treatment variable in ordered_df
 unit_var_col: column number corresponding to the unitid variable in ordered_df
 Function returns a logical vector that is used to identify t/id pairs that might have potential matched sets. 
 More specifically, the logical vector corresponds to row numbers in the ordered_df for which at time t, treatmnt = 1, but =0 at time t -1
 */
// [[Rcpp::export()]]
Rcpp::LogicalVector get_treated_indices(const Rcpp::NumericMatrix &ordered_df, const Rcpp::NumericVector &treated_indices, int treat_col_idx, int unit_var_col)
{
  
  Rcpp::LogicalVector treateds(treated_indices.size());
  
  for (int i = 0; i < treated_indices.size(); i++) //iterating throw each treated unit
  {
    if ( (treated_indices[i] > 0) &&
         ( !Rcpp::internal::Rcpp_IsNA(ordered_df(int(treated_indices[i]) - 1, treat_col_idx)) ) && 
         
         ( ordered_df(int(treated_indices[i]) - 1, treat_col_idx) == 0) ) //is the treatmentvar == 0 at time t -1 ?
    {
      if ( (!Rcpp::internal::Rcpp_IsNA(ordered_df(treated_indices[i], unit_var_col)) ) &&
           (!Rcpp::internal::Rcpp_IsNA(ordered_df(treated_indices[i] - 1, unit_var_col)) ) && //na checks just to make sure data exists
           (ordered_df(treated_indices[i] -1, unit_var_col) == ordered_df(treated_indices[i], unit_var_col) ) ) //do the unit ids matched at time t and time t - 1, mostly a safeguard against weirdly formed data and edgecases
      {
        treateds[i] = 1; //if treatedvar == 1 at t and == 0 and t -1, the data exists and doesn't have obvious formatting problems, then it should be included in the set of treated units to look at. 
        //these are included as treated units
      }
      
    }
  }
  return treateds;//then we can use this to subset our vector in R
}


/*
 compmat: matrix similar to the ordered_df described above, only contains information about t, unit id, and treatment. this should be sorted.
 ts: vector corresponding to the "t" of a particular treated unit
 ids: vector corresponding to the "ids" of a particular treated unit. This vector should be the same length as 'ts'.
 t_col: column number corresponding to the time variable 
 id_col: same as t_col, but for unit id
 treat_col: same as t_col, but for the treatment variable. 
 L: lag window size
 Function returns a list of vectors. Each vector represents the treatment history needed for a unit to be included in the matched set for the unit corresponding to the t, id pair for that entry. This list should be the
 same length as the 'ts' and 'ids' vectors passed in initially and able to be matched together by index number.
 */
// [[Rcpp::export()]]
Rcpp::List get_comparison_histories(const Rcpp::NumericMatrix &compmat, 
                                    const Rcpp::NumericVector &ts, 
                                    const Rcpp::NumericVector &ids, 
                                    int t_col, int id_col, int L, int treat_col,
                                    bool atc)
{
  
  Rcpp::List comp_hists(ts.length()); //length of its and ids should be the same
  for (int i = 0; i < ts.length(); i++)
  {
    int t = ts[i];
    int id = ids[i]; //iterating over sets of t, id pairs
    
    for (int j = 0; j < compmat.nrow(); j++) //iterate through "long" form of data
    {
      if( (compmat(j, t_col) == t) && (compmat(j, id_col) == id) ) // if time and unitid of current row matches t, id pair we are currently "investigating"...
      {
        
        Rcpp::NumericVector control_hist(L+1);
        for (int k = 0; k < L+ 1; k++)
        {
          control_hist[k] = compmat(j - L + k, treat_col); 
        } // ...read the treatment history over the window into a vector...
        if (!atc)
        {
          control_hist[control_hist.length() - 1] = 0; //... and change the last entry to give the needed treatment history of a control unit for this t,id pair. Entry here should always be 1 before we change it here
        } else { // is atc
          control_hist[control_hist.length() - 1] = 1;
        }
        
        comp_hists[i] = control_hist;
        break;
      }
      
    }
  }
  return comp_hists;
}

/*
 control_history_list: result of above get_comparison_histories C++ function
 widemat: "wide" form of data matrix with a row for each unit id and columns corresponding to every known t. Entries are 1 if unit is treated at that time, 0 if not
 t_as_col_nums: integers corresponding to the column number of the t for a (t,id) pair for which we are trying to find a matched set.
 ids: vector containing ids for treated units for which we are attempting to find matched sets
 L: lag window size
 returns a list of vectors, vectors will contain unit ids for units included in a matched set. Size should correspond to the length of t_as_col_nums vector and ids vector, although some sets might be empty.
 */
// [[Rcpp::export()]]
Rcpp:: List get_msets_helper(const Rcpp::List &control_history_list, const Rcpp::NumericMatrix &widemat, const Rcpp::NumericVector &t_as_col_nums, const Rcpp::NumericVector &ids, int L)
{
  Rcpp::List matched_sets(ids.length());
  Rcpp::NumericVector units = widemat(_, 0); //assuming that the ids are in the first column of this matrix, according to expectations
  
  for (int i = 0; i < ids.length(); ++i)
  {	
    //ids and t_as_col_nums should be the same length...parallel vectors
    //should also be the same legnth as the control_history_list
    int id = ids[i]; //just for clarity/readability, t, id represent the current t, id pair for a unit at a time we want to find a matched set for
    int t = t_as_col_nums[i];
    Rcpp::NumericVector cont_hist = control_history_list[i]; //also mostly for readability, just storing the current needed treatment history in a variable
    
    Rcpp::LogicalVector in_matched_set_idx(widemat.nrow()); // will use this to indicate which units should be included in the matched set.
    
    // last index should be of length n, where n is the number of units, should be same length as units vector created earlier
    
    for (int j = 0; j < widemat.nrow(); j++)
    {
      if (widemat(j, 0) != id) //do nothing if we are looking at the row of the unit we are looking for matched sets for
      {
        
        Rcpp::NumericVector tempcomp(L + 1);
        
        
        for (int k = 0; k < L + 1; k++)
        {
          // if(t - L < 1)
          // {
          // 	Rcpp::stop('Inadmissable Lag window valuee with given t value')
          // }
          tempcomp[k] = widemat(j, t - L + k); //retrieving treatment history for the window of interest
          //tempcomp is the actual history of a unit, cont_hist is what must be matched in order for a unit to be included in a matched set for a given t, id
        }
        
        if ( (!Rcpp::internal::Rcpp_IsNA(Rcpp::is_true(Rcpp::all(tempcomp == cont_hist)))) && //Do the actual treatment history of a unit match the needed control history? If so...
             Rcpp::is_true(Rcpp::all(tempcomp == cont_hist)) ) // checking that NOT na might be redundant, but also might prevent bug
        {
          in_matched_set_idx[j] = true; //... then that unit should be included in the matched set.
        } //otherwise, leave the default, which is false
      }
    }
    matched_sets[i] = units[in_matched_set_idx];
  }
  return matched_sets;
}

/*
 control_history_list: result of above get_comparison_histories C++ function
 widemat: "wide" form of data matrix with a row for each unit id and columns corresponding to every known t. Entries are 1 if unit is treated at that time, 0 if not
 t_as_col_nums: integers corresponding to the column number of the t for a (t,id) pair for which we are trying to find a matched set.
 ids: vector containing ids for treated units for which we are attempting to find matched sets
 L: lag window size
 returns a list of vectors, vectors will contain unit ids for units included in a matched set. Size should correspond to the length of t_as_col_nums vector and ids vector, although some sets might be empty.
 */
// [[Rcpp::export()]]
Rcpp:: List non_matching_matcher(const Rcpp::List &control_history_list, 
                                 const Rcpp::NumericMatrix &widemat, const Rcpp::NumericVector &t_as_col_nums, 
                                 const Rcpp::NumericVector &ids, int L, int missing_window) // in this case, L will always be one
{
  Rcpp::List matched_sets(ids.length());
  Rcpp::NumericVector units = widemat(_, 0); //assuming that the ids are in the first column of this matrix, according to expectations
  
  for (int i = 0; i < ids.length(); ++i)
  {	
    //ids and t_as_col_nums should be the same length...parallel vectors
    //should also be the same legnth as the control_history_list
    int id = ids[i]; //just for clarity/readability, t, id represent the current t, id pair for a unit at a time we want to find a matched set for
    int t = t_as_col_nums[i];
    Rcpp::NumericVector cont_hist = control_history_list[i]; //also mostly for readability, just storing the current needed treatment history in a variable
    
    Rcpp::LogicalVector in_matched_set_idx(widemat.nrow()); // will use this to indicate which units should be included in the matched set.
    
    // last index should be of length n, where n is the number of units, should be same length as units vector created earlier
    
    for (int j = 0; j < widemat.nrow(); j++)
    {
      if (widemat(j, 0) != id) //do nothing if we are looking at the row of the unit we are looking for matched sets for
      {
        
        Rcpp::NumericVector na_tempcomp(missing_window + 1);
        for (int k = 0; k < missing_window + 1; k++)
        {
          na_tempcomp[k] = widemat(j, t - missing_window + k); 
          
        }
        if(Rcpp::all(!Rcpp::is_na(na_tempcomp)))
        {
          Rcpp::NumericVector tempcomp(L + 1);
          
          
          for (int k = 0; k < L + 1; k++)
          {
            tempcomp[k] = widemat(j, t - L + k); //retrieving treatment history for the window of interest
        
          }
          
          
          Rcpp::NumericVector cont_hist_comp(2);
          
          for(int k = 0; k < 2; k++)
          {
            cont_hist_comp[k] = cont_hist[cont_hist.length() -2 + k];
          }
          
          if ( (!Rcpp::internal::Rcpp_IsNA(Rcpp::is_true(Rcpp::all(tempcomp == cont_hist_comp)))) && //Do the actual treatment history of a unit match the needed control history? If so...
               Rcpp::is_true(Rcpp::all(tempcomp == cont_hist_comp)) ) // checking that NOT na might be redundant, but also might prevent bug
          {
            in_matched_set_idx[j] = true; //... then that unit should be included in the matched set.
          }
        }
      }
    }
    matched_sets[i] = units[in_matched_set_idx];
  }
  return matched_sets;
}


// [[Rcpp::export]]
Rcpp::List filter_placebo_results(Rcpp::NumericMatrix expanded_data,
                                  Rcpp::NumericVector ordered_outcome_data,
                                  Rcpp::NumericVector treated_ids,
                                  Rcpp::NumericVector treated_ts,
                                  Rcpp::List sets,
                                  int lag) {
  
  //creating the mapping of key to index for easier lookups without search
  std::unordered_map<std::string, int> indexMap;
  Rcpp::List subsets(treated_ids.size());
  for(int i = 0; i < expanded_data.nrow(); i++)
  {
    int id_1 = expanded_data(i,0);
    int t_1 = expanded_data(i,1);
    
    std::string id = std::to_string(id_1);
    std::string t = std::to_string(t_1);
    
    std::string key = id + "." + t;
    
    indexMap[key] = i;
    
  }
  
  for(int i = 0; i < sets.size(); i++) //iterating over the matched sets
  {
    
    
    //check treated unit
    int id_t = treated_ids[i];
    std::string id = std::to_string(id_t);
    
    int t = treated_ts[i];
    bool check_controls = true;
    for (int j = lag; j > 0; j--)
    {
      int new_time = t - j;
      std::string xx = std::to_string(new_time);
      std::string key = id + "." + xx;

      if(Rcpp::internal::Rcpp_IsNA(ordered_outcome_data[indexMap[key]]))
      {
        check_controls = false;
      }
    }
    //check control units
    
    Rcpp::NumericVector control_ids = sets[i];
    Rcpp::LogicalVector keep(control_ids.size()); //default to false
    
    if (check_controls)
    {
      for (int y = 0; y < control_ids.size(); y++) //iterating over the controls in a particular matched set
      {
        
        int ctrl = control_ids[y];
        std::string id = std::to_string(ctrl);
        keep[y] = true;
        for (int j = lag; j > 0; j--) //iterating over the window for each control unit
        {
          int new_time = t - j;
          std::string xx = std::to_string(new_time);
          std::string key = id + "." + xx;
          
          if(Rcpp::internal::Rcpp_IsNA(ordered_outcome_data[indexMap[key]]))
          {
            keep[y] = false;
          } 
        }
      }
      subsets[i] = control_ids[keep];
    } else 
    {
      subsets[i] = control_ids[keep];
    }
    
  }
  return subsets; //then just filter empty sets when they're back.
  
}


// [[Rcpp::export]]
Rcpp::List enforce_caliper(Rcpp::NumericMatrix expanded_data,
                           Rcpp::NumericVector unique_unit_ids,
                           int idx_variable_to_check,
                           Rcpp::IntegerVector treated_ids,
                           Rcpp::IntegerVector treated_ts,
                           double caliper_level,
                           int lag) 
{
  
  //creating the mapping of key to index for easier lookups without search
  std::unordered_map<std::string, int> indexMap;
  for(int i = 0; i < expanded_data.nrow(); i++)
  {
    int id_1 = expanded_data(i,0);
    int t_1 = expanded_data(i,1);
    
    std::string id = std::to_string(id_1);
    std::string t = std::to_string(t_1);
    std::string key = id + "." + t;
    indexMap[key] = i;
    
  }
  
  //  we want to check against every plausible control unit, so we will exclude other treated units eventually.
  // here, we just create the keys which we will use to filter out later
  Rcpp::CharacterVector treated_keys(treated_ids.size());
  for(int i = 0; i < treated_ids.size(); i++) 
  {
    treated_keys[i] = std::to_string(treated_ids[i]) + "." + std::to_string(treated_ts[i]);
  }
  Rcpp::List matched_sets(treated_ids.size());
  for(int i = 0; i < treated_ids.size(); i++) //iterating over treated observations
  {
    // get treated unit history
    int id_t = treated_ids[i];
    std::string id = std::to_string(id_t);
    
    int t = treated_ts[i];
    // bool check_controls = true;
    // build treatment history to compare over from t-lag to t-1
    // todo: deal with instances where there are NAs in the treatment history
    Rcpp::NumericVector treated_history(lag);
    for (int j = lag; j > 0; j--)
    {
        int new_time = t - j;
        std::string xx = std::to_string(new_time);
        std::string key = id + "." + xx;
        treated_history[lag-j] = expanded_data(indexMap[key], idx_variable_to_check);
        
    }
    //   //check against control units
    //   // we need to check against every other non-treated unit (since we cannot match treated units to other treated units)
    // Rcpp::CharacterVector unit_ids = expanded_data_, unit_id_variable];
    // 
    // // Find the unique values in 'unit.ids'
    // Rcpp::CharacterVector unique_unit_ids = unique(unit_ids);
    
    
    int n = unique_unit_ids.size();
    Rcpp::CharacterVector potential_match_start_keys(n);
    Rcpp::LogicalVector matched_control_idx(n);
    // append the year of treatment from the treated observation. 
    // We can use this to quickly look up potential controls
    for(int w = 0; w < n; w++) {
      std::string a;
      int ap = unique_unit_ids[w];
      a = std::to_string(ap);
      potential_match_start_keys[w] = a + "." + std::to_string(t);
    }
    
    // now we need to loop over the potential matches, extract treatment history, and compare
    for(int z = 0; z < potential_match_start_keys.size(); z++)
    {
      // if the potential matched unit is one of the treated observations, filter it out before checking anything
      bool is_viable = true;
      // Rcpp::Rcout << treated_keys;
      // Rcpp::Rcout << potential_match_start_keys;
      for (int q = 0; q < treated_keys.size(); ++q) {
        if (treated_keys[q] == potential_match_start_keys[z]) 
        { // z is looping over potential starting points associated with potential control units. q is iterating over the treated observations. the condition is flipped if we find a treated observation that matches one of the potential controls
          is_viable = false;
        }
      }
      if (!is_viable)
      {
        matched_control_idx[z] = false;
      }
      else 
      {
        // at this point, the candidate is still potentially able to be matched
        Rcpp::NumericVector control_history(lag);
        std::string jl;
        jl = potential_match_start_keys[z];
        int jlx = indexMap[jl]; // this time corresponds to time t
        int startidx =  jlx - lag; // since we want t-lag to t-1, this gives us t-lag
        for (int j = 0; j < lag; j++) // again, deal with NAs
        { //extracting the treatment history for a candidate control unit
          int rowidx = (startidx + j);
          control_history[j] = expanded_data(rowidx, idx_variable_to_check);
        }
        //Rcpp::NumericVector result(control_history.size());
        bool match_control = true;
        for (int j = 0; j < control_history.size(); j++)
        {
          bool cond1 = Rcpp::NumericVector::is_na(treated_history[j]);
          bool cond2 = Rcpp::NumericVector::is_na(control_history[j]);
          
          if (cond1 || cond2)
          {
            match_control = false;
            break;
          } 
          else 
          {
            double r = abs(treated_history[j] - control_history[j]);
            if (r > caliper_level)
            {
              match_control = false;
              break;
            }
          }
          
          
        } 
        if (match_control)
        {
          matched_control_idx[z] = true;
        }
      }
        
    }
    matched_sets[i] = unique_unit_ids[matched_control_idx];
    // potential match start keys and unique unit ids should be parallel vectors, but we want to return integers without the appended year data if possible
    //matched_sets[i] = potential_match_start_keys[matched_control_idx];
  }
  return matched_sets;
}