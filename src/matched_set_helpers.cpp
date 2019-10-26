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
    if ( (treated_indices[i] > 0) &
         ( !Rcpp::internal::Rcpp_IsNA(ordered_df(int(treated_indices[i]) - 1, treat_col_idx)) ) & 
         
         ( ordered_df(int(treated_indices[i]) - 1, treat_col_idx) == 0) ) //is the treatmentvar == 0 at time t -1 ?
    {
      if ( (!Rcpp::internal::Rcpp_IsNA(ordered_df(treated_indices[i], unit_var_col)) ) &
           (!Rcpp::internal::Rcpp_IsNA(ordered_df(treated_indices[i] - 1, unit_var_col)) ) & //na checks just to make sure data exists
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
Rcpp::List get_comparison_histories(const Rcpp::NumericMatrix &compmat, const Rcpp::NumericVector &ts, const Rcpp::NumericVector &ids, int t_col, int id_col, int L, int treat_col)
{
  
  Rcpp::List comp_hists(ts.length()); //length of its and ids should be the same
  for (int i = 0; i < ts.length(); i++)
  {
    int t = ts[i];
    int id = ids[i]; //iterating over sets of t, id pairs
    
    for (int j = 0; j < compmat.nrow(); j++) //iterate through "long" form of data
    {
      if( (compmat(j, t_col) == t) & (compmat(j, id_col) == id) ) // if time and unitid of current row matches t, id pair we are currently "investigating"...
      {
        
        Rcpp::NumericVector control_hist(L+1);
        for (int k = 0; k < L+ 1; k++)
        {
          control_hist[k] = compmat(j - L + k, treat_col); 
        } // ...read the treatment history over the window into a vector...
        
        control_hist[control_hist.length() - 1] = 0; //... and change the last entry to give the needed treatment history of a control unit for this t,id pair. Entry here should always be 1 before we change it here
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
        
        if ( (!Rcpp::internal::Rcpp_IsNA(Rcpp::is_true(Rcpp::all(tempcomp == cont_hist)))) & //Do the actual treatment history of a unit match the needed control history? If so...
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
          // if(widemat(j, 0) == 4)
          // {
          //   Rcpp::Rcout << widemat(j, 0) << std::endl;  
          // }
          
          for (int k = 0; k < L + 1; k++)
          {
            tempcomp[k] = widemat(j, t - L + k); //retrieving treatment history for the window of interest
            
            //for(int s = 0; s < tempcomp.length(); s++)
            //Rcpp::Rcout << tempcomp[s]; 
            //tempcomp is the actual history of a unit, cont_hist is what must be matched in order for a unit to be included in a matched set for a given t, id
          }
          // if(widemat(j, 0) == 4)
          // {
          //   for(int s = 0; s < tempcomp.length(); s++)
          //   {
          //     Rcpp::Rcout << tempcomp[s];
          //   }
          //   Rcpp::Rcout << std::endl;
          //   
          // }
          
          Rcpp::NumericVector cont_hist_comp(2);
          
          for(int k = 0; k < 2; k++)
          {
            cont_hist_comp[k] = cont_hist[cont_hist.length() -2 + k];
          }
          
          // if(id == 3 && t == 38)
          // {
          //   for(int xx = 0; xx < cont_hist_comp.length(); xx++)
          //   {
          //     Rcpp::Rcout << cont_hist_comp[xx];
          //   }
          //   Rcpp::Rcout << std::endl;
          // }
          if ( (!Rcpp::internal::Rcpp_IsNA(Rcpp::is_true(Rcpp::all(tempcomp == cont_hist_comp)))) & //Do the actual treatment history of a unit match the needed control history? If so...
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

