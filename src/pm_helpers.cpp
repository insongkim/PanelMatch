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

