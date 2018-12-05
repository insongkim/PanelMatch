//#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List create_dmats_for_distance(NumericMatrix expanded_data, NumericVector treated_ids, NumericVector treated_ts, CharacterVector row_key, List matched_sets) 
{
  std::unordered_map<std::string, int> indexMap;
  for(int i = 0; i < row_key.size(); i++)
  {
    std::string key;
    key = row_key[i];
    indexMap[key] = i;
  } //build mapping from t/id pair to row in the expanded and sorted data matrix for easy retrieval.
  
  List df_list(matched_sets.size());
  for(int i = 0; i < matched_sets.size(); i++)
  {
    NumericVector matched_set = matched_sets[i];
    NumericVector index_vector(matched_set.size() + 1);
    int t = treated_ts[i];
    std::string tstring = std::to_string(t);
    for(int j = 0; j < matched_set.size(); j++)
    {
      int temp = matched_set[j];
      std::string tempstring = std::to_string(temp);
      std::string lookup = tempstring + "." + tstring; // might be able to omit the ., but also might keep in for consistency
      index_vector[j] = indexMap[lookup];
    }
    int t_id = treated_ids[i];
    std::string t_id_string = std::to_string(t_id);
    std::string t_lookup = t_id_string + "." + tstring;
    index_vector[index_vector.size() - 1] = indexMap[t_lookup];
    df_list[i] = index_vector + 1;

  }

  return df_list;
}

