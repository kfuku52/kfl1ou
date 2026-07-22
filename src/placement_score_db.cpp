#include <algorithm>
#include <vector>
#include <map>

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

typedef std::map<std::string, double> DataBase;
typedef std::vector<double> DoubleVector;
typedef std::tuple<double, std::string, std::string>  ScoreConfig;
typedef std::vector<ScoreConfig> ScoreConfigVec;

DataBase db;
ScoreConfigVec myConfigVec;


// [[Rcpp::export]]
void add_configuration_score_to_db(std::string str_key, double value, std::string mInfo){
    db[str_key] = value;
    auto existing = std::find_if(
        myConfigVec.begin(), myConfigVec.end(),
        [&str_key](const ScoreConfig& entry){
            return std::get<1>(entry) == str_key;
        }
    );
    if (existing == myConfigVec.end()) {
        myConfigVec.push_back(ScoreConfig(value, str_key, mInfo));
    } else {
        *existing = ScoreConfig(value, str_key, mInfo);
    }
}


struct Compare{
    bool operator()(const ScoreConfig lhs, const ScoreConfig rhs){
        return std::get<0>(lhs) < std::get<0>(rhs);
    }
} myCompare;

// [[Rcpp::export]]
Rcpp::List get_stored_config_score(){

    std::sort(myConfigVec.begin(), myConfigVec.end(), myCompare);
   
    //std::sort(myConfigVec.begin(), myConfigVec.end(),
    //        [](const ScoreConfig& lhs, const ScoreConfig& rhs) {
    //        return std::get<0>(lhs) < std::get<0>(rhs); } );

    std::vector<std::string> configVec, moreInfoVec;
    DoubleVector doubleVec;
    //for(auto &itr : myConfigVec)
    for(auto itr=myConfigVec.begin(); itr !=myConfigVec.end(); ++itr)
    {
        doubleVec.push_back   (std::get<0>(*itr));
        configVec.push_back   (std::get<1>(*itr));
        moreInfoVec.push_back (std::get<2>(*itr));
    }

    return Rcpp::List::create( 
                Rcpp::Named("scores") = doubleVec,
                Rcpp::Named("configurations") = configVec,
                Rcpp::Named("moreInfo") = moreInfoVec
                );
}

// [[Rcpp::export]]
void erase_configuration_score_db(){
    db.erase(db.begin(), db.end());
    myConfigVec.erase(myConfigVec.begin(), myConfigVec.end() );
}
 
// [[Rcpp::export]]
Rcpp::List get_score_of_configuration(std::string str_key){
    bool   valid = 0;
    double value = 0;

    auto itr = db.find(str_key);
    if(  itr != db.end() ){
        valid = 1;
        value = itr->second;
    }
    return( Rcpp::List::create( 
                Rcpp::Named("value") = value,
                Rcpp::Named("valid") = valid ) );
}


