#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

struct ScoreRecord {
    double score;
    std::string more_info;
};

typedef std::map<std::string, ScoreRecord> DataBase;
typedef std::vector<double> DoubleVector;
typedef std::tuple<double, std::string, std::string>  ScoreConfig;
typedef std::vector<ScoreConfig> ScoreConfigVec;

DataBase db;


// [[Rcpp::export]]
void add_configuration_score_to_db(std::string str_key, double value, std::string mInfo){
    if (std::isnan(value)) {
        Rcpp::stop("configuration scores cannot be NA or NaN");
    }
    db[str_key] = ScoreRecord{value, mInfo};
}


struct Compare{
    bool operator()(const ScoreConfig& lhs, const ScoreConfig& rhs) const {
        const double lhs_score = std::get<0>(lhs);
        const double rhs_score = std::get<0>(rhs);
        const bool lhs_nan = std::isnan(lhs_score);
        const bool rhs_nan = std::isnan(rhs_score);
        if (lhs_nan != rhs_nan) {
            return !lhs_nan;
        }
        if (!lhs_nan) {
            if (lhs_score < rhs_score) {
                return true;
            }
            if (rhs_score < lhs_score) {
                return false;
            }
        }
        return std::get<1>(lhs) < std::get<1>(rhs);
    }
} myCompare;

// [[Rcpp::export]]
Rcpp::List get_stored_config_score(){

    ScoreConfigVec config_entries;
    config_entries.reserve(db.size());
    std::size_t interrupt_counter = 0;
    for (const auto& entry : db) {
        config_entries.emplace_back(
            entry.second.score, entry.first, entry.second.more_info
        );
        if (((++interrupt_counter) & 0xffffU) == 0U) {
            Rcpp::checkUserInterrupt();
        }
    }
    std::sort(config_entries.begin(), config_entries.end(), myCompare);

    std::vector<std::string> configVec, moreInfoVec;
    DoubleVector doubleVec;
    doubleVec.reserve(config_entries.size());
    configVec.reserve(config_entries.size());
    moreInfoVec.reserve(config_entries.size());
    for(auto itr=config_entries.begin(); itr !=config_entries.end(); ++itr)
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
    db.clear();
}
 
// [[Rcpp::export]]
Rcpp::List get_score_of_configuration(std::string str_key){
    bool   valid = 0;
    double value = 0;

    auto itr = db.find(str_key);
    if(  itr != db.end() ){
        valid = 1;
        value = itr->second.score;
    }
    return( Rcpp::List::create( 
                Rcpp::Named("value") = value,
                Rcpp::Named("valid") = valid ) );
}
