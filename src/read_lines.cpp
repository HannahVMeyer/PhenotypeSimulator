#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <Rcpp.h>

//' Scan file for specific line numbers
//' 
//' @param filename /path/to/chromosomefile [string]
//' @param lines vector of line numbers [integer] to be read
//' @param sep [string] end-of-line delimiter
// [[Rcpp::export]]
Rcpp::CharacterVector read_lines(std::string const& filename, Rcpp::NumericVector const& lines, char sep = '\n') {
    std::ifstream file{filename.c_str()};
    Rcpp::CharacterVector ret{lines.size(), NA_STRING};
    std::string line;
    int lineno = 1;
    unsigned int next_index = 0;
    unsigned int ret_index = 0;
    while (std::getline(file, line, sep)) {
        if (lineno % 10000 == 0) Rcpp::checkUserInterrupt();
        if (next_index == lines.size()) break;
        if (lineno++ == lines[next_index]) {
            ret[ret_index++] = line;
            ++next_index;
        }
    }

    if (next_index < lines.size()) {
        Rcpp::Function warning{"warning"};
        warning("File contained fewer lines than indicated by arguments.");
    }

    return ret;
}
