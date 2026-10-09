#ifndef __MATH_FUNCTIONS__
#define __MATH_FUNCTIONS__
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>

double return_percentile_value(std::vector<double> values_table, double percentile);
void computeMedianAndMAD(std::vector<double> values, double &med, double &mad);
double computeThresholdMAD(std::vector<double> values, double facteur_mad);
#endif
