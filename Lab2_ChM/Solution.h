#pragma once
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include "Operations.h"


class Solution {
private:
	static void change_matrix_to_itereable(std::vector<std::vector<double>>& matrix);
	static double get_q_from_Yakobi(const std::vector<std::vector<double>>& matrix);
	static std::vector<double> find_roots(const std::vector<std::vector<double>>& mat, const std::vector<double>& b, bool transpose);
public:
	Solution() = default;
	static std::vector<double> Simple_iteration(const double& approx);
	static std::vector<double> Yakobi(const double& approx);
	static std::vector<double> Square_root();
};