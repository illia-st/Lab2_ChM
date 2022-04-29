#include "Solution.h"
#include <stdexcept>
#include <iostream>

template <class T>
std::ostream& operator << (std::ostream& out, const std::vector<std::vector<T>>& mat) {
	for (const auto& i : mat) {
		for (const auto& j : i) {
			out << j << " ";
		}
		out << std::endl;
	}
	return out;
}

void Solution::change_matrix_to_itereable(std::vector<std::vector<double>>& matrix)
{
	size_t n = matrix.size(), m = matrix[0].size();
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < m - 1; ++j) {
			if (i != j) {
				matrix[i][j] = -matrix[i][j];
			}
		}
	}
	return;
}

double Solution::get_q_from_Yakobi(const std::vector<std::vector<double>>& matrix)
{
	size_t n = matrix.size(), m = matrix[0].size();
	if (n != m - 1) {
		throw std::domain_error("The matrix isn't a square one");
	}
	std::pair<size_t, std::pair<double, double>> min_difference = { 0, {DBL_MAX, 0} };
	for (size_t i = 0; i < n; ++i) {
		if (matrix[i][i] == 0) {
			throw std::domain_error("One if the diagonal items is 0");
		}
		double sum = 0;
		for (size_t k = 0; k < m - 1; ++k) {
			if (k != i) {
				sum += std::fabs(matrix[i][k]);
			}
		}
		if (matrix[i][i] <= sum) {
			throw std::domain_error("Yakobi method can't solve this SLAR");
		}
		double dif = std::fabs(matrix[i][i] - sum);
		min_difference = dif < min_difference.second.first ? std::make_pair(i, std::make_pair(dif, sum)) : min_difference;
	}
	double q = min_difference.second.second / std::fabs(matrix[min_difference.first][min_difference.first]);

	return q;
}

std::vector<double> Solution::Simple_iteration(const double& approx) {
	std::vector<std::vector<double>> matrix = {// не забути змінити матрицю на нормальну
		{4, 0, 1, 0, 12},
		{0, 3, 0, 2, 19},
		{1, 0, 5, 1, 27},
		{0, 2, 1, 4, 30}
	};
	std::cout << "Your matrix is: \n";
	std::cout << matrix;
	Solution::change_matrix_to_itereable(matrix);

	std::vector<double> first_sol(4), second_sol(4);
	for (size_t i = 0; i < matrix.size(); ++i) {
		double divider = matrix[i][i];
		for (size_t j = 0; j < matrix[0].size(); ++j) {
			matrix[i][j] /= divider;
			if (j == matrix[0].size() - 1) {
				first_sol[i] = matrix[i][j];
			}
		}
	}
	int iterations = 1;
	while (true) {
		for (size_t i = 0; i < matrix.size(); ++i) {
			double sol = matrix[i][matrix[0].size() - 1];
			for (size_t j = 0; j < matrix[0].size() - 1; ++j) {
				if (i == j) {
					continue;
				}
				sol += matrix[i][j] * first_sol[j];
			}
			second_sol[i] = sol;
		}
		double last_approx = DBL_MIN;
		for (size_t i = 0; i < matrix.size(); ++i) {
			last_approx = std::max(fabs(second_sol[i] - first_sol[i]), last_approx);
		}
		//make a typical log
		if (last_approx <= approx) {
			//make a final log
			std::cout << "This method took " << iterations << " iterations" << std::endl;
			return second_sol;
		}
		first_sol = second_sol;
		++iterations;
	}
	return second_sol;
}

std::vector<double> Solution::Yakobi(const double& approx) {
	std::vector<std::vector<double>> matrix = {
		{4, 0, 1, 0, 12},
		{0, 3, 0, 2, 19},
		{1, 0, 5, 1, 27},
		{0, 2, 1, 4, 30}
	};
	std::cout << "Your matrix is: \n";
	std::cout << matrix;
	// для програми цього досить, у звіті вказана приблизна кількість ітерацій для розрахунку
	double q = Solution::get_q_from_Yakobi(matrix);
	Solution::change_matrix_to_itereable(matrix);

	std::vector<double> first_sol(4), second_sol(4);// початкове наближення - 0, 0, ... , 0
	for (size_t i = 0; i < matrix.size(); ++i) {
		double divider = matrix[i][i];
		for (size_t j = 0; j < matrix[0].size(); ++j) {
			matrix[i][j] /= divider;
		}
	}
	int iterations = 1;
	while (true) {
		for (size_t i = 0; i < matrix.size(); ++i) {
			double sol = matrix[i][matrix[0].size() - 1];
			for (size_t j = 0; j < matrix[0].size() - 1; ++j) {
				if (i == j) {
					continue;
				}
				sol += matrix[i][j] * first_sol[j];
			}
			second_sol[i] = sol;
		}
		double last_approx = DBL_MIN;
		for (size_t i = 0; i < matrix.size(); ++i) {
			last_approx = std::max(fabs(second_sol[i] - first_sol[i]), last_approx);
		}
		//make a typical log
		if (last_approx <= approx) {
			std::cout << "This method took " << iterations << " iterations" << std::endl;
			//make a final log
			return second_sol;
		}
		first_sol = second_sol;
		++iterations;
	}
	return second_sol;
}

std::vector<double> Solution::find_roots(const std::vector<std::vector<double>>& mat, const std::vector<double>& b, bool transpose) {
	std::vector<double> roots(b.size());
	if (transpose) {
		double y1 = b[0] / mat[0][0], y2 = 0, y3 = 0;
		roots[0] = y1;
		double y1a22 = mat[1][0] * y1;
		y2 = (b[1] - y1a22) / mat[1][1];
		roots[1] = y2;
		double y1y2a33 = mat[2][0] * y1 + mat[2][1] * y2;
		y3 = (b[2] - y1y2a33) / mat[2][2];
		roots[2] = y3;
	}
	else {
		double x3 = b[2] / mat[2][2], x2 = 0, x1 = 0;
		roots[2] = x3;
		double x3a22 = mat[1][2] * x3;
		x2 = (b[1] - x3a22) / mat[1][1];
		roots[1] = x2;
		double x3x2a11 = mat[0][2] * x3 + mat[0][1] * x2;
		x1 = (b[0] - x3x2a11) / mat[0][0];
		roots[0] = x1;
	}
	return roots;
}

std::vector<double> Solution::Square_root() {
	std::vector<std::vector<double>> matrix = {
		{1,2,0},
		{2,2,3},
		{0,3,2}
	};
	std::cout << "Your matrix is: \n";
	std::cout << matrix;
	size_t n = matrix.size(), m = matrix[0].size();

	std::vector<double> b = { 8,22,17 };

	std::vector<std::vector<double>> D(n, std::vector<double>(m))
		, S(n, std::vector<double>(m));

	D[0][0] = Operations::Sgn(matrix[0][0]);
	S[0][0] = sqrt(std::fabs(matrix[0][0]));
	S[0][1] = (matrix[0][1]) / (D[0][0] * S[0][0]);
	S[0][2] = (matrix[0][2]) / (D[0][0] * S[0][0]);
	D[1][1] = Operations::Sgn(matrix[1][1] - std::pow(S[0][1], 2) * D[0][0]);
	S[1][1] = sqrt(std::fabs(matrix[1][1] - std::pow(S[0][1], 2) * D[0][0]));
	S[1][2] = (matrix[1][2] - S[0][1] * D[0][0]*S[0][2]) / (D[1][1] * S[1][1]);
	S[2][2] = sqrt(std::fabs(matrix[2][2] - std::pow(S[0][2], 2) * D[0][0] - std::pow(S[1][2], 2)*D[1][1]));
	D[2][2] = Operations::Sgn(matrix[2][2] - std::pow(S[0][2], 2) * D[0][0] - std::pow(S[1][2], 2) * D[1][1]);
	
	auto trans_S = Operations::Transpose(S);
	trans_S = Operations::Multiply_Matrix(trans_S, D);
	std::vector<double> b_roots = find_roots(trans_S, b, true);
	std::vector<double> roots = find_roots(S, b_roots, false);

	return roots;
}