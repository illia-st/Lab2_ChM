#pragma once
class Operations {
public:
	template <typename T>
	static T Sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}
	template <typename T>
	static std::vector<std::vector<T>> Transpose(const std::vector<std::vector<T>>& mat) {
		std::vector<std::vector<T>> b(mat.size(), std::vector<T>(mat[0].size()));
		for (size_t i = 0; i < mat.size(); i++) {
			for (size_t j = 0; j < mat[0].size(); j++) {
				b[i][j] = mat[j][i];
			}
		}
		return b;
	}
	template <typename T>	
	static std::vector<std::vector<T>> Multiply_Matrix(const std::vector<std::vector<T>>& lhs,
		const std::vector<std::vector<T>>& rhs) {
		size_t n = lhs.size(), m = lhs[0].size();
		std::vector<std::vector<T>> H(n, std::vector<T>(m));
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				H[i][j] = 0;
				for (int t = 0; t < n; t++) {
					H[i][j] += lhs[i][t] * rhs[t][j];
				}
			}
		}
		return H;
	}
};