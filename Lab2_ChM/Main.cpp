#include <iostream>
#include "Solution.h"

enum class Method{
	Simple_iteration,
	Square_root,
	Yakobi
};

#define PRINT_VECTOR(vec)									  \
	std::cout << "The solution is:\n";						  \
	size_t counter = 1;										  \
	for(const auto& i: vec){								  \
		std::cout << "X" << counter++ << " = " << i << ";\n"; \
	}

int main() {
	try {
		double approx = 0.001;
		while (true) {
			bool change_approx;
			std::cout << "Do you want to change the approximation(1 or 0)?: ";
			std::cin >> change_approx;
			if (change_approx) {
				std::cout << "Enter a new approximation: ";
				std::cin >> approx;
			}
			std::vector<double> sol;
			int m = 0;
			Method method;
			std::cout << "Choose the method(0 - Simple_iteration, 1 - Square_root, 2 - Yakobi) > ";
			std::cin >> m;
			if (m > static_cast<int>(Method::Yakobi) || m < static_cast<int>(Method::Simple_iteration)) {
				std::cin.ignore(10'000, '\n');
				std::cin.clear();
				continue;
			}
			method = static_cast<Method>(m);

			switch (method)
			{
			case Method::Simple_iteration:
				sol = Solution::Simple_iteration(approx);
				break;
			case Method::Square_root:
				sol = Solution::Square_root();
				break;
			case Method::Yakobi:
				sol = Solution::Yakobi(approx);
				break;
			default:
				break;
			}
			PRINT_VECTOR(sol);
		}
	}
	catch (std::domain_error dm) {
		std::cerr << dm.what() << std::endl;
	}
	catch (...) {
		std::cerr << "Smt went wrong" << std::endl;
	}


	return 0;
}