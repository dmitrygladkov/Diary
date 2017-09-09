#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cstring>
#include <chrono>
#include <limits>
#include <iomanip>

const size_t TEST_ITERATIONS = 100000;
const size_t TEST_ARRAY_SIZE = 10000;

static inline bool is_number(const std::string &s)
{
	return (!s.empty() &&
	        std::find_if_not(
			s.cbegin(),
			s.cend(),
			[](char c) {
				return std::isdigit(c);
			}) == s.cend());
}

int main(int argc, char** argv) {
	std::function<int(int, char *[])> parse_num_iter = 
	[](int argc, char *argv[]) -> size_t {
		if ((argc > 1) && (is_number(argv[1]))) {
			size_t iter = std::stoi(argv[1]);
			return ((iter > 0) ? iter :
					     TEST_ITERATIONS);
		}
		return TEST_ITERATIONS;
	};
	std::function<int(int, char *[])> parse_size = 
	[](int argc, char *argv[]) -> size_t {
		if ((argc > 2) && (is_number(argv[2]))) {
			size_t iter = std::stoi(argv[2]);
			return ((iter > 0) ? iter :
					     TEST_ARRAY_SIZE);
		}
		return TEST_ARRAY_SIZE;
	};
	size_t iter = parse_num_iter(argc, argv);
	size_t arr_size = parse_size(argc, argv);
	std::vector<int> memset_v(arr_size, 0);
	std::vector<int> fill_v(arr_size, 0);
	std::vector<int> for_v(arr_size, 0);
	std::vector<int> assign_v(arr_size, 0);

	auto memset_t_start = std::chrono::high_resolution_clock::now();
	for(size_t i = 0; i < iter; ++i) {
		memset(&memset_v[0], 0, memset_v.size() * sizeof memset_v[0]);
	}
	auto memset_t_finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> memset_t_diff =
		(memset_t_finish - memset_t_start);

	auto fill_t_start = std::chrono::high_resolution_clock::now();
	for(size_t i = 0; i < iter; ++i) {
		std::fill(fill_v.begin(), fill_v.end(), 0);
	}
	auto fill_t_finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> fill_t_diff =
		(fill_t_finish - fill_t_start);

	auto for_t_start = std::chrono::high_resolution_clock::now();
	for(size_t i = 0; i < iter; ++i) {
		for (std::vector<int>::iterator it = for_v.begin(),
			end = for_v.end(); it != end; ++it)
			*it = 0;
	}
	auto for_t_finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> for_t_diff =
		(for_t_finish - for_t_start);

	auto assign_t_start = std::chrono::high_resolution_clock::now();
	for(size_t i = 0; i < iter; ++i) {
		assign_v.assign(assign_v.size(),0);
	}
	auto assign_t_finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> assign_t_diff =
		(assign_t_finish - assign_t_start);

	std::cout << "Iterations: " << iter <<
		     " Array's size: " << arr_size << std::endl;
	std::cout << "Memset: "
		<< std::setprecision(std::numeric_limits<double>::max_digits10)
		<< (memset_t_diff.count() / iter) << " ms" << std::endl;
	std::cout << "Fill: "
		<< std::setprecision(std::numeric_limits<double>::max_digits10)
		<< (fill_t_diff.count() / iter) << " ms" << std::endl;
	std::cout << "For: "
		<< std::setprecision(std::numeric_limits<double>::max_digits10)
		<< (for_t_diff.count() / iter) << " ms" << std::endl;
	std::cout << "Assign: "
		<< std::setprecision(std::numeric_limits<double>::max_digits10)
		<< (assign_t_diff.count() / iter) << " ms" << std::endl;

	return 0;
}
