#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>

const int iterations = 5;

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

static inline void fibonacci_test(int num_iter)
{
	std::vector<size_t> fib_res;
	// increases capacity of the vector storage
	fib_res.resize(num_iter);

	if (num_iter > 0)
		fib_res[0] = ((num_iter > 1) ? (fib_res[1] = 1) :
					       1);

	for_each(std::next(fib_res.begin(), 2), fib_res.end(),
		 [](size_t &val) {
			val = *(&val - 2) + *(&val - 1);
		 });

	for_each(fib_res.begin(), fib_res.end(),
		 [](size_t val) {
		     std::cout << val << " ";
		 });

	std::cout << std::endl;
}

int main(int argc, char *argv[])
{
	std::function<int(int, char *[])> parse_num_iter = 
	[](int argc, char *argv[]) -> int {
		if ((argc > 1) && (is_number(argv[1]))) {
			int iter = std::stoi(argv[1]);
			return ((iter > 0) ? iter :
					     iterations);
		}
		return iterations;
	};

	fibonacci_test(parse_num_iter(argc, argv));

	return 0;
}
