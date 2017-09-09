#include <vector>
#include <string>
#include <iostream>

struct President
{
	std::string name;
	std::string country;
	int year;
 
	President(std::string && p_name, std::string && p_country,
		  int p_year) : name(std::move(p_name)),
				country(std::move(p_country)),
				year(p_year) {
		std::cout << "I am being constructed." << std::endl;
	}

	President(President&& other) : name(std::move(other.name)),
				       country(std::move(other.country)),
				       year(other.year) {
		std::cout << "I am being moved." << std::endl;
	}

	President& operator=(const President& other) = default;
};

// the emplace_back call helps us to avoid the additional copying/moving
// operations thas are done in the push_back call
int main(int argc, char *argv[])
{
	std::vector<President> elections;
	std::cout << "emplace_back:" << std::endl;
	elections.emplace_back("Nelson Mandela", "South Africa", 1994);
 
	std::vector<President> reElections;
	std::cout << std::endl << "push_back:" << std::endl;
	reElections.push_back(President("Franklin Delano Roosevelt",
					"the USA", 1936));
 
	std::cout << std::endl << "Contents:" << std::endl;
	for (President const& president: elections)
		std::cout << president.name << " was elected president of "
			<< president.country << " in " << president.year << std::endl;
	for (President const& president: reElections)
		std::cout << president.name << " was re-elected president of "
			<< president.country << " in " << president.year << std::endl;

	return 0;
}
