#include <iostream>
#include <string>
#include <utility>
#include <cstring>
#include <memory>

class CString {
private:
	const char *s_;
public:
	CString(const char *s) : s_(strdup(s)) {
		std::cout << "Constructor(char*)" << std::endl;
	}
	CString(const std::string &s) : CString(s.c_str()) {
		std::cout << "Constructor(std::string)" << std::endl;
	}
	CString(const CString &that) : s_(strdup(that.s_)) {
		std::cout << "Copy constructor" << std::endl;
	}
	CString(CString &&that) : s_(std::move(that.s_)) {
		that.s_ = nullptr;
		std::cout << "Move constructor" << std::endl;
	}

	CString & operator= (const CString & that) {
		std::cout << "Operator= copy" << std::endl;
		if (s_)
			delete[] s_;
		this->s_ = strdup(that.get());

		return *this;
	}
	CString & operator= (CString && that) {
		std::cout << "Operator= move" << std::endl;
		if (s_)
			delete[] s_;
		s_ = std::move(that.s_);
		that.s_ = nullptr;

		return *this;
	}

	~CString() {
		std::cout << "Destructor" << std::endl;
		delete[] s_;
	}

	const char* get() const {
		return s_;
	};

	bool operator<(const CString &that) const {
		return strcmp(s_,that.s_) < 0;
	}

	bool operator==(const CString &that) const {
		return strcmp(s_,that.s_) == 0;
	}
};

static inline
CString produce(const char *c_str)
{
	CString tmp = CString(c_str);
	return tmp;
}

int main(int argc, char *argv[])
{
	std::string std_str { "Hello World!" };
	CString cstr1 { "Hello Universe!" }; // Constructor(char *)
	CString cstr2 { std_str }; // Constructor(char *); Constructor(std::string)
	CString cstr3 { "WoW!" }; // Constructor(char *)
	CString cstr4("HiHiHi"); // Constructor(char *)
	CString cstr5(std::move(cstr4)); // Move constructor
	CString cstr6 = std::move(produce("Produced String")); //Constructor(char *); Move constructor; Destructor

	std::cout << cstr2.get() << std::endl; // Hello World!
	cstr2 = cstr1; // Operator= copy
	std::cout << cstr2.get() << std::endl; // Hello Universe!
	cstr2 = std::move(cstr3); // Operator= move
	std::cout << cstr2.get() << " "
		<< (cstr3.get() == nullptr ? "nullptr" : "???") << std::endl; // WoW! nullptr
	std::cout << cstr5.get() << " "
		<< (cstr4.get() == nullptr ? "nullptr" : "???") << std::endl; // HiHiHi nullptr
	std::cout << cstr6.get() << std::endl; // Produced String
	
	return 0;
}
