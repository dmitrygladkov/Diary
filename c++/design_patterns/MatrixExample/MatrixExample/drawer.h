#pragma once 
#include <cstddef>
#include <vector>
#include <memory>
#include <iostream>

using namespace std;

class IDrawer {
public:
	virtual void DrawBorder(size_t width, size_t height) = 0;
	virtual void DrawValues(size_t width, size_t heigth, bool show_border) = 0;
	virtual void Show() = 0;
};

class IDrawable {
public:
	virtual void Draw(IDrawer* drawer) const = 0;
};

struct Item {
	Item() = default;
	Item(double value, bool show_border) {
		this->value = value;
		this->show_border = show_border;
		this->show_value = true;
	}
	double value;
	bool show_value = false;
	bool show_border = false;
};

class ConsoleDrawer : public IDrawer {
public:
	void DrawBorder(size_t width, size_t height) override {
		cout << "ConsoleDrawer::DrawBorder" << endl;
	}
	void DrawValues(size_t width, size_t height, bool show_border) override {
		cout << "ConsoleDrawer::DrawValue" << (show_border ? " with borders" : "") << endl;
	}
	void Show() override {
		cout << "ConsoleDrawer::Show" << endl;
	}
private:
	size_t width = 0;
	size_t height = 0;
	std::vector<std::vector<Item>> drawItems;
};

class HTMLDrawer : public IDrawer {
public:
	void DrawBorder(size_t width, size_t height) override {
		cout << "HTMLDrawer::DrawBorder" << endl;
	}
	void DrawValues(size_t width, size_t height, bool show_border) override {
		cout << "HTMLDrawer::DrawValue" << (show_border ? " with borders" : "") << endl;
	}
	void Show() override {
		cout << "HTMLDrawer::Show" << endl;
	}
private:
	size_t width = 0;
	size_t height = 0;
	std::vector<std::vector<Item>> drawitems;
};

class GraphicsDrawer : public IDrawer {
public:
	void DrawBorder(size_t width, size_t height) override {
		cout << "GraphicsDrawer::DrawBorder" << endl;
	}
	void DrawValues(size_t width, size_t height, bool show_border) override {
		cout << "GraphicsDrawer::DrawValue" << (show_border ? " with borders" : "") << endl;
	}
	void Show() override {
		cout << "GraphicsDrawer::Show" << endl;
	}
private:
	size_t width = 0;
	size_t height = 0;
	std::vector<std::vector<Item>> drawitems;
};
