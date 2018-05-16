#pragma once

#include <iostream>

using namespace std;

class IMatrixDrawer {
public:
	IMatrixDrawer()
	{
		cout << "constructor IMatrixDrawer" << endl;
	}
	virtual void DrawOutBorder(int weight, int height) = 0;
	virtual void DrawCellBorder(int i, int j) = 0;
	virtual  void DrawCell(int i, int j, int val) = 0;
	virtual ~IMatrixDrawer()
	{
		cout << "destructor IMatrixDrawer" << endl;
	}
};

class ConsoleMatrixDrawer : public IMatrixDrawer {
public:
	ConsoleMatrixDrawer()
	{
		cout << "constructor ConsoleMatrixDrawer" << endl;
	}
	void DrawOutBorder(int weight, int height)
	{
		cout << "ConsoleMatrixDrawer::DrawOutBoarder - " << weight << ":" << height << endl;
	}
	void DrawCellBorder(int i, int j)
	{
		cout << "ConsoleMatrixDrawer::DrawOutBoarder for " << i << ", " << j << endl;
	}
	virtual  void DrawCell(int i, int j, int val)
	{
		cout << "ConsoleMatrixDrawer::DrawCell " << i << ", " << j  <<" - " << val << endl;
	}
	virtual ~ConsoleMatrixDrawer()
	{
		cout << "destructor ConsoleMatrixDrawer" << endl;
	}
};

class GraphMatrixDrawer : public IMatrixDrawer {
public:
	GraphMatrixDrawer()
	{
		cout << "constructor GraphMatrixDrawer" << endl;
	}
	void DrawOutBorder(int weight, int height)
	{
		cout << "GraphMatrixDrawer::DrawOutBoarder - " << weight << ":" << height << endl;
	}
	void DrawCellBorder(int i, int j)
	{
		cout << "GraphMatrixDrawer::DrawOutBoarder for " << i << ", " << j << endl;
	}
	virtual  void DrawCell(int i, int j, int val)
	{
		cout << "GraphMatrixDrawer::DrawCell " << i << ", " << j << " - " << val << endl;
	}
	virtual ~GraphMatrixDrawer()
	{
		cout << "destructor GraphMatrixDrawer" << endl;
	}
};

class HTMLMatrixDrawer : public IMatrixDrawer {
public:
	HTMLMatrixDrawer()
	{
		cout << "constructor HTMLMatrixDrawer" << endl;
	}
	void DrawOutBorder(int weight, int height)
	{
		cout << "HTMLMatrixDrawer::DrawOutBoarder - " << weight << ":" << height << endl;
	}
	void DrawCellBorder(int i, int j)
	{
		cout << "HTMLMatrixDrawer::DrawOutBoarder for " << i << ", " << j << endl;
	}
	virtual  void DrawCell(int i, int j, int val)
	{
		cout << "HTMLMatrixDrawer::DrawCell " << i << ", " << j << " - " << val << endl;
	}
	virtual ~HTMLMatrixDrawer()
	{
		cout << "destructor HTMLMatrixDrawer" << endl;
	}
};
