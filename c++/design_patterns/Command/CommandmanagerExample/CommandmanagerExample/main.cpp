#include <vector>
#include <stack>
#include <memory>
#include <list>
#include <iostream>
#include <conio.h>

using namespace std;

enum CommandType {
	PUSH_BACK = 1,
	POP_BACK = 2,
	UNDO = 3,
};

void printCurrState(vector<int> *arr) {
	if (!arr->size()) {
		cout << "EMPTY" << endl;
		return;
	}

	for (size_t i = 0; i < arr->size() - 1; i++) {
		cout << (*arr)[i] << ", ";
	}
	cout << (*arr)[arr->size() - 1] << endl;
}

class ICommand {
	virtual void exec() = 0;
	virtual void undo() = 0;
public:
	virtual ~ICommand()
	{

	}
};

class PushBackCommand : public ICommand {
	vector<int> *arr;
	vector<int> prev_arr;
	int val;
	void exec() {
		arr->push_back(val);
	}
	void undo() {
		*arr = prev_arr;
	}
public:
	PushBackCommand(vector<int> *_arr, int _val) : arr(_arr), val(_val) {
		prev_arr = *arr;
		//cout << "Debug prev state: ";
		//printCurrState(&prev_arr);
		this->exec();
	}
	~PushBackCommand() {
		undo();
	}
};

class PopBackCommand : public ICommand {
	vector<int> *arr;
	vector<int> prev_arr;
	void exec() {
		return arr->pop_back();
	}
	void undo() {
		*arr = prev_arr;
	}
public:
	PopBackCommand(vector<int> *_arr) : arr(_arr) {
		prev_arr = *arr;
		//cout << "Debug prev state: ";
		//printCurrState(&prev_arr);
		this->exec();
	}
	~PopBackCommand() {
		undo();
	}
};

class EmptyCommand : public ICommand {
	vector<int> *arr;
	vector<int> prev_arr;
	void exec() {
	}
	void undo() {
		*arr = prev_arr;
	}
public:
	EmptyCommand(vector<int> *_arr) : arr(_arr) {
		prev_arr = *arr;
		//cout << "Debug prev state: ";
		//printCurrState(&prev_arr);
	}
	~EmptyCommand() {
		undo();
	}
};

class ArrayManager {
	vector<int> arr;
	stack<ICommand *> history;
	EmptyCommand *empty_com;
public:
	ArrayManager() {
		empty_com = new EmptyCommand(&arr);
		history.push(empty_com);
	}
	void executeCommand(CommandType com, int val) {
		ICommand *com_obj;

		switch (com) {
		case PUSH_BACK:
			com_obj = new PushBackCommand(&arr, val);
			history.push(com_obj);
			break;
		case POP_BACK:
			com_obj = new PopBackCommand(&arr);
			history.push(com_obj);
			break;
		case UNDO:
			com_obj = history.top();
			history.pop();
			if (com_obj == empty_com) {
				delete com_obj;
				empty_com = new EmptyCommand(&arr);
				history.push(empty_com);
			}
			else {
				delete com_obj;
			}
			break;
		}
		printCurrState(&arr);
	}
	void executeCommand(CommandType com) {
		ICommand *com_obj;

		switch (com) {
		case PUSH_BACK:
			com_obj = new PushBackCommand(&arr, 0);
			history.push(com_obj);
			break;
		case POP_BACK:
			com_obj = new PopBackCommand(&arr);
			history.push(com_obj);
			break;
		case UNDO:
			com_obj = history.top();
			history.pop();
			if (com_obj == empty_com) {
				delete com_obj;
				empty_com = new EmptyCommand(&arr);
				history.push(empty_com);
			}
			else {
				delete com_obj;
			}
			break;
		}
		printCurrState(&arr);
	}
};

int main()
{
	ArrayManager manager;

	cout << "Do several push_backs: " << endl;
	manager.executeCommand(PUSH_BACK, 5);
	manager.executeCommand(PUSH_BACK, 6);
	manager.executeCommand(PUSH_BACK, 7);
	manager.executeCommand(PUSH_BACK, 8);
	cout << endl << endl << endl;

	cout << "Do several pop_backs: " << endl;
	manager.executeCommand(POP_BACK);
	manager.executeCommand(POP_BACK);
	manager.executeCommand(POP_BACK);
	cout << endl << endl << endl;

	cout << "Test undo: " << endl;
	manager.executeCommand(UNDO);
	manager.executeCommand(UNDO);
	manager.executeCommand(UNDO);
	manager.executeCommand(UNDO);
	manager.executeCommand(UNDO);
	manager.executeCommand(UNDO);
	manager.executeCommand(UNDO);
	cout << endl << endl << endl;

	cout << "Test errorneous undo: " << endl;
	manager.executeCommand(UNDO);
	manager.executeCommand(UNDO);
	manager.executeCommand(UNDO);
	cout << endl << endl << endl;

	_getch();

	return 0;
}