#pragma once
#include <iostream>
#include <fstream>
#include <set>

using namespace std;

template <typename T> class TableConst {
private:
    set<T> table; // �������
public:
    // �����������
    TableConst() {} 
    // ����������
    ~TableConst() { table.clear(); } 

    // ���� �� �����
    bool Input(string filename)
    {
        ifstream in(filename);
        if (!in.is_open())
            return false;

        T tmp;
        while (!in.eof()) {
            in >> tmp;
            table.insert(tmp);
        };
        in.close();

        return true;
    }

    // ���������� ��������
    void addElem(T elem) {
        table.insert(elem);
    }

    // �������� �� ������������� ��������
    bool elemExists(T elem) {
        return table.find(elem) != table.end();
    }

    // ������ ������� ��������
    int getIndex(T elem) {
        return elemExists(elem) ? distance(table.begin(), table.find(elem)) : -1;
    }

    // ����� �������� �� �������
    const T* GetElemByInd(int index) {
        if (index < 0 || index >= table.size())
            return NULL;

        auto iter = table.begin();
        advance(iter, index); // ���������� iter � ������� index 
        return &(*iter);
    }
};
