#pragma once
#include <iostream>
#include <fstream>
#include <set>

using namespace std;

template <typename T> class TableConst {
private:
    set<T> table; // таблица
public:
    // конструктор
    TableConst() {} 
    // деструктор
    ~TableConst() { table.clear(); } 

    // ввод из файла
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

    // добавление элемента
    void addElem(T elem) {
        table.insert(elem);
    }

    // проверка на существование эелемнта
    bool elemExists(T elem) {
        return table.find(elem) != table.end();
    }

    // взятие индекса элемента
    int getIndex(T elem) {
        return elemExists(elem) ? distance(table.begin(), table.find(elem)) : -1;
    }

    // поиск элемента по индексу
    const T* GetElemByInd(int index) {
        if (index < 0 || index >= table.size())
            return NULL;

        auto iter = table.begin();
        advance(iter, index); // перемещает iter к позиции index 
        return &(*iter);
    }
};
