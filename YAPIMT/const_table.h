#pragma once
#include <iostream>
#include <fstream>
#include <set>

template <typename T> class const_table {
private:
    std::set<T> table;
public:
    const_table() {}
    ~const_table() {
        table.clear();
    }

    bool input(std::string filename) {
        std::ifstream in(filename.c_str());
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

    void add(T elem) {
        table.insert(elem);
    }

    bool contains(T elem) {
        return table.find(elem) != table.end();
    }

    int indexOf(T elem) {
        return contains(elem) ? distance(table.begin(), table.find(elem)) : -1;
    }

    const T* get(int index) {
        if (index < 0 && index >= table.size())
            return NULL;

        auto iter = table.begin();
        advance(iter, index);
        return &(*iter);
    }
};
