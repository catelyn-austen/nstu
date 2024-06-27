#pragma once
#include <vector>
#include <string>
#include "lexeme.h"

class variable_table {
private:
    std::vector<lexeme> table;
public:
    variable_table() {}
    ~variable_table() {
        table.clear();
    }

    bool add(std::string name) {
        if (contains(name))
            return false;
        table.push_back(lexeme(name));
        return true;
    }

    int indexOf(std::string name) {
        for (int i = 0; i < table.size(); i++)
            if (table[i].name == name)
                return i;
        return -1;
    }

    bool contains(std::string name) {
        return indexOf(name) != -1;
    }

    lexeme* get(std::string name) {
        int index = indexOf(name);
        return index != -1 ? &table[index] : NULL;
    }

    lexeme* get(int index) {
        if (index < 0 || index >= table.size())
            return NULL;
        return &table[index];
    }

    bool setType(int index, int type) {
        lexeme* l = get(index);
        if (l == NULL)
            return false;
        l->type = type;
        return true;
    }

    bool setValue(int index, bool val) {
        lexeme* l = get(index);
        if (l == NULL)
            return false;
        l->value = val;
        return true;
    }
};
