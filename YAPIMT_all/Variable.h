#pragma once
#include <vector>
#include <string>
#include "LexTok.h"

using namespace std;

class tableVariable {
private:
    vector<lexeme> table;
public:
    // конструктор
    tableVariable() {} 
    // деструктор
    ~tableVariable() { table.clear(); } 

    // добавление элемента
    bool AddElem(string name) // добавить элемент
    {
        if (elemExists(name))
            return false;
        table.push_back(lexeme(name));
        return true;
    }

    // проверка на существование элемента
    bool elemExists(string name) // проверить, существует ли элемент
    {
        return GetIndByName(name) != -1;
    }

    // взятие индекса элемента по имени
    int GetIndByName(string name) // найти индекс элемента по имени
    {
        for (int i = 0; i < table.size(); i++)
            if (table[i].name == name)
                return i;
        return -1;
    }

    // взятие элемента по имени
    lexeme* GetElemByName(string name) // взять элемент по имени
    {
        int index = GetIndByName(name);
        return index != -1 ? &table[index] : NULL;
    }

    // взятие элемента по индексу
    lexeme* GetElemByInd(int index) // взять элемент по индексу
    {
        if (index < 0 || index >= table.size())
            return NULL;
        return &table[index];
    }

    // установление типа
    bool SetType(string name, int type) // установить тип элемента
    {
        for (int i = 0; i < table.size(); i++)
        {
            if (table[i].name == name)
            {
                table[i].type = type;
                return 0;
            }
        }
        return 1;
    }

    // установление значения
    bool SetValue(string the_name, bool val) // установить значение элемента
    {
        for (int i = 0; i < table.size(); i++)
        {
            if (table[i].name == the_name)
            {
                table[i].value = val;
                return 0;
            }
        }
        return 1;
    }

    // установление вложенности
    bool SetNest(string name, int newValue) // установить вложенность
    {
        for (int i = 0; i < table.size(); i++)
        {
            if (table[i].name == name)
            {
                table[i].nest = newValue;
                return 0;
            }
        }
        return 1;
    }

    int Size()
    {
        return table.size();
    }
};
