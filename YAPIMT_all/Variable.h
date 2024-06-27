#pragma once
#include <vector>
#include <string>
#include "LexTok.h"

using namespace std;

class tableVariable {
private:
    vector<lexeme> table;
public:
    // �����������
    tableVariable() {} 
    // ����������
    ~tableVariable() { table.clear(); } 

    // ���������� ��������
    bool AddElem(string name) // �������� �������
    {
        if (elemExists(name))
            return false;
        table.push_back(lexeme(name));
        return true;
    }

    // �������� �� ������������� ��������
    bool elemExists(string name) // ���������, ���������� �� �������
    {
        return GetIndByName(name) != -1;
    }

    // ������ ������� �������� �� �����
    int GetIndByName(string name) // ����� ������ �������� �� �����
    {
        for (int i = 0; i < table.size(); i++)
            if (table[i].name == name)
                return i;
        return -1;
    }

    // ������ �������� �� �����
    lexeme* GetElemByName(string name) // ����� ������� �� �����
    {
        int index = GetIndByName(name);
        return index != -1 ? &table[index] : NULL;
    }

    // ������ �������� �� �������
    lexeme* GetElemByInd(int index) // ����� ������� �� �������
    {
        if (index < 0 || index >= table.size())
            return NULL;
        return &table[index];
    }

    // ������������ ����
    bool SetType(string name, int type) // ���������� ��� ��������
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

    // ������������ ��������
    bool SetValue(string the_name, bool val) // ���������� �������� ��������
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

    // ������������ �����������
    bool SetNest(string name, int newValue) // ���������� �����������
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
