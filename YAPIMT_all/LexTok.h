#pragma once
#include <iostream>
#include <string>

using namespace std;

struct lexeme {
    string name; // ��� �������
    int type; // ��� �������
    double value; // ���� �� ��������
    int nest;

    lexeme() {}
    lexeme(string the_name) {
        name = the_name;
        type = 0;
        value = 0; // �������� �� ���������
    }
};

enum class TypeOfTable {
    Keywords,
    Operators,
    Identifiers,
    Constants,
    Specials
};

ostream& operator <<(ostream& os, const TypeOfTable& type) // ����� ����� ������ � ����� os
{                                                          // ����� ��� ������ ������� ����
    switch (type)
    {
    case TypeOfTable::Keywords:
        os << "Keywords";
        break;
    case TypeOfTable::Operators:
        os << "Operators";
        break;
    case TypeOfTable::Identifiers:
        os << "Identifiers";
        break;
    case TypeOfTable::Constants:
        os << "Constants";
        break;
    case TypeOfTable::Specials:
        os << "Specials";
    default:
        os << "Unknown";
        break;
    }
    return os;
}

class token 
{
public:
    TypeOfTable table; // ��� �������
    int index;         // ������ �������� � �������
    int nest = 0;          // �����������

    // �����������
    token();
    token(TypeOfTable the_table, int the_index) {
        table = the_table;
        index = the_index;
        //nest = 0;
    };

    friend ostream& operator << (ostream& os, const token& out_t);
};

ostream& operator << (ostream& os, token& t) 
{
    os << "(" << t.table << ", " << t.index << ", " << t.nest << ")";
    return os;
}