#pragma once
#include <string>
#include <vector>
#include "LexTok.h"
#include <variant>

using namespace std;

// ������� ����� ��������� ���������, ��� ����� �������� ���������� � ������
struct expr_struct 
{
    vector<token> _expr;
};
// ������� ��������� �������� ��� �������� ������������
struct expr_assign
{
    token identifier;
    expr_struct expression;
};
// ������� ��������� �������� ��� ������
struct expr_cycle
{
    expr_struct expression;
};
// ������ ��������� - ������ ������������, ������� ������ �������������� � ������ variant
struct expr_empty {};
// ������� ��������� ��� �������� ���������� � ���������� ��������� ����� � �� �����������
struct expr_oper
{
    variant<expr_empty, expr_assign, expr_cycle> operation;
    int nestedness;
    
};
