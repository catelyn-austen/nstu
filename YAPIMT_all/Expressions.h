#pragma once
#include <string>
#include <vector>
#include "LexTok.h"
#include <variant>

using namespace std;

// создаем общую структуру выражения, где будет хранится информация о токене
struct expr_struct 
{
    vector<token> _expr;
};
// создаем структуру операций для операций присваивания
struct expr_assign
{
    token identifier;
    expr_struct expression;
};
// создаем структуру операций для циклов
struct expr_cycle
{
    expr_struct expression;
};
// пустая структура - важная альтернатива, которая должна присутствовать в классе variant
struct expr_empty {};
// создаем структуру для хранения информации о выражениях различных типов и их вложенности
struct expr_oper
{
    variant<expr_empty, expr_assign, expr_cycle> operation;
    int nestedness;
    
};
