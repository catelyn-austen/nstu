#pragma once
#include <iostream>
#include <fstream>
#include <stack>
#include <map>
#include <variant>
#include "LexicalScaner.h"
#include "LexTok.h"
#include "Variable.h"
#include "Constant.h"
#include "grammar.h"
#include "Expressions.h"

using namespace std;

string assem_str;
vector<expr_oper> v_oper;

class SyntaxScaner
{
    public:
    
    // Создаем указатель на объект класса лексического сканера 
    LexicalScaner* sc;
    // Конструктор класса, в параметрах которого указатель на созданный в main лекс. сканер
    // Он необходим, что передать в синтакс. сканер все то, что уже сущесвтует в лекс. сканере
    SyntaxScaner(LexicalScaner *sc)
    {
        this->sc = sc;
    }

    // потоки вывода
    ofstream out_assem;
    ofstream out_postfix;
    ofstream ferr;
    vector<string> syntax_errors;

    // Основная функция синтаксического анализатора
    bool S_Analize() 
    {
        ferr.open("out/Error_file.txt");
        out_postfix.open("out/postfix.txt");

        vector<string> errors;
        vector<expr_oper> body;
        int er_depth = 0;

        cout << "Postfix form: " << endl;

        // Пошли анализировать с начального Нетерминального символа
        int a = parsing("S", 0, 0, body, errors, er_depth);

        // Если все токены были обработаны
        bool b = true;
        if (a == sc->tokens.size()) {
            b = true;
            v_oper = body;
        }
        else
            b = false;

        // Перенос всех синтаксических ошибок в общий вектор, созданный во 2 ЛР
        if (!b)
        {
            for (int i = 0; i < errors.size(); i++) {
                sc->lex_errors.push_back(errors[i]);
                if (errors[i] == errors[i + 1])
                    break;
            }
        }

        for (int i = 0; i < sc->lex_errors.size(); i++)
        {
            ferr << sc->lex_errors[i] << endl;
            syntax_errors.push_back(sc->lex_errors[i]);
        }

        out_postfix.close();
        ferr.close();

        return b;
    }

    int parsing(string r_name, int p, int depth, vector<expr_oper> &body, vector<string>& errors, int& er_depth)
    {
        // Считыаем количество правил в строке
        int rule_count = sc->grammar.RuleCount(r_name);
        // Проверка на существование правила
        if (rule_count == 0) 
        {
            errors.push_back("Rule error on " + to_string(p) + " position");
            return -1;
        }

        string postfix;
        expr_oper op;
        bool op_found = false;

        // Разбор нетерминала
        for (int i = 0; i < rule_count; i++) 
        {
            int res, p0 = p;
            bool coinc = true;
            GrammarRule* rule = sc->grammar.GetRule(r_name, i);

            // Разбор каждого правила в нетерминале
            for (auto word : rule->rule) 
            {
                switch (word.type) 
                {
                // Если встретился еще один нетерминал, то происходит рекурсия
                case GrammarType::Nonterm:

                    res = parsing(word.value, p, depth + 1, body, errors, er_depth);
                    if (res == -1)
                        coinc = false;
                    else
                        p += res;

                    // Если попалась конструкция присваивания
                    if (coinc && r_name == "Assignment" && word.value == "Expression") 
                    {
                        // для ассемблера
                        op_found = true;
                        op.operation = expr_assign{ sc->tokens[p0], Sort(p0 + 2, p) };
                        op.nestedness = sc->tokens[p0].nest;
                        body.push_back(op);

                        for (int i = 0; i < sc->tokens[p0].index; i++)
                            postfix += "    ";

                        // Добавляем переменную, которой что-то присваивается и знак равенства
                        postfix += TokenToString(sc->tokens[p0]) + " ";
                        //postfix += TokenToString(sc->tokens[p0 + 1]) + " ";
                        for (auto t : Sort(p0 + 2, p))
                            postfix += TokenToString(t) + " ";

                        postfix += TokenToString(sc->tokens[p0 + 1]);
                    }
                    // Если попалась конструкция с циклом, сначала выводим цикл, затем условие, в котором за знаком неравенства следует постфиксная форма
                    else if (coinc && r_name == "Cycle" && word.value == "Expression") 
                    {
                        // для ассемблера
                        op_found = true;
                        op.operation = expr_cycle{ Sort(p - res, p) };
                        op.nestedness = sc->tokens[p0].nest + 1;
                        body.push_back(op);

                        for (int i = 0; i < sc->tokens[p0].nest + 1; i++)
                            postfix += "    ";
                        postfix += "while ";
                        postfix += TokenToString(sc->tokens[p0 + 2]) + " ";
                        postfix += TokenToString(sc->tokens[p0 + 3]) + " ";
                        for (auto t : Sort(p - res + 2, p))
                            postfix += TokenToString(t) + " ";
                    }
                    // Если попалось объявление с заданием начального значения
                    else if (coinc && r_name == "IDeclaration" && word.value == "Expression")
                    {
                        // для ассемблера
                        op_found = true;
                        op.operation = expr_assign{ sc->tokens[p0 + 1], Sort(p0 + 3, p) };
                        op.nestedness = sc->tokens[p0].nest;
                        body.push_back(op);

                        for (int i = 0; i < sc->tokens[p0].index; i++)
                            postfix += "    ";

                        // Добавляем переменную, которой что-то присваивается и знак равенства
                        postfix += TokenToString(sc->tokens[p0 + 1]) + " ";
                        //postfix += TokenToString(sc->tokens[p0 + 2]) + " ";
                        for (auto t : Sort(p0 + 3, p))
                            postfix += TokenToString(t) + " ";
                        postfix += TokenToString(sc->tokens[p0 + 2]);
                    }
                    break;
                // Если встретился терминал, причем в виде строки
                case GrammarType::String:
                    // Если токен есть ожидаемое значение
                    if (p < sc->tokens.size() && TokenToString(sc->tokens[p]) == word.value)
                        p++;
                    else {
                        coinc = false;
                        if (depth > er_depth) 
                        {
                            errors.clear();
                            er_depth = depth;
                        }
                        if (depth >= er_depth)
                            errors.push_back("Synt. error on " + to_string(p) + " token position: " + TokenToString(sc->tokens[p]));
                    }
                    break;
                // Если встретилось общее название группы терминалов
                case GrammarType::Term: // Ссылка на тип (таблицу), считывание токена
                    // Если токен есть ожидаемое значение
                    if (p < sc->tokens.size() && TokenTable(sc->tokens[p]) == word.value)
                        p++;
                    else {
                        coinc = false;
                        if (depth > er_depth) {
                            errors.clear();
                            er_depth = depth;
                        }
                        if (depth >= er_depth)
                            errors.push_back("Synt. error on " + to_string(p) + " token position: " + TokenToString(sc->tokens[p]));
                    }
                    break;
                default:
                    coinc = false;
                    break;
                }

                // Если есть несовпадения по правилам
                if (!coinc)
                    break;
            }

            // Если все хорошо, выводим непустой постфикс в файл
            if (coinc) 
            {
                /*if (op_found)
                {
                    body.push_back(op);
                }*/
                if (postfix != "") {
                    out_postfix << postfix << endl;
                    cout << postfix << endl;
                }
                return p - p0;
            }
            p = p0;
        }
        return -1;
    }

    bool generate_asm()
    {
        out_assem.open("out/assembler.txt");

        // объявляем секцию с данными
        assem_str += "section .data\n";
        // объявляем в ней переменные dd (4-байтовые) и приравниваем их к 0
        for (int i = 0; sc->Videntifiers.GetElemByInd(i) != NULL; i++)
        {
            lexeme* id = sc->Videntifiers.GetElemByInd(i);
            if (id->type == 1)
                assem_str += id->name + " dd 0\n";
        }
        // добавляем секции
        assem_str += "\nsection .text\nglobal _start\n\n_start:\n";

        // стек нужен для циклов
        stack<string> label_stack;
        int c = 0;
        int nest = 0;
        int i = 0;
        // проходимся по каждой занесенной операции (присваивания, циклы)
        for (expr_oper& _oper : v_oper)
        {
            
            // если в новой операции повышается вложенность идентификаторов
            if (_oper.nestedness > nest)
            {             
                assem_str += "\nL_" + to_string(c) + ":\n";
                label_stack.push("L_" + to_string(c));
                c++;
                nest = _oper.nestedness;
            }
            // если вложенность переменных уменьшилась, значит достаем метки из стека
            else if (_oper.nestedness < nest)
            {
                //label_stack.pop();
                nest = _oper.nestedness;
                assem_str += "jmp " + label_stack.top() + "\n\nend_of_cycle_" + to_string(c - 1) + ":\n";
                c--;
                label_stack.pop();
                label_stack.push("L_" + to_string(c - 1));
            }
            // если встречается вариант присваивания
            if (holds_alternative<expr_assign>(_oper.operation))
            {
                expr_assign expr_as = get<expr_assign>(_oper.operation);
                // описываем присваивание на языке ассемблера
                assem_str += generate_expr(expr_as.expression, c) + "mov dword [" + sc->Videntifiers.GetElemByInd(expr_as.identifier.index)->name + "], eax\n";
                // если идет цикл, то в конце записываем прыжок
                
                /*if (nest > 0)
                {
                    assem_str += "jmp " + label_stack.top() + "\n\nend_of_cycle_" + to_string(c - 1) + ":\n";
                    c--;
                    label_stack.pop();
                    label_stack.push("L_" + to_string(c - 1));
                }*/
            }
            // если встретился цикл в операциях
            else if (holds_alternative<expr_cycle>(_oper.operation))
            {
                expr_cycle op = get<expr_cycle>(_oper.operation);
                assem_str += generate_expr(op.expression, c);
            }
            i++;
        }
        // если в стеке что-то осталось
        if (!label_stack.empty())
        {
            assem_str += "jmp " + label_stack.top() + "\n\nend_of_cycle_" + to_string(c - 1) + ":\n";
            c--;
            label_stack.pop();
        }
        // завершение программы
        assem_str += "\n; end of program\nmov eax, 1\nmov ebx, 0\nint 0x80\n";
        
        out_assem << assem_str;
        return true;
    }

    // для подробного описания циклов и присваиваний
    string generate_expr(expr_struct expr, int c)
    {
        string code;
        for (auto token : expr._expr)
        {
            if (token.table == TypeOfTable::Identifiers)
                code += "push dword [" + sc->Videntifiers.GetElemByInd(token.index)->name + "]\n";
            else if (token.table == TypeOfTable::Constants)
                code += "push " + sc->constants.GetElemByInd(token.index)->name + "\n";
            else if (token.table == TypeOfTable::Operators)
            {
                // ecx = id2, eax = id1
                code += "pop ecx\npop eax\n";

                if (TokenToString(token) == "+")
                    code += "add eax, ecx\n";
                else if (TokenToString(token) == "-")
                    code += "sub eax, ecx\n";
                // совершаем противоположное сравнение, чтобы перескочить на метку конца цикла
                else if (TokenToString(token) == ">")
                {
                    code += "cmp eax, ecx\njnae end_of_cycle_" + to_string(c - 1) + "\n";
                }
                else if (TokenToString(token) == "<")
                {
                    code += "cmp eax, ecx\njae end_of_cycle_" + to_string(c - 1) + "\n";
                }
                else if (TokenToString(token) == ">=")
                {
                    code += "cmp eax, ecx\njb end_of_cycle_" + to_string(c - 1) + "\n";
                }
                else if (TokenToString(token) == "<=")
                {
                    code += "cmp eax, ecx\nja end_of_cycle_" + to_string(c - 1) + "\n";
                }
                code += "push eax\n";
            }
        }
        code += "pop eax\n";
        return code;
    }

    string PrintAssemCode()
    {
        return assem_str;
    }

    int Precedence(token t) 
    {
        if (sc->priority.count(t.index))
            return sc->priority[t.index];
        else
            return 0;
    }

    // Сортировка токенов в строке для вывода в постфиксной форме
    vector<token> Sort(int st, int end) 
    {
        stack<token> stack;
        vector<token> queue;

        for (int i = st; i <= end; ++i) 
        {
            token t = sc->tokens[i];

            // Если попались константы или идентификатор, то кладем их в очередь
            if (t.table == TypeOfTable::Constants || t.table == TypeOfTable::Identifiers) 
                queue.push_back(t);
            // Если попался оператор, то кладем в очередь его в зависимости от приоритета
            else if (t.table == TypeOfTable::Operators) 
            {
                while (!stack.empty() && stack.top().table == TypeOfTable::Operators && Precedence(t) <= Precedence(stack.top())) 
                {
                    queue.push_back(stack.top());
                    stack.pop();
                }
                stack.push(t);
            }
        }

        // Пока стек не станет пустым, вытаскиваем из него все
        while (!stack.empty()) 
        {
            queue.push_back(stack.top());
            stack.pop();
        }
        return queue;
    }

    // Вывод названия таблицы
    string TokenTable(token t) {
        switch (t.table) 
        {
        case TypeOfTable::Keywords:
            return "Keywords";

        case TypeOfTable::Operators:
            return "Operators";

        case TypeOfTable::Identifiers:
            return "Identifiers";

        case TypeOfTable::Constants:
            return "Constants";

        case TypeOfTable::Specials:
            return "Specials";

        default:
            return "";
        }
    }

    // Вывод имя токена
    string TokenToString(token t) {
        switch (t.table) {
        case TypeOfTable::Keywords:
            return *sc->keywords.GetElemByInd(t.index);

        case TypeOfTable::Operators:
            return *sc->operators.GetElemByInd(t.index);

        case TypeOfTable::Identifiers:
            return sc->Videntifiers.GetElemByInd(t.index)->name;

        case TypeOfTable::Constants:
            return sc->constants.GetElemByInd(t.index)->name;

        case TypeOfTable::Specials:
            return string(1, sc->specials.GetElemByInd(t.index)->data()[0]);

        default:
            return "";
        }
    }

    
};