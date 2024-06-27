#pragma once
#include <iostream>
#include <fstream>
#include <stack>
#include <map>
#include "variable_table.h"
#include "const_table.h"
#include "grammar.h"

class translator {
private:
    const_table<std::string> t_keywords; // std::set<std::string> keywords = { "void", "main", "return", "int", "do", "while" };
    const_table<std::string> t_operators; // std::set<std::string> operators = { "+=", "-=", "==", "!=", "+", "-", "=" };
    const_table<char> t_operatorChars; // std::set<char> operatorChars = { '+', '-', '=', '!' };
    const_table<char> t_delimeterChars; // std::set<char> delimeterChars = { ' ', '\t', '\n', '\r', ';', '(', ')', '{', '}', ',' };
    const_table<char> t_specialChars; // std::set<char> specialChars = { ';', '(', ')', '{', '}', ',' };
    const_table<char> t_identifierChars; /* std::set<char> identifierChars = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '_' }; */
    const_table<char> t_identifierStartChars; /* std::set<char> identifierStartChars = { 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        '_' }; */
    const_table<char> t_numberChars; // std::set<char> numberChars = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' };
    grammar_t grammar;

    std::map<int, int>  t_precedence;
    variable_table t_identifiers;
    variable_table t_consts;
    std::vector<token> v_tokens;
    std::vector<std::string> v_errors;

    std::ofstream f_postfix;
    std::ofstream f_err;

    // Чтение приоритетов операций из файла
    bool inputPrecedence(std::string filename) {
        std::ifstream file(filename);
        std::string line;
        int val = 1000;

        if (!file.is_open())
            return false;

        while (getline(file, line)) {
            std::istringstream iss(line);
            std::string t;
            while (iss >> t) {
                int index = t_operators.indexOf(t);
                t_precedence[index] = val;
            }
            val--;
        }

        file.close();
        return true;
    }

    // Проверка, является ли представленный строкой токен числом
    bool isNumber(std::string token) {
        if (token.empty()) return false;
        for (size_t i = 0; i < token.length(); ++i)
            if (!t_numberChars.contains(token[i])) return false;
        return true;
    }

    // Получение типа токена из строки и добавление его в вектор токенов
    void processToken(std::vector<token>& tokens, int& nestedness, std::string& str) {
        if (str.empty())
            return;

        if (t_keywords.contains(str)) {
            token t(table_type::Keywords, t_keywords.indexOf(str));
            if (str == "while")
                nestedness--;
            t.nestedness = nestedness;
            if (str == "do")
                nestedness++;
            tokens.push_back(t);
        }
        else if (isNumber(str)) {
            t_consts.add(str);
            t_consts.setType(t_consts.indexOf(str), 1);
            token t(table_type::Constants, t_consts.indexOf(str));
            t.nestedness = nestedness;
            tokens.push_back(t);
        }
        else if (t_identifierStartChars.contains(str[0])) {
            t_identifiers.add(str);
            token t(table_type::Identifiers, t_identifiers.indexOf(str));
            t.nestedness = nestedness;
            tokens.push_back(t);
        }
        else {
            v_errors.push_back("tokenize: Error at '" + str + "': disallowed sequence");
        }
        str.clear();
    }

    // Добавление токена оператора из строки в вектор токенов
    void processOperator(std::vector<token>& tokens, int& nestedness, std::string str) {
        token t(table_type::Operators, t_operators.indexOf(str));
        t.nestedness = nestedness;
        tokens.push_back(t);
    }

    // Добавление токена специального символа из строки в вектор токенов
    void processSpecialChar(std::vector<token>& tokens, int& nestedness, char ch) {
        token t(table_type::SpecialChars, t_specialChars.indexOf(ch));
        t.nestedness = nestedness;
        tokens.push_back(t);
    }

    // Получение строкового представления токена
    std::string getTokenString(token t) {
        switch (t.table) {
        case table_type::Keywords:
            return *t_keywords.get(t.index);
        case table_type::Operators:
            return *t_operators.get(t.index);
        case table_type::Identifiers:
            return t_identifiers.get(t.index)->name;
        case table_type::Constants:
            return t_consts.get(t.index)->name;
        case table_type::SpecialChars:
            return std::string(1, *t_specialChars.get(t.index));
        default:
            return "";
        }
    }

    // Получение типа токена в строковом виде
    std::string getTokenTable(token t) {
        switch (t.table) {
        case table_type::Keywords:
            return "Keywords";
        case table_type::Operators:
            return "Operators";
        case table_type::Identifiers:
            return "Identifiers";
        case table_type::Constants:
            return "Constants";
        case table_type::SpecialChars:
            return "SpecialChars";
        default:
            return "";
        }
    }

    // Получение приоритета токена оператора из таблицы приоритета
    int getTokenPrecedence(token t) {
        if (t_precedence.count(t.index))
            return t_precedence[t.index];
        else
            return 0;
    }

    // Разбор выражения в области вектора токенов, алгоритм сортировочной станции
    std::vector<token> shuntingYard(int start, int end) {
        std::stack<token> opstack;
        std::vector<token> queue;

        for (int i = start; i <= end; ++i) {
            token t = v_tokens[i];

            if (t.table == table_type::Constants || t.table == table_type::Identifiers) {
                queue.push_back(t);
            }
            else if (t.table == table_type::Operators) {
                while (!opstack.empty() && opstack.top().table == table_type::Operators && getTokenPrecedence(t) <= getTokenPrecedence(opstack.top())) {
                    queue.push_back(opstack.top());
                    opstack.pop();
                }
                opstack.push(t);
            }
            else if (getTokenString(t) == "(") {
                opstack.push(t);
            }
            else if (getTokenString(t) == ")") {
                while (!opstack.empty() && getTokenString(opstack.top()) != "(") {
                    queue.push_back(opstack.top());
                    opstack.pop();
                }
                if (!opstack.empty()) opstack.pop();
            }
        }

        while (!opstack.empty()) {
            queue.push_back(opstack.top());
            opstack.pop();
        }

        return queue;
    }

    // Рекурсивный разбор правила по имени
    int parseRule(std::string rule_name, int pos, int depth, std::vector<std::string>& errors, int& error_depth) {
        if (grammar.getRuleCount(rule_name) == 0) {
            errors.push_back("parse: Unknown rule '" + rule_name + " at token pos " + std::to_string(pos));
            return -1;
        }

        std::string postfix;
        for (int i = 0; i < grammar.getRuleCount(rule_name); i++) { // Для каждого правила нетерминального символа
            grammar_rule* rule = grammar.getRule(rule_name, i); // Правило

            int start = pos;
            int result;
            bool matched = true;

            for (auto ref : rule->rule) { // Для каждого элемента текущего правила
                switch (ref.type) {
                case grammar_reference_type::Rule: // Ссылка на нетерминал, рекурсивный переход
                    result = parseRule(ref.value, pos, depth + 1, errors, error_depth);
                    if (result == -1)
                        matched = false;
                    else
                        pos += result;

                    if (matched && rule_name == "ASSIGNMENT" && ref.value == "EXPRESSION") { // Успешное считывание присваивания, создание постфиксной записи значения
                        for (int i = 0; i < v_tokens[start].nestedness; i++)
                            postfix += "    ";
                        postfix += getTokenString(v_tokens[start]) + " ";
                        for (auto t : shuntingYard(start + 2, pos))
                            postfix += getTokenString(t) + " ";
                        postfix += getTokenString(v_tokens[start + 1]);
                    }
                    else if (matched && rule_name == "DOWHILE" && ref.value == "EXPRESSION") { // Успешное считывание цикла do-while, создание постфиксной записи условия
                        for (int i = 0; i < v_tokens[start].nestedness + 1; i++)
                            postfix += "    ";
                        for (auto t : shuntingYard(pos - result, pos))
                            postfix += getTokenString(t) + " ";
                        postfix += "dowhile";
                    }

                    break;
                case grammar_reference_type::String: // Ссылка на текстовое представление, считывание токена
                    if (pos < v_tokens.size() && getTokenString(v_tokens[pos]) == ref.value)
                        pos++;
                    else {
                        matched = false;
                        if (depth > error_depth) {
                            errors.clear();
                            error_depth = depth;
                        }
                        if (depth >= error_depth)
                            errors.push_back("parse: Expected '" + ref.value + "', got '" + getTokenString(v_tokens[pos]) + "' at token pos " + std::to_string(pos) + " (rule '" + rule_name + "')");
                    }
                    break;
                case grammar_reference_type::Table: // Ссылка на тип (таблицу), считывание токена
                    if (pos < v_tokens.size() && getTokenTable(v_tokens[pos]) == ref.value)
                        pos++;
                    else {
                        matched = false;
                        if (depth > error_depth) {
                            errors.clear();
                            error_depth = depth;
                        }
                        if (depth >= error_depth)
                            errors.push_back("parse: Expected [" + ref.value + "], got '" + getTokenString(v_tokens[pos]) + "' at token pos " + std::to_string(pos) + " (rule '" + rule_name + "')");
                    }
                    break;
                default:
                    matched = false;
                    break;
                }

                if (!matched)
                    break;
            }

            if (matched) { // Вывод постфиксной записи, возврат количества поглощенных токенов
                if (postfix != "")
                    f_postfix << postfix << std::endl;
                return pos - start;
            }
            pos = start;
        }

        return -1; // Ни одно правило не подошло, возврат ошибки
    }
public:
    translator() {
        t_keywords.input("const/keywords.txt");
        t_operators.input("const/operators.txt");

        t_operatorChars.input("const/chars_operator.txt");
        t_delimeterChars.input("const/chars_delimeter.txt");
        t_specialChars.input("const/chars_special.txt");
        t_identifierChars.input("const/chars_identifier.txt");
        t_identifierStartChars.input("const/chars_identifierStart.txt");
        t_numberChars.input("const/chars_number.txt");

        grammar.input("const/grammar.txt");
        inputPrecedence("const/precedence.txt");

        t_delimeterChars.add(' ');
        t_delimeterChars.add('\t');
        t_delimeterChars.add('\n');

        f_postfix.open("out/postfix.txt");
        f_err.open("out/err.txt");
    };
    ~translator() {};

    // Токенизация файла
    bool tokenize(std::string filename) {
        std::ifstream file(filename);
        std::string t;
        int errors = v_errors.size();

        char ch;
        bool inSinglelineComment = false;
        bool inMultilineComment = false;
        int nestedness = 0;

        while (file.get(ch)) {
            if (inSinglelineComment) { // Текст однострочного комментария
                if (ch == '\n')
                    inSinglelineComment = false;
            }
            else if (inMultilineComment) { // Текст многострочного комментария
                if (ch == '*' && file.peek() == '/') {
                    inMultilineComment = false;
                    file.get(ch);
                }
            }
            else if (ch == '/') { // Начало комментария
                file.get(ch);
                if (ch == '/')
                    inSinglelineComment = true;
                else if (ch == '*')
                    inMultilineComment = true;
                else
                    v_errors.push_back("tokenize: Error at pos " + std::to_string(file.tellg()) + ": incomplete comment definition");
            }
            else if (t_operatorChars.contains(ch)) { // Оператор
                processToken(v_tokens, nestedness, t);

                std::string op(1, ch);
                if (t_operatorChars.contains(file.peek()))
                    op += file.get();

                if (t_operators.contains(op))
                    processOperator(v_tokens, nestedness, op);
                else
                    v_errors.push_back("tokenize: Error at pos " + std::to_string(file.tellg()) + ": unknown operator '" + op + "'");
            }
            else if (t_identifierChars.contains(ch)) { // Символ идентификатора
                t += ch;
            }
            else if (t_delimeterChars.contains(ch)) { // Разделитель
                processToken(v_tokens, nestedness, t);
                if (t_specialChars.contains(ch))
                    processSpecialChar(v_tokens, nestedness, ch);
            }
            else { // Запрещенный символ
                v_errors.push_back("tokenize: Error at pos " + std::to_string(file.tellg()) + ": disallowed character '" + ch + "'");
            }
        }

        if (inMultilineComment)
            v_errors.push_back("tokenize: Error at EOF: unterminated multiline comment");
        else if (!t.empty())
            processToken(v_tokens, nestedness, t);

        file.close();
        return errors == v_errors.size();
    }

    // Парсинг вектора токенов, полученного на этапе токенизации
    bool parse() {
        int error_depth = 0;
        std::vector<std::string> errors;
        int consumed = parseRule("PROGRAM", 0, 0, errors, error_depth); // Начать разбор с правила PROGRAM
        bool ok = consumed == v_tokens.size(); // Если все токены были поглощены правилом, разбор успешен
        if (!ok)
            for (auto error : errors)
                v_errors.push_back(error);

        for (auto error : v_errors)
            f_err << error << std::endl;

        return ok;
    }
};
