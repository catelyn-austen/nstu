#pragma once
#include <iostream>
#include <fstream>
#include "stdlib.h"
#include "Variable.h"
#include "Constant.h"
#include "LexTok.h"
#include <vector>
#include <map>
#include "Grammar.h"

using namespace std;

class LexicalScaner 
{
public:
    tableVariable Videntifiers; // ����������
    tableVariable constants; // ���������
    TableConst<string> keywords; // �������� �����: void/main/return/int/while
    TableConst<string> operators; // ���������
    TableConst<string> specials;
    TableConst<string> chars_operator;
    TableConst<string> chars_special;
    TableConst<char> operators1; // ������ ��������� ������� ����������
    TableConst<char> delimeters; // �����
    TableConst<char> identifiers; // ��������������
    TableConst<char> identifier1; // ������ ��������� ������� ���������������
    TableConst<char> numbers; // �����
    map<int, int> priority;

    GrammarClass grammar;

    vector<token> tokens;
    vector<string> lex_errors;
    
    ofstream ferr; // ����� ������ � ����
    ofstream tokenfile; // ����� ��� ����� �������

    // �������� �� ������ ������
    bool isNumber(string str) 
    {
        if (str.empty()) 
            return false;
        for (size_t i = 0; i < str.length(); ++i) 
        {
            if (!numbers.elemExists(str[i])) 
                return false;
        }
        return true;
    }

    // ������������� �������������� "�����" ��� ��� ���� �������
    void Token_proc(string& str, int& nest) 
    {
        if (str.empty())
            return;
        // ���� �� ����� � �������� ������
        if (keywords.elemExists(str))
        {
            token t(TypeOfTable::Keywords, keywords.getIndex(str));
            if (str == "while")
            {
                t.nest = nest;
                nest++;
            }
            tokenfile << str << "\t" << t << endl;
            tokens.push_back(t);
        }
        // ���� ����� ������� �� ���� (�����)
        else if (isNumber(str)) 
        {
            constants.AddElem(str);
            constants.SetType(str, 1);
            token t(TypeOfTable::Constants, constants.GetIndByName(str));
            t.nest = nest;
            tokenfile << str << "\t" << t << endl;
            tokens.push_back(t);
        }
        // ���� ����� ���������� � ����� ��� _
        else if (identifier1.elemExists(str[0])) 
        {
            Videntifiers.AddElem(str);
            Videntifiers.SetType(str, 1);
            token t(TypeOfTable::Identifiers, Videntifiers.GetIndByName(str));
            t.nest = nest;
            tokenfile << str << "\t" << t << endl;
            tokens.push_back(t);
        }
        // ���� �� ��� �����������������, �� ������
        else 
        {
            ferr << "Error: '" << str << "'is incorrect sequence of characters" << endl; // ����� ������ � ���� ������
            lex_errors.push_back("Error: '" + str + "'is incorrect sequence of characters");
        }
        str.clear();
    }

    void Operator_proc(string str) 
    {
        token t(TypeOfTable::Operators, operators.getIndex(str));
        //cout << str << "\t" << t << endl;
        tokenfile << str << "\t" << t << endl;
        tokens.push_back(t);
    }

    // ���������� ������ ������������ ������� �� ������ � ������ �������
    void Special_Chars(string str) {
        if (specials.getIndex(str) != -1) {
            token t(TypeOfTable::Specials, specials.getIndex(str));
            tokenfile << str << "\t" << t << endl;
            tokens.push_back(t);
        }
        
    }

    bool inputPrecedence(string filename) {
        ifstream file(filename);
        string line;
        int val = 1000;

        if (!file.is_open())
            return false;

        while (getline(file, line)) {
            istringstream iss(line);
            string t;
            while (iss >> t) {
                int index = operators.getIndex(t);
                priority[index] = val;
            }
            val--;
        }

        file.close();
        return true;
    }

//public:
    // �����������
    LexicalScaner() 
    {
        keywords.Input("sym/keywords.txt");
        operators.Input("sym/operators.txt");
        operators1.Input("sym/operators1.txt");
        delimeters.Input("sym/delimeter.txt");
        identifiers.Input("sym/identifiers.txt");
        identifier1.Input("sym/identifiers1.txt");
        chars_operator.Input("sym/chars_operator.txt");
        specials.Input("sym/chars_special.txt");
        numbers.Input("sym/number.txt");
        grammar.Input("sym/Grammar1.txt");
        inputPrecedence("sym/precedence.txt");
        // ��������� ��, ��� �� ����� �� ��������
        delimeters.addElem(' ');
        delimeters.addElem('\t');
        delimeters.addElem('\n');

        ferr = ofstream("out/Error_file.txt");
        tokenfile = ofstream("out/Token_file.txt");
    };
    // ����������
    ~LexicalScaner() 
    {
        ferr.close();
        tokenfile.close();
    };

    // ������� ��������� �����
    bool Analize(string filename) 
    {
        ifstream file(filename);
        string token;
        int nest = 0;
        int s_errors = lex_errors.size();

        char word;
        bool s_comm = false;
        bool m_comm = false;

        while (file.get(word)) 
        {
            if (s_comm) // ���� � ������� ������� �� ������������ �����������
            { 
                if (word == '\n')
                    s_comm = false;
            }
            else if (m_comm) // ���� ���������� ������������� ������. �����������
            { 
                if (word == '*' && file.peek() == '/') 
                {
                    m_comm = false;
                    file.get(word);
                }
            }
            else if (word == '/') // ������ �����������
            { 
                file.get(word);
                if (word == '/') 
                {
                    s_comm = true;
                }
                else if (word == '*') 
                {
                    m_comm = true;
                }
                else 
                {
                    ferr << "Error: " << file.tellg() << " - incomplete comment" << endl;
                    lex_errors.push_back("Error: " + to_string(file.tellg()) + " - incomplete comment");
                }
            }
            else if (operators1.elemExists(word)) // ���� ����������� ������ ������-���� ���������
            { 
                Token_proc(token, nest);

                string op(1, word);
                if (operators1.elemExists(file.peek()))
                    op += file.get();

                if (operators.elemExists(op))
                    Operator_proc(op);
                else {
                    ferr << "Error: " << file.tellg() << " - incorrect operator '" << op << "'" << endl;
                    lex_errors.push_back("Error: " + to_string(file.tellg()) + " - incorrect operator '" + op + "'");
                }
            }
            else if (identifiers.elemExists(word)) // ���� ���������� ������ ��������������
            { 
                token += word;
            }
            else if (delimeters.elemExists(word)) // ���� ���������� �����������
            { 
                if (word == '}')
                    nest--;
                Token_proc(token, nest);
                if (specials.elemExists(string(1, word)))
                    Special_Chars(string(1, word));
            }
            else // ���� ����������� ���-�� ����
            { 
                ferr << "Error: " << file.tellg() << " - incorrect symbol '" << word << "'" << endl;
                lex_errors.push_back("Error: " + to_string(file.tellg()) + " - incorretc symbol '" + word + "'");
            }
        }

        if (m_comm) // ���� ����������� �� ������
        {
            ferr << "Error: incomplete multiline comment" << endl;
            lex_errors.push_back("Error: incomplete multline comment");
        }
        else if (!token.empty())
            Token_proc(token, nest);

        file.close();
        return s_errors == lex_errors.size();
    }

    // ����� ���� ���������� � ��������
    void PrintInfo()
    {
        cout << endl << "Identifiers" << endl;
        for (int i = 0; i < Videntifiers.Size(); i++)
            cout << i << ": " << Videntifiers.GetElemByInd(i)->name << endl;
        cout << endl << "Constants" << endl;
        for (int i = 0; i < constants.Size(); i++)
            cout << i << ": " << constants.GetElemByInd(i)->name << endl;
        cout << endl;
        int c = 0;
        for (int i = 0; i < tokens.size(); i++)
        {
            if (tokens[i].nest >= c)
                c = tokens[i].nest;
        }
        if (c >= 1)
            cout << "Number of nested cycles 'while': " << c + 1 << endl << endl;
    }

};

