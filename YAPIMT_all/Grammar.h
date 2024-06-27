#pragma once
#include <string>
#include <vector>
#include <sstream>

using namespace std;

// ����������� ����� ����������
// <> - �������������� �������
// " " - ������������ �������
// [] - ����������� ������ ������������ ��������
enum GrammarType 
{
    Nonterm,
    String,
    Term
};

// ��������� �������������� �������
// ������: <Program> - �������������� �����
struct GrammarUnit 
{
    GrammarType type;
    string value;
};


// ��������� ������� ��������������� �����
// ������: S = <Program> - S - �������������� ����� name, <Program> - ������ � ������ ������ rule
struct GrammarRule 
{
    string name;
    vector<GrammarUnit> rule;
};

class GrammarClass 
{
    public:

    vector<GrammarRule> rules;
    
    // �����������
    GrammarClass() {};

    // ����������
    ~GrammarClass() 
    {
        rules.clear();
    }

    // ������� ��� ����������� �������
    int RuleCount(string _name) 
    {
        int c = 0;
        for (int i = 0; i < rules.size(); i++)
        {
            if (rules[i].name == _name)
                c++;
        }
        return c;
    };

    // ������ ������� �� ������� � �����
    GrammarRule* GetRule(string _name, int index) 
    {
        int i = 0;
        int c = 0;
        for (int k = 0; k < rules.size(); k++)
        {
            if (rules[k].name == _name)
            {
                if (c == index)
                    return &rules[k];
                else
                    c++;
            }
        }
        return NULL;
    };

    // ���� ���������� �� �����
    void Input(string filename)
    {
        ifstream file(filename);
        string line;

        if (!file.is_open())
        {
            cerr << "Grammar file error: " << filename << endl;
            return;
        }

        // ���������� ������ �����
        while (getline(file, line))
        {
            GrammarRule rule;
            istringstream line_arr(line);
            string word;

            while (line_arr >> word)
            {
                // ������� � ��������� ���, ���� ��� ��� ���
                if (rule.name == "") {
                    rule.name = word;
                    line_arr >> word;
                }
                else
                {
                    GrammarUnit gr;

                    // ���� ����������
                    if (word[0] == '<' && word.back() == '>') {
                        gr.type = Nonterm;
                        gr.value = word.substr(1, word.length() - 2);
                    }
                    // ���� ��������
                    else if (word[0] == '"' && word.back() == '"') {
                        gr.type = String;
                        gr.value = word.substr(1, word.length() - 2);
                    }
                    // ���� ������ ������������ ����
                    else if (word[0] == '[' && word.back() == ']') {
                        gr.type = Term;
                        gr.value = word.substr(1, word.length() - 2);
                    }
                    rule.rule.push_back(gr);
                }
            }
            // �� ��������� ����������� ��������
            if (rule.rule.size() > 0 && rule.rule[0].type == Nonterm && rule.rule[0].value == rule.name)
                cerr << "Recurssion in gramatic rules" << endl;
            if (rule.name != "")
                rules.push_back(rule);
        }
        file.close();
    }
};
