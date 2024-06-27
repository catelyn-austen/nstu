#pragma once
#include <string>
#include <vector>
#include <sstream>

using namespace std;

// Обозначение типов грамматики
// <> - Нетерминальные символы
// " " - Терминальные символы
// [] - Обозначение группы терминальных символов
enum GrammarType 
{
    Nonterm,
    String,
    Term
};

// Структура грамматической единицы
// Пример: <Program> - нетерминальное слово
struct GrammarUnit 
{
    GrammarType type;
    string value;
};


// Структура каждого нетерминального слова
// Пример: S = <Program> - S - нетерминальное слово name, <Program> - входит в вектор правил rule
struct GrammarRule 
{
    string name;
    vector<GrammarUnit> rule;
};

class GrammarClass 
{
    public:

    vector<GrammarRule> rules;
    
    // Конструктор
    GrammarClass() {};

    // Деструктор
    ~GrammarClass() 
    {
        rules.clear();
    }

    // Сколько раз встречается правило
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

    // Взятие правила по индексу и имени
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

    // Ввод грамматики из файла
    void Input(string filename)
    {
        ifstream file(filename);
        string line;

        if (!file.is_open())
        {
            cerr << "Grammar file error: " << filename << endl;
            return;
        }

        // построчный разбор файла
        while (getline(file, line))
        {
            GrammarRule rule;
            istringstream line_arr(line);
            string word;

            while (line_arr >> word)
            {
                // Заносим в структуру имя, если его еще нет
                if (rule.name == "") {
                    rule.name = word;
                    line_arr >> word;
                }
                else
                {
                    GrammarUnit gr;

                    // Если нетерминал
                    if (word[0] == '<' && word.back() == '>') {
                        gr.type = Nonterm;
                        gr.value = word.substr(1, word.length() - 2);
                    }
                    // Если терминал
                    else if (word[0] == '"' && word.back() == '"') {
                        gr.type = String;
                        gr.value = word.substr(1, word.length() - 2);
                    }
                    // Если группа терминальных слов
                    else if (word[0] == '[' && word.back() == ']') {
                        gr.type = Term;
                        gr.value = word.substr(1, word.length() - 2);
                    }
                    rule.rule.push_back(gr);
                }
            }
            // Во избежании бесконечной рекурсии
            if (rule.rule.size() > 0 && rule.rule[0].type == Nonterm && rule.rule[0].value == rule.name)
                cerr << "Recurssion in gramatic rules" << endl;
            if (rule.name != "")
                rules.push_back(rule);
        }
        file.close();
    }
};
