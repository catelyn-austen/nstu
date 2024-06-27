#pragma once
#include <string>
#include <vector>
#include <sstream>

enum grammar_reference_type {
    Rule,   // <...>
    String, // "..."
    Table   // [...]
};

struct grammar_reference {
    grammar_reference_type type;
    std::string value;
};

struct grammar_rule {
    std::string name;
    std::vector<grammar_reference> rule;
};

class grammar_t {
private:
    std::vector<grammar_rule> rules;
public:
    grammar_t() {};
    ~grammar_t() {
        rules.clear();
    }

    void input(std::string filename) {
        std::ifstream file(filename);
        std::string line;

        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename << std::endl;
            return;
        }

        while (std::getline(file, line)) {
            if (line[0] == '#')
                continue;

            std::istringstream ss(line);
            grammar_rule rule;
            std::string segment;

            while (ss >> segment) {
                if (rule.name == "") {
                    rule.name = segment;
                    ss >> segment;
                }
                else {
                    grammar_reference ref;
                    if (segment[0] == '<' && segment.back() == '>') {
                        ref.type = Rule;
                        ref.value = segment.substr(1, segment.length() - 2);
                    }
                    else if (segment[0] == '"' && segment.back() == '"') {
                        ref.type = String;
                        ref.value = segment.substr(1, segment.length() - 2);
                    }
                    else if (segment[0] == '[' && segment.back() == ']') {
                        ref.type = Table;
                        ref.value = segment.substr(1, segment.length() - 2);
                    }
                    rule.rule.push_back(ref);
                }
            }

            if (rule.rule.size() > 0 && rule.rule[0].type == Rule && rule.rule[0].value == rule.name)
                std::cerr << "Warning: direct left recursion detected!" << std::endl;
            if (rule.name != "")
                rules.push_back(rule);
        }

        file.close();
    }

    int getRuleCount(std::string name) {
        int i = 0;
        for (auto it = rules.begin(); it != rules.end(); it++)
            if (it->name == name)
                i++;
        return i;
    };
    grammar_rule* getRule(std::string name, int index) {
        int i = 0;
        for (auto it = rules.begin(); it != rules.end(); it++)
            if (it->name == name)
                if (i == index)
                    return &(*it);
                else
                    i++;
        return NULL;
    };
};
