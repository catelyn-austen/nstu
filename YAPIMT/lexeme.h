#pragma once
#include <string>

struct lexeme {
    std::string name;
    int type;
    bool value;

    lexeme() {}
    lexeme(std::string _name) {
        name = _name;
        type = 0;
        value = false;
    }
};
