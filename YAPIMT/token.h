#pragma once
#include <iostream>
#include "table_type.h"

struct token {
    table_type table;
    int index;
    int nestedness;

    token();
    token(table_type _table, int _index) {
        table = _table;
        index = _index;
        nestedness = 0;
    };

    friend std::ostream& operator << (std::ostream& os, const token& out_t);
};

std::ostream& operator << (std::ostream& os, const token& t) {
    os << "token{t=" << t.table << ",i=" << t.index << ",n=" << t.nestedness << "}";
    return os;
}
