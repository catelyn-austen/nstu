#include "translator.h"

int main() {
    translator t;

    std::cout << (t.tokenize("in/parse.txt") ? "Tokenize: OK" : "Tokenize: ERR") << std::endl;
    std::cout << (t.parse() ? "Parse: OK" : "Parse: ERR") << std::endl;

    return 0;
}
