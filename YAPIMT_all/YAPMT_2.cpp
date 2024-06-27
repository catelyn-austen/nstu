#include "LexicalScaner.h"
#include "SyntaxScaner.h"

using namespace std;

int main() {

    // Создаем объект класса лексического сканера из ЛР 2
    LexicalScaner translator;

    // Проверяем программу на лексические ошибки
    if (translator.Analize("in/The_program.txt"))
    {
        // Выводим инфо о найденных иденификаторах, константах и циклах
        translator.PrintInfo();

        // Создаем объект класса синтаксического сканера и передаем ему все данные, полученные предыдущем сканером
        SyntaxScaner ssc(&translator);

        // Проверяем программу на синтаксические ошибки
        if (ssc.S_Analize())
        {
            cout << endl << "Success" << endl;
            if (ssc.generate_asm())
                cout << ssc.PrintAssemCode() << endl;
        }
        else
        {
            cout << "Syntax error" << endl;
            // Выводим все синтаксические ошибки
            for (auto er : ssc.syntax_errors)
                cout << er << endl;
        }
            
    }
    else
    {
        cout << "Lexical error" << endl;
        // Выводим все лексические ошибки
        for (auto er : translator.lex_errors)
            cout << er << endl;
    }
    return 0;
}
