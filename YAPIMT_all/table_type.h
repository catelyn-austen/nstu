#pragma once

using namespace std;

enum class TypeOfTable {
    Keywords,
    Operators,
    Identifiers,
    Constants
};

ostream& operator <<(ostream& os, const TypeOfTable& type) 
{
    switch (type) 
    {
        case TypeOfTable::Keywords:    
            os << "Keywords"; 
            break;
        case TypeOfTable::Operators:   
            os << "Operators"; 
            break;
        case TypeOfTable::Identifiers: 
            os << "Identifiers"; 
            break;
        case TypeOfTable::Constants:   
            os << "Constants"; 
            break;
        default:                      
            os << "Unknown"; 
            break;
    }
    return os;
}
