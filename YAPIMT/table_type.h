#pragma once

enum class table_type {
    Keywords,
    Operators,
    Identifiers,
    Constants,
    SpecialChars
};

std::ostream& operator<<(std::ostream& os, const table_type& type) {
    switch (type) {
    case table_type::Keywords:       os << "Keywords"; break;
    case table_type::Operators:      os << "Operators"; break;
    case table_type::Identifiers:    os << "Identifiers"; break;
    case table_type::Constants:      os << "Constants"; break;
    case table_type::SpecialChars:   os << "SpecialChars"; break;
    default:                         os << "Unknown"; break;
    }
    return os;
}
