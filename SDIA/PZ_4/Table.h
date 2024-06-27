#ifndef TABLE_H
#define TABLE_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

const int N = 20;
const int b = 5;

struct elem {
	std::string con[b]; //контейнер
};
struct table {
	elem el[N];
	int n = 0;
};
int findKey(table* A, std::string name); //найти ключ куда вставить (бинарный поиск)

int insertT(table* A, elem* el); //вставить в таблицу

int get_toFile(FILE* in, std::string* out); //записать строчку из файла в str
#endif