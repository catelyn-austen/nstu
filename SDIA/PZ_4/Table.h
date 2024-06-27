#ifndef TABLE_H
#define TABLE_H

#include <string>
#include <stdio.h>
#include <stdlib.h>

const int N = 20;
const int b = 5;

struct elem {
	std::string con[b]; //���������
};
struct table {
	elem el[N];
	int n = 0;
};
int findKey(table* A, std::string name); //����� ���� ���� �������� (�������� �����)

int insertT(table* A, elem* el); //�������� � �������

int get_toFile(FILE* in, std::string* out); //�������� ������� �� ����� � str
#endif