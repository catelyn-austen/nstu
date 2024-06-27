#ifndef STACK_H
#define STACK_H

#include <stdlib.h>

const int N = 5;

struct stack {
	int top;
	char el[N];
};

// ������������� �����
void create(stack* st);

// �������� �� �������
int is_empty(stack* st);

// ���������� �������� � ����
int push(stack* st, char x);

// ������ �������� �� �����
int pop(stack* st, char* x);

// ����� ����������� �����
int printst(stack* st);

#endif /* STACK_H */