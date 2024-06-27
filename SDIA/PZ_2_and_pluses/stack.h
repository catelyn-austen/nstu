#ifndef STACK_H
#define STACK_H

#include <stdlib.h>

const int N = 5;

struct stack {
	int top;
	char el[N];
};

// инициализация стека
void create(stack* st);

// проверка на пустоту
int is_empty(stack* st);

// добавление элемента в стек
int push(stack* st, char x);

// взятие элемента из стека
int pop(stack* st, char* x);

// вывод содержимого стека
int printst(stack* st);

#endif /* STACK_H */