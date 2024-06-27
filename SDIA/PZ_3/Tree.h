#ifndef TREE_H
#define TREE_H

#include <stdio.h>

struct tree {
	char val;
	tree* left = NULL;
	tree* right = NULL;
};
struct list {
	tree* data = NULL;
	list* next = NULL;
};
struct queue {
	list* beg = NULL;
	list* end = NULL;
};

struct stack {
	stack* next = NULL;
	tree* el = NULL;
};

void push(queue* A, tree* B); // ����� � ������� ���������

tree* pop(queue* A); // ���� ��������� �� �������

tree* build_tree(FILE* in); // ������ ������

void print_tree1(tree* d, FILE* out);

stack* push_stack(stack* st, tree* d);
stack* pop_stack(stack* st, tree** d);

#endif