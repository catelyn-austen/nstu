#include "Tree.h"
#include <iostream>

void push(queue* A, tree* B)
{
	list* r = new list;
	r->data = B;
	if (A->beg == NULL) {
		A->end = A->beg = r;
	}
	else {
		A->end->next = r;
		A->end = r;
	}
}
stack* push_stack(stack* st, tree* d) {
	stack* r;
	r = new stack; r->el = d; r->next = st;
	return r;
}
stack* pop_stack(stack* st, tree** d) {
	if (st != NULL) {
		stack* r;
		r = st; (*d) = st->el; st = st->next; delete r;
		return st;
	}
}
tree* pop(queue* A)
{
	list* r;
	tree* B = A->beg->data;
	r = A->beg;
	A->beg = A->beg->next;
	delete r;
	return B;
}
tree* build_tree(FILE* in)
{
	char sym;
	tree* d;
	fscanf_s(in, "%c", &sym, 1);
	switch (sym)
	{
	case '(': {d = new tree;
		fscanf_s(in, "%c", &sym, 1); 
		d->val = sym;
		d->left = build_tree(in);
		d->right = build_tree(in); fscanf_s(in, "%c", &sym, 1);
		return  d; }
	case '0':return NULL;
	case ',':d = build_tree(in); break;
	default: return NULL;
	}
}

void print_tree1(tree* d, FILE* out)
{
	stack* S = NULL;
	while (d != NULL) {
		fprintf_s(out, "%c ", d->val);
		if (d->left != NULL && d->right != NULL)
		{
			S = push_stack(S, d->right); d = d->left;
		}
		else if (d->left == NULL && d->right == NULL) {
			if (S != NULL) {
				S = pop_stack(S, &d);
			}
			else {
				d = NULL;
			}
		}
		else if (d->left != NULL) {
			d = d->left;
		}
		else d = d->right;
	}
}

int calc(tree* A, int t) {
	queue* B = new queue; tree* r = NULL;
	if (t == 0) return 1; // если уровень 0 то кол-во эл. 1
	int n = 1; int k = 0; // n - кол-во эл. на данном уровне, k - для подсчёта кол-во эл. на след. уровне
	push(B, A); // положить корень в очередь
	for (int l = 0; l < t; l++) { // пока уровень не t
		for (; n > 0; n--) { // пробегаемся по всем эл. уровня l
			r = pop(B); // берём эл с очереди
			if (r->left != NULL) { push(B, r->left); k++; } // если слева есть эл. то кладём в очередь
			if (r->right != NULL) { push(B, r->right); k++; } // если справа есть эл. то кладём в очередь
		}
		n = k; k = 0;
	}
	return n;
}

int main() {
	setlocale(LC_ALL, "ru");
	FILE* in;
	fopen_s(&in, "in.txt", "r");
	FILE* out;
	printf("Данная программа считывает дерево из файла и выводит его в порядке прямого (левого) обхода. Автор: Хамитова Е.Г.\nПожалуйста, введите число листьев");
	int t; // t - уровень введёный пользователем
	if (in != NULL)
	{
		scanf_s("%d", &t);
		fopen_s(&out, "out.txt", "w");
		tree* F = build_tree(in); //создаём дерево
		print_tree1(F, out);
	}
	else printf("Error while opening file");
}