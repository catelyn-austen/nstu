#include <stdio.h>
#include <cmath>

struct list {
	char val = 0;
	list* next = NULL;
};
typedef list dstack;

struct list1 {
	int val = 0;
	list1* next = NULL;
};
typedef list1 stack;

dstack* push(char x, dstack* st) {
	dstack* r; r = new dstack; r->val = x; r->next = st;
	return r;
}
stack* push1(int x, stack* st) {
	stack* r; r = new stack; r->val = x; r->next = st;
	return r;
}

dstack* pop(char* x, dstack* st) {
	dstack* r; r = st; *x = st->val; st = st->next; delete r;
	return st;
}
stack* pop1(int* x, stack* st) {
	stack* r; r = st; *x = st->val; st = st->next; delete r;
	return st;
}

int oper(char a, int b, int c) {
	if (a == '+') return b + c;
	else if (a == '-') return b - c;
	else if (a == '*') return b * c;
	else return round(b / c);
}

int calc(FILE* in) {
	char sv, c;
	int n, k;
	stack* A = NULL; dstack* B = NULL;
	while ((sv = fgetc(in)) != EOF) {
		if (sv == '+' || sv == '-' || sv == '*' || sv == '/') {
			B = push(sv, B);
		}
		else if (sv == ')') {
			A = pop1(&n, A); A = pop1(&k, A); B = pop(&c, B);
			n = oper(c, k, n);
			A = push1(n, A);
		}
		else if (sv != '(') A = push1(sv - '0', A);
	}
	A = pop1(&n, A);
	return n;
	
}

int main() {
	FILE* in;
	fopen_s(&in, "in.txt", "r");
	int t = calc(in);
	printf("%d", t);
}