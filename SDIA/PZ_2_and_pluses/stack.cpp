#include <stdio.h>
#include "stack.h"

/*const int N = 5;

struct stack {
	int top;
	char el[N];
};*/

void create(stack *st) {
	st->top = 0;
}

int is_empty(stack* st) {
	if (st->top == 0) return 1;
	return 0;
}

int push(stack* st, char x) {
	if (st->top != N) {
		st->el[st->top] = x;
		st->top++;
		return 1;
	}
	return 0;
}

int pop(stack* st, char* x)
{
	if (!is_empty(st))
	{
		st->top--;
		*x = st->el[st->top];
		return 1;
	}
	return 0;
}


int printst(stack* st) {
	stack a;
	create(&a);
	char x = ' ';
	while (!is_empty(st)) {
		pop(st, &x);
		printf_s("%c ", x);
		push(&a, x);
	}
	while (!is_empty(&a)) {
		pop(&a, &x);
		push(st, x);
	}
	return x;
}

/*int main() {
	stack st;
	stack* S1 = &st;
	create(S1);
	push(S1, '1');
	push(S1, '2');
	push(S1, '3');
	push(S1, '4');
	//char sym;
	//for (int i = 0; i < N; i++) {
		//scanf_s("%c", &sym, 1);
		//push(S1, sym);
	//}
	printst(S1);
	char x;
	pop(S1, &x);
	printf_s("%c\n", x);
	printst(S1);
}*/