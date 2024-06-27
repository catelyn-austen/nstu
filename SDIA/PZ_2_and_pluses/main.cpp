#include <stdio.h>
#include "stack.h"

int main() {
	stack st;
	stack* S1 = &st;
	create(S1);
	push(S1, '1');
	push(S1, '2');
	push(S1, 'a');
	push(S1, '4');
	printst(S1);
	char x;
	pop(S1, &x);
	printf_s("\n%c\n", x);
	printst(S1);
}