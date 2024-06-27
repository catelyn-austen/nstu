#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

struct list {
	list* next;
	char elem;
};

int input(FILE* in, list** r)
{
	list* b = *r;
	char a;
	if (fscanf_s(in, "%c", &a, 1) == EOF)
		return 0;
	else {
		b->elem = a;
		while (!feof(in)) {
			fscanf_s(in, "%c", &a, 1);
            b->next = new list;
			b = b->next;
			b->elem = a;
		}
		b->next = NULL;
		return 1;
	}

}

int find(list* i, list* j) {
    if (j != NULL) {
        char r = i->elem;
        char x = j->elem;
        while (j != NULL) {
            if (r == x) {
                return 1;
            }
            else {
                if (j->next != NULL) {
                    x = (j->next)->elem;
                    j = j->next;
                }
                else {
                    i = i->next;
                    j = i->next;
                    find(i, j);
                }
            }
        }
    }
    return 0;
}

int main()
{
	setlocale(LC_ALL, "Russian");
	FILE* in;
	list* r = new list;
	fopen_s(&in, "test.txt", "r");
	int m = input(in, &r);
	if (m == 0)
		printf_s("%s", "Список пуст");
	else
	{
        //list* x;
        //x = r->next;
        //m = find(r, x);
        printf_s("%c", r->elem);
	}
}

/*#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
struct node {
    char value;
    struct node* next;
};
void push(struct node** plist, char new_value) {
    struct node* p = malloc(sizeof(struct node));
    p->value = new_value;
    p->next = *plist;
    *plist = p;
}
void print_list(struct node* plist) {
    for (; plist != NULL; plist = plist->next) {
        printf("%c ", plist->value);
    }
}
bool is_empty(struct node* list) {
    if (list == NULL) {
        return true;
    }
    return false;
}
int pop(struct node** plist) {
    struct node* p = *plist;
    if (is_empty(p) == false) {
        char res = p->value;
        *plist = p->next;
        free(p);
        return res;
    }
    else {
        printf("List is empty");
        exit(0);
    }
}
int find(struct node* i, struct node* j) {
    if (j != NULL) {
        char r = i->value;
        char x = j->value;
        while (j != NULL) {
            if (r == x) {
                return 1;
            }
            else {
                if (j->next != NULL) {
                    x = (j->next)->value;
                    j = j->next;
                }
                else {
                    i = i->next;
                    j = i->next;
                    find(i, j);
                }
            }
        }
    }
    return 0;
}
int main() {
    struct node* list = NULL;
    push(&list, 'z');
    push(&list, 'c');
    push(&list, 'x');
    push(&list, 's');
    push(&list, 'v');
    push(&list, 'y');
    struct node* i = list;
    struct node* j = list->next;
    printf("%d\n", find(i, j));
    print_list(list);
    return 0;
}*/