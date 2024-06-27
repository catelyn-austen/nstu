#include <stdio.h>
#include <stdlib.h>

struct queue {
    int el = 0;
    int priority = 0;
    queue* next = NULL;
    queue* prev = NULL;
};

bool empty(queue* q) { return q == NULL; }

void push(queue*& q, int el, int priority)
{
    if (empty(q))
    {
        q->prev = q = new queue;
        q->el = el;
        q->priority = priority;
    }
    else
    {
        if (priority >= q->prev->priority)
        {
            q->prev->next = new queue;
            q->prev->next->el = el;
            q->prev->next->priority = priority;
            q->prev->next->prev = q->prev;
            q->prev = q->prev->next;
        }
        else if (priority < q->priority)
        {
            queue* t = q;
            q = new queue;
            q->el = el;
            q->priority = priority;
            q->next = t;
            q->prev = t->prev;
        }
        else
        {
            queue* m = q;
            for (; m->next->priority <= priority; m = m->next);
            queue* t = m->next;
            m->next = new queue;
            m->next->prev = m;
            m = m->next;
            m->el = el;
            m->priority = priority;
            m->next = t;
            t->prev = m;
        }
    }
}

void print(queue* q)
{
    for (; !empty(q); q = q->next)
        printf_s("%d", q->el);
    printf_s("\n");
}

int main()
{
    queue* a = NULL;
    push(a, 9, 9);
    push(a, 2, 2);
    push(a, 5, 5);
    push(a, 6, 6);
    push(a, 4, 4);
    push(a, 28, 28);
    push(a, 7, 7);
    push(a, 3, 3);
    push(a, 8, 8);
    push(a, 1, 1);
    print(a);
}