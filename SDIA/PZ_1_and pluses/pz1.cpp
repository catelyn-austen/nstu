#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct Node {
    void* value;
    struct Node* next;
    struct Node* prev;
} Node;

typedef struct DblLinkedList {
    size_t size;
    Node* head;
    Node* tail;
} DblLinkedList;

// создание двусвязного списка
DblLinkedList* createDblLinkedList() {
    DblLinkedList* tmp = new DblLinkedList;
    tmp->size = 0;
    tmp->head = tmp->tail = NULL;

    return tmp;
}

// добавление элементов в конец списка
void pushBack(DblLinkedList* list, void* value) { 
    Node* tmp = new Node;
    if (tmp == NULL) {
        exit(3);
    }
    tmp->value = value;
    tmp->next = NULL;
    tmp->prev = list->tail;
    if (list->tail) {
        list->tail->next = tmp;
    }
    list->tail = tmp;

    if (list->head == NULL) {
        list->head = tmp;
    }
    list->size++;
}

// подпрограмма для считывания данных из файла
void input(const char* filename, DblLinkedList** list)
{
    FILE* f;
    fopen_s(&f,filename, "r+");
    if (f != NULL)
    {
        do
        {
            float* a = new float;
            fscanf_s(f, "%f", a);
            pushBack(*list, a);
        } while (!feof(f));

    }
    else
    {
        printf("Error");
    }
}
// подпрограмма суммирования
float summ(DblLinkedList* list)
{
    float sum = 0;

    Node* tmp_h = list->head;
    Node* tmp_t = list->tail;
    for (int i = 0; i < int(list->size / 2); i++)
    {
        float a = 0;
        float b = 0;
        // void* --> float
        a = *((float*)tmp_h->value); 
        b = *((float*)tmp_t->value);
        sum += a * b;
        tmp_h = tmp_h->next;
        tmp_t = tmp_t->prev;

    }
    if (tmp_h == tmp_t)
    {
        float a = 0;
       a = *((float*)tmp_h->value);
        sum = 2 * sum;
        sum += a * a;
        return sum;
    }

    return 2 * sum;

}

int main() {
    DblLinkedList* list = createDblLinkedList();
    DblLinkedList* list = NULL;

    input("test.txt", &list);
 
    float s = summ(list);
    printf("%1.2f", s);

    return 0;
}