#include <iostream>
#include <stdio.h>

struct list
{
	list* prev;
	list* next;
	float elem;
};

void input(FILE* in, list** r)
{
	*r = new list;
	list* b = *r;

	float a;
	fscanf_s(in, "%f", &a);
	b->elem = a;
	while (!feof(in))
	{
		fscanf_s(in, "%f", &a);
		b->next = new list;
		b->next->prev = b;
		b = b->next;
		b->elem = a;
	}
	b->next = NULL;
	(*r)->prev = b;
}

float summ(list* r, list* d)
{
	float sum = 0;
	while (d != r && d->next != r)
	{
		sum += 2 * d->elem * r->elem;
		d = d->next;
		r = r->prev;
	}
	if (r == d)
		sum += r->elem * d->elem;
	else
		sum += 2 * r->elem * d->elem;
	return sum;
}

int main()
{
	FILE* in, * out;

	char name_read[20] = "read1.txt";
	char name_write[20] = "write.txt";
	char text_error[25] = "WARNING!!!! Dont open '";
	
	errno_t err1 = fopen_s(&in, name_read, "r"),
			err2 = fopen_s(&out, name_write, "w");

	float s;
	list* a;

	if (err1) {
		printf_s("%s%s%s", text_error, name_read, "'");
	}
	else if (err2) {
		printf_s("%s%s%s", text_error, name_write, "'");
	}
	else {
		input(in, &a);
		list* d = a->prev;
		fclose(in);
		s = summ(d, a);

		fprintf_s(out, "%.3f", s);

		fclose(out);
		
		printf_s("%.3f", s);
	}
}