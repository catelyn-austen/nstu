#include "Table.h"
#include <iostream>
using namespace std;

int findKey(table* A, std::string name) //бинарный поиск
{
	if (A->n == 0) return 0; //если таблица пуста
	else if (A->el[A->n - 1].con[1] < name) return A->n; //если в конец нужно положить
	int beg = 0; int end = A->n; int i = end / 2;
	while ((i != 0) && (A->el[i - 1].con[1] > name || A->el[i].con[1] < name) && (A->el[i].con[1] != name) && (beg != end)) { //Сложно объяснить :D 
		if (A->el[i].con[1] > name) end = i; // Если нужное место находится слева
		else beg = i; // Если справа
		i = (beg + end) / 2; //указатель на середину промежутка
	}
	return i;
}
int insertT(table* A, elem* el) { //вставить в таблицу
	if (A->n == N) return 0; // Если полная
	if (A->n == 0) { // Если пустая
		A->el[0] = *el;
		A->n++;
	}
	else { // Если не полная не пустая
		int i = findKey(A, el->con[1]); int k; // находим ключ  куда вставлять
		for (int j = A->n - 1; j >= i; j--) { // перемещаем элементы до i
			for (k = 0; k < 5; k++) A->el[j + 1].con[k] = A->el[j].con[k];
		}
		A->el[i] = *el;
		A->n++;
	}
	return 1;
}
int get_toFile(FILE* in, std::string* out) { //Взять строчку из файла и положить в str
	char val;
	fscanf_s(in, "%c", &val, 1);
	while (val != '\n' && !feof(in)) { //Если символ не равен концу строки и файл не закончился
		(*out).push_back(val); //кладём в конец строки символ
		fscanf_s(in, "%c", &val, 1);
	}
	if (feof(in)) return 0; // если мы в конце файла то сообщаем об этом
	return 1;
}

int Scan(FILE* in, table** T) {
	string all_type[14] = { "unsigned char","unsigned int","unsigned long long","unsigned long","unsigned short","long double","char","int","short","long long","long","bool","float","double" }; //массив с типами
	int all_mem[14] = { 1,4,8,4,2,10,1,4,2,8,4,1,4,8 }; //масив с размераи типов
	string str; elem* el = new elem; string n[5];
	int j, k, s, r, len, mem, num; string name, count;
	while (get_toFile(in, &str) != 0) { //берём строчку из файла
		if ((*T)->n == N) return 0;
		if (!str.empty()) { // если она не пуста
			for (int i = 0; i < 14; i++) { //находим определение переменных
				len = all_type[i].length(); //длинна строки типа
				while ((r = str.find(all_type[i])) != -1 && str[r + len] == ' ') { // находим начало
					for (j = r + len + 1; str[j] != '[' && str[j] != ';' && str[j] != ' ' && str[j] != '=' && str[j] != '(' && str[j] != ','; j++); //находим длину имени
					if (str[j] == '(') break; //убираем функции
					name = str.substr(r + len + 1, j - (r + len + 1)); // получаем имя
					if (str[j] != '[') { //если простой тип
						el->con[0] = all_type[i]; el->con[1] = name; el->con[2] = to_string(all_mem[i]); el->con[3] = el->con[4] = "";
						insertT(*T, el); // кладём в таблицу
					}
					else { //если сложный тип
						num = 1;
						for (; str[j] == '['; j++) { //находим число компонент
							for (k = j + 1; str[k] != ']'; k++);
							count = str.substr(j + 1, k - j - 1);
							num *= atoi(count.c_str()); // перевод строки в int
							j = k;
						}
						mem = all_mem[i] * num; // вычисляем память
						el->con[0] = el->con[4] = all_type[i]; el->con[1] = name; el->con[2] = to_string(mem); el->con[3] = to_string(num);
						insertT(*T, el);
					}
					for (; str[j] == ','; j++) { //если после него ","
						if (str[j + 1] == ' ') j = j + 2; // убираем возможный пробел
						for (k = j + 1; str[k] != '[' && str[k] != ';' && str[k] != ','; k++);
						name = str.substr(j, k - j); //тоже самое что и до этого ищем имя
						if (str[k] != '[') { // если простой тип
							el->con[0] = all_type[i]; el->con[1] = name; el->con[2] = to_string(all_mem[i]); el->con[3] = el->con[4] = "";
							insertT(*T, el);
						}
						else { // если сложный тип
							num = 1;
							for (; str[k] == '['; k++) { //находим число компонент
								for (s = k + 1; str[s] != ']'; s++);
								count = str.substr(k + 1, s - k - 1);
								num *= atoi(count.c_str()); // перевод строки в int
								k = s;
							}
							mem = all_mem[i] * num; //подсчёт памяти
							el->con[0] = el->con[4] = all_type[i]; el->con[1] = name; el->con[2] = to_string(mem); el->con[3] = to_string(num);
							insertT(*T, el);
						}
						j = k - 1;
					}
					str.erase(r, j); // удаляем то, что взяли
				}
			}
		}
		str.clear(); //очищаем строку
	}
	return 1;
}

void printTable(table* A) {
	for (int i = 0; i < A->n; i++) {
		for (int j = 0; j < 5; j++) {
			for (char c : A->el[i].con[j]) printf("%c", c);
			printf("\t");
		}
		printf("\n");
	}
}

int main() {
	table* A = new table; //создаём таблицу с ключами
	FILE* in;
	if (fopen_s(&in, "in.txt", "r")) {
		printf("Not opened");
		return 1;
	}
	if (Scan(in, &A) == 0) {
		printf("Переполнение таблицы");
	}
	printTable(A);
	return 1;
}