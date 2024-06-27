#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <time.h>
#include <random>
#include <windows.h>

using namespace std;
vector < vector < int >> g;
int n, m, st, en;
vector<bool> marked;

void output()
{
	ofstream out("in.txt");
	int v = 0;
	cout << "Введите кол-во вершин: ";
	cin >> v;
	out << v << " " << v * (v - 1) << endl;
	for (int i = 1; i <= v - 1; i++)
	{
		for (int j = 1; j <= v - 1; j++)
		{
			if (i != j)
				out << i << " " << j << endl;
		}
	}
}

/*void output()
{
	ofstream out("input.txt");
	n = 10000;//rand() % 10000 + 1;
	m = 1000000;//rand() % (n * (n - 1) / 2 + 1);
	out << n << " " << m << endl;
	g.resize(n + 1);
	int a, b;
	bool created = false;
	for (int i = 0; i < m; i++)
	{
		while (!created)
		{
			a = rand() % (n - 1) + 1;
			b = rand() % (n - 1) + 1;
			created = true;
			if (a == b)
				created = false;
			for (int i = 0; i < g[a].size(); i++)
			{
				if (g[a][i] == b)
					created = false;
			}
		}
		created = false;
		g[a].push_back(b);
		g[b].push_back(a);
		out << a << " " << b << endl;
	}
}*/
void output1()
{
	ofstream out("in.txt");
	n = 4000000;//rand() % 10000 + 1;
	out << n << " " << n << endl;
	g.resize(n + 1);
	for (int i = 1; i < n; i++)
	{
		out << i << " " << i + 1 << endl;
		g[i].push_back(i + 1);
		g[i + 1].push_back(i);
	}
	g[1].push_back(n);
	g[n].push_back(1);
	out << 1 << " " << n << endl;
}

void input_test() {
	int v = 0;
	cout << "Введите кол-во вершин: ";
	cin >> v;
	// << v << " " << v * (v - 1) << endl;
	g.resize(v + 1);

	for (int i = 1; i <= v - 1; i++)
	{
		for (int j = 1; j <= v - 1; j++)
		{
			if (i != j)
			{
				g[i].push_back(j);
				g[j].push_back(i);
			}
		}
	}
	cout << "ok" << endl;
}

int input()
{
	ifstream in("in.txt");
	if (!in)
	{
		cerr << "Not opened";
		return 0;
	}
	int n, m;
	in >> n >> m;
	g.resize(n + 1);
	int a, b;
	for (int i = 0; i < m; i++)
	{
		in >> a >> b;
		//cout << a << " " << b << endl;
		g[a].push_back(b);
		g[b].push_back(a);
	}
	return 1;
}

bool bfs(int start, int end)
{
	queue<int> Queue;
	Queue.push(start);

	vector<bool> visited;
	visited.resize(g.size() + 1);

	int c;

	while (!Queue.empty())
	{
		start = Queue.front();
		Queue.pop();

		for (int i = 0; i < g[start].size(); i++)
		{
			c = g[start][i];
			if (visited[c] == false)
			{
				Queue.push(c);
				visited[c] = true;
			}
			if (c == end)
				return true;
		}
	}
	return false;
}
bool dfs(int s)
{
	marked[s] = true;
	for (int i = 0; i < g[s].size(); i++)
	{
		if (g[s][i] == en)
			return 1;
		if (!marked[g[s][i]] && dfs(g[s][i]))
			return 1;
	}
	return 0;
}

int main()
{
#pragma comment(linker, "/STACK:167772166")

	setlocale(LC_ALL, "Russian");
	srand(time(NULL)); // чтобы формула функции rand() принимала разные значения

	//output();
	input_test();

	marked.resize(g.size(), false);

	cin >> st >> en;

	int start_time = clock();
	if (bfs(st, en))
		cout << "Маршрут есть" << endl;
	else
		cout << "Маршрута нет" << endl;
	int end_time = clock();
	cout << end_time - start_time << endl;


	start_time = clock();
	if (dfs(st))
		cout << "Маршрут есть" << endl;
	else
		cout << "Маршрута нет" << endl;
	end_time = clock();
	cout << end_time - start_time << endl;
}