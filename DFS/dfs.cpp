#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <time.h>
#include <chrono>

using namespace std;

vector<vector<int>> graph;
vector<int> used;

int st = 0, fn = 0;

int n, m;
void output()
{
	ofstream out("input1.txt");
	n = 1000000; //rand() % 10000 + 1;
	m = 10000000; //rand() % (n * (n - 1) / 2 + 1);
	out << n << " " << m << endl;
	graph.resize(n + 1);
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
			for (int i = 0; i < graph[a].size(); i++)
			{
				if (graph[a][i] == b)
					created = false;
			}
		}
		created = false;
		graph[a].push_back(b);
		graph[b].push_back(a);
		out << a << " " << b << endl;
	}
}
void output1()
{
	ofstream out("input1.txt");
	n = 10000000;//rand() % 10000 + 1;
	out << n << " " << n << endl;
	graph.resize(n + 1);
	for (int i = 1; i < n; i++)
	{
		out << i << " " << i + 1 << endl;
		graph[i].push_back(i + 1);
		graph[i + 1].push_back(i);
	}
	graph[1].push_back(n);
	graph[n].push_back(1);
	out << 1 << " " << n << endl;
}

int Input()
{
	ifstream in("input1.txt");
	if (!in)
	{
		cerr << "Not opened";
		return 0;
	}
	int n, m;
	in >> n >> m;
	used.resize(n + 1, 0);
	graph.resize(n + 1);
	int a, b;
	for (int i = 0; i < m; i++)
	{
		in >> a >> b;
		graph[a].push_back(b);
		graph[b].push_back(a);
	}
	return 1;
}
//bool dfs(int start)
//{
//	for (int i = 0; i < graph[start].size(); i++)
//	{
//		if (!used[graph[start][i]])
//		{
//			if (graph[start][i] == fn)
//				return 1;
//			used[graph[start][i]] = 1;
//			if (dfs(graph[start][i]))
//				return 1;
//			
//		}
//	}
//	return 0;
//}
bool dfs(int start, int finish)
{
	stack<int> s;
	s.push(start);
	vector<bool> visited;
	visited.resize(graph.size() + 1);

	int c;

	while (!s.empty())
	{
		start = s.top();
		s.pop();
		for (int i = 0; i < graph[start].size(); i++)
		{
			c = graph[start][i];
			if (used[c] == false)
			{
				s.push(c);
				used[c] = true;
			}
			if (c == finish)
				return true;
		}
	}
	return false;
}
int main()
{
#pragma comment(linker, "/STACK:167772166")
	srand(time(NULL));
	srand(time(NULL)); // чтобы формула функции rand() принимала разные значения
	output1();
	Input();
	cin >> st >> fn;
	auto start_time = chrono::system_clock::now().time_since_epoch();
	if (dfs(st, fn))
	{
		cout << "Exist" << endl;
	}
	else
	{
		cout << "Does not exist" << endl;
	}
	auto end_time = chrono::system_clock::now().time_since_epoch();
	int search_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count();
	cout << search_time << " ms";
	return 0;
}