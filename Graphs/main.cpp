#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <time.h>
#include <random>
#include <chrono>
#include <windows.h>

using namespace std;
vector < vector < int >> graph;

int st = 0, fn = 0;

vector < vector < int >> g;
int n, m;
void output()
{
	ofstream out("input.txt");
	n = 1280000;//rand() % 10000 + 1;
	m = 1280000;//rand() % (n * (n - 1) / 2 + 1);
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
}
void output1()
{
	ofstream out("in.txt");
	n = 1280000;//rand() % 10000 + 1;
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


int input()
{
	ifstream in("input.txt");
	if (!in)
	{
		cerr << "Not opened";
		return 0;
	}
	int n, m;
	in >> n >> m;
	graph.resize(n + 1);
	int a, b;
	for (int i = 0; i < m; i++)
	{
		in >> a >> b;
		//cout << a << " " << b << endl;
		graph[a].push_back(b);
		graph[b].push_back(a);
	}
	return 1;
}

bool BFS(int start, int end)
{
	queue<int> Queue;
	Queue.push(start);
	
	vector<bool> visited;
	visited.resize(graph.size() + 1);

	int c;

	while (!Queue.empty())
	{
		start = Queue.front();
		Queue.pop();
		
		for (int i = 0; i < graph[start].size(); i++)
		{
			c = graph[start][i];
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

int main()
{
#pragma comment(linker, "/STACK:167772166")
	srand(time(NULL));
	srand(time(NULL)); // чтобы формула функции rand() принимала разные значения
	//output();
	input();
	cin >> st >> fn;
	auto start_time = chrono::system_clock::now().time_since_epoch();
	if (BFS(st, fn))
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