#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <algorithm>

using namespace std;

void dfsold(int start, vector<int>& comp, vector<bool>& used, vector<vector<int>>& graph)
{
    stack<int> s;
    s.push(start);

    while (!s.empty())
    {
        start = s.top();
        s.pop();
        if (!used[start])
        {
            used[start] = true;
            comp.push_back(start);
            for (int i = 0; i < graph[start].size(); i++)
            {
                int c = graph[start][i];
                if (!used[c])
                {
                    s.push(c);
                }
            }
        }
    }
}

void find_comps(int n, vector<vector<int>>& graph)
{
    vector<int> comp;
    vector<bool> used(graph.size() + 1);
    int k = 1;
    ofstream out("out.txt");

    for (int i = 0; i < n; ++i)
        used[i] = false;
    for (int i = 0; i < n; ++i)
    {
        if (!used[i])
        {
            comp.clear();
            dfsold(i, comp, used, graph);

            out << "Component " << k << ": ";
            for (int j = 0; j < comp.size(); ++j)
                out << " " << comp[j] + 1;
            out << endl;
            k++;
        }
    }
    out.close();
}

int main()
{
    int n, m; // n - вершины, m - рёбра
    ifstream in("in.txt");
    if (!in)
    {
        cerr << "Error";
        return 0;
    }
    in >> n >> m;
    vector<vector<int>> graph(n);
    int a, b;
    for (int i = 0; i < m; i++)
    {
        in >> a >> b;
        graph[a - 1].push_back(b - 1);
        graph[b - 1].push_back(a - 1);
    }
    in.close();
    for (int i = 0; i < n; i++)
    {
        sort(graph[i].begin(), graph[i].end());
    }

        //обход в глубину
    find_comps(n, graph);
}