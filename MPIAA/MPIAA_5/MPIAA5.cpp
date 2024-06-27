#include <iostream>
#include <fstream>
#include <vector>
#include <queue>



using namespace std;

class Ver
{
public:
    int v;
    int wt;
    Ver()
    {
        v = 0;
        wt = 0;
    }
    Ver(int _v, int _wt)
    {
        v = _v;
        wt = _wt;
    }
};
class Edge
{
public:
    int a;
    int b;
    int wt;
    Edge()
    {
        a = 0;
        b = 0;
        wt = 0;
    }
    Edge(int _a, int _b, int _wt)
    {
        a = _a;
        b = _b;
        wt = _wt;
    }
};

class Graph
{
public:
    int V;
    int E;
    vector<vector<Ver>> g;
    Graph()
    {
        V = 0;
        E = 0;
    }
    void Init(string filename)
    {
        int v, w, wt;;
        ifstream in;
        in.open("in.txt");
        if (!in.is_open())
        {
            cout << "Can't open file \"in.txt\"" << endl;
            return;
        }
        in >> V >> E;
        g.resize(V);

        for (int i = 0; i < E; i++)
        {
            in >> v >> w >> wt;
            g[v-1].push_back(Ver(w-1, wt));
            //g[w-1].push_back(Ver(v-1, wt));
        }
        return;
    }

};

int operator <(Edge e1, Edge e2)
{
    return e1.wt > e2.wt;
}

void findPrimMST(Graph g, int start, int end)
{
    vector<bool> used;
    vector<int> key;
    vector<int> p;
    priority_queue<Edge> q;
    used.resize(g.V);
    key.resize(g.V, INT_MAX);
    p.resize(g.V, -1);
    used[start] = 1;
    key[start] = 0;
    p[start] = -1;

    q.push(Edge(start, start, 0));

    while (!q.empty())
    {
        Edge curr = q.top();
        q.pop();
        //if (curr.wt > key[curr.b])
          //  continue;
        used[curr.b] = 1;
        for (int i = 0; i < g.g[curr.b].size(); i++)
        {
            int w = g.g[curr.b][i].v;
            int len = g.g[curr.b][i].wt;
            if (!used[w] && key[w] > len + key[curr.b])
            {
                p[w] = curr.b;
                key[w] = len + key[curr.b];
                q.push(Edge(curr.b, w, key[w]));
            }
        }
    }

    ofstream out("out.txt");
    for (int i = 0; i < g.V; i++)
        if (p[i] != -1)
            out << p[i]+1 << "-" << i + 1<< " " << key[i] << endl;
    out.close();

}


int main()
{
    Graph g;
    g.Init("in.txt");
    findPrimMST(g, 5, 9);
}