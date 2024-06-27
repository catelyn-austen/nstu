#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

using namespace std;

int main() {
    vector<vector<int>> graph;
    ofstream out("out.txt");
    int n, ni, a; // n - количество файлов, ni - i-тый файл, a - текущая вершина

    ifstream nfile("n.txt");
    nfile >> n;
    nfile.close();

    for (int i = 1; i <= n; i++) {
        ifstream file(to_string(i) + ".txt");
        file >> ni;
        if (graph.size() < ni) graph.resize(ni);
        vector<int> v;
        for (int j = 0; j < ni; j++) { // заносим каждую вершину из i-того файла в вектор вершин
            file >> a;
            if (graph.size() < a) graph.resize(a);
            v.push_back(a - 1);
        }
        for (int k = 0; k < v.size(); k++) { // связываем вершины из i-того файла
            for (int l = k + 1; l < v.size(); l++) {
                graph[v[k]].push_back(v[l]);
                graph[v[l]].push_back(v[k]);
            }
        }
        file.close();
    }
    

    for (int i = 0; i < graph.size(); i++) {
        out << i + 1 << " : ";
        sort(graph[i].begin(), graph[i].end());
        graph[i].erase(unique(graph[i].begin(), graph[i].end()), graph[i].end());
        for (int j = 0; j < graph[i].size(); j++) {
            out << graph[i][j] + 1 << ", ";
        }
        out << endl;
    }
    out.close();
}