#include <bits/stdc++.h>
using namespace std;

// size”Å Union-Find–Ø : [0,n)
class UnionFind {
	vector<int> data;
	int root(int a) {
		return data[a] < 0 ? a : data[a] = root(data[a]);
	}
public:
	UnionFind(int n) : data(n, -1) {}
	bool Find(int a, int b) {
		return root(a) == root(b);
	}
	void Union(int a, int b) {
		a = root(a);
		b = root(b);
		if (a == b) return;
		if (data[a] < data[b]) {
			data[a] += data[b];
			data[b] = a;
		}
		else {
			data[b] += data[a];
			data[a] = b;
		}
	}
	int Size(int a) {
		return -data[root(a)];
	}
};

// rank”Å Union-Find–Ø : [0,n)
class UnionFind {
	vector<int> data;
	int root(int a) {
		return data[a] < 0 ? a : data[a] = root(data[a]);
	}
public:
	UnionFind(int n) : data(n, -1) {}
	bool Find(int a, int b) {
		return root(a) == root(b);
	}
	void Union(int a, int b) {
		a = root(a);
		b = root(b);
		if (a == b) return;
		if (data[a] < data[b]) {
			data[b] = a;
		}
		else {
			data[a] = b;
			if (data[a] == data[b]) data[b]--;
		}
	}
	int Rank(int a) {
		return -data[root(a)];
	}
};

// rank”ÅAsize•Û‘¶ Union-Find–Ø : [0,n)
class UnionFind {
	vector<int> data;
	vector<int> rank;
	int root(int a) {
		return data[a] < 0 ? a : data[a] = root(data[a]);
	}
public:
	UnionFind(int n) : data(n, -1), rank(n, 1) {}
	bool Find(int a, int b) {
		return root(a) == root(b);
	}
	void Union(int a, int b) {
		a = root(a);
		b = root(b);
		if (a == b) return;
		if (rank[a] < rank[b]) {
			data[a] += data[b];
			data[b] = a;
		}
		else {
			data[b] += data[a];
			data[a] = b;
			if (rank[a] == rank[b]) rank[b]++;
		}
	}
	int Size(int a) {
		return -data[root(a)];
	}
	int Rank(int a) {
		return rank[root(a)];
	}
};
