#include <bits/stdc++.h>
using namespace std;

// Segment–Ø ‹æŠÔ[0,n)
template <typename T>
class SegmentTree {
	using func_t = function<T(T, T)>;
	const int n;
	const T id;
	func_t merge;
	vector<T> data;
	int size(int n) {
		int res;
		for (res = 1; res < n; res <<= 1);
		return res;
	}
	T sub(int l, int r, int node, int lb, int ub) {
		if (ub <= l || r <= lb) return id;
		if (l <= lb && ub <= r) return data[node];
		return merge(sub(l, r, node * 2, lb, (lb + ub) / 2), sub(l, r, node * 2 + 1, (lb + ub) / 2, ub));
	}
public:
	SegmentTree(int n_, T id_, func_t merge_) :
		n(size(n_)), id(id_), merge(merge_), data(size(n_) * 2, id_) {}
	SegmentTree(const vector<T>& data_, T id_, func_t merge_) :
		n(size(data_.size())), id(id_), merge(merge_), data(size(data_.size()) * 2, id_) {
		for (int i = 0; i < data_.size(); i++)
			data[i + n] = data_[i];
		for (int i = n - 1; i >= 0; i--)
			data[i] = merge(data[i * 2], data[i * 2 + 1]);
	}
	void Update(int p, T val) {
		p += n;
		data[p] = val;
		while (p >>= 1) data[p] = merge(data[p * 2], data[p * 2 + 1]);
	}
	void Add(int p, T val) {
		p += n;
		data[p] += val;
		while (p >>= 1) data[p] = merge(data[p * 2], data[p * 2 + 1]);
	}
	T Find(int l, int r) {
		return sub(l, r + 1, 1, 0, n);
	}
};

// “ñŸŒ³Segment–Ø ‹æŠÔ[0,n)	–¢Š®¬
template <typename T>
class SegmentTree2 {
	using func_t = function<T(T, T)>;
	const int w, h;
	const T id;
	func_t merge;
	vector<vector<T>> data;
	int size(int n) {
		int res;
		for (res = 1; res < n; res <<= 1);
		return res;
	}
	T sub(int li, int lj, int ri, int rj, int node, int lbi, int lbj, int ubi, int ubj) {
		if (ub <= l || r <= lb) return id;
		if (l <= lb && ub <= r) return data[node];
		return merge(sub(l, r, node * 2, lb, (lb + ub) / 2), sub(l, r, node * 2 + 1, (lb + ub) / 2, ub));
	}
public:
	SegmentTree2(int w_, int h_, T id_, func_t merge_) :
		w(size(w_)), h(size(h_)), id(id_), merge(merge_), data(size(w_) * 2, vector<T>(size(h_) * 2, id_)) {}
	SegmentTree2(vector<vector<T>> data_, T id_, func_t merge_) :
		n(size(data_.size())), id(id_), merge(merge_), data(size(data_.size()) * 2, id_) {
		for (int i = 0; i < data_.size(); i++)
			for (int j = 0; j < data_[0].size(); j++)
				data[i + w][j + h] = data_[i][j];
		for (int i = n - 1; i >= 0; i--)
			data[i] = merge(data[i * 2], data[i * 2 + 1]);
	}
	void Update(int pi, int pj, T val) {
		pi += w, pj += h;
		data[pi][pj] = val;
		for (int i = pi; i != 0; i >>= 1)
			for (int j = pj; j != 0;j >>= 1)
				data[pi][pj] = merge(merge(data[pi * 2][pj * 2], data[pi * 2][pj * 2 + 1]), merge(data[pi * 2 + 1][pj * 2], data[pi * 2 + 1][pj * 2 + 1]));
	}
	void Add(int pi, int pj, T val) {
		pi += w, pj += h;
		data[pi][pj] += val;
		for (int i = pi; i != 0; i >>= 1)
			for (int j = pj; j != 0; j >>= 1)
				data[pi][pj] = merge(merge(data[pi * 2][pj * 2], data[pi * 2][pj * 2 + 1]), merge(data[pi * 2 + 1][pj * 2], data[pi * 2 + 1][pj * 2 + 1]));
	}
	T Find(int li, int lj, int ri, int rj) {
		return sub(li, lj, ri + 1, rj + 1, 1, 0, 0, n, n);
	}
};
template <typename T>
class SegmentTree2 {
	using func_t = function<T(T, T)>;
	const int H, W;
	const T id;
	vector<vector<int>> data;
	int size(int n) {
		int res;
		for (res = 1; res < n; res <<= 1);
		return res;
	}
	int find_h(int li, int lj, int ri, int rj, int si, int ti, int k) {
		if (ri <= si || ti <= li) return id;
		if (li <= si && ti <= ri) return find_w(lj, rj, 0, W, k, 0);
		return merge(find_h(li, lj, ri, rj, si, (si + ti) / 2, 2 * k + 1), find_h(li, lj, ri, rj, (si + ti) / 2, ti, 2 * k + 2));
	}
	int find_w(int lj, int rj, int sj, int tj, int i, int k) {
		if (rj <= sj || tj <= lj) return id;
		if (lj <= sj && tj <= rj) return data[i][k];
		return merge(find_w(lj, rj, sj, (sj + tj) / 2, i, 2 * k + 1), find_w(lj, rj, (sj + tj) / 2, tj, i, 2 * k + 2));
	}
public:
	SegmentTree2(const vector<vector<int>> &f, T id_, func_t merge_): H(size(f.size())), W(size(f[0].size())), id(id_), merge(merge_) {
		data.assign(2 * H - 1, vector<int>(2 * W - 1, id));
		for (int i = 0; i < (int)f.size(); i++)
			for (int j = 0; j < (int)f[0].size(); j++)
				data[i + H - 1][j + W - 1] = f[i][j];
		for (int i = 2 * H - 2; i > H - 2; i--)
			for (int j = W - 2; j >= 0; j--)
				data[i][j] = merge(data[i][2 * j + 1], data[i][2 * j + 2]);
		for (int i = H - 2; i >= 0; i--)
			for (int j = 0; j < 2 * W - 1; j++)
				data[i][j] = merge(data[2 * i + 1][j], data[2 * i + 2][j]);
	}
	int Find(int li, int lj, int ri, int rj) {
		return find_h(li, lj, ri, rj, 0, H, 0);
	}
};
