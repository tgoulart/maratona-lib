
#define MAXV 100

typedef struct {
	int u, v, cost;
} Edge;

/***********************************************
 * Union-Find.
 */

int label[MAXV], size[MAXV];

void uf_init(int n) {
	int i;
	for (i=0; i < n; i++)
		size[label[i] = i] = 1;
}

int getlabel(int v) {
	while (v != label[v]) v = label[v];
	return v;
}

int uf_find(int u, int v) {
	return (getlabel(u) == getlabel(v));
}

void uf_union(int u, int v) {
	int i = getlabel(u), j = getlabel(v);
	if (i == j) return;
	if (size[i] < size[j])
		size[label[i] = j] += size[i];
	else
		size[label[j] = i] += size[j];
}

/***********************************************
 * Arvore geradora minima.
 */

int kruskal(Edge g[], int m, int n, int mst[]) {
	int i, qtde;
	//sort(g,len);
	uf_init(n);
	for (i=qtde=0; i < m && qtde < n-1; i++) {
		if (!uf_find(g[i].u,g[i].v)) {
			uf_union(g[i].u,g[i].v);
			mst[qtde++] = i;
		}
	}
	return qtde;
}

int main() {}
