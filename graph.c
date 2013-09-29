
/***********************************************
 * Definicoes e constantes.
 *
 */

#include <string.h>

#define MAX 55
#define NO_ROUTE -1 // valor depende do caso!
#define INF 0x3f3f3f3f

typedef int Graph[MAX][MAX];

Graph g;
int n;

/***********************************************
 * Dijkstra && Prim.
 *
 */

void dijkstra_prim(int from, int dist[]) {
	char visited[MAX];
	int j, k;

	memcpy(dist,g[from],n*sizeof(int));
	memset(visited,0,sizeof(visited));
	visited[from] = 1, dist[from] = 0;

	while (1) {
		for (k=-1,j=0; j < n; j++)
			if (!visited[j] && dist[j] != NO_ROUTE && (k == -1 || dist[j] < dist[k]))
				k = j;
		if (k == -1) return;
		visited[k] = 1;
		for (j=0; j < n; j++)
			if (!visited[j] && g[k][j] != NO_ROUTE && (dist[j] == NO_ROUTE || g[k][j]+dist[k] < dist[j])) /* dijkstra */
				dist[j] = g[k][j] + dist[k];
			if (!visited[j] && g[k][j] != NO_ROUTE && (dist[j] == NO_ROUTE || g[k][j] < dist[j])) /* prim */
				dist[j] = g[k][j];
	}
}

// O(m*log(n)) -> par formado por (cost,to)

#include <vector>
#include <set>

using namespace std;

vector <vector <pair<int,int> > > g;
int n;

vector <int> dijkstra(int from) {
	vector <int> dist(n,INF);
	set <pair<int,int> > q;
	int i, j, cost, k, len;

	dist[from] = 0;
	q.insert(make_pair(0,from));

	while (!q.empty()) {
		i = q.begin()->second, len = g[i].size();
		q.erase(q.begin());
		for (k=0; k < len; k++) {
			j = g[i][k].second, cost = g[i][k].first;
			if (dist[i]+cost < dist[j]) {
				if (dist[j] != INF)
					q.erase(q.find(make_pair(dist[j],j)));
				dist[j] = dist[i] + cost;
				q.insert(make_pair(dist[j],j));
			}
		}
	}

	return dist;
}

/***********************************************
 * Menor caminho entre quaisquer vertices.
 *
 */

int min[MAX][MAX];

void floyd_warshal() {
	int i, j, k, k2;
	memcpy(min,g,sizeof(Graph));
	for (k=0; k < n; k++) {
		for (i=0; i < n; i++) {
			if (min[i][k] != NO_ROUTE) {
				for (j=0; j < n; j++) {
					if (min[k][j] != NO_ROUTE && (min[i][j] == NO_ROUTE || min[i][k] + min[k][j] < min[i][j]))
						min[i][j] = min[i][k] + min[k][j];
				}
			}
		}
	}
}

/***********************************************
 * Menor caminho com arestas negativas.
 *
 */

vector <int> bellman_ford(int from) {
	vector <int> prev(n,-1), dist(n,INF);
	int i, j, k;

	dist[from] = 0;
	for (k=1; k < n; k++) {
		for (i=0; i < n; i++) {
			for (j=0; j < g[i].size(); j++) {
				if (dist[g[i][j].second] > dist[i] + g[i][j].first) {
					dist[g[i][j].second] = dist[i] + g[i][j].first;
					prev[g[i][j].second] = i;
				}
			}
		}
	}

	for (i=0; i < n; i++) {
		for (j=0; j < g[i].size(); j++) {
			if (dist[g[i][j].second] > dist[i] + g[i][j].first)
				puts("ciclo com peso negativo!");
		}
	}

	return dist;
}

/***********************************************
 * Fluxo maximo entre 2 vertices. O(n^3)
 * -> para O(v*e), utilizar uma matriz de adjacencia com arestas bidirecionais
 *    e nao esquecer de limpar size[] e remover os memsets de g[][] e flow[][].
 */

#define MAXG 105
#define edge(i,j) g[i][j] = 1, adj[i][size[i]++] = j, adj[j][size[j]++] = i
int adj[MAXG][MAXG], size[MAXG];

#define min(a,b) ((a) < (b) ? (a) : (b))

Graph flow;
int q[MAXG], qf, qb, prev[MAXG];

int max_flow(Graph g, int n, int from, int to) {
	int total = 0, inc, i, j, u, v;

	for (i=0; i < n; i++) {
		memset(flow[i],0,sizeof(flow[i]));
	}

	while (1) {
		memset(prev,-1,sizeof(prev));

		qf = 0;
		prev[q[qf] = from] = -2;
		while (qf > -1 && prev[to] == -1) {
			//for (u=q[qf--],j=0; j < size[u]; j++) {
			//	v = adj[u][j];
			for (u=q[qf--],v=0; v < n; v++) {
				if (prev[v] == -1 && flow[u][v]-flow[v][u] < g[u][v])
					prev[q[++qf] = v] = u;
			}
		}

		if (prev[to] == -1)
			break;

		inc = INF;
		for (i=to; prev[i] >= 0; i=prev[i])
			inc = min(inc, g[prev[i]][i]-flow[prev[i]][i]+flow[i][prev[i]]);
		for (i=to; prev[i] >= 0; i=prev[i])
			flow[prev[i]][i] += inc;
		total += inc;
	}

	return total;
}

/***********************************************
 * Matching em grafo bipartido
 * -> Inicializar nl e nr.
 *
 */

int nr, nl, mate[MAX];
char visited[MAX];

int dfs(int i) {
	int j;
	if (visited[i]) return 0;
	visited[i] = 1;
	for (j=0; j < nr; j++) {
		if (g[i][j] && (mate[j] == -1 || dfs(mate[j]))) {
			mate[j] = i;
			return 1;
		}
	}
	return 0;
}

int match() {
	int i, ans = 0;
	memset(mate,-1,sizeof(mate));
	for (i=0; i < nl; i++) {
		memset(visited,0,sizeof(visited));
		ans += dfs(i);
	}
	return ans;
}

/***********************************************
 * Min-cut entre quaisquer vertices
 *
 */

int stoer_wagner(Graph g, int n) {
	int order[MAX], used[MAX], merged[MAX], i, j, k, p, cut, ans = INF;
	Graph aux;

	memset(merged,0,sizeof(merged));

	for (p=1; p < n; p++) {
		memset(used,0,sizeof(used));
		memcpy(aux,g,sizeof(Graph));
		used[order[0] = 0] = 1;
		for (k=1; k < n; k++) {
			for (i=-1,j=1; j < n; j++) {
				if (!used[j] && !merged[j] && (i == -1 || aux[0][j] > aux[0][i]))
					i = j;
			}
			if (i == -1) break;
			used[order[k] = i] = 1;
			for (j=1; j < n; j++) {
				aux[0][j] += aux[i][j];
				aux[j][0] += aux[i][j];
				aux[i][j] = aux[j][i] = 0;
			}
		}
		for (i=cut=0; i < n; i++) {
			if (order[k-2] != i) {
				g[order[k-2]][i] += g[order[k-1]][i];
				g[i][order[k-2]] += g[order[k-1]][i];
			}
			cut += g[order[k-1]][i];
			g[order[k-1]][i] = g[i][order[k-1]] = 0;
		}
		merged[order[k-1]] = 1;
		if (cut < ans) ans = cut;
	}

	return ans;
}

/***********************************************
 * Caminho por todas as arestas. (inicializar size!)
 *
 */

int path[1005], size;

void euler_tour(int i) {
	for (int j=0; j < n; j++) {
		if (g[i][j]) {
			g[i][j]--, g[j][i]--;
			euler_tour(j);
		}
	}
	path[size++] = i;
}

/***********************************************
 * Pontes e articulacoes.
 *
 */

int depth[MAX], parent[MAX], low[MAX], is_ap[MAX];
int bridges, ap;

void dfs(int i, int d) {
	int j, cont = 0;
	low[i] = depth[i] = d;
	for (j=0; j < n; j++) {
		if (g[i][j]) {
			if (depth[j] == -1) {
				parent[j] = i;
				cont++; // ap
				dfs(j,d+1);
				if (low[j] < low[i])
					low[i] = low[j];
				if (low[j] == depth[j]) /* i-j is a bridge! */
					bridges++;
				if (!is_ap[i] && d > 0 && low[j] >= depth[i]) // ap
					is_ap[i] = 1, ap++;
			}
			else if (j != parent[i] && depth[j] < low[i])
				low[i] = depth[j];
		}
	}
	if (!is_ap[i] && d == 0 && cont > 1) // ap
		is_ap[i] = 1, ap++;
}

int find(Graph g, int n) {
	int i;
	bridges = ap = 0;
	memset(is_ap,0,sizeof(is_ap)); // ap
	memset(depth,-1,sizeof(depth));
	for (i=0; i < n; i++) {
		if (depth[i] == -1) {
			parent[i] = i;
			dfs(i,0);
		}
	}
	return bridges * ap; // escolher um
}

/***********************************************
 * Componentes fortemente conexos.
 *
 */

int stack[MAX], top, instack[MAX];
int low[MAX], depth[MAX], id, comp[MAX], ncomp;

void tarjan(int i) {
	int j;
	depth[i] = low[i] = id++;
	stack[++top] = i;
	instack[i] = 1;
	for (j=0; j < n; j++) {
		if (g[i][j] != NO_ROUTE) {
			if (depth[j] == -1) {
				tarjan(j);
				if (low[j] < low[i])
					low[i] = low[j];
			}
			else if (instack[j]) {
				if (low[j] < low[i])
					low[i] = low[j];
			}
		}
	}
	if (low[i] == depth[i]) {
		do {
			comp[stack[top]] = ncomp;
			instack[stack[top]] = 0;
		} while (stack[top--] != i);
		ncomp++;
	}
}

int scc() {
	int i;
	memset(depth,-1,sizeof(depth));
	memset(instack,0,sizeof(instack));
	top = -1, ncomp = id = 0;
	for (i=0; i < n; i++)
		if (depth[i] == -1)
			tarjan(i);
	return ncomp;
}

/***********************************************
 * LCA
 * -> restriction: for each 0 <= i < n, i > parent[i]
 */

#define MAX 100000
#define MAXLOG 18

int parent[MAX];
int depth[MAX];
int table[MAX][MAXLOG];

void init(int parent[], int n, int depth[]) {
    depth[0] = 0;
    for (int i=1; i < n; ++i) {
        depth[i] = depth[parent[i]] + 1;
    }

    memset(table, -1, sizeof(table));

    for (int i=1; i < n; ++i) {
        table[i][0] = parent[i];
        for (int k=1; table[parent[i]][k-1] != -1; ++k) {
            table[i][k] = table[table[i][k-1]][k-1];
        }
    }
}

int lca(int a, int b) {

    if (depth[a] < depth[b]) {
        swap(a, b);
    }

    int subtract = depth[a] - depth[b];
    for (int i=0; (1 << i) <= subtract; ++i) {
        if (subtract & 1 << i) {
            a = table[a][i];
        }
    }

    if (a == b) return a;

    subtract = depth[b];
    for (int i=MAXLOG-1; i >= 0; --i) {
        if (table[a][i] != -1 && table[a][i] != table[b][i]) {
            a = table[a][i];
            b = table[b][i];
        }
    }

    return (a == 0) ? 0 : parent[a];
}

/***********************************************
 * Verifica se o grafo é cordal.
 * -> true se a ordem retornada por mcs (maximum cardinality search) é uma eliminacao perfeita, ou seja,
 *    uma ordenação onde para qualquer vértice v, seus vizinhos que ocorrem após v na ordenação formam um clique
 * -> a ordenação pode estar ao contrario, mas isso é inútil mesmo...
 * -> lista de adjacência bidirecional!
 */

int g[MAX][MAX], gsize[MAX], n;

char neighbor[MAX][MAX];

int is_peo(int order[]) {
	int reverse[MAX];
	int i, j, k, pos;

	for (i=0; i < n; i++) {
		reverse[order[i]] = i;
	}

	memset(neighbor,0,sizeof(neighbor));

	for (i=0; i < n; i++) {
		j = -1;
		for (pos=0; pos < gsize[order[i]]; pos++) {
			k = g[order[i]][pos];
			if (reverse[k] < i) {
				neighbor[order[i]][k] = 1;
				if (j == -1 || reverse[k] > reverse[j]) {
					j = k;
				}
			}
		}
		if (j == -1) {
			continue;
		}
		for (pos=0; pos < gsize[order[i]]; pos++) {
			k = g[order[i]][pos];
			if (reverse[k] < i && j != k && !neighbor[j][k]) {
				return 0;
			}
		}
	}

	return 1;
}

void mcs(int order[]) {
	int w[MAX], numbered[MAX];
	int i, j, k;

	memset(w,0,sizeof(w));
	memset(numbered,0,sizeof(numbered));

	for (i=n-1; i >= 0; i--) {
		k = -1;
		for (j=0; j < n; j++) {
			if (!numbered[j] && (k == -1 || w[j] > w[k])) {
				k = j;
			}
		}

		order[n-i-1] = k;
		numbered[k] = 1;

		for (j=0; j < gsize[k]; j++) {
			if (!numbered[g[k][j]]) {
				w[g[k][j]]++;
			}
		}
	}
}

int main() {

}
