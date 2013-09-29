
// ########## LIS ##########

void lis(const vector<int>& v, vector<int>& asw) {
	int n = v.size();

	if (n == 0) {
		return;
	}

	vector <int> pd(n, 0), pd_index(n), pred(n);
	int maxi = 0, x, j, ind = 0;

	for (int i=0; i < n; ++i)  {
		x = v[i];
		j = lower_bound(pd.begin(), pd.begin() + maxi, x) - pd.begin();
		pd[j] = x;
		pd_index[j] = i;
		if (j == maxi) {
			maxi++;
			ind = i;
		}
		pred[i] = j ? pd_index[j-1] : -1;
	}

	int pos = maxi - 1, k = v[ind];
	asw.resize(maxi);

	while (pos >= 0) {
		asw[pos--] = k;
		ind = pred[ind];
		k = v[ind];
	}
}

// ########## KMP ##########

#define MAXLEN 1005

int kmp(char s[], char sub[]) {
	int F[MAXLEN], i, k;

	F[0] = -1;
	for (i=0; sub[i]; i++) {
		F[i+1] = F[i] + 1;
		while (F[i+1] > 0 && sub[i] != sub[F[i+1]-1])
			F[i+1] = F[F[i+1]-1] + 1;
	}

	for (i=k=0; s[i]; i++) {
		while (1) {
			if (s[i] == sub[k]) {
				if (!sub[++k]) {
					return 1;
					k = F[k];
				}
				break;
			}
			if (!k) break;
			k = F[k];
		}
	}

	return 0;
}

// ########## LIS O(N*logN) ##########

int lis(vector <int> &v) {
	vector <int> table(v.size(), 0);
	int ans = 0;
	for (int i=0; i < v.size(); i++) {
		int j = lower_bound(table.begin(), table.begin() + ans, v[i]) - table.begin();
		table[j] = v[i];
		if (j == ans) ans++;
	}
	return ans;
}

// ########## Knight Distance ##########

long long dist(long long x1, long long y1, long long x2, long long y2) {
	long long dx = abs(x2-x1), dy = abs(y2-y1), lb = (dx+1)/2;
	lb = max( lb, (dy+1)/2 ); lb = max( lb, (dx+dy+2)/3 );
	if ( lb % 2 != (dx+dy) % 2 ) lb++;
	if ( abs(dx) == 1 && dy == 0 ) return 3;
	if ( abs(dy) == 1 && dx == 0 ) return 3;
	if ( abs(dx) == 2 && abs(dy) == 2 ) return 4;
	return lb;
}

// ########## Josephu's Problem ##########

long long survivor(long long n, long long d) {
	long long K = 1;
	while (K <= (d-1)*n) {
		K = (d * K + d - 2) / (d - 1);
	}
	return d * n + 1 - K;
}

// ########## returns the index (0-based!!) of s among all the possible anagrams ##########

int freq[256];
long long count;

long long find_index(char *s, int n) {
	long long ans;
	int i;

	if (n < 2) {
		freq[*s-'a']++;
		count = 1;
		return 0;
	}

	ans = find_index(s+1,n-1);
	count = (count * n) / ++freq[*s-'a'];

	for (i=0; i < *s-'a'; i++) {
		ans += freq[i] * count / n;
	}

	return ans;
}

// ########## returns the index k-th (1-based!!) anagram of s ##########

long long fat(int n) {
	if (n <= 1) return 1;
	return n * fat(n-1);
}

string getAnagram(string s, int k) {
	int n = s.length();
	vector <int> freq(26,0);
	vector <int> index(26);
	long long acc = 0, count;

	if (!n) return "";

	sort(s.begin(),s.end());

	for (int i=0; i < n; i++) {
		freq[s[i]-'a']++;
		index[s[i]-'a'] = i;
	}

	count = fat(n);
	for (int i=0; i < 26; i++) {
		count /= fat(freq[i]);
	}

	for (int i=0; i < 26; i++) {
		if (freq[i] && acc+count*freq[i]/n >= k) {
			s[index[i]] = s[0];
			return (char)('a'+i) + getAnagram(s.substr(1),k-acc);
		}
		acc += count * freq[i] / n;
	}

	return "";
}

// ########## Stable marriage problem ##########

#define MAX 505

int men[MAX][MAX], women[MAX][MAX];
int order[MAX][MAX];

void marry(int men[MAX][MAX], int women[MAX][MAX], int n, vector <int> &ans, vector <int> &ans_rev) {
	vector <int> proposals(n,0);
	stack <int> singles;
	int i, j, k;

	for (i=0; i < n; i++) {
		for (j=0; j < n; j++) {
			order[i][women[i][j]] = j;
		}
	}

	ans = vector <int>(n,-1);
	ans_rev = vector <int>(n,-1);

	for (i=0; i < n; i++) {
		singles.push(i);
	}

	while (!singles.empty()) {

		i = singles.top(); singles.pop();
		j = men[i][proposals[i]++];
		k = ans_rev[j];

		if (k == -1) {
			ans[i] = j;
			ans_rev[j] = i;
		}
		else if (order[j][i] < order[j][k]) {
			ans[k] = -1;
			ans[i] =  j;
			ans_rev[j] = i;
			singles.push(k);
		}
		else {
			singles.push(i);
		}
	}
}

// ########## Fenwick Tree (BIT) ##########

// 1D

#define MAX 10005

int tree[MAX];

void create(int t[], int n) {
	memset(t,0,n*sizeof(int));
}

int query(int tree[], int from, int to) {
	int sum;

	if (from != 0) {
		return query(tree, 0, to) - query(tree, 0, from-1);
	}

	for (sum=0; to >= 0; to = (to & (to + 1)) - 1)
		sum += tree[to];
	return sum;

}

void update(int tree[], int k, int inc, int n) {
	for ( ; k < n; k |= k + 1) {
		tree[k] += inc;
	}
}

// 2D

#define MAX1 1001
#define MAX2 5001

int tree[MAX1][MAX2];

void create(int n1, int n2) {
	int i;
	for (i=0; i < n1; i++) {
		memset(tree[i],0,n2*sizeof(int));
	}
}

int query2(int x, int y) {
	int sum;
	for (sum=0; y >= 0; y=(y&(y+1))-1) {
		sum += tree[x][y];
	}
	return sum;
}

int query(int x1, int y1, int x2, int y2) {
	int sum;

	if (x1 != 0) {
		return query(0,y1,x2,y2) - query(0,y1,x1-1,y2);
	}

	for (sum=0; x2 >= 0; x2=(x2&(x2+1))-1) {
		sum += query2(x2,y2);
		if (y1) {
			sum -= query2(x2,y1-1);
		}
	}

	return sum;
}

void update(int x, int y, int inc, int n1, int n2) {
	int y2;
	for ( ; x < n1; x|=x+1) {
		for (y2=y; y2 < n2; y2|=y2+1) {
			tree[x][y2] += inc;
		}
	}
}

// ########## INTERVAL TREE ##########
// Importante: todos os intervalos sao no formato [l,r)
// -> Para fazer query em intervalos, usar ideia parecida com a do count_occurrences para podar a busca

#define MAX 100005

struct Node {
	int l, r, x, count;
	int cl, cr; // apenas para count_occurrences
};

Node tree[4*MAX];

void create(int l, int r, int i = 1) {
	tree[i].l = l;
	tree[i].r = r;
	tree[i].x = (l + r) / 2;
	tree[i].count = tree[i].cl = tree[i].cr = 0;
	if (r-l == 1) return;
	if (l < tree[i].x) create(l,tree[i].x,2*i);
	if (r > tree[i].x) create(tree[i].x,r,2*i+1);
}

// inserir ou remover, w é a quantidade
int update(int l, int r, int w, int i = 1) {
	if (l <= tree[i].l && r >= tree[i].r) {
		tree[i].count += w;
	}
	else if (tree[i].r-tree[i].l > 1) {
		if (l < tree[i].x) tree[i].cl = update(l,r,w,i*2);
		if (r > tree[i].x) tree[i].cr = update(l,r,w,i*2+1);
	}
	return tree[i].count > 0 ? tree[i].r - tree[i].l : tree[i].cl + tree[i].cr;
}

// quantos intervalos cruzam o ponto x?
int query(int x, int i = 1) {
	int ans = tree[i].count;
	if (tree[i].r-tree[i].l == 1) return ans;
	if (x < tree[i].x) ans += query(x,2*i);
	if (x >= tree[i].x) ans += query(x,2*i+1);
	return ans;
}

// quantas casas no intervalo [l,r) sao cortadas por algum intervalo?
int count_occurrences(int l, int r, int i = 1) {
	int ans = 0;
	if (tree[i].count > 0) return r - l;
	if (tree[i].r - tree[i].l == 1) return 0;
	if (l == tree[i].l && r >= tree[i].x) ans += tree[i].cl;
	else if (l < tree[i].x && tree[i].cl) ans += count_occurrences(l,min(r,tree[i].x),2*i);
	if (r == tree[i].r && l < tree[i].x) ans += tree[i].cr;
	else if (r > tree[i].x && tree[i].cr) ans += count_occurrences(max(l,tree[i].x),r,2*i+1);
	return ans;
}

// ############################## BIT TODO STUFF ################################# //

// SCALE
void scale(int c){
	for (int i = 1 ; i <= MaxVal ; i++)
		tree[i] = tree[i] / c;
}

// FIND INDEX WITH GIVEN CUMULATIVE FREQUENCY
// if in tree exists more than one index with a same
// cumulative frequency, this procedure will return
// some of them (we do not know which one)
// bitMask - initialy, it is the greatest bit of MaxVal
// bitMask store interval which should be searched
int find(int cumFre){
	int idx = 0; // this var is result of function

	while ((bitMask != 0) && (idx < MaxVal)){ // nobody likes overflow :)
		int tIdx = idx + bitMask; // we make midpoint of interval
		if (cumFre == tree[tIdx]) // if it is equal, we just return idx
			return tIdx;
		else if (cumFre > tree[tIdx]){
		        // if tree frequency "can fit" into cumFre,
		        // then include it
			idx = tIdx; // update index
			cumFre -= tree[tIdx]; // set frequency for next loop
		}
		bitMask >>= 1; // half current interval
	}
	if (cumFre != 0) // maybe given cumulative frequency doesn't exist
		return -1;
	else
		return idx;
}

// FIND INDEX WITH GIVEN CUMULATIVE FREQUENCY 2
// if in tree exists more than one index with a same
// cumulative frequency, this procedure will return
// the greatest one
int findG(int cumFre){
	int idx = 0;

	while ((bitMask != 0) && (idx < MaxVal)){
		int tIdx = idx + bitMask;
		if (cumFre >= tree[tIdx]){
		        // if current cumulative frequency is equal to cumFre,
		        // we are still looking for higher index (if exists)
			idx = tIdx;
			cumFre -= tree[tIdx];
		}
		bitMask >>= 1;
	}
	if (cumFre != 0)
		return -1;
	else
		return idx;
}

