
#define MAX 1000000

int primes[80000], nprimes;
char isprime[MAX+1];

void sieve() {
	int i, j;
	memset(isprime,-1,sizeof(isprime));
	isprime[0] = isprime[1] = nprimes = 0;
	for (j=2; j <= MAX; j++) {
		if (isprime[j] == -1) {
			isprime[j] = 1;
			primes[nprimes++] = j;
			for (i=2*j; i <= MAX; i+=j)
				isprime[i] = 0;
		}
	}
}

int phi[MAX+1];

void totient() {
	int i, j;
	for (i=1; i <= MAX; i++)
		phi[i] = i;
	for (i=2; i <= MAX; i+=2)
		phi[i] >>= 1;
	for (j=3; j <= MAX; j+=2) {
		if (phi[j] == j) {
			phi[j]--;
			for (i=2*j; i <= MAX; i+=j)
				phi[i] = phi[i] / j * (j - 1);
		}
	}
}

int is_prime(int n) {
	int lim, p;
	if (n < 5 || n%2 == 0 || n%3 == 0)
		return (n == 2 || n == 3);
	lim = (int) (sqrt(n) + 1);
	for (p = 5; p < lim; p += 6)
		if (n%p == 0 || n%(p+2) == 0)
			return 0;
	return 1;
}

int big_mod(int b, int p, int MOD) {
	long long temp;
	if (b == 0 || MOD == 1) return 0;
	if (p == 0 || b == 1) return 1;
	if (p == 1) return b%MOD;
	temp = big_mod(b%MOD,p/2,MOD);
	if (p&1)
		return (((temp*temp)%MOD)*(b%MOD))%MOD;
	return (temp*temp)%MOD;
}

void euclid(int a, int b, int *x, int *y, int *mdc) {
	int x2, y2;
	if (b == 0)
		*x = 1, *y = 0, *mdc = a;
	else {
		euclid(b,a%b,&x2,&y2,mdc);
		*x = y2, *y = x2-(a/b)*y2;
	}
}

// (a * b) % MOD sem overflow. ATENCAO: mod < 2^62
uint64 mult_mod(uint64 a, uint64 b, uint64 MOD) {
	uint64 y = (float64)a*(float64)b/MOD + (float64)1/2;
	y *= MOD;
	uint64 x = a * b;
	uint64 r = x - y;
	if ((long long)r < 0) r += MOD, y--;
	return r;
}

// quantidade de fatores D em N!
int count(int n, int d) {
	return (n < d) ? 0 : n/d + count(n/d,d);
}

int gcd(int a, int b) {
	return b ? gcd(b,a%b) : a;
}

int lcm(int a, int b) {
	return a / gcd(a,b) * b;
}

// evita overflow
#define INF 10000000000ll
long long lcm(long long a, long long b) {
	long long x = a / gcd(a,b);
	return (INF / x < b) ? INF + 1 : x * b;
}

// ax == b (mod c)
// usado tambem pra achar o inverso modular (b == 1)!
long long axbmodc(long long a, long long b, long long c) {
	return a ? (axbmodc(c % a, (a - b % a) % a, a) * c + b) / a : 0;
}

// ############### FRACTIONS ############### //

typedef long long number;

number gcd(number a, number b) {
	return b ? gcd(b,a%b) : a;
}

struct Frac {
	number n, d;
	bool ok;
	Frac(number a, number b, bool v = true) {
		ok = v && b;
		if (ok) {
			n = a/gcd(a,b), d = b/gcd(a,b);
			if (d < 0) n *= -1, d *= -1;
		}
	}
	Frac operator +(const Frac &f) const {return Frac(n*f.d+f.n*d,d*f.d,ok && f.ok);}
	Frac operator -(const Frac &f) const {return Frac(n*f.d-f.n*d,d*f.d,ok && f.ok);}
	Frac operator *(const Frac &f) const {return Frac(n*f.n,d*f.d,ok && f.ok);}
	Frac operator /(const Frac &f) const {return Frac(n*f.d,d*f.n,ok && f.ok);}
	void print() const {if (!ok) puts("INVALID");else if (d == 1) printf("%lld\n",n);else printf("%lld|%lld\n",n,d);}
};

// ############### expression evaluation ############### //

Frac eval(char s[], int from, int to) {
	for (int k=0; k < 3; k++) {
		for (int i=to,d=0; i >= from; i--) {
			if (s[i] == ')') d++;
			else if (s[i] == '(') d--;
			if (d > 0) continue;
			assert(d == 0);
			if (k == 0 && s[i] == '+') return eval(s,from,i-1) + eval(s,i+1,to);
			if (k == 0 && s[i] == '-') return eval(s,from,i-1) - eval(s,i+1,to);
			if (k == 1 && s[i] == '*') return eval(s,from,i-1) * eval(s,i+1,to);
			if (k == 1 && s[i] == '/') return eval(s,from,i-1) / eval(s,i+1,to);
			if (k == 2 && s[i] == '|') return eval(s,from,i-1) / eval(s,i+1,to);
		}
	}

	if (s[from] == '(' && s[to] == ')') return eval(s,from+1,to-1);

	number n = 0;
	while (from <= to) assert(s[from] >= '0' && s[from] <= '9'), n = n*10 + s[from++] - '0';
	return Frac(n,1);
}

// retorna a quantidade de numeros <= n que sao multiplos de algum elemento de v[]
// V1: genérica. DIZ A LENDA que percorrer do maior ao menor é mais rápido, mas acho que não.
// v = vetor de divisores, index = chamar com tamanho do vetor, n = valor a ser testado, LCM = chamar com 1
long long count_v1(vector <long long> &v, int index, long long n, long long LCM) {
	long long req, ans = 0;
	for (int i=0; i < index; ++i) {
		req = lcm(LCM,v[i]);
		if (req <= n) {
			ans += n / req - cnt(v,i,n,req);
		}
	}
	return ans;
}

// V2: muito mais rápida, mas só funciona quando gcd(v[i],v[j]) == 1, para qualquer i != j.
// v = vetor de divisores, index = inicializar com tamanho do vetor, n = valor a ser testado,
int count_v2(int v[], int index, int n) {
	int ans = 0, i;
	for (i=0; i < index && v[i] <= n; i++) {
		ans += n / v[i] - count(v,i,n/v[i]);
	}
	return ans;
}

