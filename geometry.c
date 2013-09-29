#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX 105
#define EPS 1e-7

typedef struct {
	int x, y;
} Point;

typedef Point Polygon[MAX];

typedef struct {
	double x, y;
} Point_d;

/***********************************************
 * Funcoes basicas
 *
 */

void copy(Polygon to, Polygon from, int n) {
	memcpy(to,from,sizeof(Point)*n);
}

int equals(double x, double y) {
	return x-EPS < y && x+EPS > y;
}

int direction(int x) {
	if (x == 0) return 0;
	return (x < 0) ? -1 : 1;
}

int prod_vet(int x1, int y1, int x2, int y2) {
	return x1*y2 - x2*y1;
}

double dist_pp(double x1, double y1, double x2, double y2) {
	return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

double dist_pr(double x, double y, double a, double b, int inf) {
	double tgaux, baux;
	if (inf) return fabs(b-x);
	if (equals(a,0.0)) return fabs(b-y);
	tgaux = -1.0 / a;
	baux = y - tgaux * x;
	return dist_pp(x,y,(baux-b)/(a-tgaux),a*(baux-b)/(a-tgaux)+b);
}

/***********************************************
 * Convex-hull de um poligono anti-horario.
 *
 * -> destroi o poligono original!
 * -> retorna o numero de pontos do hull.
 */

int xmin, ymax;

int cmp_sort(Point *a, Point *b) {
	if (a->x == xmin && b->x == xmin)
		return b->y - a->y;
	if (a->x == xmin)
		return -1;
	if (b->x == xmin)
		return 1;
	if ((a->y == ymax && b->y == ymax) || prod_vet(a->x-xmin,a->y-ymax,b->x-xmin,b->y-ymax) == 0)
		return a->x - b->x;
	return -prod_vet(a->x-xmin,a->y-ymax,b->x-xmin,b->y-ymax);
}

void sort(Polygon p, int n) {
	int i, j, min = 0;
	Point aux;

	for (i=1; i < n; i++) {
		if (p[i].x <  p[min].x || (p[i].x == p[min].x && p[i].y > p[min].y))
			min = i;
	}
	if (min != 0)
		aux = p[0], p[0] = p[min], p[min] = aux;
	xmin = p[0].x, ymax = p[0].y;

	qsort(&p[1],n-1,sizeof(Point),(void*)cmp_sort);

	for (i=n-1; i > 0 && prod_vet(p[i].x-p[0].x,p[i].y-p[0].y,p[i-1].x-p[0].x,p[i-1].y-p[0].y) == 0; i--);
	for (j=0; j < (n-i)/2; j++)
		aux = p[i+j], p[i+j] = p[n-j-1], p[n-j-1] = aux;
}

int convex_hull(Polygon p, int n, Polygon hull) {
	int qtdepontos = 2, i;
	double x1, y1, x2, y2;

	sort(p,n);

	hull[0] = p[0], hull[1] = p[1];
	for (i=2; i <= n; i++) {
		do {
			x1 = p[i%n].x - hull[qtdepontos-1].x;
			y1 = p[i%n].y - hull[qtdepontos-1].y;
			x2 = hull[qtdepontos-2].x - hull[qtdepontos-1].x;
			y2 = hull[qtdepontos-2].y - hull[qtdepontos-1].y;
			if (prod_vet(x1,y1,x2,y2) <= 0 && i != n) // verificar a segunda condicao, pontos colineares??
				qtdepontos--;
			else
				break;
		}
		while (qtdepontos > 1);
		if (i != n)
			hull[qtdepontos++] = p[i];
	}

	return qtdepontos;
}

/***********************************************
 * Interseccao de retas e segmentos.
 *
 */

int intersect(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double *x, double *y) {
	double a1, b1, a2, b2;
	if (equals(x0,x1) && equals(x2,x3)) return 0; // tratar retas verticais?
	if (!equals(x0,x1)) a1 = (y1-y0)/(x1-x0), b1 = y0 - a1*x0;
	if (!equals(x2,x3)) a2 = (y3-y2)/(x3-x2), b2 = y2 - a2*x2;
	if (equals(x0,x1))
		*x = x0, *y = a2*x0 + b2;
	else if (equals(x2,x3))
		*x = x2, *y = a1*x2 + b1;
	else {
		if (equals(a1,a2)) return 0; // tratar retas colineares?
		*x = (b2 - b1) / (a1 - a2), *y = a1*(*x) + b1;
	}
	if (!equals(dist_pp(*x,*y,x0,y0)+dist_pp(*x,*y,x1,y1),dist_pp(x0,y0,x1,y1)))
		return 0; // se (p0,p1) eh segmento
	if (!equals(dist_pp(*x,*y,x2,y2)+dist_pp(*x,*y,x3,y3),dist_pp(x2,y2,x3,y3)))
		return 0; // se (p2,p3) eh segmento
	return 1;
}

/***********************************************
 * Interseccao de segmentos (v2).
 *
 * -> trocar para <= e >= se nao podem encostar
 */

int segment_intersect(int x1, int y1, int x2, int y2, int x3, int y3, int x4, int y4) {
	return ((prod_vet(x1-x2,y1-y2,x3-x2,y3-y2) < 0 && prod_vet(x1-x2,y1-y2,x4-x2,y4-y2) > 0) ||
			(prod_vet(x1-x2,y1-y2,x3-x2,y3-y2) > 0 && prod_vet(x1-x2,y1-y2,x4-x2,y4-y2) < 0)) &&
		   ((prod_vet(x4-x3,y4-y3,x1-x3,y1-y3) < 0 && prod_vet(x4-x3,y4-y3,x2-x3,y2-y3) > 0) ||
			(prod_vet(x4-x3,y4-y3,x1-x3,y1-y3) > 0 && prod_vet(x4-x3,y4-y3,x2-x3,y2-y3) < 0));
}

/***********************************************
 * Ponto dentro de um poligono anti-horario qualquer.
 *
 */

int point_in_pol(int x, int y, Polygon p, int n) {
	int i, cuts = 0;
	double xaux;

	for (i=1; i <= n; i++) {
		if ((p[i-1].x == x && p[i-1].y == y) || (prod_vet(p[i-1].x-x,p[i-1].y-y,p[i%n].x-x,p[i%n].y-y) == 0 &&
			equals(dist_pp(x,y,p[i-1].x,p[i-1].y)+dist_pp(x,y,p[i%n].x,p[i%n].y),dist_pp(p[i%n].x,p[i%n].y,p[i-1].x,p[i-1].y)))) {
			return 1;
		}
		if (p[i-1].y == y && p[i-1].y == p[i%n].y) {
			if (p[i-1].x > x && p[i%n].x > x && ((p[(n+i-2)%n].y < y && p[(i+1)%n].y > y) || (p[(n+i-2)%n].y > y && p[(i+1)%n].y < y)))
				cuts++;
			continue;
		}
		if (p[i%n].y == y && p[i%n].x > x && ((p[i-1].y < y && p[(i+1)%n].y > y) || (p[i-1].y > y && p[(i+1)%n].y < y))) {
			cuts++;
			continue;
		}
		if (p[i%n].x == p[i-1].x)
			xaux = p[i%n].x;
		else
			xaux = (y - p[i%n].y + (p[i%n].y-p[i-1].y)/(double)(p[i%n].x-p[i-1].x)*p[i%n].x) /
																((p[i%n].y-p[i-1].y)/(double)(p[i%n].x-p[i-1].x));
		if (xaux+EPS > x && ((y > p[i%n].y && y < p[i-1].y) || (y < p[i%n].y && y > p[i-1].y)))
			cuts++;
	}

	return cuts & 1;
}

/***********************************************
 * Ponto dentro de um poligono anti-horario e convexo.
 *
 */

int point_in_pol_2(int x, int y, Polygon p, int n) {
	int i;
	for (i=1; i <= n; i++) {
		if (prod_vet(p[i-1].x-x, p[i-1].y-y, p[i%n].x-x, p[i%n].y-y) < 0)
			return 0;
	}
	return 1;
}

/***********************************************
 * Verifica se o poligono eh convexo.
 *
 */

int is_convex(Polygon p, int n) {
	int i, dir= 0, aux;
	for (i=1; i <= n; i++) {
		aux = direction(prod_vet(p[i-1].x-p[i%n].x,p[i-1].y-p[i%n].y,p[(i+1)%n].x-p[i%n].x,p[(i+1)%n].y-p[i%n].y));
		if (!dir) dir = aux;
		if (dir && aux && aux != dir) return 0;
	}
	return 1;
}

/***********************************************
 * Verifica se o poligono eh anti-horario.
 *
 */

int ccw(Polygon p, int n) {
	int e, d, i, aux;
	for (e=d=0,i=1; i < n-1; i++) {
		aux = prod_vet(p[i-1].x-p[i].x,p[i-1].y-p[i].y,p[i+1].x-p[i].x,p[i+1].y-p[i].y);
		if (aux > 0) d++;
		else if (aux < 0) e++;
	}
	return e > d;
}

/***********************************************
 * Area do poligono. Negativa caso poligono horario.
 *
 */

double area(Polygon p, int n) {
	int i;
	double ans = 0.0;
	for (i=1; i <= n; i++)
		ans += prod_vet(p[i-1].x,p[i-1].y, p[i%n].x, p[i%n].y);
	return ans/2.0;
}

/***********************************************
 * Perimetro.
 *
 */

double perimeter(Polygon p, int n) {
	int i;
	double ans = 0.0;
	for (i=1; i <= n; i++)
		ans += sqrt((p[i%n].x-p[i-1].x)*(p[i%n].x-p[i-1].x) + (p[i%n].y-p[i-1].y)*(p[i%n].y-p[i-1].y));
	return ans;
}

/***********************************************
 * Centro de massa.
 *
 */

void centroid(Polygon p, int n, Point_d *r) {
	int i, x = 0, y = 0;
	double s = area(p,n);
	for (i=1; i <= n; i++) {
		x += (p[i-1].x + p[i].x)*prod_vet(p[i-1].x, p[i-1].y, p[i%n].x, p[i%n].y);
		y += (p[i-1].y + p[i].y)*prod_vet(p[i-1].x, p[i-1].y, p[i%n].x, p[i%n].y);
	}
	r->x = x/(6.0*s);
	r->y = y/(6.0*s);
}

/***********************************************
 * Circulo definido por 3 pontos.
 *
 * -> resultado nas variaveis globais!
 */

double x, y, r;

void find_circle(Point *p1, Point *p2, Point *p3) {
	double x1, y1, x2, y2, x3, y3, temp;
	x1 = p1->x, y1 = p1->y;
	x2 = p2->x, y2 = p2->y;
	x3 = p3->x, y3 = p3->y;
	if (equals(x1,x2)) {
		temp = y3, y3 = y2, y2 = temp;
		temp = x3, x3 = x2, x2 = temp;
	}
	y = ((x1-x2)*(x1*x1+y1*y1-x3*x3-y3*y3)-(x1-x3)*(x1*x1+y1*y1-x2*x2-y2*y2)) / (2*((y2-y1)*(x1-x3)-(y3-y1)*(x1-x2)));
	x = (x1*x1+y1*y1-x2*x2-y2*y2+2*y*(y2-y1)) / (2*(x1-x2));
	r = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1));
}

/***********************************************
 * Menor circulo que engloba um conjunto de pontos.
 *
 * -> utiliza find_circle()
 * -> soh funciona para pontos que facam parte do convex hull do conjunto!!!!!!!
 */

void min_circle(Polygon p, int n) {
	int i = 0, j = 1, k, posmin;
	double PI = acos(-1.0), min, a;
	while (1) {
		min = PI;
		for (k=0; k < n; k++) {
			if (k != i && k != j) {
				a = acos(((p[i].x-p[k].x)*(p[j].x-p[k].x)+(p[i].y-p[k].y)*(p[j].y-p[k].y)) / (sqrt((p[i].x-p[k].x)*(p[i].x-p[k].x)+(p[i].y-p[k].y)*(p[i].y-p[k].y))*sqrt((p[j].x-p[k].x)*(p[j].x-p[k].x)+(p[j].y-p[k].y)*(p[j].y-p[k].y))));
				if (a < min)
					min = a, posmin = k;
			}
		}
		if (min+EPS > PI/2) {
			x = (p[i].x + p[j].x) / 2.0;
			y = (p[i].y + p[j].y) / 2.0;
			r = sqrt((p[i].x-p[j].x)*(p[i].x-p[j].x) + (p[i].y-p[j].y)*(p[i].y-p[j].y)) / 2;
			return;
		}
		if (acos(((p[posmin].x-p[i].x)*(p[j].x-p[i].x)+(p[posmin].y-p[i].y)*(p[j].y-p[i].y)) / (sqrt((p[posmin].x-p[i].x)*(p[posmin].x-p[i].x)+(p[posmin].y-p[i].y)*(p[posmin].y-p[i].y))*sqrt((p[j].x-p[i].x)*(p[j].x-p[i].x)+(p[j].y-p[i].y)*(p[j].y-p[i].y))))-EPS > PI/2)
			i = posmin;
		else if (acos(((p[i].x-p[j].x)*(p[posmin].x-p[j].x)+(p[i].y-p[j].y)*(p[posmin].y-p[j].y)) / (sqrt((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y))*sqrt((p[posmin].x-p[j].x)*(p[posmin].x-p[j].x)+(p[posmin].y-p[j].y)*(p[posmin].y-p[j].y))))-EPS > PI/2)
			j = posmin;
		else {
			find_circle(&p[i],&p[j],&p[posmin]);
			return;
		}
	}
}

/***********************************************
 * Distancia sobre uma esfera.
 *
 */

double spherical_distance(double p_lat, double p_long, double q_lat, double q_long) {
	return acos(sin(p_lat)*sin(q_lat) +	cos(p_lat)*cos(q_lat)*cos(p_long-q_long))*6378.0;
}

/***********************************************
 * Eliminacao de Gauss
 *
 */

void gaussian_elimination(double a[MAX][MAX], double b[MAX], int n) {
	int i, j, k;
	int l, maxi;
	double f, aux;

	for (i=0; i < n; i++) {
		maxi = i;
		for (l=i; l < n; l++) {
			if (fabs(a[l][i]) > fabs(a[maxi][i])) {
				maxi = l;
			}
		}
		for (l=0; l < n; l++) {
			aux = a[i][l], a[i][l] = a[maxi][l], a[maxi][l] = aux;
		}
		aux = b[i], b[i] = b[maxi], b[maxi] = aux;

		for (k=i+1; k < n; k++) {
			f = a[k][i] / a[i][i];
			for (j=i; j < n; j++) {
				a[k][j] -= a[i][j] * f;
			}
			b[k] -= b[i] * f;
		}
	}

	for (i=n-1; i >= 0; i--) {
		b[i] = b[i] / a[i][i];
		a[i][i] = 1.0;

		for (j=i-1; j >= 0; j--) {
			b[j] -= a[j][i] * b[i];
			a[j][i] = 0.0;
		}
	}
}

/***********************************************
 * Menor distancia entre quaisquer 2 pontos: pseudo O(NlogN)
 * -> ordenar o vetor em relacao a coordenada X
 */

double closest_pair(Point p[], int n) {
	int nl = n/2, nr = n - nl;
	double dl, dr, d, ans;
	int i, j;

	dl = (nl > 1) ? closest_pair(p,nl) : INF;
	dr = (nr > 1) ? closest_pair(&p[nl],nr) : INF;
	ans = min(dl,dr);

	for (i=0; i < nl; i++) {
		if (p[nl].x-p[i].x > ans) continue;
		for (j=0; j < nr; j++) {
			if (fabs(p[nl+j].y-p[i].y) > ans) continue;
			d = dist_pp(&p[i],&p[nl+j]);
			ans = min(ans,d);
		}
	}

	return ans;
}

/***********************************************
 * Encontra retas que passam pela origem e tangenciam uma circunferencia
 * -> Retorna tudo como referencia
 */

void find_angles(double xc, double yc, double r, double *angi, double *angs, int *infi, int *infs, double *pxi, double *pyi, double *pxs, double *pys) {
	double tgi, tgs;

	if (equals(xc*xc+yc*yc,r*r)) {
		*pxi = *pyi = *pxs = *pys = 0;
		if (equals(xc,0.0)) {
			*infi = *infs = 0;
			if (yc < 0)
				*angi = M_PI, *angs = 0;
			else
				*angi = 0, *angs = M_PI;
		}
		else if (equals(yc,0.0)) {
			*infi = *infs = 1;
			if (xc < 0)
				*angi = M_PI/2, *angs = 3*M_PI/2;
			else
				*angi = 3*M_PI/2, *angs = M_PI/2;
		}
		else {
			*infi = *infs = 0;
			tgi = -xc/yc;
			if (tgi < 0) {
				if (yc-tgi*xc > 0)
					*angi = atan(tgi) + 2*M_PI, *angs = atan(tgi) + M_PI;
				else
					*angi = atan(tgi) + M_PI, *angs = atan(tgi) + 2*M_PI;
			}
			else {
				if (yc-tgi*xc > 0)
					*angi = atan(tgi), *angs = atan(tgi) + M_PI;
				else
					*angi = atan(tgi) + M_PI, *angs = atan(tgi);
			}
		}
		return;
	}

	if (equals(r,xc)) {
		if (equals(yc,0.0)) {
			*infi = *infs = 1;
			*angi = 1.5*M_PI;
			*angs = 0.5*M_PI;
			*pxi = *pyi = *pxs = *pys = 0;
		}
		else if (yc+EPS > 0) {
			*infi = 0;
			tgi = (yc*yc-r*r)/(2*yc*xc);
			if (equals(yc,r))
				*pxi = xc;
			else
				*pxi = (yc - r)/tgi;
			*pyi = tgi*(*pxi);
			*angi = atan(tgi);
			if (r-EPS > yc)
				*angi += 2*M_PI;
			*infs = 1;
			*angs = 0.5*M_PI;
			*pxs = 0;
			*pys = yc;
		}
		else {
			*infi = 1;
			*angi = 1.5*M_PI;
			*pxi = 0;
			*pyi = yc;
			*infs = 0;
			tgs = (yc*yc-r*r)/(2*yc*xc);
			if (equals(-yc,r))
				*pxs = xc;
			else
				*pxs = (-yc - r)/tgs;
			*pys = tgs*(*pxs);
			*angs = atan(tgs);
			if (r+EPS < -yc)
				*angs += 2*M_PI;
		}
	}
	else if (equals(r,-xc)) {
		if (equals(yc,0.0)) {
			*infi = *infs = 1;
			*angi = 0.5*M_PI;
			*angs = 1.5*M_PI;
			*pxi = *pyi = *pxs = *pys = 0;
		}
		else if (yc+EPS > 0) {
			*infi = 1;
			*angi = 0.5*M_PI;
			*pxi = 0;
			*pyi = yc;
			*infs = 0;
			tgs = (yc*yc-r*r)/(2*yc*xc);
			if (equals(yc,r))
				*pxs = xc;
			else
				*pxs = (-yc + r)/tgs;
			*pys = tgs*(*pxs);
			*angs = atan(tgs) + M_PI;
		}
		else {
			*infi = 0;
			tgi = (yc*yc-r*r)/(2*yc*xc);
			if (equals(-yc,r))
				*pxi = xc;
			else
				*pxi = (-yc - r)/tgi;
			*pyi = tgi*(*pxi);
			*angi = atan(tgi) + M_PI;
			*infs = 1;
			*angs = 1.5*M_PI;
			*pxs = 0;
			*pys = yc;
		}
	}
	else {
		*infi = *infs = 0;
		tgi = (xc*yc - r*sqrt(xc*xc+yc*yc-r*r)) / (xc*xc - r*r);
		*angi = atan(tgi);
		*pxi = (xc+tgi*yc)/(1+tgi*tgi);
		*pyi = tgi*(*pxi);
		if ((*pxi)+EPS < 0)
			*angi += M_PI;
		else if ((*pyi)+EPS < 0)
			*angi += 2*M_PI;

		tgs = (xc*yc + r*sqrt(xc*xc+yc*yc-r*r)) / (xc*xc - r*r);
		*angs = atan(tgs);
		*pxs = (xc+tgs*yc)/(1+tgs*tgs);
		*pys = tgs*(*pxs);
		if ((*pxs)+EPS < 0)
			*angs += M_PI;
		else if ((*pys)+EPS < 0)
			*angs += 2*M_PI;
	}
}

int main() {

}
