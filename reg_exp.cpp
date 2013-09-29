// Regular Expression evaluation
//
// O problema:
// # - numero impar de caracteres
// * - numero par de caracteres
// alfabeto eh apenas digitos (trocar isdigit por isalpha caso necessario).
//
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

vector <int> qtdes;

string reduce(string s) {
	string ans = "";
	int acc = 0, iswildcard = 0;

	for (int i=0; s[i]; i++) {
		if (isdigit(s[i])) {
			if (iswildcard) {
				iswildcard = 0;
				qtdes[ans.size()] = acc;
				ans += acc % 2 ? '#' : '*';
				acc = 0;
			}
			ans += s[i];
		}
		else {
			iswildcard = 1;
			acc += (s[i] == '#');
		}
	}

	if (iswildcard) {
		qtdes[ans.size()] = acc;
		ans += acc % 2 ? '#' : '*';
	}

	return ans;
}

int kmp(string s, string sub) {
	int F[sub.size()+1], m, i, k;

	F[0] = -1;
	for (m=0; sub[m]; m++) {
		F[m+1] = F[m] + 1;
		while (F[m+1] > 0 && sub[m] != sub[F[m+1]-1])
			F[m+1] = F[F[m+1]-1] + 1;
	}

	for (i=k=0; s[i]; i++) {
		while (1) {
			if (s[i] == sub[k]) {
				if (++k == m) {
					return i - k + 1;
					k = F[k];
				}
				break;
			}
			if (!k) break;
			k = F[k];
		}
	}

	return -1;
}


int main() {
	string s, pat;
	int t;

	for (t=1; getline(cin,pat) && pat != "QUIT"; t++) {
		int c = 1, has_wildcard, i, j;

		for (i=0; pat[i] && isdigit(pat[i]); i++);
		has_wildcard = !!pat[i];

		qtdes = vector <int>(pat.length(),0);

		pat = reduce(pat);

		while (getline(cin,s) && s != "END" && s != "QUIT") {
			string S, PAT;
			int ans, pos, erased = 0, acc;

			if (!has_wildcard) {
				ans = (pat == s);
				goto END;
			}

			for (i=0; s[i] && isdigit(pat[i]) && pat[i] == s[i]; i++);
			for (j=0; j < s.length() && isdigit(pat[pat.length()-1-j]) && s[s.length()-1-j] == pat[pat.length()-1-j]; j++);

			erased += i;

			if (isdigit(pat[i]) || isdigit(pat[pat.length()-1-j])) {
				ans = 0;
				goto END;
			}

			S = s.substr(i,s.length()-i-j);
			PAT = pat.substr(i,pat.length()-i-j);

			while (PAT.length() > 1) {
				string aux = "";
				for (i=1; isdigit(PAT[i]); i++) {
					aux += PAT[i];
				}

				acc = 0;

				while (true) {
					pos = kmp(S,aux);

					if (pos == -1) {
						ans = 0;
						goto END;
					}

					if (acc+pos >= qtdes[erased] && (acc+pos)%2 == (PAT[0] == '#')) {
						break;
					}

					S = S.substr(pos+1);
					acc += pos + 1;
				}

				PAT = PAT.substr(1+aux.length());
				erased += 1 + aux.length();

				S = S.substr(pos+aux.length());
			}

			ans = S.length() >= qtdes[erased] && (S.length()%2 == (PAT[0] == '#'));

			END:
			printf("%d.%d. %s\n",t,c++,ans ? "match" : "not");
		}

		if (s == "QUIT") {
			break;
		}
	}

	return 0;
}
