#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <time.h>
#include <assert.h>
#include <algorithm>
#include <utility>
using namespace std;

struct List {
	int id;
	List *next;
} *_list, *_list_idle, ***list_head_aa; int _n_list;

#define a2i(x) ((x)=='-'? 26: (x)-'A')

int main (int argc, char **argv)
{
	char *infilename = argv[1], *outfilename = argv[2];
	double fthreshold = atof(argv[3]);
	clock_t _beg;

	FILE *fin = fopen(infilename, "r"); assert (fin);
	struct stat st; fstat(fileno(fin), &st); int filesize = st.st_size;
#define BUFFERSIZE 8192
	char buffer[BUFFERSIZE];
	int seql = strlen (fgets (buffer, BUFFERSIZE, fin)) -1;
	int nseq = filesize/(seql+1);
	printf ("nseq = %d, seql = %d\n", nseq, seql);
	assert (nseq*(seql+1) == filesize);
	_n_list = nseq*seql;
	char **AA = new char* [nseq];
	_list = _list_idle = new List [_n_list];
	list_head_aa = new List** [seql]; for (int i = 0; i < seql; i++) memset (list_head_aa[i] = new List* [27], 0, sizeof(List*)*27);

	int **cntAA = new int* [seql]; for (int i = 0; i < seql; i++) memset (cntAA[i] = new int [27], 0, sizeof(int)*27);
	_beg = clock ();
	for (int i = 0; i < nseq; i++) {
		char *&AAi = AA[i];
		if (i) assert (strlen(fgets (AAi = new char [seql+2], seql+2, fin)) == seql+1 && AAi[seql] == '\n'), AAi[seql] = 0;
		else strcpy (AAi = new char [seql+2], buffer), AAi[seql] = 0;
		for (int j = 0; j < seql; j++) {
			cntAA[j][AAi[j] = a2i(AAi[j])]++;
			List *x = _list_idle++; x->id = i;
			List *&y = list_head_aa[j][AAi[j]];
			x->next = y; y = x;
		}
	}
	// should be optimized to O(NL) from O(N*2 + NL), by maintaining a map from distance to linked list of seqs
	fclose (fin);
	for (int i = 0; i < seql; i++) {
		int c = 0, *ci = cntAA[i]; for (int j = 0; j < 27; j++) c += ci[j];
		assert (c == nseq);
	}
	printf ("%d\n", clock () - _beg);


	int threshold = fthreshold * seql;

	int *S = new int [nseq]; memset (S, 0, sizeof(int)*nseq);
	_beg = clock ();
	int *nsame = new int [nseq];
	for (int i = 0; i < nseq; i++) {
		int offset = 0;
		memset (nsame, 0, sizeof(int)*nseq);
		char *AAi = AA[i], AAij;
		for (int j = 0; j < seql; j++) if ((cntAA[j][AAij = AAi[j]]<<1) <= nseq) {
			for (List *x = list_head_aa[j][AAij]; x && x->id > i; x = x->next) nsame[x->id]++;
		} else {
			List **list_head_aaj = list_head_aa[j];
			for (int k = 0; k < 27; k++) if (AAij != k) for (List *x = list_head_aaj[k]; x && x->id > i; x = x->next) nsame[x->id]--;
			offset++;
		}
		for (int j = i+1; j < nseq; j++) if (nsame[j] + offset > threshold) S[i]++, S[j]++; S[i]++;
		if (!((i+1)%100)) {clock_t _t = clock () - _beg; printf ("%d\t%.4lf\t%.4lf\t%.4lf\n", i+1, _t*1e-6,  _t*1.e-6/(i+1), _t*1.e-6/(i+1) * (nseq-i-1));}
	}

//	for (int i = 0; i < nseq; i++) printf ("%d, ", S[i]); printf ("\n");

	FILE *fout = fopen (outfilename, "w");
	double fS = 0; for (int i = 0; i < nseq; i++) fS += 1. / S[i];
	printf ("%.6lf\n", fS);
	fprintf (fout, "%.6lf", fS);
	fclose (fout);
}
