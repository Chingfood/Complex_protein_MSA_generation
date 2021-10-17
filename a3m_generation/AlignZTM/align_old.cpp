#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sys/stat.h>
#include <time.h>
#include <assert.h>
#include <algorithm>
#include <utility>

//https://groups.google.com/forum/#!topic/comp.std.c/d-6Mj5Lko_s
#define PP_NARG(...) \
         PP_NARG_(__VA_ARGS__,PP_RSEQ_N())
#define PP_NARG_(...) \
         PP_ARG_N(__VA_ARGS__)
#define PP_ARG_N( \
          _1, _2, _3, _4, _5, _6, _7, _8, _9,_10, \
         _11,_12,_13,_14,_15,_16,_17,_18,_19,_20, \
         _21,_22,_23,_24,_25,_26,_27,_28,_29,_30, \
         _31,_32,_33,_34,_35,_36,_37,_38,_39,_40, \
         _41,_42,_43,_44,_45,_46,_47,_48,_49,_50, \
         _51,_52,_53,_54,_55,_56,_57,_58,_59,_60, \
         _61,_62,_63,N,...) N
#define PP_RSEQ_N() \
         63,62,61,60,                   \
         59,58,57,56,55,54,53,52,51,50, \
         49,48,47,46,45,44,43,42,41,40, \
         39,38,37,36,35,34,33,32,31,30, \
         29,28,27,26,25,24,23,22,21,20, \
         19,18,17,16,15,14,13,12,11,10, \
         9,8,7,6,5,4,3,2,1,0
#define _VFUNC_(name, n) name##n
#define _VFUNC(name, n) _VFUNC_(name, n)
#define VFUNC(func, ...) _VFUNC(func, PP_NARG(__VA_ARGS__)) (__VA_ARGS__)
template<typename T> struct argument_type;
template<typename T, typename U> struct argument_type<T(U)> { typedef U type; };
//#define STRV(...) __VA_ARGS__



#define UNIPROT_ID_LEN 16
#define NSEQ 262140
struct MSA {
	char uniprot_id[NSEQ][UNIPROT_ID_LEN];
	char *seq[NSEQ];
	int nseq;
	int contig_index[NSEQ+1];
	int seql;
} msa1, msa2;


#define CONTIG_ID_LEN 32
#define NCONTIGLOC (NSEQ<<3)
struct Contig {
	char id[CONTIG_ID_LEN];
	int loc, msa_ind;
};
struct Contigs {
	Contig contigs[NCONTIGLOC];
	Contig *ptr_contigs[NCONTIGLOC];
	int n;
} contig1, contig2;


void loadMSA (char *infilename, MSA &msa) {
//#define BUFFERSIZE 8192
#define BUFFERSIZE 16000
	char AA[26] = {'A','X','C','D','E','F','G','H','I','X','K','L','M','N','O','P','Q','R','S','T','X','V','W','X','Y','X'};
	char buffer[BUFFERSIZE]; msa.seql = -1;
	FILE *fin = fopen (infilename, "r"); assert (fin);
	fgets (buffer, BUFFERSIZE, fin); for (char *d = msa.uniprot_id[msa.nseq = 0], *s = buffer+1; (*s == '\n'/*' '*/? *d = 0: *(d++) = *(s++));) ;
	printf ("%s\n", msa.uniprot_id[0]);
	msa.seql = strlen (fgets (buffer, BUFFERSIZE, fin))-1; memcpy (msa.seq[0] = new char [msa.seql+1], buffer, sizeof(char)*(msa.seql)); msa.seq[0][msa.seql] = 0;
	for (msa.nseq = 1; fgets (buffer, BUFFERSIZE, fin); msa.nseq++) {
		printf("buffer is %s\n",buffer);	
		int lv, rv; for (lv = 0; lv < BUFFERSIZE && buffer[lv] != '|'; lv++) ; assert (lv < BUFFERSIZE); for (rv = ++lv; rv < BUFFERSIZE && buffer[rv] != '|'; rv++) ; assert (rv < BUFFERSIZE); buffer[rv++] = 0;
		assert (msa.nseq < NSEQ);
		memcpy (msa.uniprot_id[msa.nseq], buffer+lv, sizeof(char)*(rv-lv));
//		fread (msa.seq[msa.nseq] = new char [msa.seql+1], sizeof(char), msa.seql, fin); msa.seq[msa.nseq][msa.seql] = 0; assert (fgetc (fin) == '\n');
		assert (fgets (buffer, BUFFERSIZE, fin)); char *s = buffer, *t = msa.seq[msa.nseq] = new char[msa.seql+1], c, *tt = t+msa.seql;
		for (; (c=*s) != '\n' && t <= tt; s++) if ('A' <= c && c <= 'Z') *(t++) = AA[c-'A']; else if (c == '-') *(t++) = '-'; else if ('a' <= c && c <= 'z') continue; else assert (false);
		assert (t == tt); assert (c == '\n'); *t = 0; for (t = msa.seq[msa.nseq]; t < tt; assert (*(t++))) ;
	}
	fclose (fin);
// to skip one line
//		if (strlen (buffer) == BUFFERSIZE && buffer[BUFFERSIZE] != '\n') while (fgets (buffer, BUFFERSIZE, fin) && strlen(buffer) == BUFFERSIZE && buffer[BUFFERSIZE-1] != '\n') ;
//		while (fgets (buffer, BUFFERSIZE, fin) && strlen(buffer) == BUFFERSIZE && buffer[BUFFERSIZE-1] != '\n') ;
#undef BUFFERSIZE
}

bool Cmp_ContigPtr (const Contig *i, const Contig *j) {
	assert (i && j);
	return strcmp (i->id, j->id) < 0;
}

void loadContig (MSA &msa, const char *databasefilename, Contigs &contig) {
#define BUFFERSIZE 32
	char buffer[BUFFERSIZE];
	FILE *fin = fopen (databasefilename, "r"); assert (fin);
	struct stat st; fstat(fileno(fin), &st); unsigned long long int filesize = st.st_size;
	printf ("file size = %lld\n", filesize);
	unsigned long long int l, r, m, s;
	contig.n = 0;
	for (int i = 1; i < msa.nseq; i++) {
		l = 0; r = filesize; msa.contig_index[i] = contig.n;
		while (l < r) {
			for (fseek (fin, s = std::max ((m=(l+r>>1))-BUFFERSIZE, l), SEEK_SET); fgets (buffer, BUFFERSIZE, fin) && ftell (fin) <= m; s = ftell (fin)) ;
			for (int _ = 0; (buffer[_]=='\t'? buffer[_]=0: true); _++) ;
			if (strcmp (buffer, msa.uniprot_id[i]) >= 0) r = s; else l = ftell (fin);
		}
		for (fseek (fin, l, SEEK_SET); fgets (buffer, BUFFERSIZE, fin); contig.n++) {
			for (l = 0; buffer[l] != '\t'; l++) ; buffer[l++] = 0; for (r = l; buffer[r] != '\t'; r++) ; buffer[r++] = 0;
			if (strcmp (buffer, msa.uniprot_id[i])) break;
			assert (contig.n < NCONTIGLOC);
			memcpy (contig.contigs[contig.n].id, buffer+l, sizeof(char)*(r-l));
			contig.contigs[contig.n].loc = atoi(buffer+r);
//			int &t = contig.contigs[contig.n].loc; t = 0; for (char *p = buffer+r, c; (c=*p) != '\n'; p++) {t *= 62; if (c <= '9') t += c-'0'; else if (c <= 'Z') t += c-'A'+36; else t += c-'a'+10;}
			contig.contigs[contig.n].msa_ind = i;
		}
	}
	fclose (fin);
#undef BUFFERSIZE
	msa.contig_index[msa.nseq] = contig.n;
	for (int i = 0; i < contig.n; i++) contig.ptr_contigs[i] = contig.contigs+i;
	clock_t _beg = clock ();
	std::sort (contig.ptr_contigs, contig.ptr_contigs + contig.n, Cmp_ContigPtr);
	printf ("sorted %d\n", clock () - _beg);
}


namespace joint_Hungarian {
#define NVERTEX NCONTIGLOC
#define NEDGE ((NCONTIGLOC)<<3)
struct Edge {
	int v;
	int w;
	Edge *next;
} _edge[NEDGE], *e_ind, *e_head[NVERTEX];
int matching[NVERTEX], mark[NVERTEX];
#undef NEDGE
#undef NVERTEX

void init () {e_ind = _edge; memset (matching, -1, sizeof(matching)); memset (mark, -1, sizeof(mark)); memset (e_head, 0, sizeof(e_head));}

void genGraph (const Contigs &c1, const Contigs &c2, const int min_dgene, const int max_dgene) {
	int cnt = 0;
	for (int i = 0, j = 0, t; i < c1.n && j < c2.n; i++) {
		char *iid = c1.ptr_contigs[i]->id;
		while (j < c2.n && (t = strcmp (iid, c2.ptr_contigs[j]->id)) > 0) j++; if (t) continue;
		int ii = c1.ptr_contigs[i]-c1.contigs;
		for (int k = j; k < c2.n && !strcmp (iid, c2.ptr_contigs[k]->id); k++) if (min_dgene <= (t = abs (c1.ptr_contigs[i]->loc-c2.ptr_contigs[k]->loc)) && t <= max_dgene)
		{Edge *x = e_ind++; x->v = c2.ptr_contigs[k]-c2.contigs; x->w = t; x->next = e_head[ii]; e_head[ii] = x; cnt++;}
	}
//	for (int i = 0, t; i < c1.n; i++) for (int j = 0; j < c2.n; j++) if (!strcmp (c1.contigs[i].id, c2.contigs[j].id) && min_dgene <= (t = abs (c1.contigs[i].loc-c2.contigs[j].loc)) && t <= max_dgene)
//	{Edge *x = e_ind++; x->v = j; x->w = t; x->next = e_head[i]; e_head[i] = x; cnt++;}
	printf ("#edges = %d\n", cnt);
}

int _hungarian_cur;

bool dfs (int i) {
	for (Edge *x = e_head[i]; x; x = x->next) if (mark[x->v] != _hungarian_cur) {
		mark[x->v] = _hungarian_cur;
		assert (!strcmp (contig1.contigs[i].id, contig2.contigs[x->v].id));
		if (matching[x->v] == -1 || dfs (matching[x->v])) {matching[x->v] = i; return true;}
	}
	return false;
}

void match (int n1, int n2) {
	int cnt = 0;
	for (int i = 0; i < n1; i++) cnt += dfs(_hungarian_cur = i);
	printf ("#match = %d\n", cnt);
}

void output (MSA &m1, Contigs &c1, MSA &m2, Contigs &c2, FILE *fout) {
	int cnt = 0;
	for (int j = 0; j < c2.n; j++) if (matching[j] != -1) {
		Contig *a = c1.contigs+matching[j], *b = c2.contigs+j;
		int x = a->msa_ind, y = b->msa_ind;
		fprintf (fout, ">%s_%s|%s %d %s %d\n%s%s\n",
				m1.uniprot_id[x], m2.uniprot_id[y],
				a->id, a->loc, b->id, b->loc,
				m1.seq[x], m2.seq[y]
				);
		cnt++;
	}
	printf ("#joint alignments = %d\n", cnt);
}

void run (MSA &m1, Contigs &c1, MSA &m2, Contigs &c2, int min_dgene, int max_dgene, FILE *fout) {
	printf ("Hungarian\n");
	init ();
	genGraph (c1, c2, min_dgene, max_dgene);
	match (c1.n, c2.n);
	output (m1, c1, m2, c2, fout);
}
}


namespace joint_Closest {

int *dist1, *dist2;
int *target1, *target2;

void run (MSA &m1, Contigs &c1, MSA &m2, Contigs &c2, int min_dgene, int max_dgene, FILE *fout) {
	dist1 = new int [m1.nseq]; dist2 = new int [m2.nseq];
	target1 = new int [m1.nseq]; target2 = new int [m2.nseq];
	memset (dist1, 127, sizeof(int)*m1.nseq); memset (dist2, 127, sizeof(int)*m2.nseq);
	memset (target1, -1, sizeof(int)*m1.nseq); memset (target2, -1, sizeof(int)*m2.nseq);
	int cnt = 0;
	for (int li = 0, ri, lj = 0, rj, t; li < c1.n && lj < c2.n; li = ri, lj = rj) {
		while (li < c1.n && lj < c2.n && (t = strcmp (c1.ptr_contigs[li]->id, c2.ptr_contigs[lj]->id))) if (t < 0) li++; else lj++;
		if (li >= c1.n || lj >= c2.n) break;
		for (ri = li+1; ri < c1.n && !strcmp (c1.ptr_contigs[li]->id, c1.ptr_contigs[ri]->id); ri++) ;
		for (rj = lj+1; rj < c2.n && !strcmp (c2.ptr_contigs[lj]->id, c2.ptr_contigs[rj]->id); rj++) ;
		for (int i = li; i < ri; i++) for (int j = lj; j < rj; j++) if ((t = abs(c1.ptr_contigs[i]->loc-c2.ptr_contigs[j]->loc)) <= max_dgene && min_dgene <= t) {
			assert (!strcmp(c1.ptr_contigs[i]->id, c2.ptr_contigs[j]->id));
			int x = c1.ptr_contigs[i]->msa_ind, y = c2.ptr_contigs[j]->msa_ind;
			if (t < dist1[x]) dist1[x] = t, target1[x] = y;
			if (t < dist2[y]) dist2[y] = t, target2[y] = x;
		}
	}
	for (int i = 0, j; i < m1.nseq; i++) if ((j = target1[i]) != -1 && target2[j] == i) {
		fprintf (fout, ">%s_%s_%d\n%s%s\n",
			m1.uniprot_id[i], m2.uniprot_id[j], dist1[i],
			m1.seq[i], m2.seq[j]
			);
		cnt++;
	}
	printf ("%d\n", cnt);
}

}

namespace joint_KM {
struct Edge {
	int v, w;
	Edge *next;
} *_edge, **e_head, *e_idle, *e_end;

std::pair<int, int> *heap; int hsz, *ih;

int *pl, *pr;
int *ml, *mr;
int *vl, *vr;
int *ql, qlf, qlr;
int *qr, qrr, *from;

inline void init (int nl, int nr) {
#define init(...) VFUNC(init, __VA_ARGS__)
#define init3(a,t,s) assert (a = new argument_type<void(t)>::type [s]);
//#define init3(a,t,s) assert (a = new t [s]);
#define init4(a,t,s,v) init3(a,t,s) memset (a, v, sizeof(t)*(s));
	init (e_idle = _edge, Edge, nl+nr<<3)	init (e_head, Edge*, nl, 0)	e_end = _edge+(nl+nr<<3);
	init (heap, (std::pair<int, int>), nr)	init (ih, int, nr, -1)	hsz = 0;
//	init (heap, STRV (std::pair<int, int>), nr)	init (ih, int, nr, -1)	hsz = 0;
	init (pl, int, nl, 127)	init (pr, int, nr, 0)
	init (ml, int, nl, -1)	init (mr, int, nr, -1)
	init (vl, int, nl, -1)	init (vr, int, nr, -1)
	init (ql, int, nl)	init (qr, int, nr)	init (from, int, nr)	qlf = qlr = qrr = 0;
#undef init
#undef init3
#undef init4
}

inline void destructor () {
	delete [] _edge, e_head;
	delete [] heap, ih;
	delete [] pl, pr, ml, mr, vl, vr;
	delete [] ql, qr, from;
}

inline void genGraph (const Contigs &c1, const Contigs &c2, const int min_dgene, const int max_dgene) {
	int cnt = 0;
	for (int i = 0, j = 0, t; i < c1.n && j < c2.n; i++) {
		char *iid = c1.ptr_contigs[i]->id;
		while (j < c2.n && (t = strcmp (iid, c2.ptr_contigs[j]->id)) > 0) j++; if (t) continue;
		int ii = c1.ptr_contigs[i]-c1.contigs;
		for (int k = j; k < c2.n && !strcmp (iid, c2.ptr_contigs[k]->id); k++) if (min_dgene <= (t = abs (c1.ptr_contigs[i]->loc-c2.ptr_contigs[k]->loc)) && t <= max_dgene)
		{Edge *x = e_idle++; assert (x < e_end); x->v = c2.ptr_contigs[k]-c2.contigs; x->w = t; x->next = e_head[ii]; e_head[ii] = x; cnt++; if (pl[ii] > t) pl[ii] = t;}
	}
	printf ("#edges = %d\n", cnt);
}

inline void match (int nl)
{
#define checkih {for (int i = 0; i < hsz; i++) assert (ih[heap[i].first] == i);}
#define checkheap {for (int i = hsz-1; i >= 0; i--) assert (heap[i].second >= heap[(i-1)>>1].second);}
#define check {checkih checkheap}
//#define check
#define rh(i,ii) heap[ih[ii.first] = i] = ii
#define	uh(i,x) for (int ii; i && heap[ii = (i-1)>>1].second > x.second; i = ii) rh (i, heap[ii]);
#define dh(i,x) for (int ii; (ii = (i<<1)+1) < hsz && heap[ii += (ii+1 < hsz && heap[ii].second > heap[ii+1].second)].second < x.second; i = ii) rh (i, heap[ii]);
#define push(x) {int i = hsz++; std::pair<int, int> _x = (x); uh (i, x) rh (i, _x); check}
#define pop {int i = 0; ih[heap[0].first] = -1; if (--hsz) {const std::pair<int, int> &x = heap[hsz]; dh (i, x) rh (i, x);} check}
#define update(i_) {int i = (i_); std::pair<int, int> x = heap[i]; uh (i, x) rh (i, x); check}

#define checkinvariant {printf ("%d\n", __LINE__); for (int k = 0; k < N; k++) for (Edge *x = e_head[k]; x; x = x->next) assert (pl[k] + pr[x->v] <= x->w);}
#define checkinvariant2 {printf ("%d\n", __LINE__); for (int k = 0; k < N; k++) for (Edge *x = e_head[k]; x; x = x->next) assert ((vl[k]==i? offset_price: 0) + pl[k] + pr[x->v] + (vr[x->v]==i? -offset_price: 0) <= x->w);}

	for (int i = 0, t, w; i < nl; i++) {
		vl[i] = ql[hsz = qlf = qrr = 0] = i; qlr = 1;
		int target = -1, offset_price = 0, offset_heap = 0;
		while (target == -1) {
			while (qlf < qlr && target == -1) {
				int cur = ql[qlf++];
				for (Edge *x = e_head[cur]; x; x = x->next) if (vr[x->v] < i) {
					if ((w = x->w - pl[cur] - pr[x->v] - offset_price) > 0) {
						w -= offset_heap;
						if ((t = ih[x->v]) == -1) {from[x->v] = cur; push (std::make_pair (x->v, w))}
						else if (heap[t].second > w) {from[x->v] = cur; heap[t].second = w; update (t)}
					}
					else if ((t = mr[x->v]) == -1) {from[target = x->v] = cur; break;}
					else if (vl[t] < i) {vl[ql[qlr++] = t] = vr[qr[qrr++] = x->v] = i; pl[t] -= offset_price; pr[x->v] += offset_price; from[x->v] = cur;}
				}
			}
			if (target != -1) break;
			while (hsz && vr[heap[0].first] == i) pop;
			if (!hsz) {int tl = i; for (int j = 1, t; j < qlr; j++) if (pl[t = ql[j]] > pl[tl]) tl = t; target = ml[tl]; ml[tl] = -1; break;}
			assert (heap[0].second + offset_heap > 0);
			offset_price += (w = heap[0].second) + offset_heap; offset_heap = -w;
			for (; hsz && heap[0].second == w; ) {
				int j = heap[0].first; pop
				if (vr[j] == i) continue;
				if (mr[j] == -1) {target = j; break;}
				vl[t = ql[qlr++] = mr[j]] = vr[qr[qrr++] = j] = i; pl[t] -= offset_price; pr[j] += offset_price;
			}
			if (target != -1) break;
		}
		if (target != -1) {for (int t; (t = mr[target] = from[target]) != i; ) target ^= ml[t], ml[t] ^= target, target ^= ml[t]; ml[i] = target;}
		for (std::pair<int, int> *h = heap+hsz-1; h >= heap; h--) ih[h->first] = -1; hsz = 0;
		for (int j = 0; j < qlr; j++) pl[ql[j]] += offset_price;
		for (int j = 0; j < qrr; j++) pr[qr[j]] -= offset_price;
	}
//	int cost = 0; for (int i = 0; i < N; i++) cost += pl[i] + pr[ml[i]];
#undef check
}

inline void check (int nl, int nr) {
	assert (!hsz); for (int i = 0; i < nr; i++) assert (ih[i] == -1);
	for (int i = 0; i < nl; i++) assert (ml[i] == -1 || mr[ml[i]] == i);
	for (int i = 0; i < nr; i++) assert (mr[i] == -1 || ml[mr[i]] == i);
	for (int i = 0; i < nl; i++) for (Edge *x = e_head[i]; x; x = x->next) if (!(ml[i] != x->v || pl[i] + pr[x->v] == x->w)) {printf ("%d\t%d\t%d\t%d\t%d\n", i, pl[i], x->v, pr[x->v], x->w);}
	for (int i = 0; i < nl; i++) for (Edge *x = e_head[i]; x; x = x->next) assert ((ml[i] != x->v || pl[i] + pr[x->v] == x->w) && pl[i] + pr[x->v] <= x->w);
}

inline void output (MSA &m1, Contigs &c1, MSA &m2, Contigs &c2, FILE *fout) {
	int cnt = 0;
	for (int j = 0; j < c2.n; j++) if (mr[j] != -1) {
		Contig *a = c1.contigs+mr[j], *b = c2.contigs+j;
		int x = a->msa_ind, y = b->msa_ind;
		fprintf (fout, ">%s_%s|%s %d %s %d\n%s%s\n",
				m1.uniprot_id[x], m2.uniprot_id[y],
				a->id, a->loc, b->id, b->loc,
				m1.seq[x], m2.seq[y]
				);
		cnt++;
	}
	printf ("#joint alignments = %d\n", cnt);
}

void run (MSA &m1, Contigs &c1, MSA &m2, Contigs &c2, int min_dgene, int max_dgene, FILE *fout) {
	init (c1.n, c2.n);
	genGraph (c1, c2, min_dgene, max_dgene);
	match (c1.n);
	check (c1.n, c2.n);
	output (m1, c1, m2, c2, fout);
	destructor ();
}
}



int main (int argc, char **argv) {
	printf ("align cpp\n");
	printf ("%d\n", argc);
	for (int i = 0; i < argc; i++) printf ("%s\n", argv[i]);

	char *inMSA1 = argv[1], *inMSA2 = argv[2];
	int min_dgene = atoi(argv[3]), max_dgene = atoi(argv[4]);
	char *database = argv[5];
	char *joint_method = argv[6];
	char *outfile = argv[7];


	clock_t _beg;

	_beg = clock ();
	loadMSA (inMSA1, msa1); printf ("loadMSA1 done %d\n", clock () - _beg);
	loadMSA (inMSA2, msa2); printf ("loadMSA2 done %d\n", clock () - _beg);
	loadContig (msa1, database, contig1); printf ("loadContig1 done %d\n", clock () - _beg);
	loadContig (msa2, database, contig2); printf ("loadContig2 done %d\n", clock () - _beg);
	printf ("%d\n", clock () - _beg);

	printf ("infile1: nseq = %d, nrecords = %d\n", msa1.nseq, contig1.n);
	printf ("infile2: nseq = %d, nrecords = %d\n", msa2.nseq, contig2.n);

	_beg = clock ();
	FILE *fout = fopen (outfile, "w");
//	FILE *fout = 0;
	fprintf (fout, ">%s_%s\n%s%s\n", msa1.uniprot_id[0], msa2.uniprot_id[0], msa1.seq[0], msa2.seq[0]);
#define is(x) if (!strcmp (joint_method, #x)) joint_##x::run(msa1, contig1, msa2, contig2, min_dgene, max_dgene, fout); else
	is (Closest)
	is (Hungarian)
	is (KM)
	printf ("unsupported joint method"), exit (-1);
	fclose (fout);
	printf ("%d\n", clock () - _beg);

	return 0;
}
