#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
using namespace std;

const int isIDLine0 = *(int*) "ID  ",
	isLocFirstLine0 = *(int*) "FT  ",
	isLocFirstLine1 = *(int*) " CDS",
	isLocFirstLine2 = *(int*) "    ",
	isLocFirstLine3 = *(int*) "    ",
	isLocFirstLine4 = *(int*) "    ",
	isLocRestLine0 = *(int*) "FT  ",
	isLocRestLine1 = *(int*) "    ",
	isLocRestLine2 = *(int*) "    ",
	isLocRestLine3 = *(int*) "    ",
	isLocRestLine4 = *(int*) "    ",
	isXrefLine0 = *(int*) "FT  ",
	isXrefLine1 = *(int*) "    ",
	isXrefLine2 = *(int*) "    ",
	isXrefLine3 = *(int*) "    ",
	isXrefLine4 = *(int*) "    ",
	isXrefLine5 = *(int*) " /db",
	isXrefLine6 = *(int*) "_xre",
	isXrefLine7 = *(int*) "f=\"U",
	isXrefLine8 = *(int*) "niPr",
	isXrefLine9 = *(int*) "otKB",
	keywordcomplement0 = *(int*) "comp",
	keywordcomplement1 = *(int*) "leme",
	keywordjoin0 = *(int*) "join";

void parseFile (char *infilename, char *outfilename) {
#define lineL 128
#define contigL 31
#define uniprotL 31
#define Nuniprot 32
#define locL 65535
	char line[lineL];
	char contig[contigL+1];
	char uniprot[Nuniprot][uniprotL+1]; int nuniprot;
	char loc[locL+1];

	char *l, *r;
	bool strand;	// 1 iff no complement iff +
	unsigned long long _s, _e, start, end, len;
	bool errflag;

	FILE *fin = fopen (infilename, "r");
	FILE *fout = fopen (outfilename, "w");
#define getLine fgets (line, lineL, fin)
#define strcmpint(str,offset,tarint) (*(((int*)(str))+offset) == tarint ## offset)
//#define isXXLine(i,XX) ((*(int*)(line+i)) == is ## XX ## Line ## i)
#define isXXLine(XX,i) strcmpint (line, i, is ## XX ## Line)
#define isIDLine (isXXLine (ID, 0)/* && line[4] == ' '*/)
#define loadContig if (!feof (fin) && isIDLine) {char *d = contig, *s = line+5; while (*s != ';') *(d++) = *(s++); *d = 0; getLine;} else *contig = 0;
#define isLocFirstLine (isXXLine (LocFirst, 0) && isXXLine (LocFirst, 1)/* && isXXLine (LocFirst, 2) && isXXLine (LocFirst, 3) && isXXLine (LocFirst, 4) && line[20] == ' '*/)
#define isLocRestLine (isXXLine (LocRest, 0) && isXXLine (LocRest, 1)/* && isXXLine (LocRest, 2) && isXXLine (LocRest, 3) && isXXLine (LocRest, 4) && line[20] == ' '*/ && line[21] != '/')
#define loadLoc if (!feof (fin) && isLocFirstLine) {char *d = loc, *s; for (s = line+21; /**s && */*s != '\n'; *(d++) = *(s++)) ; while (getLine && isLocRestLine) for (s = line+21; /**s && */*s != '\n'; *(d++) = *(s++)) ; *d = 0; r = d-1;} else *loc = 0;
#define isXrefLine (isXXLine (Xref, 0)/* && isXXLine (Xref, 1) && isXXLine (Xref, 2) && isXXLine (Xref, 3) && isXXLine (Xref, 4)*/ && isXXLine (Xref, 5)/* && isXXLine (Xref, 6) &&  isXXLine (Xref, 7) &&  isXXLine (Xref, 8) &&  isXXLine (Xref, 9) && line[40] == '/'*/)
#define loadXref nuniprot = 0; if (!feof (fin)) while (isXrefLine) {char *d = uniprot[nuniprot++], *s = line+41; while (*(s++) != ':') ; while (/**s && */*s != '\"') *(d++) = *(s++); *d = 0; if (!getLine) break;}

#define isKeyword(kw,offset) strcmpint (l, offset, keyword ## kw)
#define parseNumber(num) num = 0; if (*l == '<' || *l == '>') l++; while (l <= r && '0' <= *l && *l <= '9') num *= 10, num += *l-'0', l++; errflag |= !num;
#define parseLocDescriptor parseNumber (_s) if (start > _s) start = _s; if (l <= r && *l == '.' && l[1] == '.') {l += 2; parseNumber (_e) len += _e-_s+1; if (end < _e) end = _e; errflag |= end<start; if (end < start) printf ("ERROR ill span\t%s\t%d\t%d\n", contig, start, end), exit (-1);} else {len += 1; if (end < _s) end = _s;}
#define parseLocDescriptors if (l <= r) parseLocDescriptor while (l <= r && *l == ',' ) {l++; parseLocDescriptor}
#define parseComplement (isKeyword (complement, 0)/* && isKeyword (complement, 1) && l[8] == 'n' && l[9] == 't'*/ && l[10] == '(' && *r == ')')? l += 11, r--, strand = 0: strand = 1;
#define parseJoin if (isKeyword (join, 0) && l[4] == '(' && *r == ')') {l += 5, r--; parseLocDescriptors} else {parseLocDescriptor}
#define parseLoc parseComplement parseJoin errflag |= l<=r;

	getLine;
	int cnt = 0;
	while (!feof (fin)) {
		loadContig if (!*contig) printf ("ERROR null contig\t%s\n", line), exit (-1);
//		printf ("%s\n", contig);
		while (!feof (fin)) {
			loadLoc if (!*loc) break;
//			printf ("%s\n", loc);
			loadXref if (!nuniprot) {printf ("WARNING no xref\n"); continue;} else if (nuniprot > 1) {printf ("WARNING multiple xref\t%d\n", nuniprot); continue;}
//			for (int i = 0; i < nuniprot; i++) printf ("%s\n", uniprot[i]);
			l = loc; start = -1; end = 0; len = 0; errflag = 0;
			parseLoc
			if (!errflag) {fprintf (fout, "%s\t%s%c\t%llu\t%llu\t%llu\n", uniprot, contig, strand? '+': '-', start, end, len);}
			else printf ("ERROR ill loc\t%s\n", loc);
		}
//		if (cnt++ >= 20) break;
	}

	fclose (fin);
	fclose (fout);
}

int main () {
	char filename[64];
	char outfilename[8];
	char dirname[ ] = "CDS_list";
	DIR *dir;
	struct dirent *ent;
	dir = opendir (dirname);
	for (int i = 0; (ent = readdir (dir)) != NULL; i++) if (ent->d_name[0] != '.') {
		printf ("%s\n", ent->d_name);
		strcpy (filename, dirname);
		strcat (filename, "/");
		strcat (filename, ent->d_name);
		sprintf (outfilename, "loc_list_cpp/%d", i);
		parseFile (filename, outfilename);
	}
	closedir (dir);
	return 0;
}
