
import sys
from enaBrowserTools.python3.utils import is_wgs_sequence, is_sequence

infile, outfile_contig, outfile_wgs = sys.argv[1:]
contigs = []
wgss = []
with open(infile, 'r') as f:
	for ID in f:
		ID = ID[:-1]
		if is_sequence(ID): contigs.append(ID)
		elif is_wgs_sequence(ID): wgss.append(ID[:6])
		else: continue #sys.exit()

contigs = sorted(list(set(contigs)))
wgss = sorted(list(set(wgss)))

print("length of contigs\n")
print(len(contigs))
print("length of wgss\n")
print(len(wgss))


with open(outfile_contig, 'w') as f: f.write('\n'.join(contigs))
with open(outfile_wgs, 'w') as f: f.write('\n'.join(wgss))
