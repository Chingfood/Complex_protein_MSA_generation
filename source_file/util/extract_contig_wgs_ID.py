
## This script seperate wgs ID and contig ID using enaBrowserTools.
## wgs ID only counts up to 6th digit because: with the first 4 letters representing a project code, the first two numbers representing the assembly version, and the last 6 numbers providing unique identifiers for each contig.

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
		else: 
			print(ID)
			continue #sys.exit()

print("length of contigs before redundancy removed\n")
print(len(contigs))
print("length of wgss before redundancy removed\n")
print(len(wgss))

contigs = sorted(list(set(contigs)))
wgss = sorted(list(set(wgss)))

print("length of contigs\n")
print(len(contigs))
print("length of wgss\n")
print(len(wgss))


with open(outfile_contig, 'w') as f: f.write('\n'.join(contigs))
with open(outfile_wgs, 'w') as f: f.write('\n'.join(wgss))
