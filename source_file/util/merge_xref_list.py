import os, sys
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

infolder, outfile = sys.argv[1:]

def func(infile):
	if infile.endswith('.log'): return
	# if infile in ['.', '..']: continue
	# print infile

	start, end = map(int, infile.split('_', 2)[:2])
	assert start < end

	source = infile.split('_', 2)[2]
	infile = infolder + '/' + infile
	with open(infile, 'r') as f: c = f.read().strip()
	if c: line = [source + '\t' + _ for _ in c.split('\n')]
	else: line = []

	return line, start, end
pool = ThreadPool(8)
ret = pool.map(func, [_ for _ in os.listdir(infolder) if not _.endswith('.log') and not _.startswith('_')])
pool.close()
pool.join()
print 'merging lines'
lines = [__ for _ in ret for __ in _[0]]
print 'merging file info'
file_start, file_end = [[_[__] for _ in ret] for __ in [1, 2]]

print 'sorting'
lines = sorted(lines)

# print 'assert nored line'
# assert all(_ != __ for _, __ in zip(lines[:-1], lines[1:]))

# print 'assert unique protein_id'
# protein_ids = sorted([_.split()[3] for _ in lines])
# for _, __ in zip(protein_ids[:-1], protein_ids[1:]):
# 	if _ == __:
# 		print _
# assert all(_ != __ for _, __ in zip(protein_ids[:-1], protein_ids[1:]))

print 'writing to file'
with open(outfile, 'w') as f:
	for line in lines:
		f.write(line + '\n')

print 'coverage'
file_start, file_end = map(set, [file_start, file_end])
intersection = file_start & file_end
file_start -= intersection
file_end -= intersection
file_start, file_end = map(sorted, map(list, [file_start, file_end]))
for b in zip(file_start, file_end):
	print '\t'.join(map(str, b))