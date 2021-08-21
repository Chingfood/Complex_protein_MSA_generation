import os, sys, time
import numpy as np
from numpy.lib import recfunctions as rfn
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool

infolder, outfile = sys.argv[1:]
infiles = os.listdir(infolder)

"""
def f(infile):
	# uniprot, contig, start, end = [], [], [], []
	ret = [[], [], [], []]
	with open(infile, 'r') as f:
		for line in f:
			for arr, e in zip(ret, line.rstrip().split()):
				arr.append(e)
	return ret
pool = ThreadPool(8)
results = pool.map(f, infiles)
pool.close()
pool.join()
uniprot, contig, start, end = [
	[_[index] for result in results for _ in result]
	for index in range(4)
]
"""

content_dtype = {'names': ('uniprot', 'contig', 'start', 'end', 'length'), 'formats': ('S10', 'S20', 'i4', 'i4', 'i4')}
def func(argss):
# for nf, infile in enumerate(infiles):
	nf, infile = argss
	print('%d\t%s' % (nf, infile))
	return np.loadtxt(infolder + '/' + infile, delimiter='\t', dtype=content_dtype, ndmin=1)
# pool = ThreadPool(4)
pool = Pool(10)
content = pool.map(func, enumerate(infiles))
pool.close()
pool.join()
content = [_ if _.size==1 else np.asarray(_) for _ in content if _.size]
content = np.concatenate(content)

print(len(content))
content.sort(order=['contig', 'start', 'end'])
print('sorted')

# overlap = [float(min(a['end'], b['end'])-b['start']+1)/min(a['length'], b['length']) for a, b in zip(content[:-1], content[1:]) if a['contig'] == b['contig'] and a['end'] >= b['start']]
# print len(overlap)
# with open('overlap.log', 'w') as f: f.write('\n'.join([str(_) for _ in overlap]) + '\n')
# exit()

cnt = 0
last_contig = ''
content = rfn.append_fields(content, names='rank', dtypes='i4', data=[-1]*len(content), usemask=False)

start_time = time.time()
for i, c in enumerate(content):
	if c['contig'] != last_contig: cnt, last_contig = 0, c['contig']
	for b in content[i-1::-1] if i else []:
		# print b
		if last_contig != b['contig']: break
		if b['rank'] != cnt: break
		if b['end'] >= c['start'] and b['end']-c['start']+1 >= min(b['length'], c['length']): cnt -= 1; break
	cnt += 1
	c['rank'] = cnt

	i += 1
	if i % 1000000 == 0:
		duration = time.time() - start_time
		avg_time = duration / float(i)
		left_time = avg_time * (len(content) - i)
		print("%d\t%.4f\t%.4f\t%.4f" % (i, duration, avg_time, left_time))
		sys.stdout.flush()

content.sort(order=['uniprot', 'contig', 'rank'])

with open(outfile, 'w') as f:
	f.write('\n'.join(['\t'.join([c['uniprot'].decode(), c['contig'].decode(), str(c['rank'])]) for c in content]) + '\n')