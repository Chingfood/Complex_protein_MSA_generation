import sys, os, re, time
import numpy as np
# from multiprocessing.pool import ThreadPool

listfile, infolder, outfolder = sys.argv[1:]

start = time.time()
with open(listfile, 'r') as f: c = f.read().strip().split()[:10]
print '%.6f' % (time.time() - start)
uniprot, protein_id = [c[_::2] for _ in [0, 1]]
print '%.6f' % (time.time() - start)
mapping = dict(zip(protein_id, [[] for _ in xrange(len(protein_id))]))
print '%.6f' % (time.time() - start)
for k, v in zip(protein_id, uniprot): mapping[k].append(v)
print '%.6f' % (time.time() - start)

re_pattern_complement = re.compile('complement')
re_pattern_join = re.compile('join')
re_pattern_order = re.compile('order')
re_pattern_number = re.compile('\d+')
re_pattern_segment = re.compile('\d+\.\.\d+')
re_pattern_singlesite = re.compile('(?:[^\d\.]|^)(\d+)(?=[^\d\.]|$)')

infiles = list([_ for _ in os.listdir(infolder) if not _.endswith('.log')])

def checkLoc(loc):
	"""
	This function uses a very lazy criteria
	In the future, I'd NOT like to update it, because it needs to PARSE the string
	:param loc:
	:return:
	"""

	global checkLocFlag

	pos_complement = list(re_pattern_complement.finditer(loc))
	pos_join = list(re_pattern_join.finditer(loc))
	pos_order = list(re_pattern_order.finditer(loc))
	# loc = loc.replace('>', '').replace('<', '')
	# content_number = np.array(map(int, re_pattern_number.findall(loc)))

	errList = [
		pos_order,
		len(pos_complement) > 1,
		len(pos_join) > 1,
		loc.count(':'),
		pos_join and pos_complement and pos_join[0].start() < pos_complement[0].start(),
		# len(content_number) % 2 == 0,
		# any(content_number[1:] > content_number[:-1]),
		# any(content_number[1:] == content_number[:-1]),
		# any(content_number[::2] > content_number[1::2]),
		# any(content_number[::2] == content_number[1::2]),
	]
	errbit = '1'
	passbit = ' '
	checkLocFlag = ''.join([errbit if _ else passbit for _ in errList])

	return not checkLocFlag.count(errbit)


def parseLoc(loc):
	if re_pattern_complement.search(loc):
		strand = '-'
	else:
		strand = '+'

	loc = loc.replace('>', '').replace('<', '')
	pos = map(int, re_pattern_number.findall(loc))
	segments = re_pattern_segment.findall(loc)
	singlesites = map(int, re_pattern_singlesite.findall(loc))
	flag = len(pos) == len(segments) * 2 + len(singlesites)
	if not flag: print >> sys.stderr, loc
	assert flag

	length = len(singlesites) + len(segments)
	segments = [__ for _ in segments for __ in _.split('..')]
	segments = map(int, segments)
	length += sum(e - s for s, e in zip(segments[::2], segments[1::2]))

	start = min(pos)
	end = max(pos)

	return map(str, [strand, start, end, length])

# maxloc = 0
for nf, infile in enumerate(infiles[120:]):
# def func(argss):
	print '%d\t%s' % (nf, infile)

	contig = ''
	loc = ''
	db_xref_uniprot = []
	protein_id = []

	with open(infolder + '/' + infile, 'r') as f: content = f.read().strip().split('\n')
	i = 0
	def loadContig():
		global i, contig
		if i < len(content) and content[i].startswith('ID   '): contig = content[i][5:].split(';')[0]; i += 1
		else: contig = ''
		return contig
	def loadLoc():
		global i, loc
		loc = ''
		if i >= len(content): return loc
		if not content[i].startswith('FT   CDS             '): return loc
		else: loc = content[i][21:]; i = i+1
		while i < len(content) and content[i].startswith('FT                   ') and content[i][21] != '/': loc += content[i][21:]; i += 1
		return loc
	def loadQualifier():
		global i, db_xref_uniprot, protein_id
		db_xref_uniprot = []
		protein_id = []
		while i < len(content):
			c = content[i]
			if c.startswith('FT                   /db_xref="UniProtKB/'):
				db_xref_uniprot.append(c[c.find(':')+1:c.rfind('"')])
			elif c.startswith('FT                   /protein_id="'):
				c = c[c.find('"')+1:-1]
				ic = c.find('.')
				if ic: c = c[:ic]
				protein_id.append(c)
			else: break
			i += 1
		return db_xref_uniprot, protein_id

	with open(outfolder + '/' + infile, 'w') as f:
		while i < len(content):
			if not loadContig(): print >> sys.stderr, 'ERROR null contig\t%s\t%d' % (infile, i); print content[i]; i += 1; sys.exit(); continue
			while loadLoc():
				# if len(loc) > maxloc: maxloc = len(loc); print '%d\t%s' % (maxloc, loc)
				if not loadQualifier()[1]: print >> sys.stderr, 'ERROR no protein_id\t%s\t%d' % (infile, i); continue
				if len(db_xref_uniprot) > 1: print >> sys.stderr, 'WARNING multiple xref\t%s\t%d\t%d' % (infile, i, len(db_xref_uniprot)); continue
				if len(protein_id) > 1: print >> sys.stderr, 'WARNING multiple protein_id\t%s\t%d\t%d' % (infile, i, len(protein_id)); continue
				if not checkLoc(loc): print >> sys.stderr, 'WARNING ill loc\t%s\t%d\t%s\t%s' % (infile, i, checkLocFlag, loc); continue
				strand, start, end, length = parseLoc(loc)
				for pro_id in protein_id:
					f.write('\n'.join('\t'.join([_, contig + strand, start, end, length]) for _ in mapping.get(pro_id, [])) + '\n')
# print maxloc


print 'coverage'
for prefix in ['contig', 'wgs']:
	print prefix
	file_start, file_end = [], []
	for infile in infiles:
		t = map(int, infile.split('_'))
		if t[0] != prefix: continue
		file_start.append(t[1])
		file_end.append(t[2])
	file_start, file_end = map(set, [file_start, file_end])
	intersection = file_start & file_end
	file_start -= intersection
	file_end -= intersection
	file_start, file_end = map(sorted, map(list, [file_start, file_end]))
	for b in zip(file_start, file_end):
		print '\t'.join(map(str, b))