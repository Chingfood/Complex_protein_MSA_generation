import requests, sys, time, subprocess, StringIO, gzip, re
from enaBrowserTools.python.utils import is_wgs_set, is_sequence

# python download_contig.py contig_list contig_CDS 0 10000

infile, outfile = sys.argv[1:3]
lindex, rindex = map(int, sys.argv[3:5])
print 'loading queries'
with open(infile, 'r') as f: query_list = f.read().strip().split('\n', rindex)[:rindex]
query_list = query_list[lindex:]

with requests.Session() as s, open(outfile, 'w') as f:
	start_time = time.time()
	print 'begin'
	sys.stdout.flush()

	# for i, (data_class, ID) in enumerate(query_list):
	for i, ID in enumerate(query_list):
		print '%d\t%s' % (lindex+i, ID)

		if is_wgs_set(ID):
			for c in ['wgs/public', 'wgs/suppressed', 'tsa/public', 'tsa/suppressed']:
				r = s.get(
					'http://ftp.ebi.ac.uk/pub/databases/ena/' + c + '/' + ID[:2].lower() + '/' + ID + '.dat.gz',
					stream=False,
				).content
				if not re.findall('<title>404 Not Found</title>', r): break
			assert not re.findall('<title>404 Not Found</title>', r)
		elif is_sequence(ID):
			r = s.get(
				'http://www.ebi.ac.uk/ena/data/view/' + ID + '&display=text&download=gzip',
				stream=False,
			).content
		else: sys.exit()
		# sys.exit()
		assert r.strip() != '<h2>Too Many Requests</h2>'
		compressedFile = StringIO.StringIO()
		compressedFile.write(r)
		compressedFile.seek(0)
		decompressedFile = gzip.GzipFile(fileobj=compressedFile, mode='rb')
		r = decompressedFile.read()

		key = ''
		qualifier = ''
		for line in r.strip().split('\n'):
			if len(line) < 2: print 'ERROR short line\t' + ID + '-' + line + '-'
			if line[:2] == 'ID': f.write(line + '\n')
			if line[:2] != 'FT': continue
			if line[5] != ' ' and line[21] == '/': print 'ERROR key and qualifier\t' + ID
			if line[5] != ' ': key = line[5:20].rstrip(); qualifier = ''
			if line[21] == '/': qualifier = line[21:]
			if key == 'CDS' and (
					not qualifier or
					qualifier.startswith('/db_xref="UniProtKB/') or
					qualifier.startswith('/protein_id="')
			): f.write(line + '\n')

		i += 1
		if i % (10 if is_wgs_set(ID) else 100) == 0:
			duration = time.time() - start_time
			avg_time = duration / float(i)
			left_time = avg_time * (len(query_list) - i)
			print "%d\t%.4f\t%.4f\t%.4f" % (i, duration, avg_time, left_time)
			sys.stdout.flush()

print '%f' % (time.time() - start_time)
