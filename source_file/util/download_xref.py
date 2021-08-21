import requests, sys, time

infile, source, outfile = sys.argv[1:4]
lindex, rindex = map(int, sys.argv[4:6])
#assert source in ['UniProtKB/TrEMBL', 'UniProtKB/Swiss-Prot']

# sys.exit()

print("loading queries")
with open(infile, 'r') as f: query_list = f.read().strip().split('\n', rindex)[:rindex]
query_list = query_list[lindex:]

length = len(query_list)
print(f"length of the list is {length}")
step=200
print('begin')

start_time = time.time()

#in this way, the session won't be too long to be cancelled out.


for st in range(0,length,step):
	ed = min(st+step,length)

	with requests.Session() as s, open(outfile, 'a') as f:
	# with requests.Session() as s:
	
		sys.stdout.flush()

		for i, ID in enumerate(query_list[st:ed]):

			source='UniProtKB/TrEMBL'
			try:
				r = s.get(
				'https://www.ebi.ac.uk/ena/xref/rest/tsv/search?',
				params={
					'source':				source,
					'source_accession':		ID,
					'target':				'coding',
					'download':				'txt',
					},
					stream=False,
				)
			except:
				continue
			result = r.content
			if len(result.strip().split(b'\n'))>1:
				# print result	

				header, elements = result.strip().split(b'\n',1)
				if header == b'Source\tSource primary accession\tSource secondary accession\tSource url\tTarget\tTarget primary accession\tTarget secondary accession\tTarget url':	

					for element in elements.split(b'\n'):
						element = element.split()
						assert any(_ == source.encode() for _ in element[:6])
						assert any(_ == ID.encode() for _ in element[1:6])
						assert any(_ == 'coding'.encode() for _ in element[3:6])
						element[0] = b'UniProtKB_TrEMBL'
						result = '\t'.encode().join([element[k] for k in [0, 1, 2, 5, 6]]).decode()
						# print result
						f.write(result + '\n')
				else:
					print("header is wrong!\n")
					print(result)
					
			else:
				source='UniProtKB/Swiss-Prot'
				try:
					r = s.get(
						'https://www.ebi.ac.uk/ena/xref/rest/tsv/search?',
						params={
						'source':				source,
						'source_accession':		ID,
						'target':				'coding',
						'download':				'txt',
						},
						stream=False,
					)
				except:
					continue	
				result = r.content
				if len(result.strip().split(b'\n'))>1:
					# print result
					header, elements = result.strip().split(b'\n',1)
					if header == b'Source\tSource primary accession\tSource secondary accession\tSource url\tTarget\tTarget primary accession\tTarget secondary accession\tTarget url':
						for element in elements.split(b'\n'):
							element = element.split()
							assert any(_ == source.encode() for _ in element[:6])
							assert any(_ == ID.encode() for _ in element[1:6])
							assert any(_ == 'coding'.encode() for _ in element[3:6])
							element[0] = b"UniProtKB_Swiss-Prot"
							result = '\t'.encode().join([element[k] for k in [0, 1, 2, 5, 6]]).decode()
							# print result
							f.write(result + '\n')	

					else:
						print("header is wrong!\n")
						print(result)
						

			
			if (st+i) % 1000 == 0 and (st+i) != 0:
				duration = time.time() - start_time
				avg_time = duration / float(st+i)
				left_time = avg_time * (len(query_list)-st-i)
				print("%d\t%.4f\t%.4f\t%.4f" % ((st+i), duration, avg_time, left_time))
				time.sleep(2.3)
				sys.stdout.flush()


print("%f" % (time.time() - start_time))
print("done")
