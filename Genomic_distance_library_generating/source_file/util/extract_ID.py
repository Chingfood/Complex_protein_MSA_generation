import sys

infile, outfile = sys.argv[1:]

with open(infile, 'r') as f: IDs = [_.split('|')[1] for _ in f if len(_)>1 and _[0] == '>']
# with open(infile, 'r') as f: d = [_.split('|')[1] for _ in f if len(_)>1 and _.lstrip('\x00')[0] == '#']
print('#entries = %d' % len(IDs))

IDs = sorted(set(IDs))
print('#entries = %d' % len(IDs))

with open(outfile, 'w') as f: f.write('\n'.join(IDs) + '\n')
