#!/usr/bin/env python

from a3m import A3M_Container
from a3m import A3MFormatError
import sys

def LoadFASTAFile(seqfile):
  with open(seqfile, 'r') as fh:
     lines = fh.readlines()
  #print(lines)
  if len(lines) != 2:
     sys.stderr.write("incorrect sequence file" + seqfile)
     exit(1)

  seqname = lines[0].split()[0]
  sequence = lines[1].strip()
  return seqname, sequence

def ReformatA3M(filename, seqname, sequence):
  a3m = A3M_Container()

  if(filename.lower() == "stdin"):
    fh = sys.stdin
  else:
    fh = open(filename, "r")

  try:
    a3m.read_a3m(fh)
  except A3MFormatError as e:
    sys.stderr.write(str(e))
    exit(1)

  if not a3m.sequences[0][1].startswith(sequence):
    sys.stderr.write("inconsistent sequence for " + seqname)
    #print(a3m.sequences[0][1])
    #print(sequence)
    exit(1)

  ##seqLen = len(sequence)
  
  print(seqname)
  print(sequence)
  for seq in a3m.sequences[1:]:
     print(seq[0])
     print(seq[1])


def main():
  #print(sys.argv)
  filename = sys.argv[1]
  seqfile = sys.argv[2]
  #print(seqfile)
  seqname, sequence = LoadFASTAFile(seqfile)
  #print(seqname)
  #print(sequence)
  ReformatA3M(filename, seqname, sequence)


if __name__ == "__main__":
  main()
