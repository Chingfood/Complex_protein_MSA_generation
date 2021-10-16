#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys

infile, outfile = sys.argv[1:]


# infile = "/home/chingyuenliu/complex_contact/test.ffdata"

# In[1]:


import re


# In[2]:


"""
Find the entry names.
for example:
input: #UniRef100_A0A5A9P0L4 Peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE\n
output: A0A5A9P0L4
Mechanism: Find the word after UniRef100
"""
pattern = re.compile(r'(?:UniRef\d+)_([a-zA-Z0-9]+)')


# In[4]:


with open(infile, 'r') as f:
    IDs = [pattern.findall(_) for _ in f]
IDs = sorted(set([item for sublist in IDs for item in sublist]))


# In[6]:


print('#entries = %d' % len(IDs))


# In[ ]:


with open(outfile, 'w') as f: f.write('\n'.join(IDs) + '\n')

