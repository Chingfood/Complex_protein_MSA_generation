#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys

infile, outfile = sys.argv[1:]


# infile = "/home/chingyuenliu/complex_contact/test.ffdata"

# outfile = "/home/chingyuenliu/complex_contact/tax_id_list"

# In[6]:


import re


# In[4]:


"""
Find the entry names.
for example:
input: #UniRef100_A0A5A9P0L4 Peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE\n
output: A0A5A9P0L4
Mechanism: Find the word after UniRef100
"""
pattern = re.compile(r'(?:Tax=)(.+)(?: TaxID=)')


# In[25]:


with open(infile, 'r') as f:
    IDs = [pattern.findall(_) for _ in f]
IDs = sorted(set([item for sublist in IDs for item in sublist]))


# print('#entries = %d' % len(IDs))

# In[7]:


with open(outfile, 'w') as f: f.write('\n'.join(IDs) + '\n')


# In[ ]:





# In[ ]:





# with open(infile, 'r') as f:
#     sentence = f.read()

# pattern = re.compile(r'(?:Tax=)(.+)(?: TaxID=)')

# pattern.findall(sentence)

# sentence

# import sys
# 
# infile, outfile = sys.argv[1:]
# 

# with open(infile, 'r') as f:
#     IDs = [_.strip() for _ in f]

# import urllib.parse
# import urllib.request

# url = 'https://www.uniprot.org/uploadlists/'

# params = {
# 'from': 'UPARC',
# 'to': 'ACC',
# 'format': 'tab',
# 'query': ''
# }
# 
# for ID in IDs:
#     if ID.startswith("UPI"):
#         params["query"] = ID
#     data = urllib.parse.urlencode(params)
#     data = data.encode('utf-8')
#     req = urllib.request.Request(url, data)
#     with urllib.request.urlopen(req) as f:
#         response = f.read()
#     print(response.decode('utf-8'))

# uni_list = [ID for ID in IDs if not ID.startswith("UPI")]
# upi_list = [ID for ID in IDs if ID.startswith("UPI")]
# query = ' '.join(upi_list)

# params = {
# 'from': 'UPARC',
# 'to': 'ACC',
# 'format': 'tab',
# 'query': ''
# }
# 
# 
# params["query"] = query
# data = urllib.parse.urlencode(params)
# data = data.encode('utf-8')
# req = urllib.request.Request(url, data)
# with urllib.request.urlopen(req) as f:
#     response = f.read()
# #print(response.decode('utf-8'))

# new_response = response[8:].decode()

# p = re.compile(r'\t([a-zA-Z0-9]+)')
# uni_list_add = p.findall(new_response)

# uni_list = uni_list + uni_list_add

# uni_list = sorted(set(uni_list))

# with open(outfile, 'w') as f: f.write('\n'.join(uni_list) + '\n')
