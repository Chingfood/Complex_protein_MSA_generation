#!/usr/bin/env python
# coding: utf-8

# %config Completer.use_jedi = False

# This script converts UnirefID that starts with UPI to regular UniprotID
# This code use the cross reference service maintained by uniprot.org to convert UnirefID that starts with UPI to regular UniprotID
# This python script works conjointly with source_file/UPI_Uniprot_mapping.sh
# This script has a part precede the main function where it only loads part of the input file.  

# In[9]:

# number of entry per request. Too much entry may result broken pipe etc.
http_limit = 5000


# In[8]:


import sys


infile, outfile = sys.argv[1:3]
lindex, rindex = map(int, sys.argv[3:5])

# In[3]:

import os

# only loads part of the infile from left-index to right-index
with open(infile, 'r') as f: upi_list = f.read().strip().split('\n', rindex)[:rindex]
upi_list = upi_list[lindex:]

##with open(infile, 'r') as f:
##    upi_list = [_.strip() for _ in f]


# In[4]:


import urllib.parse
import urllib.request


# In[5]:


url = 'https://www.uniprot.org/uploadlists/'


# In[6]:

# The list that stored the converted value
uni_list = []
#uni_list = [ID for ID in IDs if not ID.startswith("UPI")]
#upi_list = [ID for ID in IDs if ID.startswith("UPI")]


# In[12]:

# For all the query ID, seperate them into chunks that based on the http_limit.
query_list = []
for i in range(0, len(upi_list), http_limit):
    query_list.append(' '.join(upi_list[i:i+http_limit]))


# In[15]:


import re


# In[16]:


params = {
'from': 'UPARC',
'to': 'ACC',
'format': 'tab',
'query': ''
}

p = re.compile(r'\t([a-zA-Z0-9]+)')


# In[17]:

for query in query_list:
    params["query"] = query
    
    data = urllib.parse.urlencode(params)
        
    data = data.encode('utf-8')
    
    count = 1
    while count <5: 
        try:
            req = urllib.request.Request(url, data)
            with urllib.request.urlopen(req) as f:
                response = f.read()
        except:
            count += 1
        else:
            break
    else:
        continue
    #print(response.decode('utf-8'))

    new_response = response[8:].decode()
    uni_list_add = p.findall(new_response)
    uni_list = uni_list + uni_list_add
    
uni_list = sorted(set(uni_list))
print('#entries = %d' % len(uni_list))


# In[18]:


with open(outfile, 'w') as f: f.write('\n'.join(uni_list) + '\n')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# new_response = response[8:].decode()

# p = re.compile(r'\t([a-zA-Z0-9]+)')
# uni_list_add = p.findall(new_response)

# uni_list = uni_list + uni_list_add

# uni_list = sorted(set(uni_list))

# print('#entries = %d' % len(uni_list))

# with open(outfile, 'w') as f: f.write('\n'.join(uni_list) + '\n')

# In[ ]:




