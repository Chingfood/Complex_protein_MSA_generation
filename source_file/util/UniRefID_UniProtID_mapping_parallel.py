#!/usr/bin/env python
# coding: utf-8

# %config Completer.use_jedi = False

# In[16]:


http_limit = 5000


import sys


infile, outfile = sys.argv[1:3]
lindex, rindex = map(int, sys.argv[3:5])

# In[8]:




# In[28]:


import os
if os.path.exists(outfile):
    os.remove(outfile)





# In[9]:


with open(infile, 'r') as f: upi_list = f.read().strip().split('\n', rindex)[:rindex]
upi_list = upi_list[lindex:]


# In[11]:


import urllib.parse
import urllib.request


# In[12]:


url = 'https://www.uniprot.org/uploadlists/'


# In[17]:


uni_list = []

query_list = []
for i in range(0, len(upi_list), http_limit):
    query_list.append(' '.join(upi_list[i:i+http_limit]))


# In[20]:



# In[16]:


params = {
'from': 'UPARC',
'to': 'ACC',
'format': 'tab',
'query': ''
}


# In[29]:



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
        
    new_response = response[8:].decode()
    with open(outfile,'a') as f:
        f.write(new_response)
        

