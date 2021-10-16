#!/usr/bin/env python
# coding: utf-8

# %config Completer.use_jedi = False

# infile = "/home/chingyuenliu/complex_contact/uniref_id_list"
# outfile = "/home/chingyuenliu/complex_contact/uniprot_id_list_uniref"

# In[9]:


http_limit = 10000


# In[8]:


import sys

infile, uniref_id_UPI, uniref_id_UNI = sys.argv[1:4]


# In[3]:


with open(infile, 'r') as f:
    IDs = [_.strip() for _ in f]


# In[4]:


import urllib.parse
import urllib.request


# In[5]:


url = 'https://www.uniprot.org/uploadlists/'


# In[6]:


uni_list = [ID for ID in IDs if not ID.startswith("UPI")]
upi_list = [ID for ID in IDs if ID.startswith("UPI")]

import os
dirname = os.path.dirname(infile)

with open(uniref_id_UPI, 'w') as f: f.write('\n'.join(upi_list) + '\n')
with open(uniref_id_UNI, 'w') as f: f.write('\n'.join(uni_list) + '\n')

# In[12]:


#query_list = []
#for i in range(0, len(upi_list), http_limit):
#    query_list.append(' '.join(upi_list[i:i+http_limit]))
#
#
## In[15]:
#
#
#import re
#
#
## In[16]:
#
#
#params = {
#'from': 'UPARC',
#'to': 'ACC',
#'format': 'tab',
#'query': ''
#}
#
#p = re.compile(r'\t([a-zA-Z0-9]+)')
#
#
## In[17]:
#
#
#for query in query_list:
#    params["query"] = query
#    data = urllib.parse.urlencode(params)
#    data = data.encode('utf-8')
#    req = urllib.request.Request(url, data)
#    with urllib.request.urlopen(req) as f:
#        response = f.read()
#    #print(response.decode('utf-8'))
#    new_response = response[8:].decode()
#    uni_list_add = p.findall(new_response)
#    uni_list = uni_list + uni_list_add
##    print("this size work!\n")
##    break
#
#uni_list = sorted(set(uni_list))
#print('#entries = %d' % len(uni_list))
#
#
## In[18]:
#
#
#with open(outfile, 'w') as f: f.write('\n'.join(uni_list) + '\n')
#

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




