#!/usr/bin/env python
# coding: utf-8

# %config Completer.use_jedi = False

## This script using uniprot.org cross reference service to find EMBL ID from UniprotID
## The server link is find through Uniprot.org database identifier mapping ('Retrieve/ID mapping') service programmatically access. Link: https://www.uniprot.org/help/api_idmapping
## Detailed usage can be checked from the linked mentioned above.
## For preventing the server shut down the connection from too much requests. Pause the process every often could solve this problem.

# In[7]:


step = 5000


# In[8]:


import sys
import time


infile, source, outfile = sys.argv[1:4]
lindex, rindex = map(int, sys.argv[4:6])


# In[5]:


with open(infile, 'r') as f: query_list = f.read().strip().split('\n', rindex)[:rindex]
query_list = query_list[lindex:]

# In[16]:


import urllib.parse
import urllib.request


# In[18]:


url = 'https://www.uniprot.org/uploadlists/'


# In[6]:





# In[8]:


query_joined_list = []
for i in range(0, len(query_list), step):
    query_joined_list.append(' '.join(query_list[i:i+step]))


# In[14]:


EMBL_list = []


# In[19]:



params_2 = {
'from': 'ACC+ID',
'to': 'EMBL_ID',
'format': 'tab',
'query': ''
}

params = {
'from': 'ACC+ID',
'to': 'EMBL',
'format': 'tab',
'query': ''
}

start_time = time.time()

for ind, query_j in enumerate(query_joined_list):
    params["query"] = query_j
    params_2["query"] = query_j
    
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()
    embl_first = response.decode()[8:].split('\n')
    
    data_2 = urllib.parse.urlencode(params_2)
    data_2 = data_2.encode('utf-8')
    req_2 = urllib.request.Request(url, data_2)
    with urllib.request.urlopen(req_2) as f:
        response_2 = f.read()
    embl_second = response_2.decode()[8:].split('\n')

    for fis in embl_first:
        if fis in embl_second:
            embl_second.remove(fis)
    
    EMBL_list += embl_second

    if ((ind+1)*step) % 1000 == 0 and (ind) != 0:
        duration = time.time() - start_time
        avg_time = duration / float((ind+1)*step)
        left_time = avg_time * (len(query_list)-(ind+1)*step)
        print("%d\t%.4f\t%.4f\t%.4f" % ((st+i), duration, avg_time, left_time))
        time.sleep(0.3)
        sys.stdout.flush()
# In[23]:


with open(outfile, 'w') as f: f.write('\n'.join(EMBL_list) + '\n')


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




