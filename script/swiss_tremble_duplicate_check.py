#!/usr/bin/env python
# coding: utf-8

# In[1]:


merged_list_path = "/home/qingyliu/ENA/merged_xref_list_2017_10"


# In[2]:


track_dict = dict()


# In[ ]:





# In[16]:


with open(merged_list_path,'rb') as f:
    for line in f:
        elements = line.strip().replace(b'\n',b'').split(sep=None,maxsplit=2)
        if elements[1] not in track_dict:
            if elements[0] == b"UniProtKB_Swiss-Prot":
                track_dict[elements[1]] = [1,0]
            else:
                track_dict[elements[1]] = [0,1]
        else:
            if elements[0] == b"UniProtKB_Swiss-Prot":
                track_dict[elements[1]][0] +=1
            else:
                track_dict[elements[1]][1] += 1
            


# In[35]:


tracking_list=[]
for key,value in track_dict.items():
    if all(value):
        tracking_list.append(key)


# In[36]:


if len(tracking_list) > 0:
    print("entries in both library")
    print(tracking_list)
else:
    print("there is no duplicate")
