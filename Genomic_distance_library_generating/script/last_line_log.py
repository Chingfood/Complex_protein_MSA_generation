#!/usr/bin/env python
# coding: utf-8

# In[1]:

## This helper script shows the last line of all the log file in the folder.

import subprocess


# In[ ]:


import sys


# In[ ]:


import os


# In[2]:


directory_path = sys.argv[1]


# In[ ]:



file_ext = sys.argv[2]


# In[3]:


def operation_x(dir_path,filename):
    fb = subprocess.run(["tail","-1", os.path.join(dir_path,filename)], stdout=subprocess.PIPE)
    if (fb.stdout[:3] == b'req'):
        fb = subprocess.run(["tail","-12", os.path.join(dir_path,filename)], stdout=subprocess.PIPE) 
    print(fb.stdout)


# In[4]:




directory = os.fsencode(directory_path)

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(file_ext):
        operation_x(directory_path, filename)

