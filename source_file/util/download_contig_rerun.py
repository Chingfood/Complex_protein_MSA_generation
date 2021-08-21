#!/usr/bin/env python
# coding: utf-8

# In[4]:


#get_ipython().run_line_magic('config', 'Completer.use_jedi = False')


# In[38]:


import os
import re
import subprocess


# In[5]:


#chunk_size = 543


# In[6]:

#the number of sections for unfinished section to split up

new_chunk_ratio = 3

category = "contig"
#category = "wgs"

# In[46]:

#The temp folder in contig dowload folder.
temp_folder = "/home/qingyliu/ENA/CDS_wgs_folder_2018_08/temp"


# In[40]:

#contig_download python script
ex_path = "/home/qingyliu/ENA/download_contig.py"
#input list file
input_list_path = "/home/qingyliu/ENA/wgs_list_2018_08"


# In[31]:





# In[30]:


def num_seperate(filename):
    temp = re.compile("([a-zA-Z]+[_][a-zA-Z]+)([0-9]+)")
    res = temp.match(filename).groups()
    return res[1]


# In[34]:


def rm_unsuccessful_chunk(start_index, directory):
    with os.scandir(directory) as it:
        for entry in it:
            if entry.name.startswith(category+"_"+str(start_index)+'_') and not entry.name.endswith(".log"):
                str_split = entry.name.split('_')
                st_idx = int(str_split[1])
                ed_idx = int(str_split[2])
                os.remove(entry.path)
                os.remove(entry.path+'.log')
                return ed_idx
            else:
                continue
    return -1


# In[55]:


def run_download_contig(ex_path, input_list_path, output_file_path, start_index, end_index ):
    with open(output_file_path+".log", 'w') as f:
        subprocess.Popen(["python3",ex_path, input_list_path, output_file_path, str(start_index), str(end_index)],stdout = f, stderr = f)


# In[42]:





# In[56]:


def run_download(chunk_size, new_chunk_ratio, ex_path, input_list_path, temp_folder, start_index):
    
    for i in range(new_chunk_ratio):
        start_idx = start_index + chunk_size//new_chunk_ratio* i
        if i != new_chunk_ratio-1:
            end_idx = start_index + chunk_size//new_chunk_ratio * (i+1)
        else:
            end_idx = start_index + chunk_size
        filename = category+"_"+str(start_idx)+'_'+str(end_idx)
        output_file_path = os.path.join(os.path.dirname(temp_folder),filename)
        run_download_contig(ex_path,input_list_path,output_file_path,start_idx,end_idx)


# In[57]:


with os.scandir(temp_folder) as it:
    for entry in it:

        if not entry.name.startswith('.') and entry.is_dir():
            
            idx = int(num_seperate(entry.name))
            
            end_index = rm_unsuccessful_chunk(idx,os.path.dirname(temp_folder))
            if end_index == -1:
                print("error")
                break
            chunk_size = end_index - idx
            run_download(chunk_size,new_chunk_ratio,ex_path,input_list_path,temp_folder, idx)

