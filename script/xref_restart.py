#!/usr/bin/env python
# coding: utf-8

# In[1]:


import subprocess


# In[ ]:


import sys


# In[ ]:


import os


# In[2]:


directory_path = sys.argv[1]


# In[ ]:



file_ext = sys.argv[2]

ID_list_path = "/home/qingyliu/ENA/temp/data/2020_06/uniprot_ID_list_2020_06"
contig_list_folder_path= "/home/qingyliu/ENA/temp/data/2020_06/xref_list_test"
program_path = "/home/qingyliu/ENA/source_file/download_xref.sh"

# In[3]:


#def operation_x(dir_path,filename):
#    fb = subprocess.run(["tail","-1", os.path.join(dir_path,filename)], stdout=subprocess.PIPE)
#    if (fb.stdout[:3] == b'req'):
#        fb = subprocess.run(["tail","-12", os.path.join(dir_path,filename)], stdout=subprocess.PIPE) 
#    print(fb.stdout)
def operation_x(dir_path,filename):
    fb = subprocess.run(["tail","-1", os.path.join(dir_path,filename)], stdout=subprocess.PIPE)
    if (fb.stdout[:4] != b'done'):
        start,end, _ = filename.strip().split('_',2)
        start = int(start)
        end = int(end)
        for i in range(10):
            fb = subprocess.run(["tail",f"-{1+i}", os.path.join(dir_path,filename)], stdout=subprocess.PIPE)
            lst = fb.stdout.split(b'\n')[0].split()
            try:
                st=int(lst[0])
            except:
                continue
            else:
                new_filename = f"{start}_{start+st}_UniProtKB_TrEMBL"
                subprocess.run(["mv", os.path.join(dir_path,filename)[:-4], os.path.join(dir_path,new_filename)])
                subprocess.run(["mv", os.path.join(dir_path,filename), os.path.join(dir_path,new_filename+".log")])
                
                with open(os.path.join(dir_path,new_filename+".log"), 'a') as f:
                    f.write("done\n")
                subprocess.run(["mv", os.path.join(dir_path,filename)[:-4]+".err", os.path.join(dir_path,new_filename+".err")])  
                ID_list = ID_list_path
                contig_list_folder = contig_list_folder_path 
                program = program_path
                subprocess.run(["bash",program, ID_list,"UniProtKB/TrEMBL", contig_list_folder, str(start+st), str(end), str(1)])              
                break
# In[4]:




directory = os.fsencode(directory_path)

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(file_ext):
        operation_x(directory_path, filename)

