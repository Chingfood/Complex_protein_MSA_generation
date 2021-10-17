#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, time, subprocess, io, gzip, re
from enaBrowserTools.python3.utils import is_wgs_set, is_sequence
import os
import subprocess
# In[2]:

enaDataGet = "/home/chingyuenliu/complex_contact/enaBrowserTools/python3/enaDataGet.py"



# In[3]:





# In[4]:


def dir_remove(dir_path, include_base_dir = False):
    # Delete everything reachable from the directory named in 'top',
    # assuming there are no symbolic links.
    # CAUTION:  This is dangerous!  For example, if top == '/', it
    # could delete all your disk files.
    if not os.path.isdir(dir_path):
        return
    top = dir_path
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    if include_base_dir == True:
        os.rmdir(dir_path)


# In[5]:





# In[6]:

def execute_download(path,executable_path,ID,wgs = False):
    if wgs:
        try:
            completed = subprocess.run(["python3", executable_path,"-w", "-f", "embl", ID, "-d", path ],timeout = 18000)
            return completed, 0
        except TimeoutExpired:
            return None, -1
        except:
            return None, -10
    else:
        try:
            completed = subprocess.run(["python3", executable_path, "-f", "embl", ID, "-d", path ],timeout = 18000)
            return completed, 0
        except TimeoutExpired:
                return None, -1
        except:
            return None, -10
        


def download_embl(path,executable_path,ID,wgs = False):
    counter = 0,
    flag = -10
    while (flag != 0 and counter<3):
        completed, flag = execute_download(path,executable_path,ID,wgs)
        counter += 1
    if completed == None:
        return False
    elif completed.returncode != 0:
        return False

    a = os.listdir(path)
    if len(a) == 0:
        return False

    completed2_returncode = 0
    with os.scandir(path) as directory:
        for entry in directory:
            if entry.name.endswith('.gz'):
                completed_2 = subprocess.run(["gunzip", entry.path])
                completed2_returncode = completed_2.returncode
    if completed.returncode != 0:
        return False
    elif completed2_returncode != 0:
        return False
    else:
        return True


# In[ ]:


# python download_contig.py contig_list contig_CDS 0 10000

infile, outfile = sys.argv[1:3]
lindex, rindex = list(map(int, sys.argv[3:5]))


# infile = "/home/chingyuenliu/complex_contact/contig_list_test"
# outfile = "/home/chingyuenliu/complex_contact/contig_CDS_0"
# lindex = 1
# rindex = 20

# In[8]:


print('loading queries')
with open(infile, 'r') as f: query_list = f.read().strip().split('\n', rindex)[:rindex]
query_list = query_list[lindex:]


# In[9]:


sec = infile.split('/')[-1].split('_')[0]


# In[10]:


download_path = os.path.join(os.path.dirname(outfile),"temp",sec+"_temp"+str(lindex))


# In[11]:


dir_remove(download_path,include_base_dir=True)


# In[12]:


os.makedirs(download_path)


# In[13]:


f = open(outfile, 'w')


# In[14]:


start_time = time.time()
print('begin')
sys.stdout.flush()

# for i, (data_class, ID) in enumerate(query_list):
for i, ID in enumerate(query_list):
    print('%d\t%s' % (lindex+i, ID))

    if is_wgs_set(ID):
        return_code = download_embl(download_path,enaDataGet,ID,wgs = True)
    elif is_sequence(ID):
        return_code = download_embl(download_path,enaDataGet,ID,wgs = False)
    if not return_code:
        continue
    
    with os.scandir(download_path) as directory:
        for entry in directory:
            if entry.is_file() and entry.name.endswith(".dat"): 
                with open(entry.path) as f1:
                    r = f1.read()
                key = ''
                qualifier = ''
                for line in r.strip().split('\n'):
                    if len(line) < 2: print('ERROR short line\t' + ID + '-' + line + '-')
                    if line[:2] == 'ID': f.write(line + '\n')
                    if line[:2] != 'FT': continue
                    if line[5] != ' ' and line[21] == '/': print('ERROR key and qualifier\t' + ID)
                    if line[5] != ' ': key = line[5:20].rstrip(); qualifier = ''
                    if line[21] == '/': qualifier = line[21:]
                    if key == 'CDS' and (
                            not qualifier or
                            qualifier.startswith('/db_xref="UniProtKB/') or
                            qualifier.startswith('/protein_id="')
                    ): f.write(line + '\n')

    i += 1
    if i % (10 if is_wgs_set(ID) else 100) == 0:
        duration = time.time() - start_time
        avg_time = duration / float(i)
        left_time = avg_time * (len(query_list) - i)
        print("%d\t%.4f\t%.4f\t%.4f" % (i, duration, avg_time, left_time))
        sys.stdout.flush()
    
    dir_remove(download_path)

print('%f' % (time.time() - start_time))


# In[15]:


f.close()


# In[16]:


dir_remove(download_path,include_base_dir = True)

