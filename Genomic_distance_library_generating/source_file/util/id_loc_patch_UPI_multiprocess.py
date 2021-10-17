#!/usr/bin/env python
# coding: utf-8

# %config Completer.use_jedi = False

# In[2]:




import time
import sys
import multiprocessing

begin_time = time.time()

idloc, UPI_UNI_mapping = sys.argv[1:3]
outfile = sys.argv[3]

#idloc = "/home/chingyuenliu/complex_contact/id_loc_test_2"
#UPI_UNI_mapping = "/home/chingyuenliu/complex_contact/uniref_id_UPI_UNI_multiprocess_test"
#outfile = "/home/chingyuenliu/Trash/test_outfile"

# In[30]:


idloc_uniprot = []
idloc_contig = []
idloc_rank = []

with open(idloc,'r') as f:
    for line in f:
        line_sep = line.strip().split()
        idloc_uniprot.append(line_sep[0])
        idloc_contig.append(line_sep[1])
        idloc_rank.append(line_sep[2])


# In[31]:


UPI_list = []
UNI_list = []
with open(UPI_UNI_mapping,'r') as f:
    for line in f:
        line_sep = line.strip().split()
        UPI_list.append(line_sep[0])
        UNI_list.append(line_sep[1])   


# In[32]:

start_time = time.time()
print(f"time to load the library: {start_time-begin_time}")
sys.stdout.flush()


def func_1 (argss):
    i, UPI, UNI = argss
    if (i%100 == 0 and i!=0):
        current_time = time.time()
        print(f"time_passed: {current_time-start_time} time_per_entry: {(current_time-start_time)/i} remaining_time: {(current_time-start_time)/i*(len(UNI_list)-i)}")
        sys.stdout.flush() 
    if UNI in idloc_uniprot:
        index = idloc_uniprot.index(UNI)
        return '\t'.join([UPI, idloc_contig[index], idloc_rank[index]])


input_list = zip( range(len(UNI_list)), UPI_list, UNI_list )


# creating a pool object
p = multiprocessing.Pool()
# map list to target function
results = p.map(func_1, input_list)

  
# print(results)

filtered_results = [value for value in results if value!=None]
# print(filtered_results)

# outfile_list = []
# for i, UNI in enumerate(UNI_list):
#     if UNI in idloc_uniprot:
#         index = idloc_uniprot.index(UNI)
#         outfile_list.append('\t'.join([UPI_list[i], idloc_contig[index], idloc_rank[index]]))
# ##         idloc_uniprot.append(UPI_list[i])
# ##         idloc_contig.append(idloc_contig[index])
# ##         idloc_rank.append(idloc_rank[index])
#     if (i%100 == 0 and i!=0):
#         current_time = time.time()
#         print(f"time_passed: {current_time-start_time} time_per_entry: {(current_time-start_time)/i} remaining_time: {(current_time-start_time)/i*(len(UNI_list)-i)}")
#         sys.stdout.flush()

# In[35]:


with open(outfile,'w') as f:
    for result in filtered_results:
        f.write(result)
        f.write('\n')

