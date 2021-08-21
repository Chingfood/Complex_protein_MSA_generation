#!/bin/bash

grep -P '^100\t' *.log | wc
tail *.log -n1
grep '^[0-9]*\.[0-9]*$' *.log | wc
for f in 4*; do if [[ `printf ${f%%_*} | wc -c` == 8 ]]; then echo $f; fi; done
while true; do for i in 4 5 6; do date; x$i "cd /home/tianmingzhou/Complex2/ENA/contig_list_2017_10; printf \$HOSTNAME'\t'; grep rror *.log | wc -l;"; sleep 60; done; done
while true; do for i in 4 5 6; do date; x$i "cd /home/tianmingzhou/Complex2/ENA/contig_list_2017_10; printf \$HOSTNAME'\t'; grep raise *.log | wc -l;"; sleep 60; done; done
while true; do for i in 4 5 6; do date; x$i "cd /home/tianmingzhou/Complex2/ENA/contig_list_2017_10; printf \$HOSTNAME'\t'; grep 'Please try again later.' *.log | wc -l;"; sleep 60; done; done
while true; do for i in 4 5 6; do date; x$i "cd /home/tianmingzhou/Complex2/ENA/CDS_folder_2017_10; echo \$HOSTNAME; grep raise *.log"; sleep 60; done; done
while true; do for i in 4 5 6; do date; x$i "cd /home/tianmingzhou/Complex2/ENA/CDS_folder_2017_10; echo \$HOSTNAME; grep '^[0-9]*\.[0-9]*$' *.log | wc; ls *.log | wc"; sleep 60; done; done
mv * ../contig_list_2017_10_done/
