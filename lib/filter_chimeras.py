'''
Created on Jan 22, 2011

@author: mkiyer
'''

#awk '$8 >= 2' ./out2/chimeras.bedpe | sort -k19,19gr -k20,20gr | cut -f7-14,18-21 | grep -v "Overlapping_Same"