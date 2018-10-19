#vcf version
import os,re
maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect2'

sample=["PE17","PE18","PE20","PE24","PE32"]
for s in sample:
	result = open("%s/%s_annovar.out3.exonic_variant_function"%(maindir,s),'r') 
	result_lines = result.readlines()
	filtered_result = open("%s/%s.txt"%(maindir,s),'w')
	for i in result_lines:
		CHR = i.split('\t')[3]
		POS = i.split('\t')[4]
		PASS = i.split('\t')[10]
#		AF = i.split('\t')[13].split(':')[4]
		if PASS == "PASS":
			filtered_result.write("%s\t%s\n"%(CHR,POS))
	result.close()
	filtered_result.close()

M = open("/storage2/Project/CSC/WES/03_SNV/MuTect/PE32.txt",'r')
M_list = set([line.rstrip("\n") for line in M.readlines()])
M.close()

M2 = open("/storage2/Project/CSC/WES/03_SNV/MuTect2/PE32.txt",'r')
M2_list = set([line.rstrip("\n") for line in M2.readlines()])
M2.close()

>>> len(M_list)
5848
>>> len(M2_list)
2030

M_list&M2_list
>>> len(M_list&M2_list)
1924

#그냥 숫자 입력하면 벤다이어그램 그려주는 웹서비스
#http://eulerr.co/
#
cut -f3-8,11,13-15 PE17_annovar.out3.exonic_variant_function > PE17_cut.txt
cut -f10 PE17_cut.txt | sed 's/:/\t/g' | cut -f5 > PE17_AF.txt
paste PE17_cut.txt PE17_AF.txt > PE17_result.txt
#
#
import os,re
maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect'

result = open("%s/PE17_cut.txt"%maindir,'r') 
result_lines = result.readlines()
filtered_result = open("%s/1.txt"%maindir,'w')
for i in result_lines:
	BC_AF = float(i.split('\t')[8].split(':')[4])
	PE_AF = float(i.split('\t')[9].split(':')[4])
	if (BC_AF != 1.0) & (PE_AF == 1.0) :
		filtered_result.write("%s"%i)

result.close()
filtered_result.close()
#
## PE_AF-BC_AF 하고 histogram 그리는 code in python
import pandas as pandas
import matplotlib; matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import os, re
import numpy as np
import seaborn as sb

sample='PE32'
AF_hist(sample)

def AF_hist(sample):
	maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect'
	result = open("%s/%s_cut.txt"%(maindir,sample),'r') 
	result_lines = result.readlines()
	AF_array = np.array([])
	for i in result_lines:
		BC_AF = float(i.split('\t')[8].split(':')[4])
		PE_AF = float(i.split('\t')[9].split(':')[4])
		dif = abs(float(PE_AF - BC_AF))
		AF_array = np.append(AF_array,dif)
	result.close()
	file_name = '%s_AF_hist.png'%sample
	save_file = os.path.join(maindir, file_name) 
	fig = plt.figure() 
	sb.distplot(AF_array, rug=True)
	fig.savefig(save_file)

#################################################
