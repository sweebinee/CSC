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
cut -f3-8,11,13-15 PE17_annovar.out3.exonic_variant_function > PE17_cut.txt #M1 M2 공통
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
import pandas as pd
import matplotlib; matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt
import os, re
import numpy as np

sample='PE32'
AF_hist(sample)

def AF_hist(sample):
	AF_dif_array(sample)
	file_name = '%s_AF_hist.png'%sample
	maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect'
	save_file = os.path.join(maindir, file_name) 
	fig = plt.figure() 
	sb.distplot(AF_array, rug=True)
	fig.savefig(save_file)

#################################################
def AF_dif_array(sample):
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
	return(AF_array)

def PASS_dif_depth(sample):
	maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect'
	result = open("%s/%s_cut.txt"%(maindir,sample),'r') 
	result_lines = result.readlines()
	for i in result_lines:
		PASS = i.split('\t')[6]
		if PASS == "PASS":
			BC_AF = float(i.split('\t')[8].split(':')[4])
			PE_AF = float(i.split('\t')[9].split(':')[4])
			dif = abs(float(PE_AF - BC_AF))
			BC_depth = sum([int(j) for j in i.split('\t')[8].split(':')[1].split(',')[:]])
			PE_depth = sum([int(j) for j in i.split('\t')[9].split(':')[1].split(',')[:]])
			print "%f\t%i\t%i"%(dif,BC_depth,PE_depth)
	result.close()

def grep_PASSM1(sample):
	maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect'
	result = open("%s/%s_cut.txt"%(maindir,sample),'r') 
	result_lines = result.readlines()
	snv_call = pd.DataFrame(columns=("Total", "Ref", "Alt", "AF"))
	for i in result_lines:
		PASS = i.split('\t')[6]
		if PASS == "PASS":
			CHR = i.split('\t')[1]
			POS = i.split('\t')[2]
			ref = i.split('\t')[4]
			alt = i.split('\t')[5]
			MUT = CHR+":"+POS+":"+ref+":"+alt
			total = sum([int(j) for j in i.split('\t')[9].split(':')[1].split(',')[:]])
			Ref = int(i.split('\t')[9].split(':')[1].split(',')[0])
			Alt = int(i.split('\t')[9].split(':')[1].split(',')[1])
			AF = float(i.split('\t')[9].split(':')[4])
			snv_call.loc[MUT] = [total, Ref, Alt, AF]
	return snv_call

def grep_PASSM2(sample):
	maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect2'
	result = open("%s/%s_cut.txt"%(maindir,sample),'r') 
	result_lines = result.readlines()
	snv_call = pd.DataFrame(columns=("Total", "Ref", "Alt", "AF"))
	for i in result_lines:
		PASS = i.split('\t')[6]
		if PASS == "PASS":
			CHR = i.split('\t')[1]
			POS = i.split('\t')[2]
			ref = i.split('\t')[4]
			alt = i.split('\t')[5]
			MUT = CHR+":"+POS+":"+ref+":"+alt
			total = sum([int(j) for j in i.split('\t')[9].split(':')[1].split(',')[:]])
			Ref = int(i.split('\t')[9].split(':')[1].split(',')[0])
			Alt = int(i.split('\t')[9].split(':')[1].split(',')[1])
			AF = float(i.split('\t')[9].split(':')[2])
			snv_call.loc[MUT] = [total, Ref, Alt, AF]
	return snv_call


def grep_PASSstrelka(sample):
	main_dir = '/storage2/Project/CSC/WES/03_SNV/Strelka'
	os.system('zcat %s/%s/results/variants/somatic.snvs.vcf.gz > %s/%s_snv.txt'%(main_dir,sample,main_dir,sample))
	result = open('%s/%s_snv.txt'%(main_dir,sample),'r')
	result_line = result.readlines()
	snv = set()
	for i in result_line:
		if i[0] == "#": continue 
		PASS = i.split('\t')[6]
		if PASS == 'PASS':
			CHR = i.split('\t')[0].replace('chr','')
			POS = i.split('\t')[1]
			ref = i.split('\t')[3]
			alt = i.split('\t')[4]
			mut = CHR + ':' + POS+":"+ref+":"+alt
			snv.add(mut)
	return snv


def snv_summary(sample):
	maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect'
	result = open("%s/%s_cut.txt"%(maindir,sample),'r') 
	result_lines = result.readlines()
	snv_call = pd.DataFrame(columns=("Total", "Ref", "Alt", "AF"))
	for i in result_lines:
		CHR = i.split('\t')[1]
		POS = i.split('\t')[2]
		MUT = CHR+":"+POS
		if MUT in union:
			total = sum([int(j) for j in i.split('\t')[9].split(':')[1].split(',')[:]])
			Ref = int(i.split('\t')[9].split(':')[1].split(',')[0])
			Alt = int(i.split('\t')[9].split(':')[1].split(',')[1])
			AF = float(i.split('\t')[9].split(':')[4])
			snv_call.loc[MUT] = [total, Ref, Alt, AF]
	return snv_call

def grep_union():
	PE17i = set(grep_PASS('PE17').index)
	PE18i = set(grep_PASS('PE18').index)
	PE20i = set(grep_PASS('PE20').index)
	PE24i = set(grep_PASS('PE24').index)
	PE32i = set(grep_PASS('PE32').index)
	union = PE17i|PE18i|PE20i|PE24i|PE32i
	maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect'
	PE17 = snv_summary("PE17")
	PE18 = snv_summary("PE18")
	PE20 = snv_summary("PE20")
	PE24 = snv_summary("PE24")
	PE32 = snv_summary("PE32")
	for i in union:
		result = open("%s/union_%s.txt"%(maindir,i),'w')
		result.write("Sample\tTotal\tRef\tAlt\tAF\n")
		if i in PE17.index:
			j = list(PE17.loc[i])
			result.write("PE17\t%s\t%s\t%s\t%s\n"%(j[0],j[1],j[2],j[3]))
		if i in PE18.index:
			j = list(PE18.loc[i])
			result.write("PE18\t%s\t%s\t%s\t%s\n"%(j[0],j[1],j[2],j[3]))
		if i in PE20.index:
			j = list(PE20.loc[i])
			result.write("PE20\t%s\t%s\t%s\t%s\n"%(j[0],j[1],j[2],j[3]))
		if i in PE24.index:
			j = list(PE24.loc[i])
			result.write("PE24\t%s\t%s\t%s\t%s\n"%(j[0],j[1],j[2],j[3]))
		if i in PE32.index:
			j = list(PE32.loc[i])
			result.write("PE32\t%s\t%s\t%s\t%s\n"%(j[0],j[1],j[2],j[3]))
		result.close()


M1 = set(grep_PASSM1('PE32').index)
M2 = set(grep_PASSM2('PE32').index)
S = grep_PASSstrelka('PE32')

maindir = '/storage2/Project/CSC/WES/03_SNV/MuTect'
file_name = 'M1M2ST_union_SNV_venn.png'
save_file = os.path.join(maindir, file_name) 
fig = plt.figure() 

labels = get_labels([PE17,PE18,PE20,PE24,PE32], fill=['number'])
fig, ax  = venn5(labels, names=['PE17', 'PE18', 'PE20', 'PE24', 'PE32'])

fig.savefig(save_file)
plt.close()

sample='PE17'
M1 = set(grep_PASSM1(sample).index)
M2 = set(grep_PASSM2(sample).index)
S = grep_PASSstrelka(sample)
PE17 = M1|M2|S

###################################################################################################
###################################################################################################
##from sangok
#/storage2/storage/Project/GC/GC/3/src/snv_support.py 
# in rocks2

import pysam , sys 
import os 

def Readcount(pos,sam):
 pos_cols=pos.split(":")
 CHR = ''.join(['chr',pos_cols[0]])
 Ref=pos_cols[2] ;  Alt=pos_cols[3]
 Total_read_count = 0  ; Refcount=0; Altcount=0; mpileup_allele=[]
 for read in sam.pileup(CHR,int(pos_cols[1])-101,int(pos_cols[1])+101):
  if read.pos == int(pos_cols[1])-1 :
   Total_read_count=read.n
   for mpileup in read.pileups:
    if not mpileup.is_del : mpileup_allele.append(mpileup.alignment.seq[mpileup.query_position].upper())
   Refcount=mpileup_allele.count(Ref);Altcount=mpileup_allele.count(Alt)
   global AF
   AF = round(float(Altcount) / float(Total_read_count),2)
 return [str(Total_read_count),str(Refcount), str(Altcount), str(AF)]


sample='PE32'
M1 = set(grep_PASSM1(sample).index)
M2 = set(grep_PASSM2(sample).index)
S = grep_PASSstrelka(sample)
PE32 = M1|M2|S
union = PE17|PE18|PE20|PE24|PE32 #1181

sample_list=['PE17','PE18','PE20','PE24','PE32','BC24']
for sample in sample_list:
	main_dir="/storage2/Project/CSC/WES/03_SNV/MuTect"
	input_dir='/storage2/Project/CSC/WES/02_Preprocessing'
	snvoutput=open("%s/%s.M1M2ST_union_readcount_table"%(main_dir,sample),"w")
	Exome="%s/%s_RECAL.bam"%(input_dir,sample)
	Exome_samfile=pysam.Samfile("%s"%Exome,"rb")
	#pysam.index(Exome)
	print "Start", sample
	for pos in union:
		POS = pos.replace(':','\t').split("\t")
		result = POS + Readcount(pos,Exome_samfile)
		#if int(result[7])>0:
		snvoutput.write("\t".join(result)); snvoutput.write("\n")
	snvoutput.close()

snv = '11:130911872:C:A'
def AFtable(snv):
	CHR = snv.split(":")[0]
	POS = snv.split(":")[1]
	Ref = snv.split(":")[2]
	Alt = snv.split(":")[3]
	main_dir="/storage2/Project/CSC/WES/03_SNV"
	sample_list = ['BC24','PE17','PE18','PE20','PE24','PE32']
	AFtable = open("%s/%s_AFtable"%(main_dir,snv),"w")
	for sample in sample_list:
		result = open("%s/%s.M1M2ST_union_readcount_table"%(main_dir,sample),"r")
		for i in result.readlines():
			s_CHR = i.split('\t')[0]
			s_POS = i.split('\t')[1]
			s_Ref = i.split('\t')[2]
			s_Alt = i.split('\t')[3]
			if (CHR==s_CHR) & (POS==s_POS) & (Ref==s_Ref) & (Alt==s_Alt):
				print sample
				AFtable.write("%s\t%s"%(sample,i))
		result.close()
	AFtable.close()
