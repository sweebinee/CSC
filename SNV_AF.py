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

def grep_PASS(sample,tools):
	if tools == 'M1': maindir = '/storage2/Project/CSC/WES/03_SNV/test_mergedBAM/MuTect'
	elif tools == 'M2' : maindir = '/storage2/Project/CSC/WES/03_SNV/test_mergedBAM/MuTect2'
	elif tools == 'S' : maindir = '/storage2/Project/CSC/WES/03_SNV/test_mergedBAM/Strelka'
	result = open("%s/%s_cut.txt"%(maindir,sample),'r') 
	result_lines = result.readlines()
	snv = set()
	#snv_call = pd.DataFrame(columns=("Total", "Ref", "Alt", "AF"))
	for i in result_lines:
		PASS = i.split('\t')[6]
		if PASS == "PASS":
			CHR = i.split('\t')[1]
			POS = i.split('\t')[2]
			ref = i.split('\t')[4]
			alt = i.split('\t')[5]
			MUT = CHR+":"+POS+":"+ref+":"+alt
			snv.add(MUT)
			#total = sum([int(j) for j in i.split('\t')[9].split(':')[1].split(',')[:]])
			#Ref = int(i.split('\t')[9].split(':')[1].split(',')[0])
			#Alt = int(i.split('\t')[9].split(':')[1].split(',')[1])
			#AF = float(i.split('\t')[9].split(':')[2])
			#snv_call.loc[MUT] = [total, Ref, Alt, AF]
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




maindir = '/storage2/Project/CSC/WES/03_SNV/'
file_name = 'M1M2ST_union_SNV_venn.png'
save_file = os.path.join(maindir, file_name) 
fig = plt.figure() 

labels = get_labels([PE17,PE18,PE20,PE24,PE32], fill=['number'])
fig, ax  = venn5(labels, names=['PE17', 'PE18', 'PE20', 'PE24', 'PE32'])

fig.savefig(save_file)
plt.close()

sample='PE32'
M1 = set(grep_PASS(sample,'M1'))
M2 = set(grep_PASS(sample,'M2'))
S =  set(grep_PASS(sample,'S'))
PE32 = M1|M2|S

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

def Readcount_merged(pos,sam):
	pos_cols=pos.split(":")
	CHR = ''.join(['chr',pos_cols[0]])
	Ref=pos_cols[2];  Alt=pos_cols[3]
	Total_read_count = 0; Refcount=0; Altcount=0; mpileup_allele=[]
	PE17={'Ref':0,'Alt':0}; PE18={'Ref':0,'Alt':0}; PE20={'Ref':0,'Alt':0}; PE24={'Ref':0,'Alt':0}; PE32={'Ref':0,'Alt':0};
	for read in sam.pileup(CHR,int(pos_cols[1])-101,int(pos_cols[1])+101):
		if read.pos == int(pos_cols[1])-1 :
			Total_read_count=read.n
			for mpileup in read.pileups:
				if not mpileup.is_del : 
					mpileup_allele.append(mpileup.alignment.seq[mpileup.query_position].upper())
					sample=mpileup.alignment.get_tag("RG")
					if mpileup.alignment.seq[mpileup.query_position].upper() == Ref:
						if sample == 'PE17': PE17['Ref'] += 1; 
						elif sample == 'PE18': PE18['Ref'] += 1
						elif sample == 'PE20': PE20['Ref'] += 1
						elif sample == 'PE24': PE24['Ref'] += 1
						elif sample == 'PE32': PE32['Ref'] += 1
					elif mpileup.alignment.seq[mpileup.query_position].upper() == Alt:
						if sample == 'PE17': PE17['Alt'] += 1; 
						elif sample == 'PE18': PE18['Alt'] += 1
						elif sample == 'PE20': PE20['Alt'] += 1
						elif sample == 'PE24': PE24['Alt'] += 1
						elif sample == 'PE32': PE32['Alt'] += 1
				Refcount=mpileup_allele.count(Ref);Altcount=mpileup_allele.count(Alt)
			global AF
			AF = round(float(Altcount) / float(Total_read_count),2)
	return [str(Total_read_count),str(Refcount), str(Altcount), str(AF), 
	':'.join([str(PE17['Ref']+PE17['Alt']),str(PE17['Ref']),str(PE17['Alt'])]), 
	':'.join([str(PE18['Ref']+PE18['Alt']),str(PE18['Ref']),str(PE18['Alt'])]),
	':'.join([str(PE20['Ref']+PE20['Alt']),str(PE20['Ref']),str(PE20['Alt'])]),
	':'.join([str(PE24['Ref']+PE24['Alt']),str(PE24['Ref']),str(PE24['Alt'])]),
	':'.join([str(PE32['Ref']+PE32['Alt']),str(PE32['Ref']),str(PE32['Alt'])]),]


sample='PE17'
M1 = set(grep_PASSM1(sample).index)
M2 = set(grep_PASSM2(sample).index)
S = grep_PASSstrelka(sample)
PE32 = M1|M2|S
union = PE17|PE18|PE20|PE24|PE32 #1181

sample_list=['PE17','PE18','PE20','PE24','PE32','BC24']
#sample_list=['Patient01_merged']; M1=grep_PASS('Patient01_merged','M1'); S=grep_PASS('Patient01_merged','S'); union = M1|S
for sample in sample_list:
	main_dir="/storage2/Project/CSC/WES/03_SNV/HaplotypeCaller"
	input_dir='/storage2/Project/CSC/WES/02_Preprocessing'
	snvoutput=open("%s/%s_readcount_table"%(main_dir,sample),"w")
	Exome="%s/%s_RECAL.bam"%(input_dir,sample)
	#Exome="%s/%s.bam"%(input_dir,sample)
	Exome_samfile=pysam.Samfile("%s"%Exome,"rb")
	#pysam.index(Exome)
	print "Start", sample
	for pos in union:
		POS = pos.replace(':','\t').split("\t")
		result = POS + Readcount(pos,Exome_samfile)
		#result = POS + Readcount_merged(pos,Exome_samfile)
		#if int(result[7])>0:
		snvoutput.write("\t".join(result)); snvoutput.write("\n")
	snvoutput.close()

snv = '19:40881751:A:T '
def AFtable(snv):
	CHR = snv.split(":")[0]
	POS = snv.split(":")[1]
	Ref = snv.split(":")[2]
	Alt = snv.split(":")[3]
	main_dir="/storage2/Project/CSC/WES/03_SNV/HaplotypeCaller"
	sample_list = ['BC24','PE17','PE18','PE20','PE24','PE32']
	#sample_list=['Patient01_merged']
	AFtable = open("%s/%s_AFtable"%(main_dir,snv),"w")
	for sample in sample_list:
		result = open("%s/%s_readcount_table"%(main_dir,sample),"r")
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

##################################################################################


def HC_TminusN():
	sample_list=['PE17','PE18','PE20','PE24','PE32']
	main_dir = '/storage2/Project/CSC/WES/03_SNV/HaplotypeCaller'
	N = set()
	normal = open("%s/BC24_cut.txt"%main_dir,'r')
	for n in normal.readlines():
		n_chr = n.split('\t')[1]
		n_pos = n.split('\t')[2]
		n_ref = n.split('\t')[4]
		n_alt = n.split('\t')[5]
		N.add(':'.join([n_chr,n_pos,n_ref,n_alt]))
	for s in sample_list:
		print s
		T = set()
		tumor = open("%s/%s_cut.txt"%(main_dir,s),'r')
		TN = open("%s/%s_TminusN"%(main_dir,s),'w')
		for t in tumor.readlines():
			t_chr = t.split('\t')[1]
			t_pos = t.split('\t')[2]
			t_ref = t.split('\t')[4]
			t_alt = t.split('\t')[5]
			T.add(':'.join([t_chr,t_pos,t_ref,t_alt]))
		#TN.write("%s\t%s - %s = %s\n"%(s,str(len(T)),str(len(N)),str(len(T-N))))
		for i in T-N:
			TN.write("%s\n"%i)
		tumor.close()
		TN.close()

PE17 = set([':'.join([i.split('\t')[3],i.split('\t')[4],i.split('\t')[6],i.split('\t')[7]]) for i in open("%s/PE17_nsyn.txt"%(main_dir,),'r').readlines()])
PE18 = set([':'.join([i.split('\t')[3],i.split('\t')[4],i.split('\t')[6],i.split('\t')[7]]) for i in open("%s/PE18_nsyn.txt"%(main_dir,),'r').readlines()])
PE20 = set([':'.join([i.split('\t')[3],i.split('\t')[4],i.split('\t')[6],i.split('\t')[7]]) for i in open("%s/PE20_nsyn.txt"%(main_dir,),'r').readlines()])
PE24 = set([':'.join([i.split('\t')[3],i.split('\t')[4],i.split('\t')[6],i.split('\t')[7]]) for i in open("%s/PE24_nsyn.txt"%(main_dir,),'r').readlines()])
PE32 = set([':'.join([i.split('\t')[3],i.split('\t')[4],i.split('\t')[6],i.split('\t')[7]]) for i in open("%s/PE32_nsyn.txt"%(main_dir,),'r').readlines()])

maindir = '/storage2/Project/CSC/WES/03_SNV/HaplotypeCaller'
file_name = 'SNV_venn.png'
save_file = os.path.join(maindir, file_name) 
fig = plt.figure() 

labels = get_labels([PE17,PE18,PE20,PE24,PE32], fill=['number'])
fig, ax  = venn5(labels, names=['PE17', 'PE18', 'PE20', 'PE24', 'PE32'])

fig.savefig(save_file)
plt.close()

def print_set_HC(target,sample,file): 
	main_dir = '/storage2/Project/CSC/WES/03_SNV/HaplotypeCaller'
	result = open("%s/%s_annovar.out3.exonic_variant_function"%(main_dir,sample),'r')
	result_lines = set([i.strip() for i in result.readlines()])
	f=open("%s/%s.txt"%(main_dir,file),'w')
#	f_syn = open("%s/%s_syn.txt"%(main_dir,file),'w')
#	f_nsyn = open("%s/%s_nsyn.txt"%(main_dir,file),'w')
	for t in target: 
		print t
		CHR = t.split(':')[0]
		POS = t.split(':')[1]
		for line in result_lines:
			if line.split('\t')[3] == CHR and line.split('\t')[4] == POS :
				f.write("%s\n"%line)
#				if line.split('\t')[1] == 'synonymous SNV':
#					f_syn.write("%s\n"%line)
#				elif line.split('\t')[1] == 'nonsynonymous SNV':
#					f_nsyn.write("%s\n"%line)
	result.close()
	f.close()
#	f_syn.close()
#	f_nsyn.close()
