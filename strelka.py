import os, re

def strelka(sample):
	main_dir = '/storage2/Project/CSC/WES/03_SNV/Strelka'
	os.system('zcat %s/%s/results/variants/somatic.snvs.vcf.gz > %s/%s_snv.txt'%(main_dir,sample,main_dir,sample))
	result = open('%s/%s_snv.txt'%(main_dir,sample),'r')
	result_line = result.readlines()
	for i in result_line:
		if i[0] == "#": continue 
		PASS = i.split('\t')