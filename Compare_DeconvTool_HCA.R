setwd("/storage2/Project/CSC/RNA/test_Deconv_tools/HCA_Raw")

library(Seurat)

filename="/storage2/Project/CSC/RNA/test_Deconv_tools/HCA_Raw/ica_cord_blood_h5.h5"
HCA <- Read10X_h5("/storage2/Project/CSC/RNA/test_Deconv_tools/HCA_Raw/ica_cord_blood_h5.h5")

library(rhdf5)
h5ls("/storage2/Project/CSC/RNA/test_Deconv_tools/HCA_Raw/ica_cord_blood_h5.h5")
    group       name       otype  dclass       dim
0       /     GRCh38   H5I_GROUP                  
1 /GRCh38   barcodes H5I_DATASET  STRING    384000
2 /GRCh38       data H5I_DATASET INTEGER 260473471
3 /GRCh38 gene_names H5I_DATASET  STRING     33694
4 /GRCh38      genes H5I_DATASET  STRING     33694
5 /GRCh38    indices H5I_DATASET INTEGER 260473471
6 /GRCh38     indptr H5I_DATASET INTEGER    384001
7 /GRCh38      shape H5I_DATASET INTEGER         2

HCA <- h5read("/storage2/Project/CSC/RNA/test_Deconv_tools/HCA_Raw/ica_cord_blood_h5.h5","GRCh38")
HCA$barcodes    HCA$gene_names  HCA$indices     HCA$shape       
HCA$data        HCA$genes       HCA$indptr      

