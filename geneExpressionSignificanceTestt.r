#pair-wise comparition of gene expression of each disease gene set with background genes in each tissue
#protein coding genes as background for example
#10-03-2015

rm(list=ls())
path_files=""; ## set the path of the files
file_Inflammatory=paste(path_files,"/Inflammatory_results.txt",sep=""); # the file saving gene expression data
file_Cardiovascular=paste(path_files,"/Cardiovascular_results.txt",sep="");
file_Cancer=paste(path_files,"/Cancer_results.txt",sep="");
file_Metabolic=paste(path_files,"/Metabolic_results.txt",sep="");
file_Neurodegenerative=paste(path_files,"/Neurodegenerative_results.txt",sep="");
file_Psychological=paste(path_files,"/Psychological_results.txt",sep="");
file_Infectious=paste(path_files,"/Infectious_results.txt",sep="");
file_background=paste(path_files,"/proteinCoding_results.txt",sep="");

#input the expression data
data_Inflammatory=read.table(file=file_Inflammatory,header=FALSE)
data_Cardiovascular=read.table(file=file_Cardiovascular,header=FALSE)
data_Cancer=read.table(file=file_Cancer,header=FALSE)
data_Metabolic=read.table(file=file_Metabolic,header=FALSE)
data_Neurodegenerative=read.table(file=file_Neurodegenerative,header=FALSE)
data_Psychological=read.table(file=file_Psychological,header=FALSE)
data_Infectious=read.table(file=file_Infectious,header=FALSE)
data_background=read.table(file=file_background,header=FALSE)

#the gene expression data are shown in the following order of tissues
tissues=c("ADIPOSE","ADRENAL GLAND","ANIMAL OVARY","APPENDIX","BLADDER","BONE MARROW","CEREBRAL CORTEX","COLON","DUODENUM","ENDOMETRIUM","ESOPHAGUS","FALLOPIAN TUBE","GALL BLADDER","HEART","KIDNEY","LIVER","LUNG","LYMPH NODE","PANCREAS","PLACENTA","PROSTATE","RECTUM","SALIVARY GLAND","SKELETAL MUSCLE","SKIN","SMALL INTESTINE","SMOOTH MUSCLE","SPLEEN","STOMACH","TESTIS","THYROID","TONSIL")



#compare between disease genes and background genes
pvalue_Inflammatory=rep(0,32)
for(i in 1:32)
{
    result=wilcox.test(log2(data_Inflammatory[,i+2]+0.25),log2(data_background[,i+2]+0.25))
    pvalue_Inflammatory[i]=result$p.value
}

pvalue_Cardiovascular=rep(0,32)
for(i in 1:32)
{
    result=wilcox.test(log2(data_Cardiovascular[,i+2]+0.25),log2(data_background[,i+2]+0.25))
    pvalue_Cardiovascular[i]=result$p.value
}


pvalue_Cancer=rep(0,32)
for(i in 1:32)
{
    result=wilcox.test(log2(data_Cancer[,i+2]+0.25),log2(data_background[,i+2]+0.25))
    pvalue_Cancer[i]=result$p.value
}

pvalue_Metabolic=rep(0,32)
for(i in 1:32)
{
    result=wilcox.test(log2(data_Metabolic[,i+2]+0.25),log2(data_background[,i+2]+0.25))
    pvalue_Metabolic[i]=result$p.value
}

pvalue_Neurodegenerative=rep(0,32)
for(i in 1:32)
{
    result=wilcox.test(log2(data_Neurodegenerative[,i+2]+0.25),log2(data_background[,i+2]+0.25))
    pvalue_Neurodegenerative[i]=result$p.value
}


pvalue_Psychological=rep(0,32)
for(i in 1:32)
{
    result=wilcox.test(log2(data_Psychological[,i+2]+0.25),log2(data_background[,i+2]+0.25))
    pvalue_Psychological[i]=result$p.value
}


pvalue_Infectious=rep(0,32)
for(i in 1:32)
{
    result=wilcox.test(log2(data_Infectious[,i+2]+0.25),log2(data_background[,i+2]+0.25))
    pvalue_Infectious[i]=result$p.value
}

#get adjusted pvalue by "fdr", adjusted in each tissue.
pvalues=cbind(tissues,pvalue_Inflammatory,pvalue_Cardiovascular,pvalue_Cancer,pvalue_Metabolic,pvalue_Neurodegenerative,pvalue_Psychological,pvalue_Infectious)
pvalues_corrected=pvalues
for(i in 1:32)
{
    pvalues_corrected[i,2:8]=p.adjust(pvalues[i,2:8],method="BH")
}


## please change to another address if want to save the result into another adress
file_pvalues="./pvalue_proteinCoding.txt"
file_pvalues_corrected="./pvalue_proteinCoding_corrected_fdr.txt"

write.table(pvalues,file=file_pvalues,quote = FALSE,sep="\t",row.names =FALSE)
write.table(pvalues_corrected,file=file_pvalues_corrected,quote = FALSE,sep="\t",row.names =FALSE)














