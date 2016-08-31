#generate boxplot comparing expression of gene sets in each tissue
#using protein coding genes as background for example
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

#####
#!! change the path if want to save it under a different address
pdf(file="./boxplot_geneExpression.pdf",width=28,height=56)
par(mfrow=c(8,4),mar= c(10,5,4,2) + 0.1,oma=c(4,2,2,1))
for(i in 1:32)
{
    data_all=c(log2(data_Inflammatory[,i+2]+0.25),log2(data_Cardiovascular[,i+2]+0.25),log2(data_Cancer[,i+2]+0.25),log2(data_Metabolic[,i+2]+0.25),log2(data_Neurodegenerative[,i+2]+0.25),log2(data_Psychological[,i+2]+0.25),log2(data_Infectious[,i+2]+0.25),log2(data_background[,i+2]+0.25))
    type=c(rep("Inflammatory",length(data_Inflammatory[,i+2])),rep("Cardiovascular",length(data_Cardiovascular[,i+2])),rep("Cancer",length(data_Cancer[,i+2])),rep("Metabolic",length(data_Metabolic[,i+2])),rep("Neurodegenerative",length(data_Neurodegenerative[,i+2])),rep("Psychological",length(data_Psychological[,i+2])),rep("Infectious",length(data_Infectious[,i+2])),rep("Background",length(data_background[,i+2])))
    data<-data.frame(data_all,type)
    #change the order of boxplot
    data$type<-factor(data$type, levels=c("Inflammatory", "Cardiovascular","Cancer","Metabolic","Neurodegenerative","Psychological","Infectious","Background"))
    
    boxplot(data_all~type,data,ylab=expression("log"[2]*"(FPKM+0.25))"),cex.lab=1.8,cex.axis=1.8,outpch = NA,las=2)
    stripchart(data_all~type,data,vertical = TRUE, method = "jitter", pch = 20,cex=0.5,bg = "bisque",add = TRUE)
    title(tissues[i], line = -2)
}

dev.off()













