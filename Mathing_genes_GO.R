
file1<-read.table ("match.txt", header=T)
file2<-read.table ("ella.txt", header=T, sep='\t')
vila<-merge(file1, file2, by="protein_id")
write.csv(vila, file="final_list_of_genes.csv")

match<-read.csv("final_list_of_genes.csv", h=T)
downregM<-read.csv("sigDownReg_larvae_new.csv", h=T)
upreg<-read.csv("sigUpReg_larvae_new.csv", h=T)
villa<-merge(match, upreg, by="locus_id")

vila<-merge(larvae, match, by="locus_id")

write.csv(villa, file="sigUpReg_larvae_final.csv")

larvae<-read.table("featurecount_results_22.txt", sep="\t",header=T,check.names=F)

