KZ_prot<-read.csv("Data/KZ_work.csv")
NZ_prot<-read.csv("Data/NZ_work.csv")
library(dplyr)
library(stringr)

NZ_uniq_prot<-anti_join(NZ_prot, KZ_prot, by = "Accession")

uni_ass<-str_split(NZ_uniq_prot$Accession, "\\|" )
NZ_uniq_prot$uniprot<-sapply(uni_ass, "[[",2)

NZ_uniq_prot<-group_by(NZ_uniq_prot, Protein.Group) %>% slice(1)
GO<-read.table("Data/Uniprot_NZ_prot.tab", sep = "\t", header = T)

colnames(GO)[1]<-"uniprot"
res<-left_join(NZ_uniq_prot, GO, by = "uniprot") %>% arrange(Protein.Group)

GO_ind<-str_split(res$Gene.ontology..GO., "; ")
vis_table_N<-as.data.frame(table(unlist(GO_ind)))

a<-str_split(vis_table_N$Var1, " \\[")
vis_table_N$Var2<-sapply(a,"[[",1)
library(scales)
library(ggplot2)
vis_table_N %>% filter(Var2!="") %>% arrange(desc(Freq)) %>% slice_head(n = 30) %>% ggplot(aes(reorder(Var2, Freq), Freq))+geom_bar(stat = "identity")+labs(title = "Топ 30 GO аннотаций по уникальным белкам NZ")+ylab(label = "Аннотации")+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3))+scale_x_discrete(labels = wrap_format(15))


#####
KZ_prot<-read.csv("Data/KZ_work.csv")
NZ_prot<-read.csv("Data/NZ_work.csv")
library(dplyr)
library(stringr)

KZ_uniq_prot<-anti_join(KZ_prot, NZ_prot, by = "Accession")

uni_ass<-str_split(KZ_uniq_prot$Accession, "\\|" )
KZ_uniq_prot$uniprot<-sapply(uni_ass, "[[",2)
write.csv(KZ_uniq_prot$uniprot, "KZ_uniprot_uniq",quote = F)
KZ_uniq_prot<-group_by(KZ_uniq_prot, Protein.Group) %>% slice(1)

GO<-read.table("Data/Uniprot_KZ_prot.tab", sep = "\t", header = T)

colnames(GO)[1]<-"uniprot"
res<-left_join(KZ_uniq_prot, GO, by = "uniprot") %>% arrange(Protein.Group) 

GO_ind<-str_split(res$Gene.ontology..GO., "; ")
vis_table_K<-as.data.frame(table(unlist(GO_ind))) 

a<-str_split(vis_table_K$Var1, " \\[")
vis_table_K$Var2<-sapply(a,"[[",1)
vis_table_K$Site<-"Caucasus"
p1<-vis_table_K %>% filter(Var2!="") %>% arrange(desc(Freq)) %>% slice_head(n = 15) %>% ggplot(aes(reorder(Var2, Freq), Freq, fill = Site))+geom_bar(stat = "identity")+theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3), axis.title.x=element_blank(),legend.title = element_blank(), axis.text=element_text(size=12), axis.title = element_text(size=12))+scale_x_discrete(labels = wrap_format(15))+scale_fill_brewer(palette="Set2")+ylab(label = "Number of annotations")

vis_table_N$Site<-"Novaya \nZemlya"

p2<-vis_table_N %>% filter(Var2!="") %>% arrange(desc(Freq)) %>% slice_head(n = 15) %>% ggplot(aes(reorder(Var2, Freq), Freq, fill = Site))+geom_bar(stat = "identity")+theme_bw()+theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.3),legend.title = element_blank(), axis.text=element_text(size=12), axis.title = element_text(size=12))+scale_x_discrete(labels = wrap_format(15))+ylab(label = "Number of annotations")+xlab(label = "GO annotations")




library(cowplot)

plot_grid(p1, p2, nrow = 2, label_x = "GO annotations")
