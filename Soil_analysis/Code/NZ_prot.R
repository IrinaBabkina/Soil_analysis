#### Dataset
NZ_prot<-read.csv("Data/NZ_work.csv")

#### Library
library(dplyr)
library(stringr)
library(ggplot2)

uni_ass<-str_split(NZ_prot$Accession, "\\|" )
NZ_prot$uniprot<-sapply(uni_ass, "[[",2)

#### Uniprot information
GO_NZ<-read.table("Data/Uniprot_NZ_prot.tab", sep = "\t", header = T)
GO_comm<-read.table("Data/Uniprot_common_prot.tab", sep = "\t", header = T , quote = "")
colnames(GO_comm)[1]<-"Zapros"
colnames(GO_NZ)[1]<-"Zapros"
GO<-rbind(GO_comm, GO_NZ)
colnames(GO)[1]<-"uniprot"

#### GO annotation
res<-left_join(NZ_prot, GO, by = "uniprot") %>% arrange(Protein.Group)
GO_ind<-str_split(res$Gene.ontology..GO., "; ")
vis_table<-as.data.frame(table(unlist(GO_ind)))

#### GO annotation graph
vis<-vis_table %>% filter(Var1!="") %>% arrange(desc(Freq)) %>% slice_head(n = 30) %>% ggplot(aes(Freq,reorder(Var1, Freq)))+geom_bar(stat = "identity")+labs(title = "Top 30 GO annotation Novaya_Zemlya")+ylab(label = "Annotation")
vis

#### Tree of life
codes_euk<-read.csv("Data/taxonomy-filtered-ancestor__Eukaryota+[2759]_.tab", sep = "\t", quote = "")#Eukaryota codes

codes_bact<-read.csv("Data/taxonomy-ancestor_Bact.tab", sep = "\t", quote = "")#Bacteria codes

codes_arc<-read.csv("Data/taxonomy-ancestor_Archaea.tab", sep = "\t", quote = "")#Archaea codes

sp_codes<-rbind(codes_arc, codes_bact, codes_euk)

ox<-sub(".*OX=", "", NZ_prot$Description)
NZ_prot$codes<-as.numeric(sub(" .*", "", ox))
filter_by_sp<-left_join(x = NZ_prot, y = sp_codes, by = c("codes"="Taxon"))
filter_by_sp<-filter_by_sp %>% filter(Lineage!="Bacteria"&Lineage!="Archaea"&!is.na(Lineage)) 


rank<-str_split(filter_by_sp$Scientific.name, " ")

line<-str_split(filter_by_sp$Lineage, " ")
rank_all<-data.frame(S_sp=sapply(rank, "[[",1),L_sp=sapply(line, "[[",2))

#### Graph protein taxonomy NZ
rank_vis_N<-rank_all %>% group_by(S_sp, L_sp) %>% summarise(n = n()) %>% arrange(desc(n)) 
rank_vis_N$L_sp<-str_remove(rank_vis_N$L_sp, pattern = ";")
ggplot(rank_vis_N[1:30,], aes(n, reorder(S_sp,n), fill = L_sp))+geom_bar(stat="identity")+labs(title = "Taxonomy for Novaya_Zemlya protein. Top 30")+ylab(label = "Taxonomy")
