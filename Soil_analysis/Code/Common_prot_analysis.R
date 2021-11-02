#### Library
library(dplyr)
library(stringr)
library(ggplot2)

#### Dataset
KZ_prot<-read.csv("Data/KZ_work.csv")
NZ_prot<-read.csv("Data/NZ_work.csv")

common_prot<-inner_join(KZ_prot, NZ_prot, by = "Accession") %>% select(!c("X.x",  "n.x", "codes.x" , "X.y","n.y","Protein.Group.y","Description.y","codes.y" ))
common_prot<-group_by(common_prot, Protein_group) %>% slice(1)
#нашла общие белки для 2 точек
 colnames(common_prot)[2:3]<-c("Protein_group", "Description")
common_prot$Accession
uni_ass<-str_split(common_prot$Accession, "\\|" )
common_prot$uniprot<-sapply(uni_ass, "[[",2)
write.csv(common_prot[,15], "Data/common_prot_KZ_NZ.csv",quote = F)

#### Uniprot information
GO<-read.table("Data/Uniprot_common_prot.tab", sep = "\t", header = T, quote = "")
colnames(GO)[1]<-"uniprot"
res<-left_join(common_prot, GO, by = "uniprot") %>% arrange(Protein_group)

#### GO annotation
GO_ind<-str_split(res$Gene.ontology.IDs, "; ")
unique(unlist(GO_ind))
GO_ind<-str_split(res$Gene.ontology..GO., "; ")
vis_table<-as.data.frame(table(unlist(GO_ind)))

#### GO annottation graph
vis<-vis_table %>% filter(Var1!="") %>% arrange(desc(Freq)) %>% slice_head(n = 30) %>% ggplot(aes(Freq,reorder(Var1, Freq)))+geom_bar(stat = "identity")+labs(title = "Топ 30 GO annotaion Caucasus|Novaya_Zemlya common protein")+ylab(label = "Annotation")


#### Tree of life
codes_euk<-read.csv("Data/taxonomy-filtered-ancestor__Eukaryota+[2759]_.tab", sep = "\t", quote = "")#Eukaryota codes

codes_bact<-read.csv("Data/taxonomy-ancestor_Bact.tab", sep = "\t", quote = "")#Bacteria codes

codes_arc<-read.csv("Data/taxonomy-ancestor_Archaea.tab", sep = "\t", quote = "")#Archaea codes

sp_codes<-rbind(codes_arc, codes_bact, codes_euk)

ox<-sub(".*OX=", "", common_prot$Description)
common_prot$codes<-as.numeric(sub(" .*", "", ox))

#### Common protein from 2 location
filter_by_sp<-left_join(x = common_prot, y = sp_codes, by = c("codes"="Taxon"))
filter_by_sp<-filter_by_sp %>% filter(Lineage!="Bacteria"&!is.na(Lineage)) 


rank<-str_split(filter_by_sp$Scientific.name, " ")

line<-str_split(filter_by_sp$Lineage, " ")
rank_all<-data.frame(S_sp=sapply(rank, "[[",1),L_sp=sapply(line, "[[",2))


#### Common protein for 2 location by taxonomy
rank_vis<-rank_all %>% group_by(S_sp, L_sp) %>% summarise(n = n()) %>% arrange(desc(n)) 
rank_vis$L_sp<-factor(rank_vis$L_sp)
rank_vis[rank_vis$L_sp=="Cyanobacteria", 2]<-"Cyanobacteria;"
ggplot(rank_vis, aes(n, reorder(L_sp,n), fill = L_sp))+geom_bar(stat="identity")+labs(title = "Common protein for 2 location by taxonomy. Top 30")+ylab(label = "Taxonomy")



