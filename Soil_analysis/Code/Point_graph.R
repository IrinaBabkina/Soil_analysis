KZ_prot<-read.csv("Data/KZ_work.csv")
NZ_prot<-read.csv("Data/NZ_work.csv")
library(dplyr)
library(ggplot2)
library(stringr)

#### Select protein group for taxonomy identification
NZ_uniq_prot<-group_by(NZ_prot, Protein.Group) %>% slice(1)
KZ_uniq_prot<-group_by(KZ_prot, Protein.Group) %>% slice(1)

##### Taxonomy base
codes_euk<-read.csv("Data/taxonomy-filtered-ancestor__Eukaryota+[2759]_.tab", sep = "\t", quote = "")#Eukaryota codes

codes_bact<-read.csv("Data/taxonomy-ancestor_Bact.tab", sep = "\t", quote = "")#Bacteria codes

codes_arc<-read.csv("Data/taxonomy-ancestor_Archaea.tab", sep = "\t", quote = "")#Archaea codes

sp_codes<-rbind(codes_arc, codes_bact, codes_euk)#Tree of life codes

#### Protein for KZ location with taxonomy
ox<-sub(".*OX=", "", KZ_uniq_prot$Description)
KZ_uniq_prot$codes<-as.numeric(sub(" .*", "", ox))
filter_by_sp<-left_join(x = KZ_uniq_prot, y = sp_codes, by = c("codes"="Taxon"))
filter_by_sp<-filter_by_sp %>% filter(Lineage!="Bacteria"&!is.na(Lineage)) 


rank<-str_split(filter_by_sp$Scientific.name, " ")

line<-str_split(filter_by_sp$Lineage, " ")
rank_all<-data.frame(S_sp=sapply(rank, "[[",1),L_sp=sapply(line, "[[",2))
rank_all$L_sp<-factor(rank_all$L_sp)

rank_vis_K<-rank_all %>% group_by( L_sp) %>% summarise(n = n()) %>% arrange(desc(n)) 


#### Protein for NZ location with taxonomy
ox<-sub(".*OX=", "", NZ_uniq_prot$Description)
NZ_uniq_prot$codes<-as.numeric(sub(" .*", "", ox))
filter_by_sp<-left_join(x = NZ_uniq_prot, y = sp_codes, by = c("codes"="Taxon"))
filter_by_sp<-filter_by_sp %>% filter(Lineage!="Bacteria"&!is.na(Lineage)) 

#### Arrange taxonomy group for identification
rank<-str_split(filter_by_sp$Scientific.name, " ")

line<-str_split(filter_by_sp$Lineage, " ")
rank_all<-data.frame(S_sp=sapply(rank, "[[",1),L_sp=sapply(line, "[[",2))
rank_all$L_sp<-factor(rank_all$L_sp)

rank_vis_N<-rank_all %>% group_by( L_sp) %>% summarise(n = n()) %>% arrange(desc(n))
rank_vis_K$Site<-c("Caucasus")
rank_vis_N$Site<-c("Novaya Zemlya")
plot_rank<-rbind(rank_vis_K, rank_vis_N)
plot_rank$L_sp<-str_remove(plot_rank$L_sp, pattern = ";")

#### Point plot
plot_rank %>% group_by(L_sp, Site) %>% summarise(n_n = sum(n)) %>% ggplot()+geom_point(aes(Site,reorder(L_sp, n_n), size = n_n, color = L_sp))+geom_text(aes(Site, L_sp, label = n_n), position = position_nudge(x = 0.1))+labs(title = "A")+ylab(label = "Microorganisms")+xlab(label = "number of identified protein")+ scale_fill_brewer(palette="Set3")+theme_bw()+theme(legend.position = "none", axis.text=element_text(size=12), axis.title = element_text(size=12))


