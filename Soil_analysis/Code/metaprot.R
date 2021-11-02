#### Library
library(tidyr)
library(dplyr)

#### Export dataset
soil_prot_1<-read.csv("Data/proteins_2_unique.csv")
str(soil_prot_1)


#### Quantification by spectrum
need_data<-soil_prot_1 %>% dplyr::select("Protein.Group",Accession,Description,starts_with("X.Spec")) 
sum(rowSums(need_data[,4:17])==1)
need_long<-need_data %>% pivot_longer(cols = 4:17, names_to = "probe_name", values_to = "n_spec")

#### Sample names
need_long$probe_name<-factor(need_long$probe_name)
levels(need_long$probe_name)
legend<-c("Au","K3", "K4", "K5", "Kp", "mor","B1","B2","B3", "B4", "B5", "B6", "K1", "K2" )
levels(need_long$probe_name)<-legend

#### Remove 3 unquality samples
need_long<-need_long %>% filter(probe_name!="Au"&probe_name!="Kp"&probe_name!="mor")

#### Taxonomy filter
codes<-read.csv("Data/codes_maml.csv", sep = "\t")
library(stringr)
ox<-sub(".*OX=", "", need_long$Description)
need_long$codes<-as.numeric(sub(" .*", "", ox))

#### Filter by Chordata
filter_by_maml<-anti_join(x = need_long, y = codes, by = c("codes"="Taxon"))
# write.csv(x = filter_by_maml, file = "maml_1.csv")
codes_chord<-read.csv("Data/taxonomy-ancestor_chordata1.tab", sep = "\t", quote = "")
filter_by_chord<-anti_join(x = filter_by_maml, y = codes_chord, by = c("codes"="Taxon"))
# write.csv(x = filter_by_chord, file = "chord_1.csv")

#### Filter by Streptophyta
codes_strepto<-read.csv("Data/taxonomy-ancestor_streptophyta.tab", sep = "\t")
filter_by_strepto<-anti_join(x = filter_by_chord, y = codes_strepto, by = c("codes"="Taxon"))
write.csv(x = filter_by_strepto, file = "Data/multi_filter_1.csv")
