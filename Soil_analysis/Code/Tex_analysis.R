#### Data

data<-read.csv("Data/multi_filter_1.csv")
data$probe_name<-as.factor(data$probe_name)
levels(data$probe_name)
data$group<-factor(data$probe_name, labels = c( "KZ", "KZ", "KZ", "KZ", "KZ", "KZ", "NZ", "NZ", "NZ", "NZ", "NZ"))

data$Accession<-as.factor(data$Accession)

####Select row with 2 and more spectrum in sample
data_spec<-data %>% filter(n_spec>1)

#### how many sample have this peptide
observ_table<-data_spec %>% group_by(Accession) %>% summarise(n=n()) 
#### peptide in 2 or more sample

obil_prot<-observ_table %>% filter(n>1)
filter_prot<-inner_join(obil_prot, data_spec, by = "Accession")

length(unique(filter_prot$Accession))

#### Working tables: technical and biological filtered
prot_wide<-filter_prot  %>% select(!X&!group) %>% pivot_wider( names_from = probe_name, values_from = n_spec)
write.csv(prot_wide, "Data/all_good_prot.csv")

NZ_wide<-filter_prot %>% filter(group=="NZ") %>% select(!X&!group) %>% pivot_wider( names_from = probe_name, values_from = n_spec)
write.csv(NZ_wide, "Data/NZ_work.csv")

KZ_wide<-filter_prot %>% filter(group=="KZ") %>% select(!X&!group) %>% pivot_wider( names_from = probe_name, values_from = n_spec)
write.csv(KZ_wide, "Data/KZ_work.csv")