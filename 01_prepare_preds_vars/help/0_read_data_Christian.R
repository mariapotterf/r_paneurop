
# from Christian


library(readxl)
library(terra)
library(dplyr)
library(ggplot2)
library(coin)
library(stringr)

setwd('C:/Users/ge45lep/Documents/2023_PanEuropean/rawData')

harz <- terra::vect("harz_final.gpkg")
austria <- terra::vect("austria_final.gpkg")
dresden <- terra::vect("dresden_final.gpkg")
cz <- terra::vect("cz_final.gpkg")
slk <- terra::vect("slk_final.gpkg")
w_ger <- terra::vect("w_ger_final.gpkg")
lpz <- terra::vect("lpz_final.gpkg")
pdb <- terra::vect("pdb_final.gpkg")
frkfrt <- terra::vect("frkfrt_final.gpkg")
pl <- terra::vect("pl_final.gpkg")
it <- terra::vect("it_final.gpkg")
ch <- terra::vect("sw_final.gpkg")
#fr <- terra::vect("fr_final.gpkg")
#n_ger <- terra::vect("n_ger_final.gpkg")
#ber <- terra::vect("ber_final.gpkg")
#slv_cr <- terra::vect("slv_cr_final.gpkg")


all <- rbind(harz,austria,dresden,cz,slk,w_ger,lpz,pdb,frkfrt,pl,it,ch)#,fr,n_ger,ber,slv_cr
terra::plot(all)


austria_att <- read_excel("austria/att_tbl/austria_att.xlsx")
harz_att <- read_excel("harz/att_tbl/harz_att.xlsx")
dresden_att <- read_excel("dresden/att_tbl/dresden_att.xlsx")
cz_att <- read_excel("cz/att_tbl/cz_att.xlsx")
slk_att <- read_excel("slk/att_tbl/slk_att.xlsx")
w_ger_att <- read_excel("w_ger/att_tbl/w_ger_att.xlsx")
lpz_att <- read_excel("lpz/att_tbl/lpz_att.xlsx")
pdb_att <- read_excel("pdb/att_tbl/pdb_att.xlsx")
frkfrt_att <- read_excel("frkfrt/att_tbl/frkfrt_att.xlsx")
pl_att <- read_excel("pl/att_tbl/pl_att.xlsx")
it_att <- read_excel("it/att_tbls/it_att.xlsx")
ch_att <- read_excel("sw/att_tbl/sw_att_tbl.xlsx")
#fr_att <- read_excel("pl/att_tbl/fr_att.xlsx")
#n_ger_att <- read_excel("pl/att_tbl/n_ger_att.xlsx")
#ber_att <- read_excel("C:/working_data_rapid_assessment/pl/att_tbl/ber_att.xlsx")
#slv_cr_att <- read_excel("C:/working_data_rapid_assessment/pl/att_tbl/slv_cr_att.xlsx")

# get attributes
cz_df <- as.data.frame(cz)


austria_att <- austria_att%>%
  dplyr::select(-fid_1,-layer,-path)
harz_att <- harz_att%>%
  dplyr::select(-fid_1,-layer,-path)
dresden_att <- dresden_att%>%
  dplyr::select(-fid_1,-layer,-path)
cz_att <- cz_att%>%
  dplyr::select(-fid_1,-layer,-path)
slk_att <- slk_att%>%
  dplyr::select(-fid_1,-layer,-path)
w_ger_att <- w_ger_att%>%
  dplyr::select(-fid_1,-layer,-path)

all_att <- rbind(harz_att,austria_att,dresden_att,cz_att,slk_att,w_ger_att,lpz_att,pdb_att,frkfrt_att,pl_att,it_att,ch_att) #,,fr_att,n_ger_att,ber_att,slv_cr_att

View(all_att)
#View(att_table)
att_table <- all_att %>%
  dplyr::select(-r_tree_spc, -ar_tree_spc, -s_tree_spc, -comments) %>% #.$stump_otsp %>% unique()
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(region, group, point, sep="_")) %>%  #dplyr::select(5:323) %>% names()
  tidyr::gather(5:322, key = "key", value="value") %>% 
  dplyr::mutate(VegType = dplyr::case_when(grepl("^(r_)", key) ~ "Regeneration",  # !! add definition of those?
                                           grepl("^(ar_)", key) ~ "advRegeneration",
                                           grepl("^(s_)", key) ~ "Survivor",
                                           TRUE ~ "Else")) %>% 
  #dplyr::filter(VegType != "Else") %>% #.$VegType %>% unique()
  dplyr::mutate(Variable = ifelse(VegType == "Else", ifelse(grepl("stump", key), "stump", key), gsub("(.*?)\\s*_", "", key))) %>% #.$key %>% unique()
  dplyr::mutate(Species = ifelse(VegType == "Else", ifelse(grepl("stump", key), paste0(ifelse(!is.na(value), value, "")), NA),
                                 ifelse(grepl("name", key), value, 
                                        ifelse(VegType == "advRegeneration", substr(key, 4, 7), substr(key, 3, 6))))) %>%
  dplyr::select(group:Species, key) %>%
  #dplyr::filter(VegType == "advRegeneration", Variable %in% c("n")) %>%
  dplyr::mutate(region = dplyr::case_when(region == "11" ~ "Harz",
                                          region == "12" ~ "WGer",
                                          region == "13" ~ "Austria",
                                          region == "14" ~ "Dresden",
                                          region == "15" ~ "Czech",
                                          region == "16" ~ "Slovakia",
                                          region == "17" ~ "Poland",
                                          region == "18" ~ "MidGer",
                                          region == "19" ~ "Leipzig",
                                          region == "20" ~ "Frankfurt",
                                          region == "21" ~ "Italy",
                                          region == "22" ~ "Switzerland",
                                          TRUE ~ "NO")) %>% #.$VegType %>% unique()
  dplyr::mutate(n = dplyr::case_when(VegType == "Survivor" ~ as.double(value),
                                     VegType == "Else" ~ NA,
                                     TRUE ~ ifelse(as.double(value) == 1, 17, as.double(value) - 1))) %>% 
  dplyr::select(-key)


#dplyr::mutate(Variable = ifelse(VegType == "Else", ifelse(grepl("stump", key), "stump", key), gsub("(.*?)\\s*_", "", key))) %>% #.$key %>% unique()
#dplyr::mutate(Species = ifelse(VegType == "Else", ifelse(grepl("stump", key), paste0(ifelse(!is.na(value), value, "")), NA)))
View(att_table)
att_table

### decoding all r_**_n and ar_**_n


r_ar_n <- all_att[grep("r_**_n", all_att)] 

names(all_att)
set.seed(123)
all_att %>% #dplyr::sample_n(1000) %>%
  dplyr::select(-r_tree_spc, -ar_tree_spc, -s_tree_spc, -comments) %>% #.$stump_otsp %>% unique()
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(region, group, point, sep="_")) %>%  #dplyr::select(6:323) %>% names()
  tidyr::gather(5:322, key = "key", value="value") %>% 
  dplyr::mutate(VegType = dplyr::case_when(grepl("^(r_)", key) ~ "Regeneration",
                                           grepl("^(ar_)", key) ~ "advRegeneration",
                                           grepl("^(s_)", key) ~ "Survivor",
                                           TRUE ~ "Else")) %>% 
  #dplyr::filter(VegType != "Else") %>% #.$VegType %>% unique()
  dplyr::mutate(Variable = ifelse(VegType == "Else", ifelse(grepl("stump", key), "stump", key), gsub("(.*?)\\s*_", "", key))) %>% #.$key %>% unique()
  dplyr::mutate(Species = ifelse(VegType == "Else", ifelse(grepl("stump", key), paste0(ifelse(!is.na(value), value, "")), NA),
                                 ifelse(grepl("name", key), value, 
                                        ifelse(VegType == "advRegeneration", substr(key, 4, 7), substr(key, 3, 6))))) %>% #gsub("_\\s*(.*?)", "", 
  dplyr::select(group:Species, -key) %>%
  dplyr::filter(VegType == "Regeneration", Variable %in% c("n")) %>% 
  dplyr::mutate(region = dplyr::case_when(region == "11" ~ "Harz",
                                          region == "12" ~ "WGer",
                                          region == "13" ~ "Austria",
                                          region == "14" ~ "Dresden",
                                          region == "15" ~ "Czech",
                                          region == "16" ~ "Slovakia",
                                          region == "17" ~ "Poland",
                                          region == "18" ~ "MidGer",
                                          region == "19" ~ "Leipzig",
                                          region == "20" ~ "Frankfurt",
                                          region == "21" ~ "Italy",
                                          region == "22" ~ "Switzerland",
                                          TRUE ~ "NO")) %>% 
  ggplot(aes(x=region, y=as.double(value), col=as.factor(region)))+ 
  geom_boxplot()+
  geom_jitter(alpha=0.2)+ 
  facet_wrap(~Species) +
  theme_bw()+
  scale_x_discrete(labels = NULL)+
  theme()


all_att %>% #dplyr::sample_n(1000) %>%
  dplyr::select(-r_tree_spc, -ar_tree_spc, -s_tree_spc, -comments) %>% #.$stump_otsp %>% unique()
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(region, group, point, sep="_")) %>%  #dplyr::select(6:323) %>% names()
  tidyr::gather(5:322, key = "key", value="value") %>% 
  dplyr::mutate(VegType = dplyr::case_when(grepl("^(r_)", key) ~ "Regeneration",
                                           grepl("^(ar_)", key) ~ "advRegeneration",
                                           grepl("^(s_)", key) ~ "Survivor",
                                           TRUE ~ "Else")) %>% 
  #dplyr::filter(VegType != "Else") %>% #.$VegType %>% unique()
  dplyr::mutate(Variable = ifelse(VegType == "Else", ifelse(grepl("stump", key), "stump", key), gsub("(.*?)\\s*_", "", key))) %>% #.$key %>% unique()
  dplyr::mutate(Species = ifelse(VegType == "Else", ifelse(grepl("stump", key), paste0(ifelse(!is.na(value), value, "")), NA),
                                 ifelse(grepl("name", key), value, 
                                        ifelse(VegType == "advRegeneration", substr(key, 4, 7), substr(key, 3, 6))))) %>% #gsub("_\\s*(.*?)", "", 
  dplyr::select(group:Species, -key) %>%
  dplyr::filter(VegType == "Regeneration", Variable %in% c("n")) %>%
  dplyr::mutate(region = dplyr::case_when(region == "11" ~ "Harz",
                                          region == "12" ~ "WGer",
                                          region == "13" ~ "Austria",
                                          region == "14" ~ "Dresden",
                                          region == "15" ~ "Czech",
                                          region == "16" ~ "Slovakia",
                                          region == "17" ~ "Poland",
                                          region == "18" ~ "MidGer",
                                          region == "19" ~ "Leipzig",
                                          region == "20" ~ "Frankfurt",
                                          region == "21" ~ "Italy",
                                          region == "22" ~ "Switzerland",
                                          TRUE ~ "NO")) %>% 
  dplyr::mutate(n = ifelse(as.double(value) == 1, 17, as.double(value) - 1)) %>%
  dplyr::group_by(region, VegType, Variable, group, point, Species) %>%
  dplyr::summarise(sum.n = sum(as.double(n), na.rm = TRUE)) %>%
  #dplyr::group_by(region, VegType, Variable, group) %>%
  #dplyr::summarise(mean.n = mean(as.double(sum.n), na.rm = TRUE)) %>%
  ggplot(aes(x=region, y=sum.n, col=as.factor(region), fill=as.factor(region)))+ 
  geom_boxplot(alpha=0.3)+
  #geom_jitter(alpha=0.2)+
  facet_wrap(~Species) +
  theme_bw()+
  scale_x_discrete(labels = NULL)+
  theme()


all_att %>% #dplyr::sample_n(1000) %>%
  dplyr::select(-r_tree_spc, -ar_tree_spc, -s_tree_spc, -comments) %>% #.$stump_otsp %>% unique()
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(region, group, point, sep="_")) %>%  #dplyr::select(6:323) %>% names()
  tidyr::gather(5:322, key = "key", value="value") %>% 
  dplyr::mutate(VegType = dplyr::case_when(grepl("^(r_)", key) ~ "Regeneration",
                                           grepl("^(ar_)", key) ~ "advRegeneration",
                                           grepl("^(s_)", key) ~ "Survivor",
                                           TRUE ~ "Else")) %>% 
  #dplyr::filter(VegType != "Else") %>% #.$VegType %>% unique()
  dplyr::mutate(Variable = ifelse(VegType == "Else", ifelse(grepl("stump", key), "stump", key), gsub("(.*?)\\s*_", "", key))) %>% #.$key %>% unique()
  dplyr::mutate(Species = ifelse(VegType == "Else", ifelse(grepl("stump", key), paste0(ifelse(!is.na(value), value, "")), NA),
                                 ifelse(grepl("name", key), value, 
                                        ifelse(VegType == "advRegeneration", substr(key, 4, 7), substr(key, 3, 6))))) %>% #gsub("_\\s*(.*?)", "", 
  dplyr::select(group:Species, -key) %>%
  dplyr::filter(VegType == "advRegeneration", Variable %in% c("hgt")) %>%
  dplyr::mutate(region = dplyr::case_when(region == "11" ~ "Harz",
                                          region == "12" ~ "WGer",
                                          region == "13" ~ "Austria",
                                          region == "14" ~ "Dresden",
                                          region == "15" ~ "Czech",
                                          region == "16" ~ "Slovakia",
                                          region == "17" ~ "Poland",
                                          region == "18" ~ "MidGer",
                                          region == "19" ~ "Leipzig",
                                          region == "20" ~ "Frankfurt",
                                          region == "21" ~ "Italy",
                                          region == "22" ~ "Switzerland",
                                          TRUE ~ "NO")) %>% 
  dplyr::group_by(region, VegType, Variable, group, point) %>%
  dplyr::summarise(mean.height = mean(as.double(value), na.rm = TRUE)) %>%
  dplyr::group_by(region, VegType, Variable, group) %>%
  dplyr::summarise(mean.height = mean(as.double(mean.height), na.rm = TRUE)) %>%
  ggplot(aes(x=region, y=mean.height, col=as.factor(region)))+ 
  geom_boxplot()+
  geom_jitter(alpha=0.2)+
  #facet_wrap(~Species) +
  theme_bw()+
  scale_x_discrete(labels = NULL)+
  theme()



plot.df <- all_att %>% #dplyr::sample_n(1000) %>%
  dplyr::select(-r_tree_spc, -ar_tree_spc, -s_tree_spc, -comments) %>% #.$stump_otsp %>% unique()
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(region, group, point, sep="_")) %>%  #dplyr::select(6:323) %>% names()
  tidyr::gather(5:322, key = "key", value="value") %>% 
  dplyr::mutate(VegType = dplyr::case_when(grepl("^(r_)", key) ~ "Regeneration",
                                           grepl("^(ar_)", key) ~ "advRegeneration",
                                           grepl("^(s_)", key) ~ "Survivor",
                                           TRUE ~ "Else")) %>% 
  #dplyr::filter(VegType != "Else") %>% #.$VegType %>% unique()
  dplyr::mutate(Variable = ifelse(VegType == "Else", ifelse(grepl("stump", key), "stump", key), gsub("(.*?)\\s*_", "", key))) %>% #.$key %>% unique()
  dplyr::mutate(Species = ifelse(VegType == "Else", ifelse(grepl("stump", key), paste0(ifelse(!is.na(value), value, "")), NA),
                                 ifelse(grepl("name", key), value, 
                                        ifelse(VegType == "advRegeneration", substr(key, 4, 7), substr(key, 3, 6))))) %>% #gsub("_\\s*(.*?)", "", 
  dplyr::select(group:Species, -key) %>%
  dplyr::filter(Species=="piab"| Species=="soau"| Species=="besp" | Species=="pisy")
plot.df$VegType <- factor(plot.df$VegType, levels=c("Regeneration","advRegeneration","Survivor"))

ggplot(plot.df, aes(x=as.double(value), fill=as.factor(Species)))+
  geom_density(alpha=0.5)+
  facet_wrap(~region)
#geom_jitter(alpha=0.05, fill="black", size=2)+
#facet_grid(Species ~ VegType)

df.plott <- all_att %>%
  dplyr::select(-r_tree_spc, -ar_tree_spc, -s_tree_spc, -comments) %>% #.$stump_otsp %>% unique()
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(region, group, point, sep="_")) %>%  #dplyr::select(6:323) %>% names()
  tidyr::gather(5:322, key = "key", value="value") %>% 
  dplyr::mutate(VegType = dplyr::case_when(grepl("^(r_)", key) ~ "Regeneration",
                                           grepl("^(ar_)", key) ~ "advRegeneration",
                                           grepl("^(s_)", key) ~ "Survivor",
                                           TRUE ~ "Else")) %>% 
  #dplyr::filter(VegType != "Else") %>% #.$VegType %>% unique()
  dplyr::mutate(Variable = ifelse(VegType == "Else", ifelse(grepl("stump", key), "stump", key), gsub("(.*?)\\s*_", "", key))) %>% #.$key %>% unique()
  dplyr::mutate(Species = ifelse(VegType == "Else", ifelse(grepl("stump", key), paste0(ifelse(!is.na(value), value, "")), NA),
                                 ifelse(grepl("name", key), value, 
                                        ifelse(VegType == "advRegeneration", substr(key, 4, 7), substr(key, 3, 6))))) %>%
  dplyr::select(group:Species, -key) %>%
  #dplyr::filter(VegType == "advRegeneration", Variable %in% c("n")) %>%
  dplyr::mutate(region = dplyr::case_when(region == "11" ~ "Harz",
                                          region == "12" ~ "WGer",
                                          region == "13" ~ "Austria",
                                          region == "14" ~ "Dresden",
                                          region == "15" ~ "Czech",
                                          region == "16" ~ "Slovakia",
                                          region == "17" ~ "Poland",
                                          region == "18" ~ "MidGer",
                                          region == "19" ~ "Leipzig",
                                          region == "20" ~ "Frankfurt",
                                          region == "21" ~ "Italy",
                                          region == "22" ~ "Switzerland",
                                          TRUE ~ "NO")) %>% 
  dplyr::mutate(n = ifelse(as.double(value) == 1, 17, as.double(value) - 1)) %>%
  dplyr::mutate(n = n * 2500) %>%
  dplyr::filter(VegType=="Regeneration"|VegType=="advRegeneration", Variable=="n", Species=="piab"| Species=="besp"| Species=="soau"| Species=="pisy") 
#dplyr::mutate(n.sum = sum(as.double(value), na.rm = TRUE)) %>%
df.plott$VegType <- factor(df.plott$VegType, levels = c("Regeneration","advRegeneration"))
ggplot(df.plott, aes(x=as.factor(region), y=as.double(value)/4, fill=as.factor(VegType)))+
  geom_boxplot()+
  facet_wrap(~Species)+
  labs(y="plants n/m²", fill="Vegtation Type")

regen_plot <-att_table %>%
  dplyr::mutate(Species_str = stringr::str_length(Species))%>%
  dplyr::filter(Species!="otsp", Species_str==4, VegType!="Else")
#dplyr::filter(VegType=="Survivor")%>%
regen_plot$VegType <- factor(regen_plot$VegType, levels = c("Regeneration","advRegeneration","Survivor"))

ggplot(regen_plot, aes(x=as.factor(VegType),y=n, fill=as.factor(VegType)))+
  geom_violin(alpha=.5, width=.5)+
  #geom_jitter(col="grey", size=2 ,alpha=.2)+
  facet_wrap(~Species)

#sum(str_length(att_table$Species) == 4, na.rm=TRUE)

#att_table$Species[nchar(att_table)]

bar_plot <- att_table %>%
  dplyr::right_join(att_table %>%
                      group_by(Species, VegType) %>%
                      dplyr::summarise(temp.var = sum(n, na.rm=TRUE)) %>%
                      dplyr::mutate(temp.var = temp.var > 0) %>%
                      dplyr::group_by(Species) %>%
                      dplyr::summarise(inside = sum(temp.var, na.rm=TRUE))) %>%
  dplyr::mutate(Species_str = str_length(Species))%>%
  dplyr::filter(Species!="otsp", Species_str==4, VegType!="Else")%>%
  dplyr::filter(is.na(n)==FALSE)%>%
  dplyr::filter(inside == 3)
bar_plot$VegType <- factor(bar_plot$VegType, levels = c("Regeneration","advRegeneration","Survivor"))
ggplot(bar_plot, aes(x=n, fill=VegType))+
  geom_bar(position = position_dodge2(width=.5, preserve = "single"), alpha=.5, col="black")+
  #geom_density(alpha=.4)
  facet_wrap(~Species)


att_table %>% .$Variable %>% unique()

# mngmnt vs nat.regen

tmp.data <- all_att %>% #dplyr::sample_n(1000) %>%
  dplyr::select(-r_tree_spc, -ar_tree_spc, -s_tree_spc, -comments) %>% #.$stump_otsp %>% unique()
  dplyr::mutate(group = group + 100) %>% 
  dplyr::mutate(ID = paste(region, group, point, sep="_")) %>%  #dplyr::select(6:323) %>% names()
  tidyr::gather(5:322, key = "key", value="value") %>% 
  dplyr::mutate(VegType = dplyr::case_when(grepl("^(r_)", key) ~ "Regeneration",
                                           grepl("^(ar_)", key) ~ "advRegeneration",
                                           grepl("^(s_)", key) ~ "Survivor",
                                           TRUE ~ "Else")) %>%
  #dplyr::filter(VegType != "Else") %>% #.$VegType %>% unique()
  dplyr::mutate(Variable = ifelse(VegType == "Else", ifelse(grepl("stump", key), "stump", key), gsub("(.*?)\\s*_", "", key))) %>% #.$key %>% unique()
  dplyr::mutate(Species = ifelse(VegType == "Else", ifelse(grepl("stump", key), paste0(ifelse(!is.na(value), value, "")), NA),
                                 ifelse(grepl("name", key), value, 
                                        ifelse(VegType == "advRegeneration", substr(key, 4, 7), substr(key, 3, 6))))) %>% #gsub("_\\s*(.*?)", "", 
  dplyr::select(group:Species, -key) %>%
  dplyr::mutate(value = dplyr::case_when(Variable %in% c("planting", "anti_browsing") & value == 1 ~ NA,
                                         Variable %in% c("planting", "anti_browsing") & value == 2 ~ "true",
                                         Variable %in% c("planting", "anti_browsing") & value == 3 ~ "false",
                                         TRUE ~ value))
#dplyr::filter(Variable != "name", Variable != "stump", Variable != "pic_plot_n",Variable != "pic_plot_e",Variable != "pic_sur_n", Variable != "pic_sur_e", Variable != "pic_sur_s", Variable != "pic_sur_w", Variable != "deadwood", Variable != "windthrow") %>%

#  .$Variable %>% unique()


mngmnt_plot <- tmp.data %>%
  dplyr::filter(Variable == "n") %>% 
  dplyr::right_join(tmp.data %>%
                      dplyr::filter(Variable %in% c("clear", "grndwrk", "logging_trail", "planting", "anti_browsing")) %>% 
                      dplyr::group_by(region, group) %>%
                      dplyr::summarise(mngmnt = ifelse("true" %in% value, 1, 0)), by=c("region", "group")) %>%
  dplyr::mutate(region = dplyr::case_when(region == "11" ~ "Harz",
                                          region == "12" ~ "WGer",
                                          region == "13" ~ "Austria",
                                          region == "14" ~ "Dresden",
                                          region == "15" ~ "Czech",
                                          region == "16" ~ "Slovakia",
                                          region == "17" ~ "Poland",
                                          region == "18" ~ "MidGer",
                                          region == "19" ~ "Leipzig",
                                          region == "20" ~ "Frankfurt",
                                          region == "21" ~ "Italy",
                                          region == "22" ~ "Switzerland",
                                          TRUE ~ "NO")) %>% 
  dplyr::mutate(n = ifelse(VegType != "Survivor" , ifelse(as.double(value) == 1, 17, as.double(value) - 1), as.double(value))) %>% 
  #dplyr::filter(VegType =="Regeneration")%>%
  dplyr::group_by(region, group, mngmnt, Species, VegType) %>%
  dplyr::summarise(mean.n = mean(as.double(n), na.rm=TRUE))

#table(mngmnt_plot$mngmnt)

ggplot(mngmnt_plot, aes(x= mngmnt, y = mean.n, fill=as.factor(mngmnt)))+
  geom_violin(draw_quantiles = 0.5)+
  facet_grid(VegType~region)





# preprocess tehdata : 11/27/2023

# correct management indicatin per plot:
crosssum <- function(x){
  result <- 0
  while (x>0){
    result <- result + (x %% 10)
    x <- floor(x/10)
  }
  return(result)
}


management <- cleaned_df_tidy %>%
  dplyr::mutate(logging_trail = ifelse(is.na(logging_trail), 0, ifelse(logging_trail == "true", 1, 0))) %>%
  dplyr::mutate(clear = ifelse(is.na(clear), 0, ifelse(clear == "true", 1, 0))) %>%
  dplyr::mutate(grndwrk = ifelse(is.na(grndwrk), 0, ifelse(grndwrk == "true", 1, 0))) %>%
  dplyr::mutate(planting = ifelse(is.na(planting), 0, ifelse(planting != 2, 0, 1))) %>%
  dplyr::mutate(anti_browsing = ifelse(is.na(anti_browsing), 0, ifelse(anti_browsing != 2, 0, 1))) %>%
  dplyr::mutate(management = ifelse(logging_trail == 1, 1, 
                                    ifelse(clear == 1, 1, 
                                           ifelse(grndwrk == 1, 1, 
                                                  ifelse(planting == 1, 1, 
                                                         ifelse(anti_browsing == 1, 1, 0)))))) %>%
  dplyr::mutate(mngmnt_cat = 0) %>%
  dplyr::mutate(mngmnt_cat = ifelse(logging_trail == 1, mngmnt_cat + 1, mngmnt_cat)) %>%
  dplyr::mutate(mngmnt_cat = ifelse(clear == 1, mngmnt_cat + 10, mngmnt_cat)) %>%
  dplyr::mutate(mngmnt_cat = ifelse(grndwrk == 1, mngmnt_cat + 100, mngmnt_cat)) %>%
  dplyr::mutate(mngmnt_cat = ifelse(planting == 1, mngmnt_cat + 1000, mngmnt_cat)) %>%
  dplyr::mutate(mngmnt_cat = ifelse(anti_browsing == 1, mngmnt_cat + 10000, mngmnt_cat)) %>% #.$mngmnt_cat %>% unique()
  dplyr::rowwise() %>%
  dplyr::mutate(crosssum = crosssum(mngmnt_cat)) #%>% View(.)

