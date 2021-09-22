# reading in and processing data
#A unified framework for probabilistic forecasting subnational populations
#(C) anonnymised (2020) 

rm(list=ls())

#required libraries ####
library(dplyr)
library(tidyr)
library(forcats)
library(readr)
library(rstan)

#labels ####
agelab<-c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")
Regn=c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT")

# internation migration & population-at-risk of internal migration data ####
multiregional_pop <- read_csv("data/revised internal migrartion, population and death 1981-86 to 2006-11 30082017.csv", 
                              col_types = cols(pop1 = col_double(),
                                               pop2 = col_double()))
#pop 2016
pop16 <- read_csv("data/2016 Census pop by age_sex_state_Nan.csv", col_types = cols(age = col_factor(levels = c("0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85")), sex = col_factor(levels = c("F", "M")), state_code = col_character()))
pop16 = pop16 %>%
  rename(A=age, S=sex, Reg=state_code, x=total) %>%
  mutate(S = plyr::mapvalues(S,c(levels(S)),c("Females","Males")), year=2011) %>%
  select(-state_name)

# cleaning pop data
xpop = multiregional_pop %>%
  group_by(A,S,year,NSD) %>%
  summarise(total_cob = sum(pop1)) %>%
  ungroup() %>%
  mutate(NSD_ID=NSD) %>%
  separate(NSD,c("Reg","SD"),sep="0") %>%
  group_by(A,S,Reg,year) %>%
  summarise(x = sum(total_cob)) %>%
  ungroup() %>%
  mutate(A = factor(A, ordered = F), S = factor(S,labels = c("Males","Females"))) %>%
  bind_rows(pop16)

#Reading in data for 2011-16
od1116 <- read_csv("data/2011-2016 interstate migration flows by age and sex_Nan.csv", col_types = cols(age = col_factor(levels = c("0", "5", "10", "15", "20", "25", "30", "35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85")), sex = col_factor(levels = c("F", "M")),state_code = col_character()), skip = 3)

od1116 =od1116 %>% rename(Reg_O=state_code, S=sex, A=age) %>%
  gather(Reg_D,Flow,D1:D8) %>%
  mutate(Reg_D = substr(Reg_D,start = 2,stop = 2),
         Corridor = paste0(Reg_O,"->",Reg_D),
         year=as.integer(2011),
         S = plyr::mapvalues(S,c(levels(S)),c("Females","Males"))) %>%
  filter(Reg_O != Reg_D, A != "0") %>%
  droplevels() %>%
  select(-stayers, - state_name)

#multiregional migration ####
xinout = multiregional_pop %>%
  group_by(A,S,year,NSD) %>%
  summarise_at(vars(D101:O801),sum) %>%
  ungroup() %>%
  mutate(NSD_ID=NSD) %>%
  separate(NSD,c("Reg","SD"),sep="0") %>%
  group_by(A,S,Reg,year) %>%
  summarise_at(vars(D101:O801),sum) %>%
  ungroup() %>%
  mutate(
    D_1=rowSums(.[grep("D1",names(.))]),
    D_2=rowSums(.[grep("D2",names(.))]),
    D_3=rowSums(.[grep("D3",names(.))]),
    D_4=rowSums(.[grep("D4",names(.))]),
    D_5=rowSums(.[grep("D5",names(.))]),
    D_6=rowSums(.[grep("D6",names(.))]),
    D_7=rowSums(.[grep("D7",names(.))]),
    D_8=D801,
    O_1=rowSums(.[grep("O1",names(.))]),
    O_2=rowSums(.[grep("O2",names(.))]),
    O_3=rowSums(.[grep("O3",names(.))]),
    O_4=rowSums(.[grep("O4",names(.))]),
    O_5=rowSums(.[grep("O5",names(.))]),
    O_6=rowSums(.[grep("O6",names(.))]),
    O_7=rowSums(.[grep("O7",names(.))]),
    O_8=O801  
  ) %>% 
  select(A:year,D_1:D_8) %>%
  rename(Reg_O = Reg) %>%
  # names(xinout)[3] <- "Reg_O"
  # xinout = xinout %>%
  # select(1:12) %>%
  gather(Reg_D,Flow,D_1:D_8) %>%
  mutate(Reg_D = gsub(".*_","",Reg_D),
         Corridor = paste0(Reg_O,"->",Reg_D),
         A = factor(A,ordered=F),
         S = factor(S, labels = c("Males","Females"))) %>%
  filter(Reg_O != Reg_D, A != 0) %>%
  droplevels() %>%
  select(Reg_O,Reg_D,Corridor,A,S,year,Flow) %>%
  bind_rows(od1116) %>%
  arrange(A,S,year,Reg_O,Reg_D) %>%
  mutate(unicor=(Reg_O<Reg_D)*1) %>%
  mutate(uniflow=case_when(unicor==1~paste0(Reg_O,Reg_D), unicor==0 ~ paste0(Reg_D,Reg_O)))



# pop for pop-at-risk ####
xpop1 = xpop %>%
  filter(A!=85) %>%#,year!=2011
  droplevels() %>%
  mutate(A = plyr::mapvalues(A,c(levels(A)),levels(xinout$A))) %>%
  rename(Reg_O = Reg) 

input_au = left_join(xinout,xpop1)
rm(xinout,xpop1,xpop,multiregional_pop)
corridor=as.numeric(input_au$Reg_D[1:56])
#creating index for unidirectional flows
uniflow=as.numeric(input_au$uniflow)
unif_ind=1
for (i in 2:56) unif_ind[i]=if(prod(uniflow[i]!= uniflow[1:(i-1)])>0) max(unif_ind)+1 else unif_ind[which(uniflow[1:(i-1)]==uniflow[i])]
rm(i)


# ERP population at risk and baselines for event data ####
ERP_pop <- read_csv("data/ERP_1981-2016.csv",
                    col_types = cols(age = col_factor(levels = c("0",
                                                                 "5", "10", "15", "20", "25", "30",
                                                                 "35", "40", "45", "50", "55", "60",
                                                                 "65", "70", "75", "80", "85")), 
                                     sex = col_factor(levels = c("F", "M")), 
                                     state = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8"))))


ERP_pop = ERP_pop %>%
  mutate(A = plyr::mapvalues(age,c(levels(age)),agelab)) %>%
  rename(S = sex, Reg = state) 

#deaths data ####
load("data/death full.RData")
deaths9311 = dat %>%
  group_by(A,S,year,NSD) %>%
  summarise(total_cob = sum(x)) %>%
  ungroup() %>%
  mutate(NSD_ID=NSD) %>%
  separate(NSD,c("Reg","SD"),sep="0") %>%
  group_by(A,S,Reg,year) %>%
  summarise(x = sum(total_cob)) %>%
  ungroup() %>%
  filter(A!="TT") %>%
  droplevels() %>%
  mutate(A = plyr::mapvalues(A,c(levels(A)),agelab)) %>%
  mutate(Reg = as_factor(Reg)) %>% 
  left_join(ERP_pop) %>%
  arrange(A,S,year,Reg)
rm(dat,SD)

# births ####
births_8116 <- read_csv("data/births_state_1981_2016.csv", col_types = cols(age = col_factor(levels = c("15", "20", "25", "30", "35", "40", "45")), state = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8"))))

births_8116 = births_8116 %>%
  mutate(A = plyr::mapvalues(age,c(levels(age)),agelab[4:10])) %>%
  rename(Reg = state) %>%
  left_join(filter(ERP_pop,A%in%agelab[4:10], S=="F") %>% droplevels()) %>%
  droplevels() %>%
  arrange(A,S,year,Reg) 

# international migration ####
intmigr_8116 <- read_csv("data/immigration and emigration_state_1981_2016.csv", 
                         col_types = cols(age = col_factor(levels = c("0", "5", "10", "15", "20", "25", "30","35", "40", "45", "50", "55", "60", "65", "70", "75", "80", "85")), sex = col_factor(levels = c("F", "M")), state = col_factor(levels = c("1", "2", "3", "4", "5", "6", "7", "8"))))

intmigr_8116 = intmigr_8116 %>%
  mutate(A = plyr::mapvalues(age,c(levels(age)),agelab)) %>%
  rename(S = sex, Reg = state) %>%
  left_join(ERP_pop) %>%
  # mutate(Reg=plyr::mapvalues(Reg,levels(Reg),Regn)) %>% #comment for multireg
  arrange(A,S,year,Reg)
# intmigr_8116 = filter(intmigr_8116,year<2012)

