# The following script is adapted from the code repository accompanying:
# Tebaldi et al. (2021) Nature Climate Change
# https://github.com/DOE-ICoM/tebaldi-etal_2021_natclimchange

# Original authors: Tebaldi et al. (2021)
# Modifications: Khin Nawarat (2026, @IHE Delft)

remove(list=ls())
filedir<-"D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/" # change to your local directory
source(paste0(filedir,"Rcode/allfunctions.r"))


rl1k<-read.csv(filedir%&%"CSV/Kirezci/Kirezci_ESLs.csv")
rl1r<-read.csv(filedir%&%"CSV/Rasmussen/Rasmussen_ESLs.csv")
rl1v<-read.csv(filedir%&%"CSV/Vousdoukas/Vousdoukas_ESLs.csv")

for(whichdecade in seq(2020,2100,by=10)){
  proj1k<-array(dim=c(11,6,8086,2)) 
  proj1r<-array(dim=c(11,6,200,2))
  proj1v<-array(dim=c(11,6,7472,2))


dimnames(proj1k)<-list(c("ID" ,    "year" ,  "q0.01" , "q0.05",  "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99",  "q0.995", "q0.999"),
                       rasmussenScenarios,NULL,c("rasmussen","bvw"))
dimnames(proj1r)<-list(c("ID" ,    "year" ,  "q0.01" , "q0.05",  "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99",  "q0.995", "q0.999"),
                       rasmussenScenarios,NULL,c("rasmussen","bvw"))
dimnames(proj1v)<-list(c("ID" ,    "year" ,  "q0.01" , "q0.05",  "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99",  "q0.995", "q0.999"),
                       rasmussenScenarios,NULL,c("rasmussen","bvw"))



for(scenario in rasmussenScenarios){
  for(proj in c("rasmussen","bvw")){
    na=filedir%&%"Input_data/Kirezci/Kirezci_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv"
    print(na)
    temp<-read.csv(file=filedir%&%"Input_data/Kirezci/Kirezci_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv")
    #temp = subset(temp, select = -c(index) )
    proj1k[,scenario,,proj]<-t(temp)
    na=filedir%&%"Input_data/Kirezci/Kirezci_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv"
    print(na)
    temp<-read.csv(file=filedir%&%"Input_data/Rasmussen/Rasmussen_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv")
    #temp = subset(temp, select = -c(index) )
    proj1r[,scenario,,proj]<-t(temp)
    na=filedir%&%"Input_data/Kirezci/Kirezci_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv"
    print(na)
    temp<-read.csv(file=filedir%&%"Input_data/Vousdoukas/Vousdoukas_SLR_"%&%scenario%&%"_"%&%proj%&%"_"%&%whichdecade%&%".csv")
    #temp = subset(temp, select = -c(index) )
    proj1v[,scenario,,proj]<-t(temp)
  }
}

save(list=c("rl1k","rl1r","rl1v","proj1k","proj1r","proj1v"),file=filedir%&%"Rdatasets/RData_esl_slr_"%&%whichdecade)
}
