# The following script is adapted from the code repository accompanying:
# Tebaldi et al. (2021) Nature Climate Change
# https://github.com/DOE-ICoM/tebaldi-etal_2021_natclimchange

# Modifications: Khin Nawarat (2026, @IHE Delft)

remove(list=ls())
filedir<-"D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/" # change to your local directory
source(paste0(filedir,"Rcode/allfunctions.r"))

whichdecade<-2100
for (RL in c(10,20,50,100,200)){

load(filedir%&%"Rdatasets/RData_esl_slr_2100")
load(filedir%&%"Rdatasets/RData_convolutions_fisherinfo_"%&%RL%&%"yr_"%&%whichdecade)


eliminate.rasmussen<-numeric(0)
for(i in 1:200)if(all(is.na(rasmussen.shape.scale.TWL.samples[i,,])))eliminate.rasmussen<-c(eliminate.rasmussen,i)
if(length(eliminate.rasmussen)>0)rasmussensubsetID<-seq(length(dimnames(rasmussen.shape.scale.TWL.samples)[[1]]))[-eliminate.rasmussen]

eliminate.kirezci<-numeric(0)
for(i in 1:8086)if(all(is.na(kirezci.shape.scale.TWL.samples[i,,])))eliminate.kirezci<-c(eliminate.kirezci,i)
if(length(eliminate.kirezci)>0)kirezcisubsetID<-seq(length(dimnames(kirezci.shape.scale.TWL.samples)[[1]]))[-eliminate.kirezci]

eliminate.vousdoukas<-numeric(0)
for(i in 1:7472)if(all(is.na(vousdoukas.shape.scale.TWL.samples[i,,])))eliminate.vousdoukas<-c(eliminate.vousdoukas,i)
if(length(eliminate.vousdoukas)>0)vousdoukassubsetID<-seq(length(dimnames(vousdoukas.shape.scale.TWL.samples)[[1]]))[-eliminate.vousdoukas]
if(length(eliminate.vousdoukas)==0)vousdoukassubsetID<-seq(length(dimnames(vousdoukas.shape.scale.TWL.samples)[[1]]))



rmatched<-rl1r[rasmussensubsetID,]
kmatched<-rl1k[kirezcisubsetID,]
vmatched<-rl1v[vousdoukassubsetID,]

rfreq<-rasmussen.freq.probproj[,,rasmussensubsetID,]
kfreq<-kirezci.freq.probproj[,,kirezcisubsetID,]
vfreq<-vousdoukas.freq.probproj[,,vousdoukassubsetID,]

rasmussenScenarios<-c("1p5degree","2p0degree","2p5degree","3p0degree","4p0degree","5p0degree")


rasmussentokirezci<-find.match(rmatched,kmatched)
rasmussentovousdoukas<-find.match(rmatched,vmatched)
dim(rasmussentokirezci)  #187 2
dim(rasmussentovousdoukas)  #183 2

tt<-table(c(rasmussentokirezci[,1],rasmussentovousdoukas[,1]))
common<-as.numeric(names(tt)[tt==2])
rasmussentokirezciandvousdoukas<-cbind(common,rasmussentokirezci[,2][match(common,rasmussentokirezci[,1])],rasmussentovousdoukas[,2][match(common,rasmussentovousdoukas[,1])])



commonlocs.kirezci<-rasmussentokirezciandvousdoukas[,2]
commonlocs.vousdoukas<-rasmussentokirezciandvousdoukas[,3]
commonlocs.rasmussen<-rasmussentokirezciandvousdoukas[,1]

qq<-"q0.5"
whenone<-array(0,dim=c(length(commonlocs.kirezci),3,2)) # assigned value, rows of common locs, 3 columns, 2 matrices
dimnames(whenone)<-list(NULL,c("RP1r","RP1k","RP1v"),c("rasmussen","bvw"))

for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locr<-commonlocs.rasmussen[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){

        tempr<-rfreq[qq,,locr,slrm]
        tempk<-kfreq[qq,,lock,slrm]
        tempv<-vfreq[qq,,locv,slrm]


        whenone[ii,"RP1r",slrm]<-seq(length(tempr)-1)[tempr[-1]==1][1]
        whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
        whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]
    }

}


whenone[is.na(whenone)]<-7


label.100to1.q0.5<-vote3.label.100to1.q0.5<-numeric(length(commonlocs.kirezci))



for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]
    tt<-table(temp)
    whichtt<-as.numeric(names(tt)[cumsum(tt)>3])
    vote3.label.100to1.q0.5[ii]<-label.100to1.q0.5[ii]<-ifelse(length(whichtt)>0,whichtt,median(temp))

}

###q0.05

qq<-"q0.05"
whenone<-array(0,dim=c(length(commonlocs.kirezci),3,2))
dimnames(whenone)<-list(NULL,c("RP1r","RP1k","RP1v"),c("rasmussen","bvw"))





for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locr<-commonlocs.rasmussen[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){

        tempr<-rfreq[qq,,locr,slrm]
        tempk<-kfreq[qq,,lock,slrm]
        tempv<-vfreq[qq,,locv,slrm]

        whenone[ii,"RP1r",slrm]<-seq(length(tempr)-1)[tempr[-1]==1][1]
        whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
        whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]

         }
}

whenone[is.na(whenone)]<-7



vote3.label.100to1.q0.05<-label.100to1.q0.05<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]

    vote3.label.100to1.q0.05[ii]<-label.100to1.q0.05[ii]<-min(temp)

}



##q0.95
qq<-"q0.95"
whenone<-array(0,dim=c(length(commonlocs.kirezci),3,2))
dimnames(whenone)<-list(NULL,c("RP1r","RP1k","RP1v"),c("rasmussen","bvw"))

for(slrm in c("rasmussen","bvw")){
  
  tempr<-rfreq[qq,,locr,slrm]
  tempk<-kfreq[qq,,lock,slrm]
  tempv<-vfreq[qq,,locv,slrm]
  
  whenone[ii,"RP1r",slrm]<-seq(length(tempr)-1)[tempr[-1]==1][1]
  whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
  whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]
  
}


for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locr<-commonlocs.rasmussen[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){

        tempr<-rfreq[qq,,locr,slrm]
        tempk<-kfreq[qq,,lock,slrm]
        tempv<-vfreq[qq,,locv,slrm]

        whenone[ii,"RP1r",slrm]<-seq(length(tempr)-1)[tempr[-1]==1][1]
        whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
        whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]

         }
}

whenone[is.na(whenone)]<-7


tempr<-rfreq["q0.95",,1,"bvw"]




vote3.label.100to1.q0.95<-label.100to1.q0.95<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]

    vote3.label.100to1.q0.95[ii]<-label.100to1.q0.95[ii]<-max(temp)

}

gridcoo.commonlocs<-rmatched[,c("longitude","latitude")][commonlocs.rasmussen,]
outfile <- filedir%&%"Rdatasets/common3grids.csv"

if (!file.exists(outfile)) {
  write.csv(gridcoo.commonlocs, outfile, row.names = FALSE)
}


ll<-length(label.100to1.q0.5)
round(table(label.100to1.q0.5)/ll,dig=4)
#100 to 1, q=0.5
#1     2     3     4     6     7     8     9
#0.5419 0.1117 0.0279 0.0056 0.0112 0.0056 0.0950 0.2011

round(table(label.100to1.q0.05)/ll,dig=4)
#100 to 1, q=0.05

#    1      2      7
#0.9888 0.0056 0.0056

round(table(label.100to1.q0.95)/ll,dig=4)
#100 to 1, q=0.95

#     1      2      3      6      7      8      9
#0.0223 0.0112 0.0503 0.0223 0.0279 0.0615 0.8045


####now only match vousdoukas and kirezci

#only consider rows that have been sampled successfully

vousdoukastokirezci<-find.match(vmatched,kmatched)

commonlocs.kirezci<-vousdoukastokirezci[,2]
commonlocs.vousdoukas<-vousdoukastokirezci[,1]

qq<-"q0.5"
whenone<-array(0,dim=c(length(commonlocs.kirezci),2,2))
dimnames(whenone)<-list(NULL,c("RP1k","RP1v"),c("rasmussen","bvw"))



for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){


        tempk<-kfreq[qq,,lock,slrm]
        tempv<-vfreq[qq,,locv,slrm]


        whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]



        whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]

}
}

whenone[is.na(whenone)]<-7


vote2.label.100to1.q0.5<-label.100to1.q0.5<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]
    tt<-table(temp)
    whichtt<-as.numeric(names(tt)[cumsum(tt)>2])
    vote2.label.100to1.q0.5[ii]<-label.100to1.q0.5[ii]<-ifelse(length(whichtt)>0,whichtt,median(temp))


}

###q0.05

qq<-"q0.05"
whenone<-array(0,dim=c(length(commonlocs.kirezci),2,2))
dimnames(whenone)<-list(NULL,c("RP1k","RP1v"),c("rasmussen","bvw"))



for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]

    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){


                tempk<-kfreq[qq,,lock,slrm]
                tempv<-vfreq[qq,,locv,slrm]


                whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]



                whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]

            }
        }



whenone[is.na(whenone)]<-7



vote2.label.100to1.q0.05<-label.100to1.q0.05<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]

    vote2.label.100to1.q0.05[ii]<-label.100to1.q0.05[ii]<-min(temp)
    }



##q0.95


qq<-"q0.95"
whenone<-array(0,dim=c(length(commonlocs.kirezci),2,2))
dimnames(whenone)<-list(NULL,c("RP1k","RP1v"),c("rasmussen","bvw"))


for(ii in seq(length(commonlocs.kirezci))){

    lock<-commonlocs.kirezci[ii]
    locv<-commonlocs.vousdoukas[ii]


    for(slrm in c("rasmussen","bvw")){


                tempk<-kfreq[qq,,lock,slrm]
                tempv<-vfreq[qq,,locv,slrm]


                    whenone[ii,"RP1k",slrm]<-seq(length(tempk)-1)[tempk[-1]==1][1]
                    whenone[ii,"RP1v",slrm]<-seq(length(tempv)-1)[tempv[-1]==1][1]



            }
        }


whenone[is.na(whenone)]<-7


label.100to1.q0.95<-vote2.label.100to1.q0.95<-numeric(length(commonlocs.kirezci))


for(ii in seq(length(commonlocs.kirezci))){

    temp<-whenone[ii,,]

    vote2.label.100to1.q0.95[ii]<-label.100to1.q0.95[ii]<-max(temp)

}


gridcoo.commonlocs<-kmatched[,c("longitude","latitude")][commonlocs.kirezci,]
outfile <- filedir%&%"Rdatasets/common2grids.csv"

if (!file.exists(outfile)) {
  write.csv(gridcoo.commonlocs, outfile, row.names = FALSE)
}


ll<-length(label.100to1.q0.5)
round(table(label.100to1.q0.5)/ll,dig=4)
#100 to 1, q=0.5
#     1      2      3      4      5      6      7      8      9
#0.4333 0.1026 0.0383 0.0043 0.0008 0.0373 0.0319 0.1178 0.2337


round(table(label.100to1.q0.05)/ll,dig=4)
#100 to 1, q=0.05
#     1      2      3      4      6      7      8
# 0.9842 0.0100 0.0019 0.0001 0.0008 0.0023 0.0005


round(table(label.100to1.q0.95)/ll,dig=4)
#100 to 1, q=0.95
#    1      2      3      4      5      6      7      8      9
# 0.0736 0.0154 0.0685 0.0019 0.0003 0.0144 0.0519 0.0943 0.6797

## results based on VS method are collected at this point, now collect the results from FC
source(paste0(filedir,"Rcode/convolve_alltheway_fordatapub_allRPs.R")) ## this line go run the script based on FC method: "convolve_alltheway_fordatapub_allRPs.R"

#now compare distributions after running the full convolution code and producing all3/all2
###one single plot:
res <- convolve_alltheway_fordatapub(
  RL = RL,
  whichdecade = whichdecade,
  filedir = filedir
)

assign(paste("V3.",RL,"to1.q0.5", sep = ""), vote3.label.100to1.q0.5)
assign(paste("all3.",RL,"to1.q0.5", sep = ""), all3.label.100to1.q0.5)
assign(paste("V3.",RL,"to1.q0.05", sep = ""), vote3.label.100to1.q0.05)
assign(paste("all3.",RL,"to1.q0.05", sep = ""), all3.label.100to1.q0.05)
assign(paste("V3.",RL,"to1.q0.95", sep = ""), vote3.label.100to1.q0.95)
assign(paste("all3.",RL,"to1.q0.95", sep = ""), all3.label.100to1.q0.95)

assign(paste("V2.",RL,"to1.q0.5", sep = ""), vote2.label.100to1.q0.5)
assign(paste("all2.",RL,"to1.q0.5", sep = ""), all2.label.100to1.q0.5)
assign(paste("V2.",RL,"to1.q0.05", sep = ""), vote2.label.100to1.q0.05)
assign(paste("all2.",RL,"to1.q0.05", sep = ""), all2.label.100to1.q0.05)
assign(paste("V2.",RL,"to1.q0.95", sep = ""), vote2.label.100to1.q0.95)
assign(paste("all2.",RL,"to1.q0.95", sep = ""), all2.label.100to1.q0.95)

assign(paste0("all3.", RL, "to1.q0.05"), res$all3$q0.05)
assign(paste0("all3.", RL, "to1.q0.5"),  res$all3$q0.5)
assign(paste0("all3.", RL, "to1.q0.95"), res$all3$q0.95)

assign(paste0("all2.", RL, "to1.q0.05"), res$all2$q0.05)
assign(paste0("all2.", RL, "to1.q0.5"),  res$all2$q0.5)
assign(paste0("all2.", RL, "to1.q0.95"), res$all2$q0.95)


save(list=c(paste("V3.",RL,"to1.q0.5",sep = ""),paste("V3.",RL,"to1.q0.05",sep = ""),paste("V3.",RL,"to1.q0.95",sep = ""),paste("V2.",RL,"to1.q0.5",sep = ""),paste("V2.",RL,"to1.q0.05",sep = ""),paste("V2.",RL,"to1.q0.95",sep = ""),
            paste("all3.",RL,"to1.q0.5",sep = ""),paste("all3.",RL,"to1.q0.05",sep = ""),paste("all3.",RL,"to1.q0.95",sep = ""),paste("all2.",RL,"to1.q0.5",sep = ""),paste("all2.",RL,"to1.q0.05",sep = ""),paste("all2.",RL,"to1.q0.95",sep = "")),
     file=filedir%&%"Rdatasets/RData_labels_"%&%RL%&%"to1_FC_VS_both_fisherinfo_2100")
}
#######################################################################################
