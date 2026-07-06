# The following script is adapted from the code repository accompanying:
# Tebaldi et al. (2021) Nature Climate Change
# https://github.com/DOE-ICoM/tebaldi-etal_2021_natclimchange

# Modifications: Khin Nawarat (2026, @IHE Delft)

# this function is to be called from "votingsystem_fordatapub_allRPs.R"

convolve_alltheway_fordatapub <- function (RL, whichdecade, filedir) {

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

rasmussentokirezci<-find.match(rmatched,kmatched)
rasmussentovousdoukas<-find.match(rmatched,vmatched)
dim(rasmussentokirezci)  #187 2
dim(rasmussentovousdoukas)  #183 2

tt<-table(c(rasmussentokirezci[,1],rasmussentovousdoukas[,1]))
common<-as.numeric(names(tt)[tt==2])
rasmussentokirezciandvousdoukas<-cbind(common,rasmussentokirezci[,2][match(common,rasmussentokirezci[,1])],rasmussentovousdoukas[,2][match(common,rasmussentovousdoukas[,1])])
###############################################################################################################################################################################################


qs<-c(0.010, 0.050, 0.167, 0.500, 0.833, 0.950, 0.990, 0.995, 0.999)
rasmussenScenarios<-c("1p5degree","2p0degree","2p5degree","3p0degree","4p0degree","5p0degree")
vMC<-vousdoukas.shape.scale.TWL.samples[vousdoukassubsetID,,]
kMC<-kirezci.shape.scale.TWL.samples[kirezcisubsetID,,]
rMC<-rasmussen.shape.scale.TWL.samples[rasmussensubsetID,,]


load(paste0(filedir%&%"Rdatasets/RData_esl_slr_",whichdecade))
subproj1v<-proj1v[,,vousdoukassubsetID,]
subproj1k<-proj1k[,,kirezcisubsetID,]  #get these "subsetID" from convolve_fisherinfo.r
subproj1r<-proj1r[,,rasmussensubsetID,]


commonlocs.kirezci<-rasmussentokirezciandvousdoukas[,2]   #get this from convolve_fisherinfo.r
commonlocs.vousdoukas<-rasmussentokirezciandvousdoukas[,3]
commonlocs.rasmussen<-rasmussentokirezciandvousdoukas[,1]


all3f.int.probproj<-all3f.freq.probproj<-array(dim=c(11,7,length(commonlocs.kirezci)))
dimnames(all3f.int.probproj)<-dimnames(all3f.freq.probproj)<-list(c("mean","sd","q0.01","q0.05" , "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99" , "q0.995", "q0.999"),
                                    c("current",rasmussenScenarios),NULL)

locv<-commonlocs.vousdoukas[1]
eslv<-vmatched[locv,]    #these "matched" are from convolve_fisherinfo.r
eslsamplev<-100*vMC[locv,,"TWL100"]  #randomized sample from shape/scale generation
scalev<-vMC[locv,,"scale"]
shapev<-vMC[locv,,"shape"]
rpscurrentv<-RP(eslsamplev/100,shapev,scalev,eslv$thresh,eslv$lambda)

set.seed(123)
for(ii in seq(length(commonlocs.kirezci))){
    lock<-commonlocs.kirezci[ii]
    locr<-commonlocs.rasmussen[ii]
    locv<-commonlocs.vousdoukas[ii]
    print(ii)

    eslv<-vmatched[locv,]    #these "matched" are from convolve_fisherinfo.r
    eslsamplev<-100*vMC[locv,,"TWL100"]  #randomized sample from shape/scale generation
    scalev<-vMC[locv,,"scale"]
    shapev<-vMC[locv,,"shape"]
    rpscurrentv<-RP(eslsamplev/100,shapev,scalev,eslv$thresh,eslv$lambda)

    eslk<-kmatched[lock,]
    eslsamplek<-100*kMC[lock,,"TWL100"]  #randomized sample from shape/scale generation
    scalek<-kMC[lock,,"scale"]
    shapek<-kMC[lock,,"shape"]
    rpscurrentk<-RP(eslsamplek/100,shapek,scalek,eslk$thresh,eslk$lambda)

    eslr<-rmatched[locr,]
    eslsampler<-100*rMC[locr,,"TWL100"]  #randomized sample from a normal approximation to 100-yr event
    scaler<-rMC[locr,,"scale"]
    shaper<-rMC[locr,,"shape"]
    rpscurrentr<-RP(eslsampler/100,shaper,scaler,eslr$thresh,eslr$lambda)


    eslsample<-c(eslsamplev,eslsamplek,eslsampler)

    rpscurrent<-c(rpscurrentv,rpscurrentk,rpscurrentr)

        all3f.int.probproj[,"current",ii]<-c(mean(eslsample),sqrt(var(eslsample)),
                                                           quantile(eslsample,prob=qs-qs[1],na.rm=TRUE))
        all3f.freq.probproj[,"current",ii]<-c(mean(rpscurrent),sqrt(var(rpscurrent)),
                                                            quantile(rpscurrent,prob=qs-qs[1],na.rm=TRUE))
            for(scenario in rasmussenScenarios){

                for(slrmethod in c("rasmussen", "bvw")){
                slr.qs<-subproj1v[-c(1,2),scenario,locv,slrmethod] #-c(1,2) means removing first two rows. Getting quantiles and their SLR values
                if(!all(is.na(slr.qs))){
                    slrsample<-unfold.quantiles(slr.qs)

                    ll<-length(slrsample)
                    slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                    eslsamplev<-eslsamplev[1:ll]


                    sssample<--slrsample+eslsamplev
                    rpssamplev<-RP(sssample/100,shapev[1:ll],scalev[1:ll],eslv$thresh,eslv$lambda)
                    sssamplev<-slrsample+eslsamplev
                    }
                else{
                    rpssamplev<-sssamplev<-rep(NA,ll)
                }


                slr.qs<-subproj1k[-c(1,2),scenario,lock,slrmethod]
                if(!all(is.na(slr.qs))){
                    slrsample<-unfold.quantiles(slr.qs)

                    ll<-length(slrsample)
                    slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                    eslsamplek<-eslsamplek[1:ll]


                    sssample<--slrsample+eslsamplek
                    rpssamplek<-RP(sssample/100,shapek[1:ll],scalek[1:ll],eslk$thresh,eslk$lambda)
                    sssamplek<-slrsample+eslsamplek
                    }
                else{
                    rpssamplek<-sssamplek<-rep(NA,ll)
                }

                slr.qs<-subproj1r[-c(1,2),scenario,locr,slrmethod]
                if(!all(is.na(slr.qs))){
                    slrsample<-unfold.quantiles(slr.qs)

                    ll<-length(slrsample)
                    slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                    eslsampler<-eslsampler[1:ll]


                    sssample<--slrsample+eslsampler
                    rpssampler<-RP(sssample/100,shaper[1:ll],scaler[1:ll],eslr$thresh,eslr$lambda)
                    sssampler<-slrsample+eslsampler
                    }
                else{
                    rpssampler<-sssampler<-rep(NA,ll)
                }

                if(slrmethod=="rasmussen"){
                    rpssample<-c(rpssamplev,rpssamplek,rpssampler)
                sssample<-c(sssamplev,sssamplek,sssampler)}

                else{
                    rpssample<-c(rpssample,rpssamplev,rpssamplek,rpssampler)
                    sssample<-c(sssample,sssamplev,sssamplek,sssampler)}
                }

                all3f.freq.probproj[,scenario,ii]<-c(mean(rpssample),sqrt(var(rpssample)),
                                                              quantile(rpssample,prob=qs-qs[1],na.rm=TRUE))

                all3f.int.probproj[,scenario,ii]<-c(mean(sssample),sqrt(var(sssample)),
                                                              quantile(sssample,prob=qs-qs[1],na.rm=TRUE))
            }
}


vousdoukastokirezci<-find.match(vmatched,kmatched)

commonlocs.kirezci<-vousdoukastokirezci[,2]
commonlocs.vousdoukas<-vousdoukastokirezci[,1]


all2f.int.probproj<-all2f.freq.probproj<-array(dim=c(11,7,length(commonlocs.kirezci)))
dimnames(all2f.int.probproj)<-dimnames(all2f.freq.probproj)<-list(c("mean","sd","q0.01","q0.05" , "q0.167", "q0.5",   "q0.833", "q0.95",  "q0.99" , "q0.995", "q0.999"),
                                                                c("current",rasmussenScenarios),NULL)

set.seed(123)
for(ii in seq(length(commonlocs.kirezci))){
    lock<-commonlocs.kirezci[ii]
    locv<-commonlocs.vousdoukas[ii]
    print(ii)

    eslv<-vmatched[locv,]
    eslsamplev<-100*vMC[locv,,"TWL100"]  #randomized sample from shape/scale generation
    scalev<-vMC[locv,,"scale"]
    shapev<-vMC[locv,,"shape"]
    rpscurrentv<-RP(eslsamplev/100,shapev,scalev,eslv$thresh,eslv$lambda)

    eslk<-kmatched[lock,]
    eslsamplek<-100*kMC[lock,,"TWL100"]  #randomized sample from shape/scale generation
    scalek<-kMC[lock,,"scale"]
    shapek<-kMC[lock,,"shape"]
    rpscurrentk<-RP(eslsamplek/100,shapek,scalek,eslk$thresh,eslk$lambda)

    eslsample<-c(eslsamplev,eslsamplek)

    rpscurrent<-c(rpscurrentv,rpscurrentk)

    all2f.int.probproj[,"current",ii]<-c(mean(eslsample),sqrt(var(eslsample)),
                                                      quantile(eslsample,prob=qs-qs[1],na.rm=TRUE))
    all2f.freq.probproj[,"current",ii]<-c(mean(rpscurrent),sqrt(var(rpscurrent)),
                                                       quantile(rpscurrent,prob=qs-qs[1],na.rm=TRUE))

        for(scenario in rasmussenScenarios){


            for(slrmethod in c("rasmussen","bvw")){

            slr.qs<-subproj1v[-c(1,2),scenario,locv,slrmethod]
            if(!all(is.na(slr.qs))){
                slrsample<-unfold.quantiles(slr.qs)

                ll<-length(slrsample)
                slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                eslsamplev<-eslsamplev[1:ll]


                sssample<--slrsample+eslsamplev
                rpssamplev<-RP(sssample/100,shapev[1:ll],scalev[1:ll],eslv$thresh,eslv$lambda)
                sssamplev<-slrsample+eslsamplev
            }
            else{
                rpssamplev<-sssamplev<-rep(NA,ll)
            }


            slr.qs<-subproj1k[-c(1,2),scenario,lock,slrmethod]
            if(!all(is.na(slr.qs))){
                slrsample<-unfold.quantiles(slr.qs)

                ll<-length(slrsample)
                slrsample<-sample(slrsample, size=ll)  #randomized sample from the quantiles of the SLR projection
                eslsamplek<-eslsamplek[1:ll]


                sssample<--slrsample+eslsamplek
                rpssamplek<-RP(sssample/100,shapek[1:ll],scalek[1:ll],eslk$thresh,eslk$lambda)
                sssamplek<-slrsample+eslsamplek
            }
            else{
                rpssamplek<-sssamplek<-rep(NA,ll)
            }


if(slrmethod=="rasmussen"){
            rpssample<-c(rpssamplev,rpssamplek)
            sssample<-c(sssamplev,sssamplek)}
            else{
                rpssample<-c(rpssample,rpssamplev,rpssamplek)
                sssample<-c(sssamplev,sssamplek)}

}
            all2f.freq.probproj[,scenario,ii]<-c(mean(rpssample),sqrt(var(rpssample)),
                                                          quantile(rpssample,prob=qs-qs[1],na.rm=TRUE))

            all2f.int.probproj[,scenario,ii]<-c(mean(sssample),sqrt(var(sssample)),
                                                         quantile(sssample,prob=qs-qs[1],na.rm=TRUE))
        }
}



save(list=c(objects(pattern="all3f"),objects(pattern="all2f")),file=paste0(filedir,"Rdatasets/RData_fullconvolutions_fisherinfo_",RL,"yr_",whichdecade))



whenone<-matrix(0,nrow(rasmussentokirezciandvousdoukas),3)
dimnames(whenone)<-list(NULL,c("q0.05","q0.5","q0.95"))


for(ii in seq(nrow(rasmussentokirezciandvousdoukas))){

for(qq in c("q0.05","q0.5","q0.95")){


        temp<-all3f.freq.probproj[qq,,ii]

        whenone[ii,qq]<-seq(length(temp)-1)[temp[-1]==1][1]


     }
}

whenone[is.na(whenone)]<-7



all3.label.100to1.q0.05<-label.100to1.q0.05<-whenone[,"q0.05"]

all3.label.100to1.q0.5<-label.100to1.q0.5<-whenone[,"q0.5"]

all3.label.100to1.q0.95<-label.100to1.q0.95<-whenone[,"q0.95"]



ll<-length(label.100to1.q0.05)
round(table(label.100to1.q0.5)/ll,dig=4)
#100 to 1, q=0.5
#     1      2      3      4      5      6      7      9

#0.5810 0.1285 0.0838 0.0279 0.0335 0.0894 0.0279 0.0279
round(table(label.100to1.q0.05)/ll,dig=4)
#100 to 1, q=0.05
#     1      2      7

#0.9888 0.0056 0.0056
round(table(label.100to1.q0.95)/ll,dig=4)
#100 to 1, q=0.95
#     1      2      3      4      5      6      7      8      9

#0.1006 0.0279 0.0670 0.0056 0.0223 0.0894 0.0279 0.0503 0.6089


assign(paste0("all3.whenone.",whichdecade),whenone)

####now only match vousdoukas and kirezci

whenone<-matrix(0,nrow(vousdoukastokirezci),3)
dimnames(whenone)<-list(NULL,c("q0.05","q0.5","q0.95"))


for(ii in seq(nrow(vousdoukastokirezci))){

    for(qq in c("q0.05","q0.5","q0.95")){


            temp<-all2f.freq.probproj[qq,,ii]

            whenone[ii,qq]<-seq(length(temp)-1)[temp[-1]==1][1]

    }
    }
whenone[is.na(whenone)]<-7


all2.label.100to1.q0.05<-label.100to1.q0.05<-whenone[,"q0.05"]

all2.label.100to1.q0.5<-label.100to1.q0.5<-whenone[,"q0.5"]

all2.label.100to1.q0.95<-label.100to1.q0.95<-whenone[,"q0.95"]

return(list(
  all3 = list(
    q0.05 = all3.label.100to1.q0.05,
    q0.5  = all3.label.100to1.q0.5,
    q0.95 = all3.label.100to1.q0.95
  ),
  all2 = list(
    q0.05 = all2.label.100to1.q0.05,
    q0.5  = all2.label.100to1.q0.5,
    q0.95 = all2.label.100to1.q0.95
  )
))

}

