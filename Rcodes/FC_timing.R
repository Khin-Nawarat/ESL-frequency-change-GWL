# Author: Khin Nawarat (2026, @IHE Delft)
# Compiling the results for later saving as netcdf files

remove(list=ls())
filedir<-"D:/Khin Nawarat/ESLs at different GWLs/0_SD_runs/"
source(paste0(filedir,"Rcode/allfunctions.r"))


alldecades<-seq(2020,2100,by=10)
RLs <- c(10, 20, 50, 100, 200)
qs  <- c("q0.05", "q0.5", "q0.95")
gwls <- c("1p5degree","2p0degree","2p5degree","3p0degree","4p0degree","5p0degree")
alldecades <- seq(2020, 2100, by = 10)
nlocs <- c(179,7283)

for (k in 1:2) {
nloc<-nlocs[k]
allprobproj <- if (nloc > 179) {
  all3f.freq.probproj
} else {
  all2f.freq.probproj
}

first_decade <- array(NA_integer_,dim = c(nloc, length(qs), length(gwls), length(RLs)), 
                      dimnames = list(location = NULL,quantile = qs,GWL = gwls, RL = as.character(RLs)))

for (r in seq_along(RLs)) {
  
  RL <- RLs[r]
  message("Processing RL = ", RL)
  
  for (whichdecade in alldecades) {
    
    message("  decade = ", whichdecade)
    load(paste0(filedir,"Rdatasets/RData_fullconvolutions_fisherinfo_",RL, "yr_", whichdecade))
    
    for (q in qs) {
      q_idx <- match(q, dimnames(allprobproj)[[1]]) 
      
      for (g in gwls) {
        g_idx <- match(g, dimnames(allprobproj)[[2]])
        
        vals <- allprobproj[q_idx, g_idx, ]
        
        hit <- vals == 1 & is.na(first_decade[, q, g, r])
        
        first_decade[hit, q, g, r] <- whichdecade
      }
    }
  }
}

save(
  first_decade,
  file = paste0(
    filedir,
    if (nloc > 179) {
      "Data/RData_first_decade_FC_allRLs_q_wls3"
    } else {
      "Data/RData_first_decade_FC_allRLs_q_wls2"
    }
  )
)
}
                   