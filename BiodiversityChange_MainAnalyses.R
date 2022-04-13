
### CODE INFORMATION ###########################################################
# This R Code runs the main analyses for the paper:
# Crockett ETH, M Vellend, and EM Bennett. Tree biodiversity in northern forests
# shows temporal stability over 35 years at different scales, levels, and 
# dimensions of diversity.

# Since some aspects of the analyses are time or memory consuming, the code is designed
# to run separately in different chunks. There five main stages include:
#1 calculate diversity indices
#2 calculate alpha diversity change
#3 explain alpha diversity change
#4 assess changes in temporal turnover and spatial beta diversity.
#5 generate summary analyses and figures

#At the top of each stage there may be options for the user to specify which specific
#part of the analysis (and sensitivity analyses) they would like to run.
#User must ensure that the appropriate directories have been specified for file
#input and outputs

#Written by Erin Crockett
#Updated 2022-Mar-01



################################################################################
#################### -- STAGE 1 -- #############################################
################################################################################

#This stage uses the community data matrices, phylogeny, traits, and spatial scale
#information to generate values of alpha diversity and generate distance matrices
#(used later for calculations of FD and PD) at all spatial scales

rm( list=ls() )


### (1) READ IN DATA AND SET UP PARAMETERS  ####################################

#!! Set your own path here  (and again throughout the code in each section/subsection)
BaseDirectory <- "..."

library(vegan)     #alpha diversity calculations
library(entropart) #PD alpha div Hill numbers
library(picante)   #read.tree function for phylogeny
library(GUniFrac)  #for GUniFrac - phylogenetic beta diversity function
library(beepr)


##[1.1] Set Parameters for the Replicates/Loops --------------------------------

#Indicate Which Parts of the Code to Run (some sections take a long time to run)
run.aDiv.Regular    <- "Yes"   #"Yes" or "No"
run.aDiv.ForestLoss <- "Yes"
run.distMatrices    <- "Yes"

#Number of replicates
nrep <- 100

#Spatial Scales (different names for different folders and columns)
resLevels <-  c("gCell200","gCell100","gCell50","QC")  
resLevel.foldr <- c("r200", "r100", "r50","rQC")
resLevel.names <- c("200km","100km","50km","Plot")

#For Sensitivity Analysis - #Low, Med, High levels of resampling within each grid cell
nplot <- c("pMid", "pLo", "pHi")  

#Aspects (or Dimensions) of Diversity
BDaspect <- c("TD", "PD", "FD")


#For Forest Loss Analyses
forest.reps <- 100                     #how many replicates to run
t.multiplier <- 35/26                  #time multiplier- since dont have disturbance data for ALL years
percentRem <- c(25,50,75,90,95,100)    #what percent of trees are removed from disturbed plots
            
#Create directories to save files, aDivVals and DistMats
dir.create( paste( BaseDirectory, "aDivVals", sep="/") )
dir.create( paste( BaseDirectory, "DistMats", sep="/") )


##Load the functions
setwd( BaseDirectory )
source("BiodiversityChange_Functions.R")


for(rr in 1:length(resLevels)){ #all levels of spatial grid cells

   #Note the only one replicate for QC (all plots used), but multiple reps for the other scaps
   #And only pMid for QC (rather than pLo, pHi as well)
   if( resLevels[rr] == "QC" ){
      loop.reps <- 1
      nplot.rr <- "pMid"
   }else{
      loop.reps <- nrep
      nplot.rr <- nplot
   }
  
   dir.create( paste( BaseDirectory, "aDivVals", resLevel.foldr[rr], sep="/") )
   dir.create( paste( BaseDirectory, "DistMats", resLevel.foldr[rr], sep="/") )

   
   for(pp in 1:length(nplot.rr)){    # different amounts of plots within each grid cell - sensitivity analysis

      dir.create( paste( BaseDirectory, "DistMats", resLevel.foldr[rr], nplot.rr[pp], sep="/") )
   
      ### (2) DATA MATCHING WITH COMM MATRICES  ################################

      #Create a name variable - useful for reading in data

      nameString.XX <- paste(BaseDirectory, "DataFiles", "CommMatrixFiles", resLevel.foldr[rr], sep="/") 
      setwd(nameString.XX)
      #Read in comm matrix: within each folder, each RData file stores an
      #array, where dim3 is the matrix for that particular replication
      if( resLevels[rr] == "QC" ){
         load( "commList1.RData" )
         load( "commList2.RData" )
      }else{
         load( paste("commList1_", nplot.rr[pp], ".RData", sep="") )
         load( paste("commList2_", nplot.rr[pp], ".RData", sep="") )
      }

      #Create Lists to Save the Mean Outputs  (list since dont know dimensions a priori)
      aDiv.Reg.TD <- list()
      aDiv.Reg.PD <- list()
      aDiv.Reg.FD <- list()

      ###[2.1] Run nn replicates loop ------------------------------------------------
      for(nn in 1:loop.reps){    
         
         #Save out the object for easy reference
         comm1 <- commList1[[nn]]  
         comm2 <- commList2[[nn]]
                 
         ###[2.2] Match Data - CommMatrix, Phylogeny, Traits----------------------------------
         #Only keep columns which have a species present at either time interval
         commBoth <- rbind(comm1, comm2)                 #Bind them together
         commBoth <- commBoth[  ,colSums(commBoth) > 0]  #Extract species with abundance >0
         spComm <- colnames(commBoth)                    #Save a list of the species found
         comm1 <- comm1[ ,spComm]                        #Trim comm1 and comm2 matrices to these species
         comm2 <- comm2[ ,spComm]

         ##(2.2.1) PHYLOGENY - Trim & Match --------------------------------------------
         #note: in overall Tiges dataset - 3 species are NOT found in the phylo  (AUC, ERN, ORT)
         #(a) Read in the data
         #(b) Adapt phylogeny to deal with species NOT found in the phylo
         #(c) Match species across phylogeny and comm matrix
    
         ###(a) Read in Phylogeny & Traits Data
         setwd( paste( BaseDirectory, "DataFiles", sep="/") )
         QCphylo <- read.tree("QCphyloCodes.tre")	#Phylogeny, with tip labels as species codes rather than full names

         comm1.p <- comm1  #create new objects so that dont write over the existing ones (used again later)
         comm2.p <- comm2
         
         ###(b) Not ALL species are contained in the phylogeny. This section of code (a) 
         #makes accommodations when a sp. is present in comm matrix but is NOT found in the 
         #phylogeny. In the CommMatrix, it either adds the name of that species directly 
         #(when only one of two closely related sp are found) OR merges the two together. 
         #For the merge, its better to have some phylo information than none at all. 
         #The species NOT found in phylogeny are species that are quite rare across QC (minor impacts on results)

         ##Write Function
         matchSpeciesInPhylo <- function(spA, spB, comm1.p, comm2.p, my.phylo){
            #spA - species that may be found in the sites, BUT is NOT contained in the phylogeny
            #spB - the closest relative of spA that IS found in the phylogeny 
            #comm1.p, comm2.p - the community data matrices
            #my.phylo - the phylogeny to use
     
            spComm.p <- colnames(comm1.p)
            #Create object with lists if one or both of the species are present
            sp.phylo.X <- spComm.p[ spComm.p %in% c(spA,spB) ]
            ##If both are present - merge the species comm1 comm2 matrices
            #ie collapse the sp info in two columns into one column
            if( length(sp.phylo.X) == 2 ){
               comm1[ ,spB] <- comm1[ ,spA] + comm1[ ,spB]
               comm1[ ,spA] <- NULL
               comm2[ ,spB] <- comm2[ ,spA] + comm2[ ,spB]
               comm2[ ,spA] <- NULL 
            }
            ##If only one is present - change the phylogeny to the correct name to match the species that is found   
            if( length(sp.phylo.X) == 1 && sp.phylo.X == spA ){	#If comm1 only contains the species not in the phylo...
               xx <- which( my.phylo$tip.label == spB )	#Find the label of the sp that is there, and change it to the species that is found in the comm matrix
               my.phylo$tip.label[xx] <- spA
            }
            return(my.phylo)
         } 
   
         #Run Function for three pairs of species (overwriting each time)
         QCphylo <- matchSpeciesInPhylo("ORT", "ORR", comm1.p,comm2.p, QCphylo)    #ulmus rubra ORR = ORR + ORT
         QCphylo <- matchSpeciesInPhylo("AUC", "AUR", comm1.p,comm2.p, QCphylo)    #incana  AUR = AUR + AUC 
         QCphylo <- matchSpeciesInPhylo("ERN", "ERS", comm1.p,comm2.p, QCphylo) 

         ###(c) Match the Phylogeny and Comm matrix
         spPhylo <- QCphylo$tip.label  #Create a list of species within the phylogeny
         #Find which species are different between comm and phylo
         SpNotInPhylo <- setdiff(spComm, spPhylo)  #Sp found in the Comm matrices, but NOT in phylogeny
         SpToDrop.phylo <- setdiff(spPhylo, spComm)   #Sp found in the Phylogeny, but NOT in Comm matrices

         #Drop columns of species that are not in the phylogeny
         comm1.p[ ,SpNotInPhylo] <- NULL   #drop the columns
         comm2.p[ ,SpNotInPhylo] <- NULL     

         #Drop Tips not Present in the Community Matrix
         phylo <- drop.tip( phy=QCphylo, tip=SpToDrop.phylo)   

         #Arrange colnames (sp) of comm1 and comm2 to have same order as the phylo tip labels
         #because this arrangement is needed for some of the phylo diversity calculations
         spOrderPhylo <- phylo$tip.label      #Create an order object to use for traits matrix (below) as well
         comm1.p <- comm1.p[ ,spOrderPhylo]   #Arrange the columns in order
         comm2.p <- comm2.p[ ,spOrderPhylo]


         ##(2.2.2) TRAITS - Trim & Match --------------------------------------------
         setwd( paste( BaseDirectory, "DataFiles", sep="/") )
         
         ###(a) Write function to clean up the traits data - that I can use both both 5 and 17 traits
         trimTraits <- function( nTraits, spInCommMatrix ){
            #nTraits - number of traits (5 or 17) used in analysis
            #spInCommMatrix - the list of species names contained within this particular comm matrix
     
            load( paste0("traits", nTraits, ".RData") )  #Load in Traits data
            traitsPhylo <- traits17$tPhylo               #phylo object based on traits dendrogram
            
            ###(a) Traits phylogeny object (see Code section 3_Traits for how this was created)
            #Looking at subset of resulting distance matrix suggests that this doenst change distance scores obtained
            SpToDrop.traits <- setdiff(traitsPhylo$tip.label, spInCommMatrix) #Note: Using spOrderPhylo here, to make sure have matching species for PD and FD analyses
            traitsPhylo <- drop.tip( phy=traitsPhylo, tip=SpToDrop.traits)

            #Arrange colnames (sp) of comm1 and comm2 to have same order as the traitsPhylo tip labels
            #because this arrangement is needed for some of the phylo diversity calculations
            spOrderTraits <- traitsPhylo$tip.label  #Create an order object to use for traits matrix (below) as well
            comm1.t <- comm1[ ,spOrderTraits]   #Arrange the columns in order
            comm2.t <- comm2[ ,spOrderTraits]

            ###(c) Return the object
            traitsClean <- list( comm1.t, comm2.t, traitsPhylo)
            names(traitsClean) <- c("comm1","comm2","traitsPhylo")
            return(traitsClean)
         }#end function trimTraits
    
         ###(b) Use the function to clean the data
         t17 <- trimTraits(nTraits=17, spComm)    
  

         ##(2.2.3) Test that all matching (&listed same order)------------------------------------
         #Make sure all files list species in the same order - need for diversity analyses
         if( ! all.equal( colnames(comm1.p), colnames(comm2.p), phylo$tip.label ) ){
            paste0("!! Sp.Order of Comm1, Comm2, Phylo do NOT match!!:", " rr",rr," pp",pp," nn",nn )
            break 
         }
         if( ! all.equal( colnames(t17$comm1), colnames(t17$comm2), t17$traitsPhylo$tip.label ) ){
            paste0("!! Sp.Order of Comm1, Comm2, TraitsPhylo do NOT match !!: ", " rr",rr," pp",pp," nn",nn )
            break 
         }
         
         
         
         ### (3) CALCULATE ALPHA DIVERSITY #####################################
         #for each location at time1 and time 2
         #Note: about 2hrs to run 100 reps for r50,100,200, pLo,pMid,pHi
         
         ###[3.1A] ALPHA DIVERSITY - REGULAR -----------------------------------
         if(run.aDiv.Regular == "Yes"){ 
            
            for(aa in 1:length(BDaspect)){
               #Calculate BD, call different function depending on the aspect (TD,PD,FD, 5v17traits)
               if(BDaspect[aa]=="TD"){ aDiv.Reg.TD[[nn]] <- prepTDalpha.Hill(comm1.p, comm2.p) }         #Taxonomic Diversity
               if(BDaspect[aa]=="PD"){ aDiv.Reg.PD[[nn]] <- prepPDalpha.Hill(comm1.p, comm2.p, phylo) }  #Phylogenetic Diversity
               if(BDaspect[aa]=="FD"){ aDiv.Reg.FD[[nn]] <- prepPDalpha.Hill(t17$comm1, t17$comm2, t17$traitsPhylo) }  #Functional Diversity: using dendrogram/phylo (PD) metrics
            }#close for aa loop
         }#close if run.aDiv.Regular

     
  	      ###[3.1B] ALPHA DIVERSITY - ADJUST FOR LANDUSE CHANGE -----------------
         if( run.aDiv.ForestLoss == "Yes" ){
 	         
            dir.create( paste( BaseDirectory, "aDivForestLoss", sep="/") )
            dir.create( paste( BaseDirectory, "aDivForestLoss/Replicates", sep="/") )
         
            #Note: Only run for the Plot Level, not all spatial scales.
            if( resLevels[rr] == "QC" ){
          	              
               #Read Disturbance Files
 	            setwd( paste( BaseDirectory, "DataFiles", sep="/") )
 	            load("infoFirst_tiges.RData") 
               info.first$gCellBioDomain <- as.character(info.first$gCellBioDomain)
               
               load("ForestChange_byBioDomain.RData") #object: Chg.bioDomain.new
               percent.Disturb <- Chg.bioDomain.new[ ,"TotalDisturbPrc"]

               #!Note code is quite slow for PD and FD... 
               for(aa in 1:length( BDaspect) ){     
            
                  for(rrr in 1:length( percentRem ) ){

                     if( BDaspect[aa]=="TD" || BDaspect[aa]=="PD" ){ 
                        aDiv <- runForestL.reps.IndividualTrees( comm1.p, comm2.p, phylo, 
                                        BDaspect[aa], percentRem[rrr], percent.Disturb, 
                                        info.first, t.multiplier, forest.reps )
                     }
                     if( BDaspect[aa]=="FD" ){ 
                        aDiv <- runForestL.reps.IndividualTrees( t17$comm1, t17$comm2, t17$traitsPhylo, 
                                        BDaspect[aa], percentRem[rrr], percent.Disturb, 
                                        info.first, t.multiplier, forest.reps )
                     }
                        
                     setwd( paste( BaseDirectory, "aDivForestLoss", "Replicates", sep="/") )
                     save(aDiv, file= paste0("aDiv_", BDaspect[aa], "_Rem", percentRem[rrr],"_AllReps", ".RData") )  

                     #Calculate the mean across all reps
                     aDiv.XX <- apply( aDiv, MARGIN=c(1,2), mean )
                     aDiv.XX <- as.data.frame(aDiv.XX)
                     colnames(aDiv.XX) <- c("q0","q1","q2")
                     aDiv.XX$timePeriod <- rep( c("t1","t2"), each = nrow(aDiv.XX) / 2 )
                     nsites <- nrow(aDiv.XX)/2
                     aDiv.XX$SiteID <- rep( rownames(aDiv.XX)[1:nsites], 2 )
   
                     setwd( paste( BaseDirectory, "aDivForestLoss", sep="/") )
                     save(aDiv.XX, file= paste0("aDiv_", BDaspect[aa], "_Rem", percentRem[rrr], ".RData") )  
 
                  }#close rrr PercentRemoved categories
 
               }#close aa
            }#close if Plot scale
         }#close if Run Forest loss replicates
               

         ### (4) CALCUALTE BETA DIVERSITY ######################################
         #Create and SAVE Distance MAtrix to use for all later sections
         #Note: creating distance matrices can be slow and memory intensive...
         if( run.distMatrices == "Yes" ){
                        
            for( aa in 1:length(BDaspect) ){

               ##[4.1] Specify which Distance metrics, CommMatrices, Phylogeny to use ----------------
               commB <- NULL      #clear from previous run
               phyloToUse <- NULL

               #Taxonomic Beta Diversity
               if( BDaspect[aa] == "TD" ){ 
                  distMethod <- c("Sorensen","Bray-Curtis")
               #Phylogenetic Beta Diversity (or for Functional beta diversity calculated with PD metrics)
               }else{
                  distMethod <- c("UniFrac","Dpw")
               }
             
               #Functional Traits
               if( BDaspect[aa] == "FD" ){ 
                  commB <- rbind(t17$comm1, t17$comm2)  
                  phyloToUse <- t17$traitsPhylo
               #Taxonomic or Phylogenetic    
               }else{
                  commB <- rbind(comm1.p, comm2.p)  
                  phyloToUse <- phylo
               }

               ##[4.2] Create Distance Matrices --------------------------------
               #Create a presence-absence matrix (needed for TD and PD)
               commBPA <- commB
               commBPA[commBPA > 0] <- 1      
   
               for(bb in 1:length(distMethod)){
               
                  ##[4.2.1] Dist for TD metrics --------------------------------
                  if( BDaspect[aa] == "TD" ){
                  #note: need 1- for Sorensen to turn similarity to dissimilarity index
                     if(distMethod[bb]=="Sorensen"){ dMat.XX <- 1 - betadiver(x=commBPA, method="sor") }        #Presence-Absence
                     if(distMethod[bb]=="Bray-Curtis"){ dMat.XX <- vegdist(x=commB, method="bray",  binary=FALSE) } #Abundance
                     
                  ##[4.2.2] Dist for PD metrics --------------------------------   
                  }else{
                     #Note: Dpw doesnt work when only have one species in the community
                     #(since cant calc mean pairwise distance for one of two communities that are looking at)
                     #Therefore need to remove any communities which are monocultures from the commB matrix,
                     #and also from the timeMat matrix
         
                     comms.with.1sp <- which(rowSums(commBPA) < 2 )
                     #Note: function fails if there are NO sites with 1 species - hence include if-else statement
                     if( length(comms.with.1sp) >= 1 ){
                        #Note: need to remove the same community at time 1 and time2
                        ncomm1 <- nrow(commB)/2  #number of unique communities in analysis

                        comm1.2sp <- comms.with.1sp[comms.with.1sp <= ncomm1] #which sites in comm1 have at least 2 sp
                        comm2.2sp <- comms.with.1sp[comms.with.1sp > ncomm1]  #which sites in comm2 have at least 2 sp
   
                        comm2.match <- comm1.2sp + ncomm1  #the matching community at time 2
                        comm1.match <- comm2.2sp - ncomm1  #the matching community at time 1
                        #All the communities (and their matching buddies) that only have 1 species
                        rowsToRemove <- sort( unique(c(comm1.2sp,comm2.2sp,comm2.match,comm1.match)) )
                        ## The updated matrices
                        commB.trimmed <- commB[-rowsToRemove, ]   #keep all the rows, except these ones that I want to remove
                     }else{
                        commB.trimmed <- commB
                        rowsToRemove <- ""
                     }

                     #Create Distance Matrices
                     if( distMethod[bb] =="UniFrac" ){
                        #Run UniFrac function, alpha=0 just to include a value (will use the unweighted function anyways, so doesnt matter the value)
                        dMat.XX.int <- GUniFrac(otu.tab=commB, tree=phyloToUse, alpha=0 ) 
                        dMat.XX <- as.dist(dMat.XX.int$unifracs[ , ,"d_UW"])  #Extract the "unweighted" value from list/array that is returned
                        rm(dMat.XX.int)
                     }
                     if( distMethod[bb] =="Dpw" ){ dMat.XX <- Dpw.Abun(comm=commB.trimmed, phylo=phyloToUse) }
                  }#close if/else TD/PD metrics

                  ##[4.2.3] Save the distance matrix for later use -------------
                  if( resLevel.foldr[rr] == "rQC" ){
                     setwd( paste( BaseDirectory, "DistMats",  resLevel.foldr[rr],  sep="/") )
                     save(dMat.XX, file = paste(BDaspect[aa],"_",distMethod[bb],".RData", sep="") )
                  }else{
                     setwd( paste( BaseDirectory, "DistMats", resLevel.foldr[rr], nplot.rr[pp], sep="/") )
                     save(dMat.XX, file = paste(BDaspect[aa],"_",distMethod[bb],"_nn", nn,".RData", sep="") )
                  }

                  #Save the Rows Removed - for Dpw (to refer to later)
                  if(distMethod[bb]=="Dpw" ){ 
                     save(rowsToRemove,  file=paste("rowsToRemove_", BDaspect[aa],"_",distMethod[bb],"_nn",nn,".RData", sep="") )
                     save(commB.trimmed, file=paste("commBtrimmed_", BDaspect[aa],"_",distMethod[bb],"_nn",nn,".RData", sep="") )
                  
                  }
                  print( paste("--- Completed", BDaspect[aa], distMethod[bb],"---", sep=" ") )
                  
               }#close for bb distMethod
            }#close for aa BDaspects

         }#close if run.distMatrices


      }#close nn loop - loop.reps      

      ##Alpha Diversity: Calculate Mean & Save to File
      if( run.aDiv.Regular == "Yes" ){ 
         setwd( paste( BaseDirectory, "aDivVals", resLevel.foldr[rr], sep="/") )
         if(resLevels[rr] == "QC"){
            pasteAddOn <- "_"
         }else{
            pasteAddOn <- paste0("_",nplot.rr[pp],"_")
         }
         aDiv.XX <- calcMean( aDiv.Reg.TD ) 
         save(aDiv.XX, file= paste0("aDiv_TD", pasteAddOn, "mean.RData") )  
         aDiv.XX <- calcMean( aDiv.Reg.PD ) 
         save(aDiv.XX, file= paste0("aDiv_PD", pasteAddOn, "mean.RData") )  
         aDiv.XX <- calcMean( aDiv.Reg.FD ) 
         save(aDiv.XX, file= paste0("aDiv_FD", pasteAddOn, "mean.RData") ) 
      }
      
      print( paste("At pp", pp ))
   } #close pp loop - nplot
         
   print( paste("At rr", rr," -----------------------------") )
}#close rr loop - resLevels
   




################################################################################
#################### -- STAGE 2 -- #############################################
################################################################################

#This R Script runs a T-Test to determine if there has been diversity change
#between the two time periods. Runs loops for different dimensions/aspects of diversity
#and for different types of change (reg difference, percent change, log change).
#Option to specify for regular repeated measurements, or corrections for forest loss/land cover change.


rm( list=ls() )


library(dplyr)


### (1) SET PARAMETERS #########################################################
   
#!! Set your path here:
BaseDirectory <- "..."

#Indicate if alpha Diversity is absolute/regular change, percentage change, or log change
Chng.Type <- "AbsChg"   #options:  "AbsChg" or "LogChg"  

#Indicate whether to run regular Repeated Measurements or with ForestLoss Corrections
   #if No, then runs regular repeated measures
   #if Yes, then runs forest loss replicates
run.aDiv.ForestLoss <- "No"  #options:  "No"or"Yes"   



## Read in data
resLevels <-  c("gCell200","gCell100","gCell50","QC")   #Resolution - the colnames in InfoSptl file
resLevel.foldr <- c("r200", "r100", "r50","rQC")
resLevel.names <- c("200km","100km","50km","Plot")

BDaspect <- c("TD","PD","FD")

nplot <- c("pMid","pLo","pHi")  #number of plots within each grid cell to use in comm matirx

percentRem <- c(25,50,75,90,95,100)    #what percent of trees are removed from disturbed plots
 

#Load file
setwd( paste( BaseDirectory, "DataFiles", sep="/") )
load("infoFirst_tiges.RData")
#Load Functions
setwd( BaseDirectory )
source("BiodiversityChange_Functions.R")

#Function to summarize the data at each level, add a column, and convert from List to datafrmae
EC_addColumns <- function(dataList, colToAdd, loopVar, colToAddName){
   dataMat <- dplyr::bind_rows(dataList)
   nc <- ncol(dataMat)
   dataMat <- data.frame(dataMat, colToAdd= rep(colToAdd[loopVar], nrow(dataMat)) )
   colnames(dataMat)[nc+1] <- colToAddName
   dataMat <- dataMat %>% mutate_if(is.factor, as.character)
   return( dataMat)  
}


### (2) LOOP THROUGH aDIV CALCULATIONS ###############################################

#Adjust for Forest Loss --- only done at QC/Plot scale
if( run.aDiv.ForestLoss == "Yes" ){ 
   resLevels <-  "QC"   
   resLevel.foldr <- "rQC"
   resLevel.names <- "Plot"
}

aDiv.resLevels.dlist <- list() 

for(rr in 1:length(resLevels)){ #all levels of spatial grid cells

   aDiv.nplot.dlist <- list()
   
   ##[2.1] Change so that dont have to loop through QC (only pMid) ---------------
   if( resLevels[rr] == "QC"){ 
      nplot = "" 
      location.df <- info.first[ ,c("ID_PE","gCellBioDomain")]
      colnames(location.df) <- c("SiteID","BioDomain")
   }else{
   
      #Find The Most Common Biodomain for each Grid Cell   
      unique.cells <- unique( info.first[ ,resLevels[rr]] )
      location.df <- data.frame( SiteID = unique.cells, BioDomain = rep(1, length(unique.cells)),
                              nBioD = rep(1, length(unique.cells) ) )
      for(cc in 1:length(unique.cells) ){
         info.cc <- info.first[ which( info.first[ ,resLevels[rr]] == unique.cells[cc] ),  ]
         cc.table <- table( info.cc$gCellBioDomain )
         location.df[cc,"BioDomain"] <- names( which.max(cc.table) )
         location.df[cc,"nBioD"] <- length(cc.table) 
      }
   }     
   location.df$SiteID <- as.character(location.df$SiteID)
   location.df$BioDomain <- as.character(location.df$BioDomain)

   #Adjust for Forest Loss --- use % Removed instead of nPlot (pMid, Lo Hi) since its all pMid
   if( run.aDiv.ForestLoss == "Yes" ){ 
      nplot <- percentRem 
   }

   for(pp in 1:length(nplot)){    # different amounts of plots within each grid cell - sensitivity analysis

      #Create empty list to work summarize at the end
      aDiv.aspect.dlist <- list()   

      ##Loop through all metrics of diversity (4)
      for(aa in 1:length(BDaspect)){
            
         ##Load the Data (mean values across all replicates)
         #Regular Values Replicates
         if( run.aDiv.ForestLoss == "No" ){ 
            setwd( paste( BaseDirectory, "aDivVals", resLevel.foldr[rr], sep="/") )
            if( resLevels[rr] == "QC" ){
               load( paste("aDiv", BDaspect[aa], "mean.RData", sep="_") ) 
            }else{
               load( paste("aDiv", BDaspect[aa], nplot[pp], "mean.RData", sep="_") )  
            }
         }
         #Forest Loss Replicates
         if( run.aDiv.ForestLoss == "Yes" ){ 
            setwd( paste( BaseDirectory, "aDivForestLoss", sep="/") )
            load( paste0("aDiv_", BDaspect[aa], "_Rem", nplot[pp], ".RData") )    #note nplot instead of PRem
         }
            

         ##[2.2] Absolute Change (x2 - x1) ---------------------------------------
         aDiv.aspect.values <- TTest.aDivChg(aDiv.XX, Chng.Type )

         aDiv.aspect.values$BDaspect <- rep( BDaspect[aa], 3)
         aDiv.aspect.dlist[[aa]] <- aDiv.aspect.values  

      } #close aa BDaspects loop

      #Add columns in indicate the BDaspect and turn into next level of the list
      aDiv.nplot.dlist[[pp]] <- EC_addColumns(aDiv.aspect.dlist, nplot, pp, "nPlot") 

   } #close pp loop - nplot
   aDiv.resLevels.dlist[[rr]] <- EC_addColumns(aDiv.nplot.dlist, resLevel.names, rr, "resLevels") 
} #close rr loop - resLevels

#Bind list results into one big df
aDiv.Results <- dplyr::bind_rows(aDiv.resLevels.dlist)   


### (3) SAVE RESULTS ###########################################################

dir.create( paste( BaseDirectory, "Results", sep="/") )
setwd( paste( BaseDirectory, "Results", sep="/" ) )

#REGULAR - ALL SPATIAL SCALES
if( run.aDiv.ForestLoss == "No" ){
   write.csv(aDiv.Results, file= paste0("aDiv_Results_", Chng.Type, ".csv") , row.names=FALSE)
}

#FOREST LOSS REPLICATIONS
if( run.aDiv.ForestLoss =="Yes" ){
   colnames(aDiv.Results)[ which( colnames(aDiv.Results)=="nPlot" )] <- "PercentRem"
   write.csv(aDiv.Results, file= paste("aDiv_ForestLoss_", Chng.Type, ".csv", sep=""),  row.names=FALSE)
}




################################################################################
####################  -- STAGE 3 --  ###########################################
################################################################################

#This R Script uses explanatory variables to explain the change in alpha diversity
#between two time periods using multiple regression, AIC model selection,
#and averaging coefficients of models within 2 AIC units.


rm(list=ls())

library(MuMIn)
library(dplyr)
library(foreign)
library(rgdal)  #readOGR
library(gstat)  #variogram
library(rgeos)  #gCentroid


### (1) SET PARAMETERS #########################################################

#!! Set your path here
BaseDirectory <- "..."

#Indicate which time lag to use for the modelling approach
time.lag <- "_0yr"   #options:  "_0yr", "_10yr" or "_20yr"

#Indicate whether to Use Full Model Coefficients or the Subset
use.full.mod.coeffs <- "Yes"  #options:  "Yes" or "No"


#
nPlot <- c("pMid")  #note: too few points (16) at pHi 200km to run stats - only run for pMid

resLevels <-  c("gCell200","gCell100","gCell50","QC")   #Resolution - the colnames in InfoSptl file
resLevel.foldr <- c("r200", "r100", "r50","rQC")
resLevel.names <- c("200km","100km","50km","Plot")
resLevels.num <- c(200,100,50,1)
 
BDaspect <- c("TD","PD","FD") 
BDmetric <- c("q0","q1","q2")

## Read in Main Variables
setwd( paste( BaseDirectory, "DataFiles", sep="/") )
load("Explantory_Variables.RData")  #object: mod.data
#Load Functions
setwd( BaseDirectory )
source("BiodiversityChange_Functions.R")

#To summarize the data at each level, add a column, and convert from List to datafrmae
EC_addColumns <- function(dataList, colToAdd, loopVar, colToAddName){
   dataMat <- dplyr::bind_rows(dataList)
   nc <- ncol(dataMat)
   dataMat <- data.frame(dataMat, colToAdd= rep(colToAdd[loopVar], nrow(dataMat)) )
   colnames(dataMat)[nc+1] <- colToAddName
   dataMat <- dataMat %>% mutate_if(is.factor, as.character)
   return( dataMat)  
}

if( use.full.mod.coeffs == "Yes") {
   add.name <- "Full"
}else{
   add.name <- "Subset"
}

dir.create( paste0( BaseDirectory, "/Results/Explain_aDiv_", add.name ) )
dir.create( paste0( BaseDirectory, "/Results/Explain_aDiv_", add.name,"/ModDiagnostics" ) )


###(2) RUN THE MODELS ##########################################################

###[2.1] For Different Spatial Scales ------------------------------------------
resLevels.coeffs.list <- list()  #Create Empty lists to fill during the loop
resLevels.modvals.list <- list()

for(rr in 1:length(resLevels)){

   ###[2.2] Extract Explanatory Variables - for that Scale ---------------------
   ## Extract Spatial Scale
   mod.data.rr <- mod.data[ which(mod.data$resLevels == resLevel.names[rr] ), ]
   ## Extract Climate Time Lag Variables
   clim.vars <- mod.data.rr[  ,grep( time.lag, colnames( mod.data.rr ) )]

   ## Climate PCAs
   #Scale the climate data - mean 0, sd 1
   clim.vars <- scale( clim.vars )
   Grow.season.days <- paste0("tc0_days", time.lag)
   
   ## Temperature
   temp.vars <- c("tc0_temp","maxt","mint")
   temp.vars <- paste0( temp.vars, time.lag )
   temp.pca <- princomp(clim.vars[ , temp.vars])
   #Extract the same number of traits
   Temperature.Ax1 <- temp.pca$scores[ ,1] #First 2 axes provide ~93% of variation (depending on scale)
   Temperature.Ax2 <- temp.pca$scores[ ,2]
   
   ## Precipitation
   precip.vars <- c("tc0_pcp","pcp")
   precip.vars <- paste0( precip.vars, time.lag )
   precip.pca <- princomp(clim.vars[ , precip.vars])
   #Extract the first precipitation axis
   Precipitation <- precip.pca$scores[ ,1] 
  
   clim.df <- data.frame( Temperature.Ax1, Temperature.Ax2, 
                          Grow.Season.Days=clim.vars[ ,Grow.season.days], Precipitation )

   ## Transform the Proportion Variables
   proportion.vars <- mod.data.rr[ ,c("FirePrc","HarvestPrc","PA","Hunt-Fish","Private")]     #as percentages
  
   if( resLevels[rr] == "QC" ){
      proportion.df <- proportion.vars
   }else{
      proportion.vars <- proportion.vars / 100  #convert percetange to proportion data
      #Logit Transformation
      proportion.df <- apply( proportion.vars, MARGIN=2, function(x) log( (x+0.0001)/(1-(x+0.0001)) )   ) 
   }
   
   ## Recombine into one df
   mod.df <- data.frame( proportion.df, clim.df)

   ## Scale all the variables - to be on same relative scale
   mod.df <- scale(mod.df)
   mod.df <- data.frame(  SiteID=mod.data.rr$SiteID, mod.df)


   ###[2.3] XY Coordinates - For Sptl correlation -------------------------------------------------------------------
   #Read in the Location Files    
   setwd( paste( BaseDirectory, "DataFiles", sep="/") )
   if( resLevels[rr] == "QC" ){
      XX.coords <- read.dbf("Placette_UTM19_Included.dbf")
      colnames(XX.coords)[2:3] <- c("y","x")
   }else{
      XX.km.int <- readOGR( dsn=".", paste("QC_", resLevels.num[rr],"km", sep="") )
      
      #Calculate Center of Each Polygon and Bind to PolygonLayer
      XX.centroids <- gCentroid( XX.km.int, byid=TRUE )
      if( resLevels.num[rr] == 50 ){
         col.xx <- paste0( "QC", resLevels.num[rr],"_OBJID")
      }else{
         col.xx <- paste0( "QC", resLevels.num[rr],"_OBID")
      }
      XX.coords <- cbind( XX.km.int@data[ ,col.xx], as.data.frame(XX.centroids) )
      colnames(XX.coords)[1] <- "SiteID" #Rename for the merge below
   }
   #Centre the coordinate data (but dont scale it)
   XX.coords[ ,c("y","x")] <- scale( XX.coords[ ,c("y","x")] , center=TRUE, scale=FALSE)

   ###[2.4] For Different Aspects (TD,PD,FD) -------------------------------------------------
   BDaspect.coeffs.list <- list()  #Create Empty lists to fill during the loop
   BDaspect.modvals.list <- list()

   dir.create( paste0( BaseDirectory, "/Results/Explain_aDiv_", add.name,"/ModDiagnostics/", resLevel.foldr[rr] ) )
            
   
   for( aa in 1:length(BDaspect) ){

      ## Create empty matrix - to fill with summary info during the loops
      m.colnames <- c("R2_best","R2adj_best","R2_mean","R2adj_mean", "n", "BDaspect","BDmetric" )
      multireg.df <- matrix(NA, nrow=length(BDmetric), ncol=length(m.colnames) )
      colnames(multireg.df) <- m.colnames
      multireg.df <- as.data.frame(multireg.df)
      multireg.df$BDaspect <- rep( BDaspect[aa], length(BDmetric)) #Add values 
      BDmetric.ceoffs.list <- list()
   
      #Create Directory to Save Files
      dir.create( paste0( BaseDirectory, "/Results/Explain_aDiv_", add.name,"/ModDiagnostics/",
                           resLevel.foldr[rr], "/", BDaspect[aa] ) )
               
      
      #Repeat for Multiple bb values
      for(bb in 1:length(BDmetric)){
            
         #Load the Data (mean values across all replicates)
         setwd( paste( BaseDirectory, "aDivVals", resLevel.foldr[rr], sep="/") )
         if( resLevels[rr] == "QC" ){
            load( paste("aDiv", BDaspect[aa], "mean.RData", sep="_") ) 
         }else{
            load( paste("aDiv", BDaspect[aa], nPlot, "mean.RData", sep="_") )  
         }
            
         #Calculate change in diversity between the two time intervals
         aDiv.chg <- calc.aDivChg( aDiv.XX, "AbsChg" )
         aDiv.chg <- data.frame( div.val= aDiv.chg[ ,bb], SiteID=rownames(aDiv.chg))
         ##Merge Data - so that both in the same order
         stats.model.data <- merge( mod.df, aDiv.chg, by="SiteID", all.x=TRUE )
         rownames(stats.model.data) <- stats.model.data$SiteID
         stats.model.data$SiteID <- NULL  #drop for stats model

         #Extract the Sites in XY Coordinates file actually used in the analysis (not all grid cells)
         id.df <- data.frame( SiteID = mod.df[ ,"SiteID"] )
         sptl.locs.touse <- merge( id.df, XX.coords, by="SiteID", all.x=TRUE, all.y=FALSE)      
         #Check Spatial Locations File
         if( ! all.equal( rownames(aDiv.chg) , rownames(sptl.locs.touse) ) ){
            print("aDiv Error - Spatial Locations - Check!")
            break
         }
        
         ##[2.5] Run Regression with AIC -------------------------------------------

         ##(2.5.1) Model with ALL Combos of Variables
         #Run the Full Model (with all variables)
         full.model.lm <- lm( div.val ~ ., data=stats.model.data )  
         #Let the model fail if any NA values
         options(na.action = "na.fail")  
         #Run all combos of variables and rank by AIC
         all.mod.results <- dredge(full.model.lm, beta = "none", rank = "AIC")
       
         #Extract Model Averaged Coefficients
         good.models <- subset( all.mod.results, delta <= 2 )
         m.avg <- model.avg( good.models )      #gives both full and subset
         m.weights <- m.avg$msTable$weight      #weights for the different models
         my.coeffs <- m.avg$coefficients
         if( use.full.mod.coeffs == "Yes" ){
            conf.ints <- confint(m.avg, full=TRUE)
            coeffs.df <- data.frame( Variable=NA, Slope=my.coeffs[1, ], conf.ints )  #"full" mod with 0s if not selected
         }else{
            conf.ints <- confint(m.avg, full=FALSE)
            coeffs.df <- data.frame( Variable=NA, Slope=my.coeffs[2, ], conf.ints )  #"full" mod with 0s if not selected
         }
         colnames(coeffs.df)[3:4] <- c("CIlow","CIupp")
         coeffs.df$Variable <- rownames(conf.ints)
         coeffs.df$OverlapZero <- "z"
         coeffs.df$BDmetric <- BDmetric[bb]
   
         #Determine if model average coefficients overlap with Zero
         for(ii in 1:nrow(coeffs.df) ){
            if( coeffs.df[ii,"CIlow"] < 0 && coeffs.df[ii,"CIupp"] > 0 ){
               coeffs.df[ii,"OverlapZero"] <- "Yes"
            }else{
               coeffs.df[ii,"OverlapZero"] <- "No"
            }
         }  

         #(2.5.3) Save R2 information
         #R2 values - from the best model
         best.model <- get.models(all.mod.results, subset = 1)[[1]]
         multireg.df[bb,"R2_best" ] <- summary(best.model)$r.squared
         multireg.df[bb,"R2adj_best" ] <- summary(best.model)$adj.r.squared

         #R2 values - weighted among the top models
         equivalent.models <- get.models(all.mod.results, subset = delta <= 2)
         r2.vec <- r2.adj.vec <- rep(999, length(equivalent.models) ) #Empty vec to fill in loop
         for(gg in 1:length(equivalent.models) ){
            r2.vec[gg] <- summary( equivalent.models[[gg]] )$r.squared
            r2.adj.vec[gg] <- summary(equivalent.models[[gg]])$adj.r.squared
         }
         multireg.df[bb,"R2_mean" ] <- weighted.mean(r2.vec, m.weights)  
         multireg.df[bb,"R2adj_mean" ] <- weighted.mean(r2.adj.vec, m.weights)  
         multireg.df[bb,"n" ] <- nrow( stats.model.data )
         multireg.df[bb, "BDmetric"] <- BDmetric[bb]

         
         ## Check Diagnostics of Good Models
         #! Note: calc.sptlCorrelation function below will fail since the MoransI function is not included, hence blocked with a #
         setwd( paste0( BaseDirectory, "Results/Explain_aDiv_", add.name,"/ModDiagnostics/",
                           resLevel.foldr[rr], "/", BDaspect[aa] ) )
         for(gg in 1:length(equivalent.models) ){
            figure.name <- paste( BDaspect[aa], BDmetric[bb], gg, sep="_" )
            ## Spatial Correlation of Residuals
            residuals.gg <- residuals( equivalent.models[[gg]] )
            residuals.gg <- data.frame( SiteID = names(residuals.gg), Resid = residuals.gg)
            sptl.locs.touse.gg <- merge( sptl.locs.touse, residuals.gg, by="SiteID")
            if( resLevel.foldr[rr] == "r200" ){ 
                #my.sp <-  calc.sptlCorrelation( sptl.locs.touse.gg, fnameF2=figure.name, dClassesF=40 )
            }else{
                #my.sp <-  calc.sptlCorrelation( sptl.locs.touse.gg, fnameF2=figure.name, dClassesF=60 )
            }
            ## General Model Diagnostics
            plotDiagnostics( equivalent.models[[gg]], "mlm", fnameF2=figure.name )
         }#close gg model loop
            
         #Add the row to the list
         BDmetric.ceoffs.list[[bb]] <- coeffs.df
         
      }#close bb loop
      BDaspect.coeffs.list[[aa]]  <- EC_addColumns(BDmetric.ceoffs.list, BDaspect, aa, "BDaspect")
      BDaspect.modvals.list[[aa]] <- multireg.df
      
   }#close aa loop
   resLevels.coeffs.list[[rr]] <- EC_addColumns(BDaspect.coeffs.list, resLevel.names, rr, "resLevels")
   resLevels.modvals.list[[rr]] <- EC_addColumns(BDaspect.modvals.list, resLevel.names, rr, "resLevels")

}#close rr loop
            
## For loop above - Bind these list results into one big df
coeff.results <- dplyr::bind_rows(resLevels.coeffs.list)   
modvals.results <- dplyr::bind_rows(resLevels.modvals.list) 

##Save Files
setwd( paste0( BaseDirectory, "/Results/", "Explain_aDiv_", add.name ) )
write.csv(coeff.results,   file= paste("ModCoeffs_", time.lag, "_",nPlot, ".csv", sep="") , row.names=FALSE)
write.csv(modvals.results, file= paste("ModR2vals_", time.lag, "_",nPlot, ".csv", sep="") , row.names=FALSE)





################################################################################
#################### -- STAGE 4 -- #############################################
################################################################################

# This R Script calculates turnover and changes in beta diversity using 
# permanova and permdisp analyses.


rm(list=ls())

library(vegan)
library(dplyr)  #%>%

### (1) SET PARAMETERS ###########################################################

#!! Set your path here
BaseDirectory <- "..."

#Indicate the number of permutations for the statistical test (changes in beta diversity) 
nperm = 999  



#
resLevels <-  c("gCell200","gCell100","gCell50","QC")   #Resolution - the colnames in InfoSptl file
resLevel.foldr <- c("r200", "r100", "r50","rQC")
resLevel.names <- c("200km","100km","50km","Plot")

nPlot <- c("pMid","pLo","pHi")  #number of plots within each grid cell to use in comm matirx

BDaspect <- c("TD","PD","FD")

#Read in the functions
setwd( BaseDirectory )
source("BiodiversityChange_Functions.R") 

#Create folders within "Results" called "PERMANOVADISP" with subfolders for each resLevel.foldr
dir.create( paste( BaseDirectory, "Results", "PERMANOVADISP", sep="/") )
for(rr in 1:length(resLevel.foldr) ){
   dir.create( paste( BaseDirectory, "Results", "PERMANOVADISP", resLevel.foldr[rr], sep="/") )
}

#Function To summarize the data at each level, add a column, and convert from List to datafrmae
EC_addColumns <- function(dataList, colToAdd, loopVar, colToAddName){
   dataMat <- dplyr::bind_rows(dataList)
   nc <- ncol(dataMat)
   dataMat <- data.frame(dataMat, colToAdd= rep(colToAdd[loopVar], nrow(dataMat)) )
   colnames(dataMat)[nc+1] <- colToAddName
   dataMat <- dataMat %>% mutate_if(is.factor, as.character)
   return( dataMat)  
}


### (2) RUN ANALYSIS ###########################################################

bDiv.resLevels.dlist <- list() 
for(rr in 1:length(resLevels)){ #all levels of spatial grid cells

   if( resLevels[rr] == "QC" ){ nPlot <- "pMid" }
      
   bDiv.nPlot.dlist <- list() 
   for(pp in 1:length(nPlot)){    # different amounts of plots within each grid cell - sensitivity analysis
      
      bDiv.aspect.dlist <- list()
      for( aa in 1:length(BDaspect) ){        
               
         bDiv.distMethod.dlist <- list()

         ## Indicate Distance Type
         #Taxonomic Beta Diversity
         if( BDaspect[aa] == "TD" ){ 
            distMethod <- c("Sorensen","Bray-Curtis")
         #Phylogenetic Beta Diversity (or for Functional beta diversity calculated with PD metrics)
         }else{
            distMethod <- c("UniFrac","Dpw")
         }
        
         ### PERMANOVA /  PERMDSIP
         #Loop through calculations for all BD metrics
         for(bb in 1:length(distMethod)){

            FileDir = paste( BaseDirectory,"DistMats",sep="/")
            
            ## Calculate Mean of distance matrices across all replicates
            #The function will read in each data file
            if( resLevel.foldr[rr] != "rQC" ){
               dMat.mean <- create.MeanDist( BDaspect[aa], distMethod[bb], nPlot[pp], resLevel.foldr[rr], FileDir)
            
            }else{
               setwd( paste( FileDir, resLevel.foldr[rr], sep="/" ) )
               load( paste0( BDaspect[aa] , "_", distMethod[bb], ".RData") )  
               dMat.mean <- bDiv.vals  #no replicates for plot scale - used all the data
               rm(dMat.XX)
            }
            #Create a dataframe for the factors
            nsites <- length( labels(dMat.mean) ) / 2
            timePeriod <-  c( rep("t1",nsites), rep("t2", nsites )  ) 
            site.names.t1 <- head(labels(dMat.mean), nsites)
            site <-  c(site.names.t1, site.names.t1)
            timeMat <- data.frame(site, timePeriod)

            ## Calculate permanova and permdisp
            bDiv.results.int <- calcPermNovaDisp( dMat.mean, timeMat, npermF=nperm,
                           BDaspect[aa], distMethod[bb], nPlot[pp],
                           saveDirF = paste( BaseDirectory, "Results", "PERMANOVADISP", resLevel.foldr[rr] ,  sep="/") )

            #Add a column to note which metric it was  (eg sorensen, bray-curtis, etc)
            bDiv.results.int <- data.frame(bDiv.results.int, BDMetric=distMethod[bb])
            #Convert columns from factor to character
            intrmdte <- sapply(bDiv.results.int, is.factor)
            bDiv.results.int[intrmdte] <- lapply(bDiv.results.int[intrmdte], as.character)
            #Add the row to the list
            bDiv.distMethod.dlist[[bb]] <- bDiv.results.int  

         }#close for bb loop
         
         #Bind these list results into one big dataframe
         bDiv.metric.mat <- dplyr::bind_rows(bDiv.distMethod.dlist)

         #Add columns for "Beta Diversity" and for which BD metric it will be
         bDiv.metric.mat <- bDiv.metric.mat %>%
                  mutate( BDaspect = rep(BDaspect[aa], nrow(bDiv.metric.mat))) %>%
                  mutate( BDdimension = rep("Beta Diversity", nrow(bDiv.metric.mat))) %>%     
                  mutate_if(is.factor, as.character)   #Change factor columns to character 
         
         #Add the row to the list
         bDiv.aspect.dlist[[aa]] <- bDiv.metric.mat
         print( paste( "completed ", BDaspect[aa], distMethod[bb], sep=" ") )
      
      }#close aa loop   

      bDiv.nPlot.dlist[[pp]] <- EC_addColumns(bDiv.aspect.dlist, nPlot, pp, "nPlot")  
   
   } #close pp loop - nPlot
   bDiv.resLevels.dlist[[rr]] <- EC_addColumns(bDiv.nPlot.dlist, resLevel.names , rr, "resLevels")  
   #Print a check up though the loops
   print( paste("At rr", rr))

}#close rr loop - resLevels


#Bind these list results into one big dataframe
bDiv.Results <- dplyr::bind_rows(bDiv.resLevels.dlist)   

## Save Results
setwd( paste( BaseDirectory, "Results", sep="/") )
write.csv( bDiv.Results, "PERMANOVADISP.csv", row.names=FALSE)



################################################################################
#################### -- STAGE 5 -- #############################################
################################################################################

#The code in this section generates summary files and figures presented in the paper.
#Each section within here is designed to be a standalone code snippet based
#on the data generated in previous steps.

#5.1:  Fig3 
#5.2:  Fig4
#5.3:  Fig5, S3, S4
#5.4:  Fig6, S7
#5.5   Table 1
#5.6   Fig S5



###(5.1) FIG 3 - Create 3 Histograms as a panel for TD,FD,PD ###################
#Run as an example for Plot scale, with q=1 (shown in paper)

rm(list=ls())

## Specify Parameters

#Specify colours for TD, FD, PD
plot.cols <- c("#FFCC4D", "#8CBA70", "#746296")


## Load Data
BaseDirectory <- "..."

setwd( paste( BaseDirectory, "Results", sep="/" ) )
aDiv.AbsChg <- read.csv("aDiv_Results_AbsChg.csv")
aDiv.LogChg <- read.csv("aDiv_Results_LogChg.csv")

#Call in Function from Stats Code  (to calculate change)
setwd( BaseDirectory )
source("BiodiversityChange_Functions.R")


## Function to create an individual histogram with specified settings
individual.hist <- function(aDiv.FF, aDiv.ResultsFF, qqFF, BDaspectFF, resLevel.nameFF, nPlotFF,
                     plotColFF, titleF, Chng.TypeFF="AbsChg" ){

   #Calculate Change Over Time                      
   div.chgFF <- calc.aDivChg( aDiv.FF, Chng.TypeFF )

   if( qqFF=="q0" ){ colChoose <- 1 }
   if( qqFF=="q1" ){ colChoose <- 2 }
   if( qqFF=="q2" ){ colChoose <- 3 }
   
   if( qqFF=="q0" ){ ymax <- 3000 }
   if( qqFF=="q1" ){ ymax <- 2500 }
   if( qqFF=="q2" ){ ymax <- 2500 }
   ax.ticks <- seq(0,ymax,500) 
   
   #Extract out the Effect Size
   aDiv.hist <- aDiv.ResultsFF[ which(aDiv.ResultsFF$resLevels == resLevel.nameFF), ]
   if( resLevel.nameFF != "Plot" ){
      aDiv.hist <- aDiv.hist[ which(aDiv.hist$nPlot == nPlotFF), ]
   }
   aDiv.hist <- aDiv.hist[ which(aDiv.hist$BDaspect == BDaspectFF), ]
   effect.size <- aDiv.hist[ which(aDiv.hist$BDmetric == qqFF), "Chg"]

   #Create Histogram
   if( BDaspectFF == "TD" ){
      hist( div.chgFF[ ,colChoose],  col=plotColFF, ylim=c(0, ymax),
            main=NULL, xlab="Diversity Change", cex.lab=1.5, cex.axis=1.5 )
   }else{
      hist( div.chgFF[ ,colChoose],  col=plotColFF, ylim=c(0, ymax), 
            yaxt="n",
            main=NULL, xlab="Diversity Change", cex.lab=1.5, cex.axis=1.5 )
      axis(2, at=ax.ticks, labels= rep("", length(ax.ticks)) ) 
   }
   
   title( titleF, line=0.75, cex.main=1.75)
   abline( v= effect.size, lwd=2 )    
   if( BDaspectFF == "TD"){
      mtext( "Frequency", side=2, line=3.5, cex=1.1)
   }
   
}   


## Function to create 3-Panel Histograms with TD,FD,PD
panel.hists <- function( qqF, resLevelsF, resLevel.nameF, aDiv.ResultsF, 
                         plotColsF=plot.cols, nPlotF="pMid", Chng.TypeF="AbsChg",
                         saveDirectoryF= paste0(BaseDirectory,"/Results")  ){
   #aDiv.ResultsF: df with the stats values
   #plotColsF: ch with the colours to use in the plots

   
   #Load in the Data & Create Histograms
   setwd( paste( BaseDirectory, "aDivVals", resLevelsF, sep="/") )
   if( resLevelsF == "rQC" ){  
      fileAddOn <- "" 
   }else{
      fileAddOn <- nPlotFF
   }
      
   load( paste0("aDiv_TD_", fileAddOn, "mean.RData") )  
   aDiv.TD <- aDiv.XX
   
   load( paste0("aDiv_FD_", fileAddOn, "mean.RData") )  
   aDiv.FD <- aDiv.XX
   
   load( paste0("aDiv_PD_", fileAddOn, "mean.RData") )  
   aDiv.PD <- aDiv.XX
   
   #Save Figure to Folder
   setwd(saveDirectoryF)
   png( paste0("Fig3_Hist_", Chng.TypeF,"_", resLevel.nameF,"_",qqF,"_",nPlotF, ".png"), 
               width=20,height=8.5,units="cm", res=300)
     
      par( mfrow=c(1,3) )
      par(oma=c(0,6,0,0)) #outer margins
      par(mar=c(6.5, 0, 3, 1)) #inner margins 
   
      #par(mar=c(5, 5, 3, 1) )
      #par( mfrow=c(1,3) )

      individual.hist(aDiv.TD, aDiv.ResultsF, qqF, "TD", resLevel.nameF, nplotF, 
                      plotColsF[1], titleF="Taxonomic", Chng.TypeFF=Chng.TypeF )
      box()

      individual.hist(aDiv.FD, aDiv.ResultsF, qqF, "FD", resLevel.nameF, nplotF, 
                      plotColsF[2], titleF="Functional", Chng.TypeFF=Chng.TypeF )
      box()
      
      individual.hist(aDiv.PD, aDiv.ResultsF, qqF, "PD", resLevel.nameF, nplotF, 
                      plotColsF[3], titleF="Phylogenetic", Chng.TypeFF=Chng.TypeF )
      box()     
   dev.off()   
}#close function

#Example Plots
#(may need to change heights of y axis for viewing at other spatial scales)
panel.hists( "q0", "rQC", "Plot", aDiv.AbsChg)
panel.hists( "q1", "rQC", "Plot", aDiv.AbsChg)
panel.hists( "q2", "rQC", "Plot", aDiv.AbsChg)

panel.hists( "q1", "rQC", "Plot", aDiv.LogChg,  Chng.TypeF="LogChg")



###(5.2) FIG 4 - mean alpha diversity Across Spatial Scales ####################
#Produces Fig 4 for main text, and Fig S3 and S4 in supplemental info
#Fig 4 - mean across scales
#Fig S3 - mean across scales for Percent and Log Changes
#Fig S4 - mean across scales for Low,Med,High Threshold sampling


rm(list=ls())

##Set Parameters
BaseDirectory <- "..."

## Set Plot Styles and Colours
plot.pchs <- c(1,5,8)

td.cols <- c("#FFCC4D","#E6B845","#D1A83F","#BF9939")
fd.cols <- c("#8CBA70","#73995C","#648550","#567345")
pd.cols <- c("#746296","#685887","#55476E","#453A59")
plot.cols <- list( td.cols, fd.cols, pd.cols) 

## Function to create individual boxplot with specified settings
individual.boxplot <- function( aDiv.Results.FF, BDaspectFF, titleFF, plotColsFF, plot.pchFF, 
                                nPlotFF=nPlotF, Chng.TypeFF=Chng.TypeF, y.rngFF=NULL ){
   
   ## Prep df for plotting
   aDiv.Results.FF$BDmetric <- as.factor( aDiv.Results.FF$BDmetric )
   aDiv.Results.FF$resLevels <- factor( aDiv.Results.FF$resLevels, levels=c("Plot","50km","100km","200km") )
   #Extract out data of interest
   aDiv.trim.largeScales <- aDiv.Results.FF[ which(aDiv.Results.FF$nPlot == nPlotFF) , ]
   aDiv.trim.Plot <- aDiv.Results.FF[ which(aDiv.Results.FF$resLevels == "Plot") , ]
   aDiv.trim <- rbind(aDiv.trim.largeScales, aDiv.trim.Plot)  #add row here since Plot scale doesnt have pMid,pLo,pHi
   aDiv.trim <- aDiv.trim[ which(aDiv.trim$BDaspect == BDaspectFF) , ]
   #Order by the Metric and Spatial Scale
   aDiv.trim <- aDiv.trim[ order(aDiv.trim$BDmetric, aDiv.trim$resLevels), ]
   #Add empty row between q0/q1 nad q1/q2
   aDiv.q0 <- aDiv.trim[ which(aDiv.trim$BDmetric == "q0") , ]
   aDiv.q1 <- aDiv.trim[ which(aDiv.trim$BDmetric == "q1") , ]
   aDiv.q2 <- aDiv.trim[ which(aDiv.trim$BDmetric == "q2") , ]
   aDiv.df <- rbind( aDiv.q0, rep(NA,ncol(aDiv.trim)), aDiv.q1, rep(NA,ncol(aDiv.trim)), aDiv.q2 )
   #Add white space for the empty row
   plotColsFF <- c(plotColsFF,"white" ,plotColsFF,"white", plotColsFF)

   #Set Y Ranges:
   if( is.null(y.rngFF) ){
      y.rngFF <- c( min( aDiv.Results.FF$CIlow ), max(aDiv.Results.FF$CIupp ) )
   }

   ## Make Empty BoxPlot (since border set to white)
   if( BDaspectFF == "TD" ){
      boxplot( Chg ~ resLevels * BDmetric, data=aDiv.df ,
         ylim=y.rngFF,
         at=c(1,2,3,4, 6,7,8,9, 11,12,13,14),        #spacing 1
         border="white",
         xlab=NULL, ylab=NULL,
         names= rep( levels(aDiv.df$resLevels), 3), las=2, cex.axis=1.25
      )
   #For FD and PD  (no Y axes labels)
   }else{
      boxplot( Chg ~ resLevels * BDmetric, data=aDiv.df ,
         ylim=y.rngFF,
         yaxt="n",
         at=c(1,2,3,4, 6,7,8,9, 11,12,13,14),        #spacing 1
         border="white",
         xlab=NULL, ylab=NULL,
         names= rep( levels(aDiv.df$resLevels), 3), las=2, cex.axis=1.25
      )
      if( Chng.TypeFF == "AbsChg" ){ axis(2, at=c(0, 0.5, 1), labels=c("","","")) }
      if( Chng.TypeFF == "LogChg" ){ axis(2, at= seq(-0.1, 0.15, 0.05), labels= rep("",6) ) }
   }

   title( titleFF, line=0.75, cex.main=1.75 )
   abline(h=0, lty=2, col="grey52")          
   
   xcoords.xvals <- 1:14                                   
   #Add 95% conf intervals
   arrows(x0=xcoords.xvals, y0=aDiv.df$CIlow, x1=xcoords.xvals, y1=aDiv.df$CIupp, 
            col= plotColsFF, angle=0, lwd=3, length=0.0)
   #Make Plot with Circles
   points( x=xcoords.xvals, y=aDiv.df$Chg, pch=plot.pchFF, cex=1.2, lwd=2, col=plotColsFF ) 
   #Add Lables   
   mtext("q=0", side=1, line=5, adj=0.13, cex=1.2)  #line=4.5
   mtext("q=1", side=1, line=5, adj=0.55, cex=1.2)
   mtext("q=2", side=1, line=5, adj=0.9,  cex=1.2)
   if( BDaspectFF == "TD"){
      mtext( "Mean Alpha Diversity Change", side=2, line=3.5, cex=1.1)
   }
 
   #Add Asterisk if significant
   asterisk.positions <- which( aDiv.df$Signif == "Yes" )
   mtext("*", side=1, line=-1, at=asterisk.positions, cex=1.1)
}


## Function to crate 3-Panel Boxplots with TD,FD,PD
create.boxplots.panel <- function( nPlotF="pMid" , Chng.TypeF="AbsChg" ,
                                   plotColsF=plot.cols , plot.pchF=plot.pchs , 
                                   saveDirectoryF= paste0(BaseDirectory,"/Results")  ){
   #aDiv.ResultsF: df with the stats values
   #plotColsF: ch with the colours to use in the plots
   
   #Load the Results File
   setwd(saveDirectoryF)
   aDiv.Results.xx <- read.csv( paste0("aDiv_Results_", Chng.TypeF, ".csv") )
   
   #Save Figure to Folder
   png( paste0("Fig4_BoxPlot_Scale_",nPlotF,"_",Chng.TypeF,".png"), 
               width=20, height=8.5, units="cm", res=300)
      #par(mar=c(5, 5, 3, 1) )
      par( mfrow=c(1,3) )
      par(oma=c(0,6,0,0)) #outer margins
      par(mar=c(6.5, 0, 3, 1)) #inner margins 

      individual.boxplot( aDiv.Results.xx, "TD", "Taxonomic", plotColsF[[1]], plot.pchF[1], 
                                 nPlotF, Chng.TypeF )
      individual.boxplot( aDiv.Results.xx, "FD", "Functional", plotColsF[[2]], plot.pchF[2], 
                                 nPlotF, Chng.TypeF )
      individual.boxplot( aDiv.Results.xx, "PD", "Phylogenetic", plotColsF[[3]], plot.pchF[3], 
                                 nPlotF, Chng.TypeF )
   dev.off()
}

## Example Figures
create.boxplots.panel( nPlotF="pMid" , Chng.TypeF="AbsChg" )

create.boxplots.panel( nPlotF="pLo", Chng.TypeF="AbsChg" )
create.boxplots.panel( nPlotF="pHi", Chng.TypeF="AbsChg" )

create.boxplots.panel( nPlotF="pMid", Chng.TypeF="LogChg" )



##(5.3) FIG 5 - Explaining Diversity Change - R2 graphs ########################
# This R Script reads in the R2 or R2adj values from the previous ExplainDiversity code
# and creates a plot showing the R2 values across spatial scales, with different
# shapes and colours for metrics and dimensions of diversity.


rm(list=ls())

## Set Up Parameters -------------------------------------------------------------
BaseDirectory <- "..."


#Function to Plot R2 values
plot.R2vals <- function(  time.lagF, my.R2F="R2adj", R2.typeF="mean", mod.coeffsF="Full", nPlotF = "pMid" ){
   #time.lagF:  time lag for climate variables to run:  #"0yr,"10yr","20yr"
   #my.R2F  #"R2" or "R2adj"
   #mod.coeffsF: ch,  "Full" or "Subset"

   #Set Values for Plotting
   BDaspect <- c("TD","FD","PD")
   resLevels <- c("Plot","50km","100km","200km")

   my.cols <- c("#E6B845","#73995C","#685887")
   plot.pch <- c(15, 19, 17)
   xcoords <- matrix( c(0.5, 1.0, 1.5, 
                        3.0, 3.5, 4.0,  
                        5.5, 6.0, 6.5,
                        8.0, 8.5, 9.0), nrow=4, byrow=TRUE)

   ## Read in Data
   setwd( paste0( BaseDirectory, "/Results/", "Explain_aDiv_", mod.coeffsF ) )
   my.results <- read.csv( paste0("ModR2vals__", time.lagF, "_",nPlotF, ".csv") )

   #Make into factors
   my.results$resLevels <- factor( my.results$resLevels, levels=c("Plot","50km","100km","200km") )
   my.results$BDaspect <- factor( my.results$BDaspect, levels=c("TD","FD","PD") )


   ## Create Plot with R2 Values ---------------------------------------------------
   png(filename= paste0("Fig5_MultiReg_",time.lagF,"_",my.R2F,"_", mod.coeffsF, ".png") ,
       width=14, height=9, units="cm", res=300)
   
    par(mar=c(3, 5, 1, 1))

   if( time.lagF == "0yr" ){
       ymx.r2 <- 0.55
       ymx.r2adj <- 0.40
    }else{
       ymx.r2 <- 0.40
       ymx.r2adj <- 0.35
    }
 
    if( my.R2F =="R2" ){ 
       xlab.zz <- expression(paste("R"^"2", " Value", sep="") )  
       ymx <- ymx.r2
    }else if( my.R2F =="R2adj" ){ 
       xlab.zz <- expression(paste("Adjusted R"^"2", " Value", sep="") )  
       ymx <- ymx.r2adj
    }

    #Open Empty Plot
    plot( 1, type="n", xaxt="n", xlim=c(0, max(xcoords)+0.5), ylim=c(0, ymx), 
          xlab="", ylab=xlab.zz )
    axis(1, at= xcoords[ ,2], labels=resLevels)

    for(rr in 1:length(resLevels)){
       my.results.rr <- my.results[ which(my.results$resLevels==resLevels[rr]), ]
       for(aa in 1:length(BDaspect)){
         #Get the Appropriate X locations
         xcoords.aa <- rep(xcoords[rr,aa], 3)
         #Extract the Aspect/Dimension
         my.results.aa <- my.results.rr[ which(my.results.rr$BDaspect == BDaspect[aa]), ]
         #Plot it Out
         my.colnameF <- paste( my.R2F,  R2.typeF, sep="_" )
         points(x=xcoords.aa, y=my.results.aa[ ,my.colnameF], pch=plot.pch, col=my.cols[aa], cex=2) 
      }#lose aa
    }#close rr
    legend( x=0, y= ymx, legend=c("Taxonomic","Functional","Phylogenetic"), bty="n", 
            col=my.cols, lty=1, lwd=3, seg.len=1, cex=1, y.intersp=1)
    legend( x=4.0, y= ymx, legend=c("q=0","q=1","q=2"), bty="n", 
            col="grey45", pch=plot.pch, pt.cex=1.75, cex=1, y.intersp=1)
 
   dev.off() 
}#close function

#Example Figures
plot.R2vals( "0yr",  my.R2F="R2adj", mod.coeffsF="Full", nPlotF = "pMid" )
plot.R2vals( "10yr", my.R2F="R2adj", mod.coeffsF="Full", nPlotF = "pMid" )
plot.R2vals( "20yr", my.R2F="R2adj", mod.coeffsF="Full", nPlotF = "pMid" )



###(5.4) FIG 6 - Alpha Diversity Change after Simulating Forest Loss ###########
# This R script creates Figure 6 in the main text and Figure S7 in the supplementary info,
# showing how mean diversity change estimates change with forest loss simulations.


rm(list=ls())


library(Hmisc)  #minor.tick()

BDmetric <- c("q0","q1","q2")
BDaspect <- c("TD","FD","PD")
PercentRem <- c(0,25,50,75,90,95,100)

BaseDirectory <- "..."

## Read in Data
setwd( paste( BaseDirectory, "Results", sep="/") )
aDiv.Results.reg <- read.csv("aDiv_Results_AbsChg.csv")
aDiv.Results.forest <- read.csv("aDiv_ForestLoss_AbsChg.csv")

#Change to Factors
my.cleaning <- function(aDiv.ResultsF){
   aDiv.ResultsF$BDmetric <- factor( aDiv.ResultsF$BDmetric, levels=BDmetric )
   aDiv.ResultsF$BDaspect <- factor( aDiv.ResultsF$BDaspect, levels=BDaspect )   
   return(aDiv.ResultsF)
}
aDiv.reg <- my.cleaning(aDiv.Results.reg) 
aDiv.forest <- my.cleaning(aDiv.Results.forest) 

#Extract only Plot from the data
aDiv.reg <- aDiv.reg[ aDiv.reg$resLevels == "Plot", ]
colnames(aDiv.reg)[ which( colnames(aDiv.reg)=="nPlot" )] <- "PercentRem"
aDiv.reg$PercentRem <- 0
#Bind into one large dataframe
aDiv.all <- rbind(aDiv.reg, aDiv.forest)


## Create Function to Create an Individual Plot
ForestChangePlot <- function( aDiv.dfF, BDmetricF, addLegendF, 
                              add.Y.labF="Yes", add.X.labF="Yes", lab.sizeF=0.75 ){

   ## Set Plot Parameters
   colsF      <- c("#FFDA7D","#8CBA70","#746296")         #TD,FD,PD
   col.linesF <- c("#E0C06E","#6A8C54", "#645582")
   pchsF <- c(1,5,8)    #TD,FD,PD
   ltysF <- c(1,3,5) 
   line.lwdF  <- 2
   pts.cexF   <- 1
   add.widthF <- 1
   loess.spanF = 0.7
   plot.xlims <- c(0,100)
   plot.ylims <- c(-0.40, 0.30)

   ## Create an Empty Plot
   plot(1, type= "n", xaxt="n", yaxt= "n", 
      xlim = plot.xlims, ylim=plot.ylims, cex.axis=0.9,
      xlab="", ylab="")
   
   axis(1, at = seq(0, 100, by = 25), labels=FALSE) 
   minor.tick(nx=4, ny="", tick.ratio=0.5)
   if(add.X.labF=="Yes"){ 
      mtext("Percent of Trees Removed",    side=1, line=3, cex=lab.sizeF) 
      axis(1, at = c(0,25,50,75,90,100), tick=FALSE, cex.axis=1) 
   }
   axis(2, at = seq( plot.ylims[1], plot.ylims[2], by = 0.1)) 
   abline(h=0, lty=2) 
   if(add.Y.labF=="Yes"){ mtext("Mean Alpha Diversity Change", side=2, line=3, cex=lab.sizeF) }

   ## Add the TD, FD, PD info
   for(aa in 1:length(BDaspect)){
      aDiv.aa <- aDiv.all[ aDiv.all$BDaspect == BDaspect[aa],  ] 
      aDiv.aa.bb <- aDiv.aa[ aDiv.aa$BDmetric == BDmetricF,  ]

      #Create coordinates for the 95% confidence intervals
      xcoords <- c(aDiv.aa.bb$PercentRem, rev(aDiv.aa.bb$PercentRem) )
      xcoords[which(xcoords==100)] <- xcoords[which(xcoords==100)] + add.widthF
      xcoords[which(xcoords==0)]   <- xcoords[which(xcoords==0)]   - add.widthF
      ycoords <- c(aDiv.aa.bb$CIlow, rev(aDiv.aa.bb$CIupp) )
      polygon( xcoords, ycoords, col=adjustcolor( colsF[aa], alpha.f = 0.45), border = NA)
   
      #Add the points back
      points( Chg ~ PercentRem, data=aDiv.aa.bb , pch=pchsF[aa] , cex=pts.cexF , col=col.linesF[aa] )
      lines( Chg ~ PercentRem, data=aDiv.aa.bb , lwd=line.lwdF, lty=ltysF[aa] , col=col.linesF[aa] )
   }

   #Legend
   if( addLegendF == "Yes" ){
      legend( x=0, y= -0.17, legend=c("Taxonomic","Functional","Phylogenetic"), bty="n", 
              col=colsF, pch=pchsF, lty=ltysF, lwd=2)
   }
}


## Figure 6 Style (individual figs)
for( bb in 1:length(BDmetric) ){
   
   png( paste("Fig6_aDiv_LandChg_TTest_",BDmetric[bb],".png", sep="") , 
        width=14, height=12, units="cm", res=300)
      par( mar=c(5,5,1,1) )
      ForestChangePlot( aDiv.all , BDmetric[bb], addLegendF="Yes", lab.sizeF=1 )  
   dev.off()
}

## Figure S7 (stacked figures)
png( "FigS7_aDiv_LandChg_TTest_AllMetrics.png" , 
        width=10, height=16, units="cm", res=300)

   par( mfrow = c(3,1) )
   par( oma=c(5,5,1,1) )
   par( mar=c(0.5,0,1,0) )

   ForestChangePlot( aDiv.all , "q0", addLegendF="Yes", add.Y.labF="No", add.X.labF="No", lab.sizeF=0.75 )
   title( "q=0", line= -1.5, cex=1.3)
   
   ForestChangePlot( aDiv.all , "q1", addLegendF="Yes", add.Y.labF="Yes", add.X.labF="No", lab.sizeF=0.75 )
   title( "q=1", line= -1.5, cex=1.3)
   
   ForestChangePlot( aDiv.all , "q2", addLegendF="Yes", add.Y.labF="No", add.X.labF="Yes", lab.sizeF=0.75 ) 
   title( "q=2", line= -1.5, cex=1.3)

dev.off()


###(5.5) TABLE 1 - Values for Permanova and Permdisp Analyses ##################
#This R Script create a table similar to that of Table 1 in the manuscript.
#Reads in permanova results and makes into a table format. Some basic formatting
#was completed in excel and values from individual permdisp analyses were entered
#manually.


rm(list=ls() )

library(dplyr)


BaseDirectory <- "..."

#Set Up Results Table
td.metric.names <- c("Sorensen", "Bray-Curtis")
pd.metric.names <- c("UniFrac","Dpw")
resLevels.names <- c("Plot","50km","100km","200km")

#Read in File
setwd( paste( BaseDirectory, "Results", sep="/") )
bDiv <- read.csv("PERMANOVADISP.csv")

#Matrix Maker Function
make.matrix <- function( bDiv.fileF, stats.testF, metric.namesF, resLevelsF, BDaspectF, nPlotF="pMid"){
   #stats.testF: "Permanova" or "Permdisp"
   
   #Make Matrix to Fill
   metric.names <- rep(metric.namesF, each=length(resLevelsF))
   resLevels.namesF <- rep(resLevelsF, 2)
   
   results.mat <- as.data.frame( matrix( NA, nrow=2, ncol=length(metric.names) ) )
   rownames(results.mat)[2] <- BDaspectF  
   colnames(results.mat) <- resLevels.namesF
   #Fill with values
   results.mat[1, ] <- metric.names
   #Create a new vector to make the loop filling easier

   #Fill the Results Mat with Values
   for(ii in 1:ncol(   results.mat ) ){         
         #Extract the right value
         bDiv.value <- bDiv.fileF %>% 
            filter( StatsTest == stats.testF ) %>%
            filter( nPlot == nPlotF ) %>% 
            filter( BDaspect == BDaspectF ) %>% 
            filter( resLevels == resLevels.namesF[ii] ) %>%
            filter( BDMetric == metric.names[ii] )
         #Fill the Matrix with this new value
         results.mat[2,ii] <-  bDiv.value[1,1]
   }
   return(results.mat)
}
      
#Run Function for the Different Aspects         
td.permanova <- make.matrix(bDiv, "Permanova", td.metric.names, resLevels.names, "TD", "pMid")
fd.permanova <- make.matrix(bDiv, "Permanova", pd.metric.names, resLevels.names, "FD", "pMid")
pd.permanova <- make.matrix(bDiv, "Permanova", pd.metric.names, resLevels.names, "PD", "pMid")

td.permdisp <- make.matrix(bDiv, "Permdisp", td.metric.names, resLevels.names, "TD", "pMid" )
fd.permdisp <- make.matrix(bDiv, "Permdisp", pd.metric.names, resLevels.names, "FD", "pMid" )
pd.permdisp <- make.matrix(bDiv, "Permdisp", pd.metric.names, resLevels.names, "PD", "pMid" )

#Bind Together into One Big Table
nova.df <- rbind(td.permanova, fd.permanova, pd.permanova)
disp.df <- rbind(td.permdisp,  fd.permdisp,  pd.permdisp)
results.table <- cbind(nova.df, disp.df)

#Save All results
#disp.df just provides Pvalues, not the distance to centroid and the change
write.csv(results.table, file="TableValues_bDivResults_basic.csv")


### (5.6) FIG S5 - Explaining Diversity Coefficients ###########################
#This R Script creates a box and whisker plot for all the variables
#listed in the analysis.


rm(list=ls())

library(dplyr)


## Set Parameters
BaseDirectory <- "..."

time.lag <- "0yr"  #options:  "0yr","10yr","20yr"
  
use.full.mod.coeffs = "Yes"  #options: "Yes" or "No"


##
nPlot <- c("pMid")

resLevels <- c("Plot","50km","100km","200km")
BDaspect <- c("TD","FD","PD")
BDmetric <- c("q0","q1","q2")

#Colours for TD, FD, PD
my.cols <- c("#E6B845","#73995C","#685887")
plot.pch <- c(15, 19, 17)


if( use.full.mod.coeffs == "Yes") {
   add.name <- "Full"
}else{
   add.name <- "Subset"
}


##Read in Data
setwd( paste0( BaseDirectory, "/Results/", "Explain_aDiv_", add.name ) )
aDiv.coeffs <- read.csv( paste0("ModCoeffs__", time.lag, "_",nPlot, ".csv") )


##Set Up Variable Names

#Create Matrix with the Full Names for Axes Labels, in order of plotting
var.names.df <- as.data.frame( matrix( 
                c(1, "Disturbance - Fire", "FirePrc",
                  2, "Disturbance - Harvesting", "HarvestPrc",
                  3, "Protected Area", "PA",
                  4, "Private Land", "Private",
                  5, "Hunting-Fishing Zone", "Hunt.Fish",
                  6, "Temperature Axis 1", "Temperature.Ax1",
                  7, "Temperature Axis 2", "Temperature.Ax2",
                  8, "Precipitation", "Precipitation" ),  ncol=3, byrow=TRUE ) )
colnames(var.names.df) <- c("Num","VariableFullName","Variable")
var.full.names <- var.names.df$VariableFullName
var.names <- var.names.df$Variable

#Add these full names to the coefficients matrix
aDiv.coeffs.full <- merge( aDiv.coeffs, var.names.df, by="Variable", all.x=TRUE )

#Order the matrix for plotting
##Create a Matrix with all the combinations of Levels, Dimensions, MEtrics, Variables
#with variable names in the order that want to plot
#Note: doing this since sometimes a variable was not selected for the model (empty)
new.matrix <- expand.grid( term=var.names.df$Variable,
                           BDmetric = BDmetric,
                           BDaspect = BDaspect,
                           ResLevel = resLevels)                           
#Re-order the column names
new.matrix <- new.matrix[ ,c("ResLevel","BDaspect","BDmetric","term")]
new.matrix <- mutate_if(new.matrix, is.factor, as.character)

##Create Graph Spacing
getSpacing <- function( sF, spacing1F, spacing2F, n.linesF){
   #sF: is the *previous one*
   my.spacing <- seq( max(sF)+spacing2F, 
                      max(sF)+spacing2F+(n.linesF*spacing1F)-spacing1F,  
                      spacing1F)
   my.spacing <- rev(my.spacing)
   return(my.spacing)
}
n.lines <- 9     #3 aspects x 3 metrics: number of lines per potential driver
start.position <- 0.05
spacing1 <- 0.25
spacing2 <- 0.5
s1 <- seq(start.position, n.lines * spacing1,  spacing1)
s1 <- rev(s1)
s2 <- getSpacing( s1, spacing1, spacing2, n.lines)
s3 <- getSpacing( s2, spacing1, spacing2, n.lines)
s4 <- getSpacing( s3, spacing1, spacing2, n.lines)
s5 <- getSpacing( s4, spacing1, spacing2, n.lines)
s6 <- getSpacing( s5, spacing1, spacing2, n.lines)
s7 <- getSpacing( s6, spacing1, spacing2, n.lines)
s8 <- getSpacing( s7, spacing1, spacing2, n.lines)
#Create Spacing Mat
my.spacing.mat <- matrix( c(s8,s7,s6,s5,s4,s3,s2,s1), nrow=8, ncol=9, byrow=TRUE)
colnames(my.spacing.mat) <- c("TD q0","TD q1","TD q2",
                              "FD q0","FD q1","FD q2", 
                              "PD q0","PD q1","PD q2")
rownames(my.spacing.mat) <- var.names.df$Variable
my.ylim <- c(start.position-spacing1, max(my.spacing.mat) + (2*spacing1) )
name.position <- rowMeans(my.spacing.mat)
my.xlim <- 0.6  #1.0

##Graph the Results

#Create Function to Add Shaded Polygon for each variable
poly.yvals <- function(s1F,s2F, prev.valF){
   mm <- mean( c(max(s1F), min(s2F)) )
   return( c(prev.valF, mm, mm, prev.valF) )
}

#Save the Figure
png( filename= paste0("FigS5_MultiRegCoeffs_", add.name,"_", time.lag, ".png") , 
      width=20, height=10, units="cm", res=300)

 par( mfrow=c(1,5) )
 par( oma=c(1,13,1,1) )
 par( mar=c(3.5,0.25,1,0.25) )
 
 #Extract the ResLevel - to a specific graph for each spatial scale
 #Adding +1 to allow for legend at the right
 for(rr in 1:(length(resLevels)+1) ){
   
   #Extract out the ResLevel
   level.rr <- resLevels[rr]
   matrix.rr <- new.matrix[ new.matrix$ResLevel == level.rr, ]
   
   if(rr==5){
      plot(1, type="n", ylim=my.ylim, xlim=c(-my.xlim, my.xlim), 
         axes=FALSE, xlab="", ylab="", yaxs="i")
      legend( x= -0.55, y= 18, 
           legend=c("Taxonomic","Functional","Phylogenetic"), 
           bty="n", cex=1, col=my.cols, lty=1, lwd=3, seg.len=0.85)
      
      legend( x= -0.5, y= 13, 
           legend=c("q=0","q=1","q=2"), 
           bty="n", cex=1, col="grey45", pch=c(15,19,17), x.intersp=1.5)
      
   }else{
      #Open an Empty Plot
      plot(1, type="n", xlim=c(-my.xlim, my.xlim), ylim=my.ylim, 
         xlab="", yaxt="n", ylab="", yaxs="i")
      if(rr==1){ axis(2, at = name.position, labels=var.full.names, las=1 ) }
 
      title(resLevels[rr], line=0.3)

      #Add Shaded Areas for Polygons
      v1 <- poly.yvals(s1, s2, -5)
      v2 <- poly.yvals(s2, s3, max(v1) )
      v3 <- poly.yvals(s3, s4, max(v2) )
      v4 <- poly.yvals(s4, s5, max(v3) )
      v5 <- poly.yvals(s5, s6, max(v4) )
      v6 <- poly.yvals(s6, s7, max(v5) )
      v7 <- poly.yvals(s7, s8, max(v6) )
      v8 <- poly.yvals(s8, s9, max(v7) )
      v9 <- poly.yvals(s9, s10, max(v8) )

      #BotLeft, TopLeft, TopRight, BotRight
      polygon( x=c(-10,-10, 10,10), y=v1, col="grey95", border = NA)
      polygon( x=c(-10,-10, 10,10), y=v3, col="grey95", border = NA)
      polygon( x=c(-10,-10, 10,10), y=v5, col="grey95", border = NA)
      polygon( x=c(-10,-10, 10,10), y=v7, col="grey95", border = NA)
      polygon( x=c(-10,-10, 10,10), y=v9, col="grey95", border = NA)
      box()
      #Add 0 Line
      abline(v=0, lty=1, col="grey45")

      #Loop Through All Potential Combinations
      for(ii in 1:nrow(matrix.rr)){
     
         #See if this variable was selected by AIC (ie foud in in the AIC coefficients file)
         aDiv.coeffs.selected <- aDiv.coeffs %>%
            filter(resLevels == level.rr ) %>% #resLevels[rr]
            filter(BDaspect == matrix.rr[ii,"BDaspect"] ) %>%
            filter(BDmetric == matrix.rr[ii,"BDmetric"] ) %>%
            filter(Variable == matrix.rr[ii,"term"] )
 
         #If there was something selected
         if( nrow(aDiv.coeffs.selected) >= 1 ){     

            ##Find the Right Line to Plot On
            #Get the Line of the variable
            var.row <- which( var.names == matrix.rr[ii,"term"])
            space.vec <- my.spacing.mat[var.row, ]
         
            if(matrix.rr[ii,"BDaspect"]=="TD"){ space.vec <- space.vec[1:3] }
            if(matrix.rr[ii,"BDaspect"]=="FD"){ space.vec <- space.vec[4:6] }
            if(matrix.rr[ii,"BDaspect"]=="PD"){ space.vec <- space.vec[7:9] }
         
            if(matrix.rr[ii,"BDmetric"]=="q0"){ space.vec <- space.vec[1] }
            if(matrix.rr[ii,"BDmetric"]=="q1"){ space.vec <- space.vec[2] }
            if(matrix.rr[ii,"BDmetric"]=="q2"){ space.vec <- space.vec[3] }
         
            ##Get the Right Colours
            if(matrix.rr[ii,"BDaspect"]=="TD"){ plot.col <- "#E6B845" } #gold
            if(matrix.rr[ii,"BDaspect"]=="FD"){ plot.col <- "#73995C" } #green
            if(matrix.rr[ii,"BDaspect"]=="PD"){ plot.col <- "#685887" } #purple
            ##Get the Symbol STyle
            if(matrix.rr[ii,"BDmetric"]=="q0"){ plot.pch <- 15 } #square
            if(matrix.rr[ii,"BDmetric"]=="q1"){ plot.pch <- 19 } #circle
            if(matrix.rr[ii,"BDmetric"]=="q2"){ plot.pch <- 17 } #triangle

            ##Add the Value to the Graph
            points( x=aDiv.coeffs.selected$Slope, y=space.vec, pch=plot.pch, cex=1, col=plot.col )
            #95% Conf Intervals
            conf.upp <- aDiv.coeffs.selected$CIupp 
            conf.low <- aDiv.coeffs.selected$CIlow 
            lines(x=c(conf.low, conf.upp), y=c(space.vec, space.vec), col=plot.col, lty=1, lwd=1.5 )
         }

      }#close ii loop

   }#close else loop
 }#close rr loop
 mtext("Coefficient Value", side = 1, outer = TRUE, line = -1, cex=0.7, adj=0.39 )
 
dev.off()


###(5.7) - Post-Hoc Power Analysis ######################################
#for alpha diversity change

library(pwr)

BaseDirectory <- "..."

setwd( paste0( BaseDirectory, "/Results" ) )
aDiv.results <- read.csv( "aDiv_Results_AbsChg.csv")

aDiv.results$Power <- NA
aDiv.results$Pwr.EffectSize <- NA

for(ii in 1:nrow(aDiv.results) ){

   nn <- aDiv.results[ii, "df"] + 1
   standard.dev <- aDiv.results[ii, "ChgSE"] * sqrt(nn)
   stnd.effect.size <- abs( aDiv.results[ii, "Chg"] / standard.dev )
   
   stats.power <- pwr.t.test(n=nn, d=stnd.effect.size, sig.level=0.05, type="one.sample")
   aDiv.results$Power[ii] <- stats.power$power

   pwr.results <- power.t.test( nn, power=0.80, sig.level = 0.05, type="one.sample")
   aDiv.results$Pwr.EffectSize[ii] <- pwr.results$delta * standard.dev
}



## END

