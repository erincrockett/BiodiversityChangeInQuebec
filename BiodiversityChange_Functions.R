
### CODE INFORMATION ###########################################################

#This R Script includes 11 functions to run the main analyses
#1  calcMean - to calculate mean across all replicates
#2  prepTDalpha.Hill - for taxonomic alpha diversity 
#3  prepPDalpha.Hill - for phylogenetic alpha diversity (for functional diversity as well)
#4  Dpw.Abun - for phylogenetic beta diversity
#5  calc.aDivChg - calculate absolute, percentage, or log change between the two time periods
#6  TTest.aDivChg - run t-test to determine if alpha diversity different between two time periods
#7  plotDiagnostics - plot and save model diagnostic plots
#8  calc.sptlCorrelation - test data/residuals for spatial correlation
#9  create.MeanDist - create mean distance matrix from across all replicates
#10 calcPermNovaDisp - run permanova and permdisp analyses
#11 Forest Loss replicates - create files for all the forest loss simulations

#Written by Erin Crockett
#Updated 2022-Mar-01



##[1] Calculate Mean Alpha Diversity Across All Replicates ---------------------
calcMean <- function( aDiv.listF ){
   #Convert list to an array
   aDiv.arrayF <- array(as.numeric(unlist(aDiv.listF)), 
                                      dim=c( nrow(aDiv.listF[[1]]), ncol(aDiv.listF[[1]]), nrep) )
                  #Calculate the mean across all reps
   aDiv.meanF <- apply( aDiv.arrayF, MARGIN=c(1,2), mean )
   #Add rownames and T1/T2
   aDiv.meanF <- as.data.frame(aDiv.meanF)
   colnames(aDiv.meanF) <- c("q0","q1","q2")
   rownames(aDiv.meanF) <- rownames( aDiv.listF[[1]] )
   aDiv.meanF$timePeriod <- rep( c("t1","t2"), each = nrow(aDiv.listF[[1]]) / 2 )
   nsitesF <- nrow(aDiv.meanF)/2
   aDiv.meanF$SiteID <- rep( rownames( aDiv.listF[[1]] )[1:nsitesF], 2 )
   return(aDiv.meanF)
}


##[2] TAXONOMIC ALPHA DIVERSITY - prepTDalpha.Hill -----------------------------
#Hill numbers of taxonomic diversity
require(vegan)

prepTDalpha.Hill <- function(comm1F, comm2F){

   ## Create presence-absence community matrices, in order to calc sp.richness
   comm1FPA <- comm1F
   comm1FPA[comm1FPA > 0] <- 1
   comm2FPA <- comm2F
   comm2FPA[comm2FPA > 0] <- 1

   ## Species Metrics - time 1  (based on Hill numbers of order 0,1,2)
   aN.SR1 <- rowSums(comm1FPA)					               #Hill Order 0: Calc species richness
   aN.Shan1 <- exp( diversity(x=comm1F, index="shannon") )	#Hill Order 1: Calc shannon div - different value for each row
   aN.Simp1 <- diversity(x=comm1F, index= "invsimpson")     #Hill Order 2: Calc Simpson div (the inverse)

   ## Species Metrics - time 2
   aN.SR2 <- rowSums(comm2FPA)	
   aN.Shan2 <- exp( diversity(x=comm2F, index="shannon") )
   aN.Simp2 <- diversity(x=comm2F, index="invsimpson")
   
   ## Create df with these values (sites a rows, metrics as columns)
   TD.t1 <- data.frame(TDq0=aN.SR1, TDq1=aN.Shan1, TDq2=aN.Simp1) #, timePeriod=rep("t1",nrow(comm1F)) )
   TD.t2 <- data.frame(TDq0=aN.SR2, TDq1=aN.Shan2, TDq2=aN.Simp2) #, timePeriod=rep("t2",nrow(comm2F)) )
  
   #Give both df rownames of the first -- but note that these will change since cant have duplicate rownames
   rownames(TD.t1) <- rownames(comm1F)
   rownames(TD.t2) <- rownames(comm1F)

   #Bind into one large df
   TDa <- rbind(TD.t1, TD.t2)
   return(TDa)
}

##[3] PHYLOGENETIC ALPHA DIVERSITY - prepPDalpha.Hill --------------------------------------------------------
require( entropart ) #Hill numbers for PD (and FD)

prepPDalpha.Hill <- function(comm1F, comm2F, phyloF){
   #comm1F, comm2F:  community data matrices for time 1 and time 2
   #phyloF: the phylogeny (or phylogeny based on trait dendrogram)

   #Write a Function to Calculate PD for a Community Matrix
   #(rather than lines individually), and for q=0,q=1,q=2
   calcHillPD <- function(commF, phyloF2, timeF){
      #Using the ChaoPD function from R package entropart
      #Depends on:  entropart::ChaoPD, vegan::decostand
   
      #Calculate Relative Abundances in Each Comm
      commF.relative <- decostand(commF, method="total")
      #Create Empty Matrix to fill during the loop  - ncol=3 for q=0,1,2 
      pd.mat <- as.data.frame( matrix(NA, ncol=3, nrow=nrow(commF)) )
      rownames(pd.mat) <- rownames(commF)
      colnames(pd.mat) <- c("PDq0","PDq1","PDq2") #, "timePeriod") #,"ID_PE") #Different Hill Numbers
      #pd.mat$timePeriod <- rep(timeF, nrow(commF) )   #Create column that says "t1" or "t2"

      #Calculate Diversity for each community(site) and for q=0,1,2
      for(ppp in 1:nrow(commF)){
         #Extract community to run for PD (ChaoPD needs it as a vector, with names)
         comm.XX <- as.numeric( commF.relative[ppp, ] )
         names(comm.XX) <- colnames(commF)

         for(qqq in 0:2){          #Order q=0,=1,=2
            #Add a correction if statement - to deal with plots after forest change
            #'Regular' communities that have at least 1 individual
            if( sum(comm.XX) > 0 ){
               #Calculate Phylogenetic Diversity of Hill Number order q
               pd.value <- ChaoPD(comm.XX, q=qqq, PhyloTree=phyloF2)
            #IF there are no species in the community (eg after forest change) - diversity is 0
            }else{
               pd.value <- 0
            }
               
            #Put this value in the appropriate spot in the matrix
            pd.mat[ppp,qqq+1] <- pd.value  #add +1 to make sure get the right column (q=0,1,2)
         }
      }
      return(pd.mat)
   }#close internal function

   #Calculate PD at time 1 and time 2
   PD.t1 <- calcHillPD(comm1F, phyloF, timeF="t1")
   PD.t2 <- calcHillPD(comm2F, phyloF, timeF="t2")
   
   PDa <- rbind(PD.t1, PD.t2)
   return(PDa)
}#End Function


##[4] PHYLOGENETY BETA DIVERSITY - Dpw.Abun -------------------------------------------------------
#Based on Code outline in Nate Swenson's Book (2014) Functional & Phylogenetic Ecology in R
#(Abundance Weighted)  (Dpw'  (or prime) in Swenson's Book)
Dpw.Abun <- function(commF, phyloF){

   #Create a phylogenetic distance matrix
   p.dist.mat <- cophenetic(phyloF)
   
   #Extract species that are present
   #note: returns a df - therefore need to use as.numeric() below
   get.abuns <- function(x){
      x[x>0]
   }
   
   #Calculate relative abundance
   commF.ra <- commF / rowSums(commF)
   
   #Apply function over all rows of community data matrix (my.sample)
   #To get a list containing names of species present in each of the communities

   list.of.abuns <- apply(commF.ra, 1, get.abuns)
   #
   Dpw.apply.function <- function(x){
      tmp.function <- function(z){
         #Abundance Weighted
         sum( p.dist.mat[names(x),names(z)] * outer( as.numeric(x), as.numeric(z) ) )
      }
      lapply(list.of.abuns, FUN=tmp.function)
   }
   dpw.output <- lapply(list.of.abuns, Dpw.apply.function)
   output.matrix <- do.call(cbind, dpw.output)
   output.dist <- as.dist(output.matrix)
   
   return(output.dist)
}


##[5] CALCULATE ALPHA DIVERSITY CHANGE -----------------------------------------
#Calculates the change in alpha diversity for absolute, percentage, or log change

calc.aDivChg <- function( aDiv.FF, Chng.TypeFF ){ 
   ##Calculate Change between the two time periods
   aDiv.FF$timePeriod <- NULL
   aDiv.FF$SiteID <- NULL
   nsitesF <- nrow(aDiv.FF)/2
   x1 <- head( aDiv.FF, nsitesF )
   x2 <- tail( aDiv.FF, nsitesF )
   
   if( Chng.TypeFF == "AbsChg" ){  
      div.chgF <- x2 - x1
   }
   if( Chng.TypeFF == "LogChg" ){
      div.chgF <- log(x2/x1)
   }
   rownames(div.chgF) <- rownames(x1)
   return(div.chgF)
}


##[6] ALPHA DIVERSITY -- T-TEST for Absolute, Percentage, or Log Change --------   
#Conducts T-Test for significant changes in alpha diversity between the two time periods

TTest.aDivChg <- function( aDiv.F, Chng.TypeF ){ 
   #aDiv.F: df, aDiv.XX
   #Chng.TypeF:  "PrcChg" or "LogChg" or "AbsChg"

   #Calculate Change (use function above)
   div.chg.df <- calc.aDivChg( aDiv.F, Chng.TypeF )

   ## Create empty object to fill during the qq loop  (nrow= for q=0,1,2)
   ttest.cols <- c("Chg", "ChgSE","Pval","Signif", "CIlow", "CIupp", "testStat", "df")
   ttest.results <- as.data.frame( matrix(NA, nrow=3, ncol=length(ttest.cols)) )
   colnames(ttest.results) <- ttest.cols

   ## One Sample T-Test -- determine if Diversity Change different than zero
   #Clear this variable (if its used elsewhere in code)
   qq <- NULL  #using qq inside the function since bb used outside the function
   #For q0,q1,q1
   for(qq in 1:3){
      div.chg.vec <- div.chg.df[ ,qq]  #Select out the corresponding column
      t.results <- t.test( div.chg.vec, alternative="two.sided" )

      div.chg.SE <- ( t.results$conf.int[2] - mean(div.chg.vec) ) / 1.96
      t.stats.vals <- c(t.results$estimate, div.chg.SE, t.results$p.value, 0, t.results$conf.int[1], t.results$conf.int[2], 
                           t.results$statistic, t.results$parameter)
      ttest.results[qq, ] <- t.stats.vals
   }
   
   #Add Character Column Here for Significant to Make Plotting Easier Later
   #Writing another loop since the rest of items are numeric above
   for(qq in 1:3){
      if( ttest.results[qq,"Pval"] < 0.05 ){
         ttest.results[qq,"Signif"] <- "Yes" 
      }else{
         ttest.results[qq,"Signif"] <- "No"
      }
   }

   #Add column values to note which metrics,aspects,dimensions and stats test they are
   ttest.results$BDmetric <- c("q0","q1","q2")

   return(ttest.results)
} 


##[7] Plot Diagnostics from a Linear Model -------------------------------------
plotDiagnostics <- function( model.object, mod.type, fnameF2, mod.formula=NULL){
   #model.object is m.lm or m.best
   #mod.type: is either "mlm" or "mmixed"  (linear or mixed)
   #fnameF2: character - is the name to paste when saving variables   

   png(filename= paste(fnameF2, mod.type, ".png", sep="" ),
     width=20,height=20,units="cm", res=300)
  
      par(mfrow=c(2,2))   

      ### Normality of Residuals
      qqnorm( resid(model.object))
      qqline( resid(model.object), lwd=2, col="blue")
      if(mod.type=="mmixed"){
         title(main= mod.formula , line=-1, cex.main=0.8, adj=0.01)
      }
      
      ## Distribution of residuals
      res.m <- resid(model.object)
      hist(res.m, freq=FALSE, main="Distribution of Residuals")
      #Get values to add lines overtop
      xfit <- seq( min(res.m), max(res.m), length=40)
      yfit <- dnorm(xfit)
      lines(xfit, yfit)

      #Plot Residuals by fitted values
      res.m <- resid(model.object)
      fit.m <- fitted(model.object)
      plot(x=fit.m, y=res.m, xlab="Fitted values", ylab="Normalized Residuals")
      abline(0,0, lty=2)

  dev.off()

}#Close Function


##[8] CALCULATE SPATIAL CORRELATION --------------------------------------------
#Function is used to calculate spatial correlation for a vector (eg the residuals 
# of a model). It plots a Variogram and Morans I correlelogram.
require(gstat)

calc.sptlCorrelation <- function( residualsF, fnameF2, my.cutoffF=1000000, dClassesF=60){ 
   #residualsF: the model residuals dataframe
   #dClassesF:  number of classes 
   #my.cutoffF: max upper value
   #fnameF2: character string - name to sae the files with
   
   #Create Spatial Points Object
   colnames(residualsF)[4] <- "Yvar"
   resids.sptl <- SpatialPointsDataFrame(coords=residualsF[ ,c("y","x")], data=residualsF)

   ## Calculate variogram model
   vario.resids <- gstat::variogram(Yvar ~ 1, resids.sptl, cutoff=my.cutoffF, width=my.cutoffF/dClassesF) #, cutoff=distanceF, width=distance.blocksF)

   ## Calculate Morans I by distance classes
   #Use the function the G.Laroque Wrote: (note: somewhat slow to calculate)
   corgramI <- MoranI(cbind(resids.sptl$x, resids.sptl$y), resids.sptl$Yvar, dclasses=dClassesF ) 

   ## Graph Both Results:
   png( paste(fnameF2, "_Variogram.png", sep="") , width=15, height=15, units="cm",res=200)
      plot(vario.resids)
   dev.off()

   png( paste(fnameF2, "_MoransI.png", sep=""), width=15, height=15, units="cm",res=200)
      plot(corgramI)
   dev.off()

   return(corgramI)
}#close function


##[9] CREATE MEAN DISTANCE MATRIX ---------------------------------------------
#Across all replicates, create a mean distance matrix
#has component to deal with missing rows for Dpw

create.MeanDist <- function( BDaspectF, distMethodF, nPlotF, resLevel.foldrF, FileDirF ){
   
   #Set the Directory to Read in the Files
   setwd( paste( FileDirF, resLevel.foldrF, nPlotF, sep="/" ) )

   #Create a Vector of ALL the files within this directory, then create a smaller subset for aDiv and bDiv separately
   filesInFolder <- list.files( )
   stringLookup <- paste( BDaspectF, distMethodF, sep="_")  #Create the string to look for, based on TD,PD,FD
   bDivFiles.int <- filesInFolder[grepl(stringLookup, filesInFolder)] #Extract all the files that START wtih aDiv

   #Remove any "trimmed names"
   bDivFiles.int <- bDivFiles.int[! grepl("*trimmed", bDivFiles.int)] #Extract all the files that START wtih aDiv
   #Remove any "rowsToRemove" names
   bDivFiles <- bDivFiles.int[! grepl("*Remove", bDivFiles.int)] #Extract all the files that START wtih aDiv

   #Read in Template File if have files to remove
   #Note: use UniFrac here as the template file for Dpw since may have sites with 0s
   if( distMethodF == "Dpw" ){
      stringLookup.temp <- "PD_UniFrac"
   }else{
      stringLookup.temp <- stringLookup
   }
            
   ##Create an array to store the data
   #load an initial file to get the dimensions
   load( paste(stringLookup.temp,"_nn1.RData", sep="") )   #object: dMat.XX
   dMat.XX.template <- dMat.XX #Save the template file for later
   rm(dMat.XX)
   
   #Create the array
   #Note: for bDiv, it has all sites listed twice (t1,t2) - unlike the aDiv files
   nsites <- nrow( as.matrix( dMat.XX.template ) )  
   bDiv.array <- array(NA, dim=c(nsites, nsites, length(bDivFiles) ) )
   dimnames(bDiv.array)[[1]] <- labels(dMat.XX.template)
   dimnames(bDiv.array)[[2]] <- labels(dMat.XX.template)
   
   #Load files an add Files to correct position in array
   for(ff in 1:length(bDivFiles)){
            
      load( bDivFiles[ff] ) #object: dMat.XX
      dMat.matrix <- as.matrix(dMat.XX)

      #Check to see if the lengths are different (if they are, thats because of missing rows)
      if( length(labels(dMat.XX)) != length(labels(dMat.XX.template) )  ){
                     
         #Find any sites that are missing
         label.add <- setdiff(labels(dMat.XX.template), labels(dMat.XX))
         #Add a column with the site IDs in order to order the values later on
         dMat.df <- as.data.frame(dMat.matrix)
         dMat.df$my.label <- rownames(dMat.df)

         #Add X extra rows for the missing sites
         df.row <- as.data.frame( matrix(NA, nrow=length(label.add), ncol=ncol(dMat.df) ) )
         rownames(df.row) <- label.add
         colnames(df.row) <- colnames(dMat.df)
         df.row$my.label <- label.add
         dMat.df <- rbind(dMat.df, df.row)
  
         #Add X extra columns for the missing sites
         df.col <- matrix(NA, nrow=nrow(dMat.df), ncol=length(label.add) ) 
         colnames(df.col) <- label.add
         dMat.df <- cbind(dMat.df, df.col)

         #Order ROWS to have the same rows as the original template file
         dMat.match <- dMat.df[ match(labels(dMat.XX.template), dMat.df$my.label),  ]
         #Order COLUMNS to have same colmum order as original
         dMat.match <- dMat.match[  ,labels(dMat.XX.template)]

         ## check to make sure all rows and columns in the same order
         if( all.equal( rownames(dMat.match), colnames(dMat.match), labels(dMat.XX.template) )  ){
            #Fill the array in the replicate position
            bDiv.array[ , ,ff] <- as.matrix( dMat.match )
            #If Not then give an error message:   
         }else{
            print( paste("Error 1: Missing Labels -- ff", ff, sep=" ") )
            break
         }
                     
      #If the Lengths are Equal Already   (No missing rows)
      }else{  
         #Check to make sure that they are in the same order   
         if( all.equal( rownames(dMat.matrix), colnames(dMat.matrix), labels(dMat.XX.template) )  ){
            #Fill the array in the replicate position
            bDiv.array[ , ,ff] <- dMat.matrix
         }else{
            print( paste("Error 2: Same Lengths, Wrong Order -- ff", ff, sep=" ") )
            break
         }
      }
      #Remove any existing distance matrices
      rm(dMat.XX)
   }#close ff loop

   #Calculate the MEAN across all the replicates, and do it for all the different columns 
   bDiv.mean <- apply( bDiv.array, MARGIN=c(1,2), mean, na.rm=TRUE)
   dMat.XX.new <- as.dist(bDiv.mean)

   #Clear Objects Between Runs
   rm(bDiv.mean, bDiv.array, bDivFiles, dMat.XX.template)
            
   return( dMat.XX.new )
}

         
##[10] PERMANOVA - PERMDISP ----------------------------------------------------
#Test for differences in centroid between groups, and differences in dispersion

calcPermNovaDisp <- function(dMatF, timeMatF, npermF, BDaspectF, BDmetricF, nPlotF, saveDirF ){
   #dMatF: distance matrix (calculated earlier)
   #timeMatF: df - with SiteID and time t1 or t2 for each site
   #npermF:  number of replicates to run
   #BDaspectF: BD aspect (ie dimension)      
   #BDmetricF:  BD "metric" or "distance method"
   #nPlotF: whether its "pMid","pLow","pHi"
   #saveDirF: directory where to save the files
   
   ## PERMANOVA analysis 
   #Statistics to asses if there are differences in centroid between the two time periods
   #Note: adonis function cannot use Sorensen input directly, therefore supplying the distance matrices
   #Extract the stats variables - Pval, FStat, R2, df residuals
   permMan.XX <- adonis(dMatF ~ timePeriod, data=timeMatF, strata=timeMatF$site )
   permMan.XX.vals <- c( permMan.XX$aov.tab$Pr[1], permMan.XX$aov.tab$F.Model[1], 
                            permMan.XX$aov.tab$R2[1], permMan.XX$aov.tab$Df[2], NA, NA)  #NAs since more info in permdisp below

   ## PERMDISP analysis 
   #Statistics to asses if there are differences in dispersion between two time periods
   #Use a Distance matrix with both time periods together in order to constrain the permutations
   #by time 1 and time 2
   permDisp.XX <- betadisper(d=dMatF, group=timeMat$timePeriod, type="centroid")
   #Calc Stats by permutation
   permDisp.XX.perm <- permutest(permDisp.XX, pairwise = TRUE, permutations = npermF)
   # Extract the stats values (Pval and F)  
   #Note: adding an NA, since function doesnt provide stats equivalent to R2 of permanova
   permDisp.XX.vals <- c( permDisp.XX.perm$tab[1,6], permDisp.XX.perm$tab$F[1], 
                             NA,  permDisp.XX.perm$tab$Df[2] ) 

   ## Get Avg Distance to centroid and change over time
   nsitesF <- length( labels(dMatF) )/2
   avg.cent <- mean( permDisp.XX$distances )
   cent.chg <- mean( tail(permDisp.XX$distances, nsitesF) ) - mean( head(permDisp.XX$distances, nsitesF) )
   permDisp.XX.vals <- c(permDisp.XX.vals, avg.cent, cent.chg)
   
   ## Save Outputs of Function to a File
   setwd( saveDirF )
   bDiv.vals <- list(permMan.XX, permDisp.XX, permDisp.XX.perm )   
   names(bDiv.vals) <- c("Permanova","Permdisp","PermdispPermuts")
   save(bDiv.vals, file= paste(BDaspectF,"_",BDmetricF,"_",nPlotF, ".RData", sep="") )   
    
   ## Return Function Values 
   #Group objects together as a single dataframe, along with info about the stats test and BDmetric used
   bDiv.stats <- rbind(permMan.XX.vals, permDisp.XX.vals)
   bDiv.stats <- data.frame( bDiv.stats, c("Permanova","Permdisp"))
   colnames(bDiv.stats) <- c("Pval","Fstat","R2","Df", "CentroidAvg","CentroidChg", "StatsTest")

   return(bDiv.stats)
}#end function



##[11] RUN FOREST LOSS REPLICATES ----------------------------------------------
#This function will remove X% of trees from a selected % forest plots across Quebec,
#based on the amount of disturbance in each Biological domain.

runForestL.reps.IndividualTrees <- function( comm1F, comm2F, phyloF, BDaspectF, 
              percent.RemF, percent.DisturbF, 
              info.firstF, time.multiplierF, n.replicatesF){
  
   #comm1F:  community data matrix at T1
   #comm2F:  community data matrix at T2
   #phyloF:  phylogeny (PD) or trait dendrogram based phylo object (FD)
   #BDaspectF:  ch, "TD","FD","PD"
   #percent.RemF: num,  % of individuals/trees in community to remove
   #percent.DisturbF:  vec,  % of disturbance in each of the Biodomains
   #info.firstF:  info.first to indicate which grid cells are within each Biodomain
   #time.multiplierF: num, to account that disturbance data only spans part of the study period
   #n.replicatesF: num, # of replicates for forest loss
     
     
   #Get proportion forest loss for each biodomain -- multiply by the multiplier to account for entire study period
   forest.prop.loss <- time.multiplierF * ( percent.DisturbF / 100) #/100 tochange from % to proportion   
   forest.prop.loss[7] <- 0   	 #For Biodomain 7 - where very very few sites - not included in disturbance calcs

   #Read in Names for Each Region 
   #In this case: its Biodomain - but could use 50km or other grid cells
   unique.regions <- c("1","2","3","4","5","6","7")
  	
 	#Create Array to Store all the Replicates
 	aDiv.array <- array(NA, dim=c( 2* nrow(comm1F), 3, n.replicatesF ) )  #3 cols for q0,q1,q2

 	for(gg in 1:n.replicatesF){

 	   #Create an empty datafame that will rebuild during the loop of 50km sites
 	   comm2.forestloss <- as.data.frame(matrix(NA, nrow=0, ncol=ncol(comm2F) ) )
 	   colnames(comm2.forestloss) <- colnames(comm2F)

      #For the List of different regions (biodomains or 50km sites)
      for(hh in 1:length(unique.regions) ){
 	              
 	       comm2.hh <- NULL  #clear from previous
 	       #Extract forest plots within this region (BioDomain)
 	       sites.region.hh <- info.firstF[ info.firstF$gCellBioDomain %in% unique.regions[hh],  ]
 	       #Find these plots in the Community Data matrix & extract these rows
 	       positions.region.hh <- which( rownames(comm1F) %in% sites.region.hh$ID_PE )
 	       comm2.hh <- comm2F[positions.region.hh, ]
          #Randomly Select X% of sites in comm2 from which to remove individuals
 	       forest.rem.hh <- round( forest.prop.loss[hh] * nrow(comm2.hh), digits=0 ) #round to a whole number
 	    
 	       #Only run *IF* there was enough forest loss in that region to remove at least 1 plot (since rounded to whole number above)
 	       if( forest.rem.hh > 0 ){ 
 	          #Randomly sample these
 	          rows.hh <- sample( x=1:nrow(comm2.hh), size = forest.rem.hh, replace=FALSE)
	          #Extract these rows
 	          comm2.hh.0 <- comm2.hh[rows.hh, ]
       
 	          #Use the Short(er) way if all diversity is going to zero (removing all trees)
 	          if( percent.RemF == 100 ){
	              comm2.hh.0[ comm2.hh.0 > 0 ] <- 0   #Make BD go to zero at time 2
 	              comm2.hh.0.new <- comm2.hh.0
               
 	          #Otherwise use my long code...     
 	          }else{
 	             #Create empty matrix to grow during the loop
 	             comm2.hh.newdf <- matrix(NA, nrow=0, ncol=ncol(comm2.hh.0))
 	             colnames(comm2.hh.newdf) <- colnames(comm2.hh.0)
 	             #Loop through all rows that were selected to experience disturbance
 	             for( mmm in 1:nrow(comm2.hh.0)){
 	                #Find positions of which species are present
 	                sp.positions.mmm <- which(comm2.hh.0[mmm, ] > 0)
 	                #Find the colnames where there are species
                  sp.names.mmm <- colnames(comm2.hh.0)[sp.positions.mmm]
                  #Find the Values of these species
                  nvals.mmm.int <- comm2.hh.0[mmm, sp.positions.mmm]
                  nvals.mmm <-  as.numeric(nvals.mmm.int )
                  names(nvals.mmm) <- colnames(nvals.mmm.int)
                     
                  ##Grow a Vector With All these Names (repeated by their occurrence)
                  all.sp.names <- NULL
                  for( nnn in 1:length(nvals.mmm)){
                     all.sp.names <- c(all.sp.names, rep( sp.names.mmm[nnn], nvals.mmm[nnn]))
                  }
   
 	               ### Randomly Sample X% of these
                  how.many.to.keep <- length(all.sp.names) * (1 - (percent.RemF/100)) #Keep = 1-Remove, and need proportion rather than %
                  trees.to.keep <- sort( sample(1:length(all.sp.names), size=how.many.to.keep, replace=FALSE) )
 	               #Select out just these trees from the full vector
                  remaining.trees <- all.sp.names[trees.to.keep]
                  
                  ### Recreate the Community Matrix
                  #Get a vector of how many instances of each name are within
                  trees.table <- table(remaining.trees)
 	               trees.vec <- as.vector(trees.table)
 	               names(trees.vec) <- names(trees.table)
 	               #Find out which species were not present - and need to be added back to the vector
 	               c.nam.add <- setdiff( colnames(comm2.hh.0), names(trees.vec)   )
 	               #Add these values
 	               c.vals.add <- rep(0, length(c.nam.add))
 	               names(c.vals.add) <- c.nam.add
 	               #Bind these together into one long vector (with sp presnt and not present)
 	               final.vec.int <- c(trees.vec, c.vals.add)
 	               #Order this vector to be the same as colnames of original df
                  final.vec <- final.vec.int[ colnames(comm2.hh.0) ]

 	               #Add this to the growing vector:
 	               comm2.hh.newdf <- rbind(comm2.hh.newdf, final.vec)
                	                
 	             }#close for mmm loop
 	             #Make sure add back the rownames with the plot IDs
 	             rownames(comm2.hh.newdf) <- rownames(comm2.hh.0)

               #Rename this as the old matrix-- same as above
               comm2.hh.0.new <- comm2.hh.newdf 
            }#close else loop 
                        
 	         #Extract the Remaining Rows
 	         comm2.hh.1 <- comm2.hh[ -rows.hh, ]
 	         #Re-Bind Together
 	         comm2.hh <- rbind(comm2.hh.0.new, comm2.hh.1)
            
 	       }#close if not removing All tree
          comm2.forestloss <- rbind(comm2.forestloss, comm2.hh)
         
      }#close for hh (number of unique regions/biodomains)
 	    
   
 	   #Order the rows - to make sure that these are in the same order for aLL REPLICATES
 	   comm1F <- comm1F[ order(rownames(comm1F)),  ]
      comm2.forestloss <- comm2.forestloss[ order(rownames(comm2.forestloss)),  ]

      #CHECK to See That All are the same - break if they are not
      if( !all.equal(rownames(comm1F), rownames(comm2.forestloss) ) ){ 
         print("Check rownames comm1F and comm2.forestloss")
         break 
      }
 	   if( !all.equal( colnames(comm1F), colnames(comm2.forestloss) ) ){ 
         print("Check colnames comm1F and comm2.forestloss")
         break 
      }
      if( !all.equal(colnames(comm1F), phyloF$tip.label ) ){ 
         print("Check colnames comm1F and phylo tip labels")
         break 
      }
      
      #### Calculate BD Change
 	   #Taxonomic Diversity
      if( BDaspectF=="TD" ){ 
 	       aDiv.XX <- prepTDalpha.Hill(comm1F, comm2.forestloss)
 	       #Note: Hill q=2 gives Inf when there are no species in the community
 	       #(like in the ones I just erased) Therefore, replace these values with 0s
 	       aDiv.XX[ aDiv.XX == Inf ] <- 0

 	   #Phylogenetic of Functional Diversity   
 	   }else{    
         aDiv.XX <- prepPDalpha.Hill(comm1F, comm2.forestloss, phyloF) 
      }
      
 	   #Store Results in Array with all replicates
      aDiv.array[ , ,gg] <- as.matrix( aDiv.XX )
    
      print( paste( "--completed rep", gg) )
   }#close for gg replicates
	 
   dimnames(aDiv.array)[[1]] <- rownames(aDiv.XX)
   dimnames(aDiv.array)[[2]] <- colnames(aDiv.XX)

 	return(aDiv.array)
}


