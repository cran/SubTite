#' Gives the subgroup specific optimal dose vector.
#'  Returns a list containing the optimal doses to enroll each subgroup at and the subgroups that should have their accrual suspended temporarily.
#' @param Y Vector containing observed event or censoring times.
#' @param I Vector containing event indicators (1 if patient experiences an event for a patient).
#' @param Doses Vector containing Doses of patients in trial.
#' @param Groups Vector containing group assignment of patients, 1 is baseline group.
#' @param T1 Reference time for toxicity.
#' @param Target Target cumulative toxicity probability vector at time T1.
#' @param Dose Vector containing the standardized doses considered.
#' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
#' @param DoseTried  Matrix that contains counts in each subgroup for the number of times each dose has been assigned.
#' @param cohort Number of patients needed to be assigned at a dose level prior to escalation.
#' @param meanmu Prior mean for baseline intercept.
#' @param meanslope Prior mean for baseline slope.
#' @param MeanInts Vector of prior means for the group specific intercept parameters.
#' @param MeanSlopes Vector of prior means for the group specific slope parameters.
#' @param varint Prior variance for the intercept parameters.
#' @param varbeta Prior variance for the slope parameters.
#' @param phetero Prior probability of heterogeneous subgroups.
#' @param Borrow Parameter to specify subgroup borrowing/clustering. 0=No borrowing, 1=Borrowing but no clustering, 2=Borrowing and clustering.
#' @param B Number of Iterations to run for MCMC
#' @return Returns a list with two objects, a vector of optimal doses for each subgroup and matrix of posterior toxicity probabilities at each dose level within each subgroup.
#' @references
#' [1] Chapple and Thall (2017), Subgroup Specific Dose Finding in Phase I Clinical Trials Based on Time to Toxicity Within a Fixed Follow Up Period.
#' @examples
#' T1=6
#' ## Reference Time for Toxicity
#' ##Number of subgroups
#' nGroups=3
#' ##Number of Doses
#' nDose=4
#' ##What is the starting dose? We want to set all dose tried > cohort that's less than this
#' DoseStart=4
#' ##Sample Size
#' n=90
#' Target=rep(.3,nGroups)
#' Upper=rep(.95,nGroups)
#' Y=rep(NA,n)
#' I=rep(NA,n)
#' Groups = sample(1:nGroups,n,replace=TRUE) - 1
#' ##Group assignment of patients (MUST BE CODED 0,1,2,...)
#' Doses = sample(1:DoseStart,n,replace=TRUE)
#' Dose2=Doses ##Going to hold the numeric dose numbers
#' ##Randomly Generate Dose values
#' x=sort(runif(nDose))  ##Doses are in ascending order
#' Dose=(x-mean(x))/sd(x)
#' ##Vector of standardized doses
#' ##Randomly generate TRUE group probabilties
#' GroupProb = matrix(ncol=nDose,nrow=nGroups)
#' for(k in 1:nGroups){
#' GroupProb[k,]=sort(runif(nDose,0,Target[k]+.2))
#' }
#' ##Randomly generate patient data from a uniform TTE dist with given probabilities.
#' for(b in 1:n){
#' I[b]= rbinom(1,1   , GroupProb[(Groups[b]+1),Dose2[b]])
#' if(I[b]==0){ Y[b]=T1 }else{ Y[b]=runif(1,0,T1) }}
#' ##How many patients in each subgroup have been assigned at each dose level?
#' DoseTried=cbind(table(Groups,Doses))
#' cohort=1 ##Cohort size required for escalation
#' ##Matrix of umber of patients tried or fully evaluated at each dose level.
#' DoseTried[,1:DoseStart]=cohort
#' ##Hyperparameters
#' meanmu=-0.4467184 ##Common Intercept hypermean
#' meanslope= 0.8861634 ##Common slope hypermean
#' MeanInts = rep(-0.5205379,nGroups-1) ##Group Intercept hypermeans
#' MeanSlopes = rep(0.1888923,nGroups-1) ##Group slope hyperneabs
#' Dose2=Doses ##Numeric Doses
#' Doses=Dose[Doses]  ##Standardize dose value
#' varint=5 #Prior Variance of the intercept betas
#' varbeta=1 ##Prior Variance of slope betas
#' phetero=.9 ##Prior Probability of hetergeneity
#' Borrow=0 ##Borrowing specification, 0=none, 1=some, 2=clustering.
#' B=5000 ##Number of iterations
#' Borrow=2
#' GetSubTite(Y, I,Doses, Groups, DoseTried,cohort, T1,
#'         Target,  Upper, Dose,  meanmu, meanslope,
#'        MeanInts,  MeanSlopes ,varint,varbeta,phetero, Borrow,B)
#'@export
GetSubTite=function(Y, I,Doses, Groups, DoseTried, cohort,T1, Target,
                    Upper, Dose,  meanmu, meanslope,
                    MeanInts,  MeanSlopes ,varint,varbeta,phetero,Borrow,B){

  ERRHOLD=c(length(Target), nrow(DoseTried), length(Upper), length(MeanInts)+1, length(MeanSlopes)+1)

  HOLD=0
  ##Check for errors in dimension specification
  for(k in 1:length(ERRHOLD)){
    for(m in 1:length(ERRHOLD)){
      if(ERRHOLD[k] != ERRHOLD[m]){
        HOLD=1
      }
    }
  }

  if(HOLD==1){
    message("Dose tried matrix,target toxicity vector, toxicity threshold, or subgroup hyperparameter vector has incorrect dimensions")
  }else{

    if(min(Groups)>0){
      message("Be sure groups are labeled g=0,1,...,G-1")
    }

    ##Reset dose tried.






    ##Repackage MeanInts and MeanSlopes
    MeanInts=c(0,MeanInts)
    MeanSlopes=c(0,MeanSlopes)
    Stopped=rep(0,nrow(DoseTried))
    ##This matrix contains posterior mean toxicity probabilities at each dose for each subgroup.
    ##The last column in the matrix has whether or not each group should be stopped.
    RESULTS=MCMC( Y,I,  Doses,  Groups,  T1,  Target,  Upper, Dose, meanmu,  meanslope,
                  MeanInts,  MeanSlopes, varint,  varbeta, phetero, Stopped,  length(Y),  Borrow,B)

    CLUST=RESULTS[[2]]


    RESULTS1=RESULTS ##HOLDER FOR STOPPED GROUPS
    RESULTS=RESULTS[[1]]
    Stopped= RESULTS[,ncol(DoseTried)+1]
    ##Get optimal Dose
    OptDose= Stopped
    ##Check if ALL groups are stopped, if so run a separate trial in each.
    if(Borrow>0){
      if(sum(Stopped)==length(Stopped)){
        message("Borrowing has caused all groups to stop due to excessive toxicity.
Separate models will be fit to each group to ensure the design is not stopping too early due to borrowing.")
        RESULTS=MCMC( Y,I,  Doses,  Groups,  T1,  Target,  Upper, Dose, meanmu,  meanslope,
                      MeanInts,  MeanSlopes, varint,  varbeta, phetero,    Stopped=rep(0,nrow(DoseTried)),  length(Y),  0,B)
        RESULTS1=RESULTS ##HOLDER FOR STOPPED GROUPS
        RESULTS=RESULTS[[1]]
      }
    }



    PROBS = data.frame(RESULTS[,1:ncol(DoseTried)])  ##Used for looping over probs
    ##Can only choose a dose level if we've escalated correctly
    PROBS1=PROBS
    Y1=DoseTried<cohort ##Flag which dose levels haven't been treated enough.
    ##Cant escalate past that dose

    for(k in 1:length(Stopped)){
      j=1
      if(sum(1-Y1[k,])<ncol(DoseTried)){
        ##Checks if all doses meet cohort criteria
        while(Y1[k,j]==FALSE){
          j=j+1
        }
        ##We can go up one more
        j=j+1
        ##Reset PROB1 with -1000, so now we can't pick it
        PROBS1[k,j:ncol(DoseTried)]=-10000
      }
    }




    ##Now get optimal doses
    for(k in 1:length(Stopped)){
      if(Stopped[k]==0){
        a1 = abs(PROBS1[k,]-Target[k])
        OptDose[k]=which(a1==min(a1)) ##Minimum distance from Target probability
      }else{
        OptDose[k]=NA ##Na for stopped groups
      }
    }



    Z=as.list(c(0,0,0,0))


    Z[[1]]=OptDose
    for(k in 1:nrow(PROBS)){
      rownames(PROBS)[k]=paste0("Subgroup ", k)
    }
    ##Columns by dose levels
    colnames(PROBS)=1:ncol(DoseTried)

    Z[[2]]=PROBS
    ##Posterior probability of stopping
    Z[[3]]=RESULTS[,ncol(RESULTS)]

    names(Z)=c("Optimal Dose","Posterior Mean Toxicity Probability","Posterior Probability of Overly Toxic Subgroup","Clustering Parameters")

    if(Borrow==2){

      ###Cluster membership...
      CLUST=CLUST+1
      ##Cluster assignment
      CLUST1 = as.data.frame(matrix(nrow=(nrow(PROBS)),ncol=nrow(PROBS)))


      for(k in 1:nrow(PROBS)){
        for(j in 1:nrow(PROBS)){
          CLUST1[k,j]=mean(CLUST[,k]==j)
        }
      }

      for(k in 1:nrow(PROBS)){
        rownames(CLUST1)[k]=paste0("Subgroup ", k)
        colnames(CLUST1)[k]=paste0("Latent Subgroup",k)

      }

      NCLUST=rep(0, nrow(CLUST))

      for(j in 1:nrow(CLUST)){
        NCLUST[j] = length(unique(CLUST[j,]))
      }


      ## print(RESULTS1[[2]])

      Z2=as.list(c(0,0,0))
      Z2[[1]]=CLUST1
      Z2[[2]]=table(NCLUST)/nrow(CLUST)


      ##Finally, let's get the clustering breakdown
      ##For G groups, there are (G choose 2) + (G choose 3) + .... (G choose G-1) + 2 cluster configurations

      G=nrow(DoseTried)

      if(G==2){
        HOLD = c(0,0)
        HOLD[1]=mean(CLUST[,1]==CLUST[,2])
        HOLD[2]=mean(CLUST[,1]!=CLUST[,2])

        names(HOLD)=c("1-2","1,2")
      }else{

        if(G<=4){
          NClust=2
          for(j in 2:(G-1)){
            NClust = NClust + choose(G,j)
          }



          ##Make Matrix
          HOLD = rep(NA,NClust)
          ##First do the pairs
          if(G==3){
            ##5 clusters
            ##1-2,3
            ##2-3, 1
            ##1-3, 2
            ##1,2,3
            ##1-2-3

            HOLD[1]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,3]!=CLUST[,1]))
            HOLD[2]=mean((CLUST[,2]==CLUST[,3])*(CLUST[,2]!=CLUST[,1]))
            HOLD[3]=mean((CLUST[,1]==CLUST[,3])*(CLUST[,2]!=CLUST[,1]))
            HOLD[4]=mean(NCLUST==G)
            HOLD[5]=mean(NCLUST==1)
            names(HOLD)=c("1-2,3", "2-3,1","1-3,2","1,2,3","1-2-3")

          }



          if(G==4){
            HOLD=rep(NA,15)
            names(HOLD)=c("1-2,3-4", "1-3,2-4","1-4,2-3",
                          "1-2,3,4", "1-3,2,4", "1-4,2,3",
                          "2-3,1,4","2-4,1,3","3-4,1,3",
                          "1-2-3,4", "1-2-4,3","1-3-4,2",
                          "2-3-4,1","1,2,3,4","1-2-3-4")
            ##15 clusters
            ##1-2,3-4
            ##1-3, 2-4
            ##1-4, 2-3 *
            HOLD[1]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,3]==CLUST[,4])*(CLUST[,3] != CLUST[,1]))
            HOLD[2]=mean((CLUST[,1]==CLUST[,3])*(CLUST[,2]==CLUST[,4])*(CLUST[,4] != CLUST[,1]))
            HOLD[3]=mean((CLUST[,1]==CLUST[,4])*(CLUST[,2]==CLUST[,3])*(CLUST[,3] != CLUST[,1]))
            ##1-2, 3,4
            #1-3, 2,4
            #1-4, 2,3
            #2-3,1,4
            #2-4,1,3
            #3-4, 1, 2
            HOLD[4]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,3]!=CLUST[,1])*(CLUST[,3] != CLUST[,4])*(CLUST[,4] != CLUST[,1]))
            HOLD[5]=mean((CLUST[,1]==CLUST[,3])*(CLUST[,2]!=CLUST[,4])*(CLUST[,4] != CLUST[,1])*(CLUST[,4] != CLUST[,2]))
            HOLD[6]=mean((CLUST[,1]==CLUST[,4])*(CLUST[,2]!=CLUST[,3])*(CLUST[,3] != CLUST[,1])*(CLUST[,2] != CLUST[,1]))
            HOLD[7]=mean((CLUST[,2]==CLUST[,3])*(CLUST[,1]!=CLUST[,2])*(CLUST[,1] != CLUST[,4])*(CLUST[,4] != CLUST[,2]))
            HOLD[8]=mean((CLUST[,2]==CLUST[,4])*(CLUST[,1]!=CLUST[,4])*(CLUST[,1] != CLUST[,2])*(CLUST[,3] != CLUST[,2]))
            HOLD[9]=mean((CLUST[,2]==CLUST[,3])*(CLUST[,1]!=CLUST[,2])*(CLUST[,4] != CLUST[,2])*(CLUST[,1] != CLUST[,4]))
            ##1-2-3,4
            ##1-2-4,3
            ##1-3-4,2*
            ##2-3-4,1
            ##1,2,3,4
            ##1-2-3-4

            ##3 pairs

            ##4 3 ways
            HOLD[10]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,1]==CLUST[,3])*(CLUST[,1] != CLUST[,4]))
            HOLD[11]=mean((CLUST[,1]==CLUST[,2])*(CLUST[,1]==CLUST[,4])*(CLUST[,1] != CLUST[,3]))
            HOLD[12]=mean((CLUST[,1]==CLUST[,3])*(CLUST[,1]==CLUST[,4])*(CLUST[,1] != CLUST[,2]))
            HOLD[13]=mean((CLUST[,2]==CLUST[,3])*(CLUST[,3]==CLUST[,4])*(CLUST[,1] != CLUST[,4]))
            HOLD[14]=mean(NCLUST==G)
            HOLD[15]=mean(NCLUST==1)

          }




        }else{
          HOLD = "Cannot Return Unique Clusters for more than 4 subgroups"
        }
      }






      Z2[[3]]=HOLD
      names(Z2)= c("Latent Subgroup Posterior","Posterior # of Clusters","Posterior Probability Cluster Configuration")
      Z[[4]]=Z2

    }else{
      Z[[4]]="No Clustering"
    }


    return(Z)

  }

}







