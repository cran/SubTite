#' Gives the subgroup specific optimal dose vector.
#'
#'  Returns a list containing the optimal doses to enroll each subgroup at and the subgroups that should have their accrual suspended temporarily.
#' @param Y   Vector containing observed event or censoring times.
#' @param I   Vector containing event indicators (1 if patient experiences an event for a patient).
#' @param Doses Vector containing Doses of patients in trial.
#' @param Groups Vector containing group assignment of patients, 1 is baseline group.
#' @param Hyper List of size 4 containing, in order, the prior mean of the baseline slope, the baseline intercept, and the prior mean vectors for group specific intercepts and slopes.
#' @param T1 Reference time for toxicity.
#' @param Target Target cumulative toxicity probability vector at time T1.
#' @param Dose Vector containing the standardized doses considered.
#' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
#' @param DoseTried #Groups X #Doses Matrix that contains 1s and 0s corresponding to whether a given dose in a subgroup had been tried.
#' @return Returns a list with two objects, a vector of optimal doses and a vector of stopped groups
#' @references
#' [1] Chapple and Thall (2017), Subgroup Specific Dose Finding in Phase I Clinical Trials Based on Time to Toxicity Within a Fixed Follow Up Period.
#' [2] Package Tutorial, https://adventuresinstatistics.wordpress.com/2017/08/24/the-subtite-package-tutorial/
#' @examples
#' T1=6
#'##Reference Time for Toxicity
#' Target=.3
#' Upper=c(.9,.95,.95)
#' n=30
#' Y=rep(NA,n)
#' I=rep(NA,n)
#' Groups = sample(1:3,n,replace=TRUE) - 1
#' ##Group assignment of patients (MUST BE CODED 0,1,2,...)
#' Doses = sample(1:5,n,replace=TRUE)
#' ##Randomly Generate Dose assignment
#' x=c(1,2,3,5,7)
#' Dose=(x-mean(x))/sd(x)
#' ##Vector of standardized doses
#' ##Next we generate the toxicity times based on a true toxicity probability matrix
#' GroupProb = matrix(rep(NA,length(Dose)*3),nrow=3)
#' GroupProb[1,]=c(.18,.25,.45,.66,.74)
#'  ## These are the true toxicity probabilities for each dose and subgroup
#' GroupProb[2,]=c(.10,.15,.30,.50,.60)
#' GroupProb[3,]=c(.10,.15,.30,.50,.60) +.1
#' for(b in 1:n){
#' I[b]= rbinom(1,1   , GroupProb[(Groups[b]+1),Doses[b]])
#' if(I[b]==0){ Y[b]=T1 }else{ Y[b]=runif(1,0,T1) }}
#' DoseTried=matrix(c(1,1,0,0,0,1,1,1,1,1,1,1,0,0,0),ncol=5,byrow=TRUE)
#'  ##Doses Tried so far in trial
#' ##Hypermeans for linear terms
#' meanmu=-0.4467184
#' meanslope= 0.8861634
#' MeanInts = c(-0.5205379, -1.3752529)
#' MeanSlopes = c(0.1888923, 0.1148791)
#' ##Fill in Hyperparameter list for MCMC
#' Hyper=as.list(c(0,0,0,0))
#' Hyper[[1]]=meanmu
#' Hyper[[2]]=meanslope
#' Hyper[[3]]=MeanInts
#' Hyper[[4]]=MeanSlopes
#' Doses=Dose[Doses]
#'  ##Standardized patient doses
#' GetSubTite(Y,I,Doses,Groups,Hyper,T1,Target,Dose,Upper,DoseTried)
#' @export
GetSubTite=function(Y,I,Doses,Groups,Hyper,T1,Target,Dose,Upper,DoseTried){

  nDose=length(Dose)
  nDoses=length(Dose)
  nGroups=max(Groups)
  nGroup=max(Groups)
  Dosetried=DoseTried
  upper=Upper

  B=2000
  meanmu = Hyper[[1]]
  meanslope=Hyper[[2]]
  MeanInts=Hyper[[3]]
  MeanSlopes=Hyper[[4]]
  nDoses=max(Doses)
  target=Target





  if(length(Target)==1){

  ##First Do MCMC using ALL data to pick an optimal dose
  G1=GroupMCMCLN( Y,I,Doses,Groups,meanmu,meanslope,MeanInts,MeanSlopes, B,T1)


  B=nrow(G1[[1]])
  B1=B/2
  IntG = as.matrix(G1[[1]][(B-B1):B,])
  SlopeG = as.matrix(G1[[2]][(B-B1):B,])

  Mu = G1[[3]][(B-B1):B]
  Slope = G1[[4]][(B-B1):B]
  sigma = sqrt(G1[[5]][(B-B1):B])

  B1=length((B-B1):B)




  etaG = matrix(rep(NA,(ncol(G1[[1]])+1)*B1),nrow=B1)

  DoseProb = matrix(rep(NA,(ncol(G1[[1]])+1)*B1),ncol=(ncol(G1[[1]])+1))


  MeanDose=matrix(rep(NA,length(Dose)*(ncol(G1[[1]])+1)),ncol=length(Dose))

  for(m in 1:length(Dose)){


    for(b in 1:B1){
      etaG[b,1]=Mu[b]+Slope[b]*Dose[m]

      for(k in 1:ncol(G1[[1]])){
        etaG[b,k+1]=Mu[b]+IntG[b,k]+(Slope[b]+SlopeG[b,k])*Dose[m]
      }




etaG=exp(etaG)




      for(k in 1:(ncol(G1[[1]])+1)){

        DoseProb[b,k] = etaG[b,k]/(1+etaG[b,k])

        if(is.infinite(etaG[b,k])){
          DoseProb[b,k]=1
        }


      }


    }


    MeanDose[,m]=colMeans(DoseProb)



  }


  OptimalDose = rep(NA,ncol(G1[[1]]+1))


  for(k in 1:(ncol(G1[[1]])+1)){
    if(length(which(abs(MeanDose[k,]-target)==min(abs(MeanDose[k,]-target))))>1){
      ##Check if all doses are bigger or less
      if(MeanDose[k,which(abs(MeanDose[k,]-target)==min(abs(MeanDose[k,]-target)))[1]]>target){
        OptimalDose[k]=which(abs(MeanDose[k,]-target)==min(abs(MeanDose[k,]-target)))[1]

      }else{
        OptimalDose[k]=which(abs(MeanDose[k,]-target)==min(abs(MeanDose[k,]-target)))[length(which(abs(MeanDose[k,]-target)==min(abs(MeanDose[k,]-target))))]
      }
    }else{
      OptimalDose[k]=  which(abs(MeanDose[k,]-target)==min(abs(MeanDose[k,]-target)))
    }
  }

  ##  print(MeanDose)

  TriedDoses = as.list(rep(NA,nGroups))

  OptDose=OptimalDose

  for(k in 1:nGroups){
    dose1 = c(0)

    for(m in 1:nDose){
      if(Dosetried[k,m]==1)
        dose1 = c(dose1,m)
    }

    TriedDoses[[k]]=dose1

  }

  for(m in 1:nGroups){
  if( sum(OptDose[m] == TriedDoses[[m]])==0){
    if(OptDose[m]>TriedDoses[[m]][length(TriedDoses[[m]])]){
      OptDose[m] = TriedDoses[[m]][length(TriedDoses[[m]])] + 1
    }else{
      OptDose[m] = TriedDoses[[m]][2]-1
    }

  }
  }


  RETURN1 = as.list(c(0,0))


  RETURN1[[1]]=OptDose





  B=nrow(G1[[1]])
  B1=B/2

  IntG = as.matrix(G1[[1]][(B-B1):B,])
  SlopeG = as.matrix(G1[[2]][(B-B1):B,])

  Mu = G1[[3]][(B-B1):B]
  Slope = G1[[4]][(B-B1):B]
  sigma = sqrt(G1[[5]])[(B-B1):B]

  B1=length((B-B1):B)






  etaG = matrix(rep(NA,(ncol(G1[[1]])+1)*B1),nrow=B1)

  DoseProb = matrix(rep(NA,(ncol(G1[[1]])+1)*B1),ncol=(ncol(G1[[1]])+1))


  MeanDose=matrix(rep(NA,length(Dose)*(ncol(G1[[1]])+1)),ncol=length(Dose))

  m=1

  for(b in 1:B1){
    etaG[b,1]=Mu[b]+Slope[b]*Dose[m]

    for(k in 1:ncol(G1[[1]])){
      etaG[b,k+1]=Mu[b]+IntG[b,k]+(Slope[b]+SlopeG[b,k])*Dose[m]
    }





etaG=exp(etaG)



    for(k in 1:(ncol(G1[[1]])+1)){
      ##Calculate probability that
      DoseProb[b,k] = etaG[b,k]/(1+etaG[b,k])

      if(is.infinite(etaG[b,k])){
        DoseProb[b,k]=1
      }

    }


  }



  ##Contains Mean Survival Probabilities
  MeanDose[,m]=colMeans(DoseProb>target)



  # print(MeanDose)

  StopGroup=rep(NA,nrow(MeanDose))

  for(k in 1:nrow(MeanDose)){
    if(MeanDose[k,m]>upper[k]){
      StopGroup[k]=1
    }else{
      StopGroup[k]=0
    }
  }

  StopGroup = StopGroup


  }else{

    ##Subgroup Specific Toxicity Targets

    ##First Do MCMC using ALL data to pick an optimal dose
    G1=GroupMCMCLN( Y,I,Doses,Groups,meanmu,meanslope,MeanInts,MeanSlopes, B,T1)


    B=nrow(G1[[1]])
    B1=B/2
    IntG = as.matrix(G1[[1]][(B-B1):B,])
    SlopeG = as.matrix(G1[[2]][(B-B1):B,])

    Mu = G1[[3]][(B-B1):B]
    Slope = G1[[4]][(B-B1):B]

    B1=length((B-B1):B)




    etaG = matrix(rep(NA,(ncol(G1[[1]])+1)*B1),nrow=B1)

    DoseProb = matrix(rep(NA,(ncol(G1[[1]])+1)*B1),ncol=(ncol(G1[[1]])+1))


    MeanDose=matrix(rep(NA,length(Dose)*(ncol(G1[[1]])+1)),ncol=length(Dose))

    for(m in 1:length(Dose)){


      for(b in 1:B1){
        etaG[b,1]=Mu[b]+Slope[b]*Dose[m]

        for(k in 1:ncol(G1[[1]])){
          etaG[b,k+1]=Mu[b]+IntG[b,k]+(Slope[b]+SlopeG[b,k])*Dose[m]
        }





etaG=exp(etaG)



        for(k in 1:(ncol(G1[[1]])+1)){

          DoseProb[b,k] = etaG[b,k]/(1+etaG[b,k])

          if(is.infinite(DoseProb[b,k])){
            DoseProb[b,k]=1
          }


        }


      }


      MeanDose[,m]=colMeans(DoseProb)



    }


    OptimalDose = rep(NA,ncol(G1[[1]]+1))


    for(k in 1:(ncol(G1[[1]])+1)){
      if(length(which(abs(MeanDose[k,]-target[k])==min(abs(MeanDose[k,]-target[k]))))>1){
        ##Check if all doses are bigger or less
        if(MeanDose[k,which(abs(MeanDose[k,]-target[k])==min(abs(MeanDose[k,]-target[k])))[1]]>target[k]){
          OptimalDose[k]=which(abs(MeanDose[k,]-target[k])==min(abs(MeanDose[k,]-target[k])))[1]

        }else{
          OptimalDose[k]=which(abs(MeanDose[k,]-target[k])==min(abs(MeanDose[k,]-target[k])))[length(which(abs(MeanDose[k,]-target[k])==min(abs(MeanDose[k,]-target[k]))))]
        }
      }else{
        OptimalDose[k]=  which(abs(MeanDose[k,]-target[k])==min(abs(MeanDose[k,]-target[k])))
      }
    }

    ##  print(MeanDose)

    TriedDoses = as.list(rep(NA,nGroups))

    OptDose=OptimalDose

    for(k in 1:nGroups){
      dose1 = c(0)

      for(m in 1:nDose){
        if(Dosetried[k,m]==1)
          dose1 = c(dose1,m)
      }

      TriedDoses[[k]]=dose1

    }

    for(m in 1:nGroups){
      if( sum(OptDose[m] == TriedDoses[[m]])==0){
        if(OptDose[m]>TriedDoses[[m]][length(TriedDoses[[m]])]){
          OptDose[m] = TriedDoses[[m]][length(TriedDoses[[m]])] + 1
        }else{
          OptDose[m] = TriedDoses[[m]][2]-1
        }

      }
    }


    RETURN1 = as.list(c(0,0))


    RETURN1[[1]]=OptDose









    B=nrow(G1[[1]])
    B1=B/2

    IntG = as.matrix(G1[[1]][(B-B1):B,])
    SlopeG = as.matrix(G1[[2]][(B-B1):B,])

    Mu = G1[[3]][(B-B1):B]
    Slope = G1[[4]][(B-B1):B]

    B1=length((B-B1):B)






    etaG = matrix(rep(NA,(ncol(G1[[1]])+1)*B1),nrow=B1)

    DoseProb = matrix(rep(NA,(ncol(G1[[1]])+1)*B1),ncol=(ncol(G1[[1]])+1))


    MeanDose=matrix(rep(NA,length(Dose)*(ncol(G1[[1]])+1)),ncol=length(Dose))

    m=1

    for(b in 1:B1){
      etaG[b,1]=Mu[b]+Slope[b]*Dose[m]

      for(k in 1:ncol(G1[[1]])){
        etaG[b,k+1]=Mu[b]+IntG[b,k]+(Slope[b]+SlopeG[b,k])*Dose[m]
      }









      for(k in 1:(ncol(G1[[1]])+1)){

        DoseProb[b,k] = etaG[b,k]/(1+etaG[b,k])

        if(is.infinite(DoseProb[b,k])){
          DoseProb[b,k]=1
        }

      }


    }



    ##Contains Mean Survival Probabilities
    MeanDose[,m]=colMeans(DoseProb>target)

    for(k in 1:nrow(MeanDose)){
      MeanDose[k,1]=mean(DoseProb[,k]>target[k])
    }


    # print(MeanDose)

    StopGroup=rep(NA,nrow(MeanDose))

    for(k in 1:nrow(MeanDose)){
      if(MeanDose[k,m]>upper[k]){
        StopGroup[k]=1
      }else{
        StopGroup[k]=0
      }
    }





  }

  for(k in 1:nrow(MeanDose)){
    if(StopGroup[k]==1){
      OptDose[k]=0
    }
  }

  RETURN1[[1]]=OptDose

  RETURN1[[2]]=StopGroup

  return(RETURN1)


}








