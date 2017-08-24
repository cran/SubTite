#' Simulates a Sub-TITE trial design
#'
#' Simulates replicates from a Sub-TITE trial with user specified true toxicity time distributions for different doses and subgroups and returns average summary statistics of the trial.
#' @importFrom stats pexp pweibull pgamma plnorm rexp rbinom runif rweibull rgamma rlnorm
#' @param nSims Number of Trials to Simulate.
#' @param Nmax Maximum Number of Patients to enroll in the trial.
#' @param T1 Reference time for toxicity.
#' @param Target Target cumulative toxicity probability (or subgroup specific vector) at time T1.
#' @param Dose Standardized vector of doses to try.
#' @param DoseStart Dose (or vector of Doses) to enroll the first patient in each subgroup at.
#' @param TGroup Amount of time to suspend accrual after three patients are enrolled in a subgroup.
#' @param Accrue Expected montly patient accrual rate.
#' @param groupprob Probability vector of subgroup assignment.
#' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
#' @param Hyper List of size 4 containing the prior mean of the baseline slope, the baseline intercept, and the prior mean vectors for group specific intercepts and slopes.
#' @param Family What distribution Family to simulate from. Options include: Exponential,Gamma, Lognormal, Uniform, Weibull.
#' @param Param1 #Groups X #Doses Matrix containing the first parameter for each subgroup and dose. For the uniform distribution, this is the probability of toxicity in a given group.
#' @param Param2 #Groups X #Doses Matrix containing the second parameter for each subgroup and dose for the Weibull, Gamma and Lognormal Distributions. This argument is not used for uniform and exponential distribution families.
#' @return Returns a list with three simulation outputs: The vector of optimal doses chosen, the vector containing the doses administered in a trial and the group assignments of each patient in a simulated trial.
#' @references
#' [1] Chapple and Thall (2017), Subgroup Specific Dose Finding in Phase I Clinical Trials Based on Time to Toxicity Within a Fixed Follow Up Period.
#' [2] Package Tutorial, https://adventuresinstatistics.wordpress.com/2017/08/24/the-subtite-package-tutorial/
#' @examples
#' ##Note: nSims  should be set larger than the example below.
#' nSims=1
#' ##Specify reference toxicity time and target
#' T1=6
#' Target=.3
#' ##Number of Groups
#' ##Specify upper bound for determining if the lowest dose is too toxic in a subgroup
#' Upper=c(.9,.95)
#' TGroup=0
#' ##Maximum Sample Size
#' Nmax=60
#' ##Standardized Dose Values and starting dose index
#' Dose=sort(rnorm(4))
#' DoseStart=1
#' ##Hypermeans for linear terms
#' meanmu=2.21
#' meanslope=-.57
#' MeanInts = c(.46)
#' MeanSlopes = c(.04)
#' ##Accrual Rate
#' Accrue=2
#' groupprob=c(.5,.5)
#' ##Fill in Hyperparameter list for MCMC
#' Hyper=as.list(c(0,0,0,0))
#' Hyper[[1]]=meanmu
#' Hyper[[2]]=meanslope
#' Hyper[[3]]=MeanInts
#' Hyper[[4]]=MeanSlopes
#' Family="Uniform"
#' Param1 = matrix(c(.2,.3,.4,.5,.6,.1,.2,.3,.4,.5),byrow=TRUE,nrow=2)
#' Param2=Param1
#' SimTrial(nSims,Nmax,T1,Target,Dose,DoseStart,TGroup,
#' Upper,Accrue,groupprob,Hyper,Family,Param1,Param2)
#' @export
SimTrial = function(nSims,Nmax,T1,Target,Dose, DoseStart,TGroup,Upper,Accrue,groupprob,Hyper,Family,Param1,Param2){



  nGroups=nrow(Param1)
  nGroup=nGroups
  nDose=ncol(Param1)
  upper=Upper

  Dist = matrix(rep(NA,nSims*nrow(Param1)),nrow=nSims)


  target=Target
  B=2000





  B=2000
  meanmu = Hyper[[1]]
  meanslope=Hyper[[2]]
  MeanInts=Hyper[[3]]
  MeanSlopes=Hyper[[4]]
  nDoses=length(Dose)
  target=Target


  ##If length(Target)==1 or length(DoseStart)==1, let's make a vector containing the group specific targets

  if(length(Target)==1){
    Target1=Target
    Target=rep(NA,nrow(Param1))
    for(k in 1:nrow(Param1)){
      Target[k]=Target1
    }

  }


  if(length(DoseStart)==1){
    DoseStart1=DoseStart
    DoseStart=rep(NA,nrow(Param1))
    for(k in 1:nrow(Param1)){
      DoseStart[k]=DoseStart1
    }

  }


  if(length(Upper)==1){
    Upper1=Upper
    Upper=rep(NA,nrow(Param1))
    for(k in 1:nrow(Param1)){
      Upper[k]=Upper1
    }

  }



  if(Family=="Exponential"){
    GroupProb = pexp(T1,1/Param1)
    Fam=0
  }

  if(Family=="Uniform"){

    GroupProb=Param1
    Fam=3
  }

  if(Family=="Weibull"){
    GroupProb=pweibull(T1,Param1,Param2)
    Fam=4
  }

  if(Family=="Gamma"){
    GroupProb=pgamma(T1,Param1,Param2)
    Fam=1
  }

  if(Family=="Lognormal"){
    GroupProb=plnorm(T1,Param1,Param2)
    Fam=2
  }




  POPT=rep(NA,nGroups)
  WHICHOPT=POPT
  ##What is POPT

  for(m in 1:nGroups){
    min1 = min(abs(GroupProb[m,]-Target[m]))
    which1 =  which(abs(GroupProb[m,]-Target[m])==min1)
    POPT[m]=GroupProb[m,which1]
    WHICHOPT[m]=which1

  }



  ##Feed Everything into the trial
Results= SimTrial1( nSims, Nmax,  T1, Target,  Dose,  DoseStart, TGroup,Upper, Accrue, groupprob,
                    Fam, Param1,Param2,meanmu, meanslope,MeanInts, MeanSlopes )

DoseOpt = Results[[1]]
NTox=Results[[2]]


##Calculate \Delta_1,\Delta_2
for(b in 1:nSims){

  for(m in 1:nGroups){
    if(DoseOpt[b,m]==0){
      Dist[b,m]=POPT[m]
    }else{
      Dist[b,m]=abs(GroupProb[m,DoseOpt[b,m]]-POPT[m])
    }
  }







}



  cat("Simulations finished: Displaying results

      ")


  ##Obtain the prob of best

  OptProb=WHICHOPT


  for(m in 1:nGroups){
    OptProb[m]=mean(DoseOpt[,m]==WHICHOPT[m])
  }


  cat("Subgroup Specific Selection Probability of the Optimal Dose", OptProb )


  cat("

      Frequencies of Optimal Dose Selected for each Subgroup")

  for(m in 1:nGroups){

    cat("
        Group", m,"
        ")

    print(table(DoseOpt[,m])/nSims)
  }


  cat("
      Average true toxicity probability difference between chosen dose and the optimal dose

      ")

  print(colMeans(Dist))



  cat("
      Average number of Toxicities for each Subgroup

      ")

  print(colMeans(NTox))

TrialTimes=Results[[3]]


  cat("
      Average Trial Time

      ")
  print(mean(TrialTimes))

  List1 = as.list(rep(0,6))
  List1[[1]]=DoseOpt
  List1[[2]]=Dist
  List1[[3]]=NTox
  List1[[4]]=TrialTimes


  List1[[5]]=Results[[4]]
  List1[[6]]=Results[[5]]


  return(List1)


}

