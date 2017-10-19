#' Gives the subgroup specific optimal dose vector.
#'  Returns a list containing the optimal doses to enroll each subgroup at and the subgroups that should have their accrual suspended temporarily.
#' @param Y   Vector containing observed event or censoring times.
#' @param I   Vector containing event indicators (1 if patient experiences an event for a patient).
#' @param Doses Vector containing Doses of patients in trial.
#' @param Groups Vector containing group assignment of patients, 1 is baseline group.
#' @param T1 Reference time for toxicity.
#' @param Target Target cumulative toxicity probability vector at time T1.
#' @param Dose Vector containing the standardized doses considered.
#' @param Upper Cutoff values used to determine if accrual in a subgroup should be suspended.
#' @param DoseTried #Groups X #Doses Matrix that contains 1s and 0s corresponding to whether a given dose in a subgroup had been tried.
#' @param meanmu Prior mean for baseline intercept.
#' @param meanslope Prior mean for baseline slope.
#' @param MeanInts Vector of prior means for the group specific intercept parameters.
#' @param MeanSlopes Vector of prior means for the group specific slope parameters.
#' @param varint Prior variance for the intercept parameters.
#' @param varbeta Prior variance for the slope parameters.#' @return Returns a list with two objects, a vector of optimal doses and a vector of stopped groups
#' @references
#' [1] Chapple and Thall (2017), Subgroup Specific Dose Finding in Phase I Clinical Trials Based on Time to Toxicity Within a Fixed Follow Up Period.
#' [2] Package Tutorial, https://adventuresinstatistics.wordpress.com/2017/08/24/the-subtite-package-tutorial/
#' @examples
#' T1=6
#'##Reference Time for Toxicity
#' Target=.3
#' Upper=c(.95,.95)
#' n=30
#' Y=rep(NA,n)
#' I=rep(NA,n)
#' Groups = sample(1:2,n,replace=TRUE) - 1
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
#' GroupProb[2,]=c(.10,.15,.30,.50,.60) +.1
#' for(b in 1:n){
#' I[b]= rbinom(1,1   , GroupProb[(Groups[b]+1),Doses[b]])
#' if(I[b]==0){ Y[b]=T1 }else{ Y[b]=runif(1,0,T1) }}
#' DoseTried=matrix(c(1,1,0,0,0,1,1,1,1,1),ncol=5,byrow=TRUE)
#'  ##Doses Tried so far in trial
#' ##Hypermeans for linear terms
#' meanmu=-0.4467184
#' meanslope= 0.8861634
#' MeanInts = -0.5205379
#' MeanSlopes = 0.1888923
#' Doses=Dose[Doses]
#' varint=5
#' varbeta=1
#' GetSubTite(Y, I,Doses, Groups, DoseTried, T1,
#'  Target,  Upper, Dose,  meanmu, meanslope,
#'   MeanInts,  MeanSlopes ,varint,varbeta)
#' @export
GetSubTite=function(Y, I,Doses, Groups, DoseTried, T1, Target,
                    Upper, Dose,  meanmu, meanslope,
                    MeanInts,  MeanSlopes ,varint,varbeta){




  OptDose=GetDose1( Y, I,Doses, Groups, DoseTried, T1, Target,  Upper, Dose,
                    meanmu, meanslope,  MeanInts,  MeanSlopes ,varint,varbeta)


  if(sum(OptDose)==0){
    cat("Stop Trial if at least one cohort of both subgroup have been fully evaluated.")
  }else{
    cat(paste("Enroll the next patient in subgroup 0 at Dose",OptDose[1],"and the next patient in subgroup 1 at Dose",OptDose[2]))
  }


  return(OptDose)


}








