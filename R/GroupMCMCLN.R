#' Performs MCMC for model used in the SubTite design
#'
#' Returns a list containing the posterior samples for the logistic model parameters.
#' @importFrom Rcpp evalCpp
#' @param Y   Vector containing observed event or censoring times.
#' @param I   Vector containing event indicators (1 if patient experiences an event for a patient).
#' @param Doses Vector containing Doses of patients in trial.
#' @param Groups Vector containing group assignment of patients, 1 is baseline group.
#' @param meanmu Prior mean of common intercept parameter.
#' @param meanslope Prior mean of common slope parameter.
#' @param MeanInts Vector of prior means for the group specific intercept parameters.
#' @param MeanSlopes Vector of prior means for the group specific slope parameters.
#' @param B Number of iterations to perform in the MCMC.
#' @param T1 Reference time for Toxicity
#' @return Returns a list containing the posterior samples for the common intercept, the common slope, the vector of group specific intercepts, the vector of group specific slopes and sigma^2.
#' @references
#' [1] Chapple and Thall (2017), Subgroup-specific Dose Finding in Phase I Clinical Trials Based on Time to Toxicity.
#' @examples
#' ##Specify reference toxicity time and target
#' T1=6
#' B=2000
#' ##Randomly generate data for use.
#' n=60
#' Y=runif(n,0,100)
#' I=rbinom(n,1,.5)
#' Y[I==0]=T1
#' Groups = sample(1:3,n,replace=TRUE)
#' Doses = sample(1:5,n,replace=TRUE)
#' ##Doses Tried so far in trial
#' ##Hypermeans for linear terms
#' meanmu=2.21
#' meanslope=-.57
#' MeanInts = c(.46, .55)
#' MeanSlopes = c(.04,.10)
#' GroupMCMCLN( Y,I,Doses,Groups,meanmu,meanslope,MeanInts,MeanSlopes, B,T1)
#' @export
GroupMCMCLN <- function(Y, I, Doses, Groups, meanmu, meanslope, MeanInts, MeanSlopes, B,T1) {
  Dose=Doses
  Group=Groups
  .Call('SubTite_GroupMCMCLN', PACKAGE = 'SubTite', Y, I, Dose, Group, meanmu, meanslope, MeanInts, MeanSlopes, B,T1)
}



