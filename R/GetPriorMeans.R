#' Calibrates prior means for Dose Finding Trial
#'
#'  Uses the clinician elicited prior reference probabilities for each subgroup and dose to obtain prior means for the Bayesian logistic regression model used in the SubTite trial design.
#' @param Clinician #Groups X #Doses matrix containing the elicited prior toxicity probabilities at the reference time for each dose and subgroup.
#' @param Dose Vector containing standardized doses.
#' @return Returns a list containing in order: the prior mean of the common baseline intercept, the prior mean of the common slope, the vector of subgroup specific intercepts and the vector of subgroup specific slopes.
#' @references
#' [1] Chapple and Thall (2017), Subgroup-specific Dose Finding in Phase I Clinical Trials Based on Time to Toxicity.
#' @examples
#' ##Specify elicited reference toxicity probabilities
#' Clinician = matrix(c(.2,.3,.4,.5,.6,.1,.2,.3,.4,.5,.05,.1,.15,.2,.3),byrow=TRUE,nrow=3)
#' Dose=rnorm(5)
#' GetPriorMeans(Clinician,Dose)
#' @export
GetPriorMeans = function(Clinician,Dose){

  X=Clinician
  nGroups=nrow(X)

  Y = rep(NA,length(X))

  for(m in 1:nrow(X)){

    for(k in 1:ncol(X)){

      Y[(m-1)*ncol(X)+k] = log(X[m,k]/(1-X[m,k]))

    }


  }

  ##Now we have our Y lets make X
  ##It's going to be structured intercept, group ints, slope, group slopes
  ##Number of Groups
  G=nrow(X)
  G2=G-1
  D=length(Dose)

  COV = matrix(rep(0,length(Y)*2*G),nrow=length(Y))
  COV[,1]=1


  DOSEVEC=Y

  for(m in 1:G){
    DOSEVEC[((m-1)*length(Dose)+1):(m*length(Dose))]=Dose
  }


  ##Fill rest of Cov with Doses and we will zero out the group scale ones
  for(m in (G+1):ncol(COV)){
    COV[,m]=0
  }

  for(m in (G+2):ncol(COV)){
    k=m-G
    COV[((k-1)*D+1):(k*D),m] = Dose

    COV[((k-1)*D+1):(k*D),k] = 1


  }

  COV[,G+1]=DOSEVEC



  beta=solve(t(COV)%*%COV)%*%t(COV)%*%Y


  MeanInts = beta[2:(nGroups)]
  meanslope=beta[nGroups+1]
  MeanSlopes=beta[(nGroups+2):length(beta)]

  X1=as.list(c(0,0,0,0))

  X1[[1]]=beta[1]
  X1[[2]]=meanslope
  X1[[3]]=MeanInts
  X1[[4]]=MeanSlopes

  return(X1)




}
