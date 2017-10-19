#include "RcppArmadillo.h"
#include <time.h>
#include <Rmath.h>
#include <math.h>
#define ARMA_DONT_PRINT_ERRORS
#include <numeric>
#include <algorithm>
#include <vector>
#include <iterator>
#include <list>
#include <iostream>     // std::cout
#include <cmath>
#include <cfloat>
#include <stdio.h>
#include <float.h>
#define PI 3.14159265


// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;


int Sample2(arma::vec groupprob){
  arma::vec cumprob=groupprob;
  int m;

  for(m=1;m<groupprob.n_rows;m++){
    cumprob[m]=cumprob[m]+cumprob[m-1];
  }
  //Now we have the vector of cumulative probabilities, let's draw a random unif
  double U=as_scalar(arma::randu(1));

  int Which;


  if(U<cumprob[0]){
    Which=0;
  }else{

    for(m=0;m<(groupprob.n_rows-2);m++){
      if(U>cumprob[m] & U<cumprob[m+1]){
        Which=m+1;
      }

    }

    if(U>cumprob[groupprob.n_rows-2]){
      Which=groupprob.n_rows -1;
    }


  }


  return(Which);

}



int GetFullyFollowed(arma::vec Y, //Survival Times
                     arma::vec I, //Censoring Indicators
                     arma::vec Doses, //Doses given to patients
                     arma::vec Dose, //Standardized doses in consideration in trial
                     double T1){
  int sum1=0;

  int m=0;

  for(m=0;m<Y.n_rows;m++){
    if(Doses[m]==Dose[0]){
      if(I[m]==1){
        sum1++;
      }

      if(Y[m]==T1){
        sum1++;
      }




    }


  }





  return(sum1);


}



double MaxVec(arma::vec Y){
  int J1=Y.n_rows;
  int j;
  double max=Y[0];
  for(j=1;j<J1;j++){
    if(Y[j]>max){
      max=Y[j];
    }
  }

  return(max);

}

arma::vec ReturnStoppedY(arma::vec Y, arma::vec I,arma::vec Groups, arma::vec StoppedGroups, int i){


//How many entries are NOT Stopped
int m=0;
int count=0;
arma::vec IN(i);
IN.zeros();
int k=0;
for(m=0;m<i;m++){
  for(k=0;k<StoppedGroups.n_rows;k++){
  if(Groups[m]==k){
    if(StoppedGroups[k]==0){
      count++;
      IN[m]=1;
    }
  }

  }
}

arma::vec Y1(count);
arma::vec I1(count);
arma::vec Groups1(count);
k=0;

for(m=0;m<count;m++){
  while(IN[k]==0){
  k++;
  }

  Y1[m]=Y[k];
  I1[m]=I[k];
  Groups1[m]=Groups[k];
  k=k+1;
}








  return(Y1);


}

arma::vec ReturnStoppedI(arma::vec Y, arma::vec I,arma::vec Groups, arma::vec StoppedGroups, int i){


  //How many entries are NOT Stopped
  int m=0;
  int count=0;
  arma::vec IN(i);
  IN.zeros();
  int k=0;
  for(m=0;m<i;m++){
    for(k=0;k<StoppedGroups.n_rows;k++){
      if(Groups[m]==k){
        if(StoppedGroups[k]==0){
          count++;
          IN[m]=1;
        }
      }

    }
  }

  arma::vec Y1(count);
  arma::vec I1(count);
  arma::vec Groups1(count);
  k=0;

  for(m=0;m<count;m++){
    while(IN[k]==0){
      k++;
    }

    Y1[m]=Y[k];
    I1[m]=I[k];
    Groups1[m]=Groups[k];
    k=k+1;
  }








  return(I1);


}

arma::vec ReturnStoppedGroups(arma::vec Y, arma::vec I,arma::vec Groups, arma::vec StoppedGroups, int i){


  //How many entries are NOT Stopped
  int m=0;
  int count=0;
  arma::vec IN(i);
  IN.zeros();
  int k=0;
  for(m=0;m<i;m++){
    for(k=0;k<StoppedGroups.n_rows;k++){
      if(Groups[m]==k){
        if(StoppedGroups[k]==0){
          count++;
          IN[m]=1;
        }
      }

    }
  }

  arma::vec Y1(count);
  arma::vec I1(count);
  arma::vec Groups1(count);
  k=0;

  for(m=0;m<count;m++){
    while(IN[k]==0){
      k++;
    }

    Y1[m]=Y[k];
    I1[m]=I[k];
    Groups1[m]=Groups[k];
    k=k+1;
  }








  return(Groups1);


}








double MinVec(arma::vec Y){
  int J1=Y.n_rows;
  int j;
  double max=Y[0];
  for(j=1;j<J1;j++){
    if(Y[j]<max){
      max=Y[j];
    }
  }

  return(max);

}


double min1(double a, double b){
  double z=0;
  if(a>=b){
    z=b;
  }else{
    z=a;
  }

  return(z);
}



double max1(double a, double b){
  double z=0;
  if(a>=b){
    z=a;
  }else{
    z=b;
  }

  return(z);
}





double Like1(arma::vec Y,  //Vector of Toxicity Times
             arma::vec I,  // Vector of Censoring Indi
             arma::vec Dose, //Vector of Doses given to Patients
             arma::vec Group, //Vector of Group Membership, 0 is baseline
             double mu, // Reference Time
             double slope, //Baseline Slope for Dose Toxicity Relationship
             arma::vec a, //Vector of intercepts for non-baseline groups
             arma::vec b, //vector of slopes for non-baseline groups
             double T1, // Reference Time
             int npats //Number of patients currently
){

  // npats=npats-1;
  //This may need to be commented out. Not totally sure here.

  double LogL = 0; // Will store likelihood to return


  //First we calculate the eta term in the GLM for each patient.

  //Baseline vector
  arma::vec eta(npats);

  eta.zeros();

  //Needed for For loops
  int m;
  int k=1;




  for(m=0;m<eta.n_rows;m++){
    for(k=0;k<a.n_rows;k++){
      if(Group[m]==(k+1)){
        eta[m] = a[k]+exp(b[k]+slope)*Dose[m]+mu;
      }
    }

    if(Group[m]==0){
      eta[m]=mu+exp(slope)*Dose[m];

    }
    //For all groups add the intercept and slope*Groups


  }



  //Now we have the linear predictors in the vector \eta for all the patients
  //Now let's compute the likelihood.
  arma::vec Y1=Y/T1; //Fraction of followup time needed to get the likelihood


  for(m=0;m<eta.n_rows;m++){
    if(I[m]==1){
      //Not censored so we use the PDF
      LogL = LogL + eta[m] - log(1+exp(eta[m]));

    }else{
      LogL = LogL + log(1+(1-Y1[m])*exp(eta[m]))-log(1+exp(eta[m]));

      //      LogL = LogL + log(1-Y1[m]*exp(eta[m])/(1+exp(eta[m])));


    }

  }








  return(LogL);


}


double LikeMULTI(arma::vec Y,  //Vector of Toxicity Times
             arma::vec I,  // Vector of Censoring Indi
             arma::vec Dose, //Vector of Doses given to Patients
             arma::vec Group, //Vector of Group Membership, 0 is baseline
             double mu, // Reference Time
             double slope, //Baseline Slope for Dose Toxicity Relationship
             arma::vec a, //Vector of intercepts for non-baseline groups
             arma::vec b, //vector of slopes for non-baseline groups
             double T1, // Reference Time
             int npats, //Number of patients currently
             arma::vec GroupMem1
){

  // npats=npats-1;
  //This may need to be commented out. Not totally sure here.



  arma::vec GroupMem(GroupMem1.n_rows+1);


  double LogL = 0; // Will store likelihood to return


  //First we calculate the eta term in the GLM for each patient.

  //Baseline vector
  arma::vec eta(npats);

  eta.zeros();

  //Needed for For loops
  int m;
  int k=1;

  for(m=1;m<GroupMem.n_rows;m++){

    GroupMem[m]=GroupMem1[m-1];

  }

  GroupMem[0]=0;
  for(m=0;m<eta.n_rows;m++){
    for(k=0;k<a.n_rows;k++){
      if(GroupMem(Group[m])==(k+1)){
        eta[m] = a[k]+exp(b[k]+slope)*Dose[m]+mu;
      }
    }

    if(GroupMem(Group[m])==0){
      eta[m]=mu+exp(slope)*Dose[m];

    }
    //For all groups add the intercept and slope*Groups


  }



  //Now we have the linear predictors in the vector \eta for all the patients
  //Now let's compute the likelihood.
  arma::vec Y1=Y/T1; //Fraction of followup time needed to get the likelihood


  for(m=0;m<eta.n_rows;m++){
    if(I[m]==1){
      //Not censored so we use the PDF
      LogL = LogL + eta[m] - log(1+exp(eta[m]));

    }else{
      LogL = LogL + log(1+(1-Y1[m])*exp(eta[m]))-log(1+exp(eta[m]));

      //      LogL = LogL + log(1-Y1[m]*exp(eta[m])/(1+exp(eta[m])));


    }

  }








  return(LogL);


}


double Like3(arma::vec Y,  //Vector of Toxicity Times
             arma::vec I,  // Vector of Censoring Indi
             arma::vec Dose, //Vector of Doses given to Patients
             arma::vec Group, //Vector of Group Membership, 0 is baseline
             double mu, // Reference Time
             double slope, //Baseline Slope for Dose Toxicity Relationship
             arma::vec a, //Vector of intercepts for non-baseline groups
             arma::vec b, //vector of slopes for non-baseline groups
             double T1, // Reference Time
             int npats, //Number of patients currently
               int which1 //This is the Group that's dropped
){

  // npats=npats-1;
  //This may need to be commented out. Not totally sure here.

  double LogL = 0; // Will store likelihood to return


  //First we calculate the eta term in the GLM for each patient.

  //Baseline vector
  arma::vec eta(npats);

  eta.zeros();

  //Needed for For loops
  int m;
  int k=1;



if(which1<0){
 //The Baseline Group was NOT REMOVED

  for(m=0;m<eta.n_rows;m++){
    for(k=0;k<a.n_rows;k++){
      if(Group[m]==(k+1)){
        eta[m] = a[k]+exp(b[k]+slope)*Dose[m]+mu;
      }
    }

    if(Group[m]==0){
      eta[m]=mu+exp(slope)*Dose[m];

    }
    //For all groups add the intercept and slope*Groups


  }

}else{
  //The Baseline Group (and maybe others) WAS REMOVED

  for(m=0;m<eta.n_rows;m++){
    for(k=0;k<a.n_rows;k++){
      if(Group[m]==(k+1)){
        eta[m] = a[k]+exp(b[k]+slope-b[which1-1])*Dose[m]+mu-a[which1-1];
      }
    }

    if(Group[m]==which1){
      eta[m]=mu+exp(slope)*Dose[m];

    }
    //For all groups add the intercept and slope*Groups


  }



}

  //Now we have the linear predictors in the vector \eta for all the patients
  //Now let's compute the likelihood.
  arma::vec Y1=Y/T1; //Fraction of followup time needed to get the likelihood


  for(m=0;m<eta.n_rows;m++){
    if(I[m]==1){
      //Not censored so we use the PDF
      LogL = LogL + eta[m] - log(1+exp(eta[m]));

    }else{
      LogL = LogL + log(1+(1-Y1[m])*exp(eta[m]))-log(1+exp(eta[m]));

      //      LogL = LogL + log(1-Y1[m]*exp(eta[m])/(1+exp(eta[m])));


    }

  }








  return(LogL);


}


//This is if the trial drops to just 1 subgroup
double Like2(arma::vec Y,  //Vector of Toxicity Times
             arma::vec I,  // Vector of Censoring Indi
             arma::vec Dose, //Vector of Doses given to Patients
             double mu, // Reference Time
             double slope, //Baseline Slope for Dose Toxicity Relationship
             double T1, // Reference Time
             int npats //Number of patients currently
){

  // npats=npats-1;
  //This may need to be commented out. Not totally sure here.

  double LogL = 0; // Will store likelihood to return


  //First we calculate the eta term in the GLM for each patient.

  //Baseline vector
  arma::vec eta(npats);

  eta.zeros();

  //Needed for For loops
  int m;
  int k=1;




  for(m=0;m<eta.n_rows;m++){

      eta[m]=mu+exp(slope)*Dose[m];


    //For all groups add the intercept and slope*Groups


  }



  //Now we have the linear predictors in the vector \eta for all the patients
  //Now let's compute the likelihood.
  arma::vec Y1=Y/T1; //Fraction of followup time needed to get the likelihood


  for(m=0;m<eta.n_rows;m++){
    if(I[m]==1){
      //Not censored so we use the PDF
      LogL = LogL + eta[m] - log(1+exp(eta[m]));

    }else{
      LogL = LogL + log(1+(1-Y1[m])*exp(eta[m]))-log(1+exp(eta[m]));

      //      LogL = LogL + log(1-Y1[m]*exp(eta[m])/(1+exp(eta[m])));


    }

  }








  return(LogL);


}





























int Sample1(int J1){

  arma::vec groupprob(J1);

  double J=J1;



  int m;

  for(m=0;m<J;m++){
    groupprob[m]=1/J;
  }



  arma::vec cumprob=groupprob;




  for(m=1;m<groupprob.n_rows;m++){
    cumprob[m]=1/J+cumprob[m-1];
  }


  //Now we have the vector of cumulative probabilities, let's draw a random unif
  double U=as_scalar(arma::randu(1));

  int Which;


  if(U<cumprob[0]){
    Which=0;
  }else{

    for(m=0;m<(groupprob.n_rows-2);m++){
      if(U>cumprob[m] & U<cumprob[m+1]){
        Which=m+1;
      }

    }

    if(U>cumprob[groupprob.n_rows-2]){
      Which=groupprob.n_rows -1;
    }


  }


  return(Which);

}





int GetIn(arma::vec INVEC){
  int k=0;

  arma::vec RetVal(sum(INVEC));
  int m=0;
  if(INVEC[0]==1){
    RetVal[0]=1;
    //Different Here

    m=1;
    for(k=1;k<RetVal.n_rows;k++){
      while(INVEC[m]==0){
        m++;
      }
      RetVal[k]=m;
      m=m+1;


    }



  }else{

    m=0;
    for(k=0;k<RetVal.n_rows;k++){
      while(INVEC[m]==0){
        m++;
      }
      RetVal[k]=m;
      m=m+1;


    }

  }




  //Now Randomly Select One

  return(RetVal(Sample1(RetVal.n_rows)));



}






int GetOut(arma::vec INVEC){
  int k=0;

  arma::vec RetVal(INVEC.n_rows-sum(INVEC));
  int m=0;
  if(INVEC[0]==0){
    RetVal[0]=1;
    //Different Here

    m=1;
    for(k=1;k<RetVal.n_rows;k++){
      while(INVEC[m]==1){
        m++;
      }
      RetVal[k]=m;
      m=m+1;


    }



  }else{

    m=0;
    for(k=0;k<RetVal.n_rows;k++){
      while(INVEC[m]==1){
        m++;
      }
      RetVal[k]=m;
      m=m+1;


    }

  }




  //Now Randomly Select One

  return(RetVal(Sample1(RetVal.n_rows)));



}





int GetNewGroup(arma::vec INVEC){

  double J1 = sum(INVEC);
  int J=sum(INVEC);
  double prob = 1/(J1+1);

  int m=0;
  int k=0;
  arma::vec cumprob(J+1);

  cumprob.zeros();
  cumprob[0]=prob;


  //Fill in w probs
  for(m=1;m<cumprob.n_rows;m++){
    cumprob[m]=prob+cumprob[m-1];
  }


  //Now we have the vector of cumulative probabilities, let's draw a random unif
  double U=as_scalar(arma::randu(1));

  int Which;


  if(U<cumprob[0]){
    Which=0;
  }else{

    for(m=0;m<(cumprob.n_rows-2);m++){
      if(U>cumprob[m] & U<cumprob[m+1]){
        Which=m+1;
      }

    }

    if(U>cumprob[cumprob.n_rows-2]){
      Which=cumprob.n_rows -1;
    }


  }


  //Now Which Contains the entry that we are going to combine into:

  if(Which==0){
    return(0);
  }else{
    arma::vec entries(J);

    k=0;
    for(m=0;m<J;m++){

      while(INVEC[k]==0){
        k++;
      }
      entries(m)=k;
      k=k+1;



    }

    entries=entries+1;

    return(entries(Which-1));

  }




}












//[[Rcpp::export]]
List SimTrial1(int nSims, //Number of Simulations to Run
               int Nmax, //Max patients to enroll in trial
               double  T1, //Reference Time
               arma::vec Target, //Reference Probability Target (or vector)
               arma::vec Dose, //Standardized Vector of Doses
               arma::vec DoseStart, //Vector (or just one) Starting Dose Level
               arma::vec Upper, //Thresholds for subgroup suspension
               double Accrue, //Expected Monthly accrual rate
               arma::vec groupprob, //Subgroup Enrollment Probabilities
               int Family, //Indicator of Distribution family to simulate from. Alphabetical order
               arma::mat Param1, //Group Specific Parameter matrix for each dose and subgroup to simulate from
               arma::mat Param2,
               double meanmu, //Prior mean of intercept
               double meanslope, //Prior mean of slope
               arma::vec MeanInts, //Prior means of group intercepts
               arma::vec MeanSlopes,
               double varint,
               double varbeta){

  int g=0;

  //Change this to be an input
  arma::vec PROBIN(Param1.n_rows-1);
  for(g=0;g<PROBIN.n_rows;g++){
    PROBIN(g)=.9;
  }

  int MEANS=0;

  int IntIN=0;
  int IntOUT=0;

  //Group Slope Prior Var
  double varbeta1=varbeta;
  //Group Int Prior Var
  double varint1=varint;



  int Group=0;
  //Important For loop integer
  int m =0; //For the inner MCMC looprep
  int i=0; //For the patient index
  int rep=0; //For the simulation repetition.


  //First Fill out Target Vector
  int nDose=Dose.n_rows;


  //Important Integer Quantities
  int B=2000; //Number of iterations for MCMC


  int B1=B/2;

  //Make List objects we'll use

  //Trial Vectors
  arma::vec Y(Nmax);
  arma::vec I(Nmax);
  arma::vec Groups(Nmax);
  arma::vec Doses(Nmax);
  arma::vec Times(Nmax);
  arma::vec ACC(Nmax);


  //Innitialize parameters for MCMC
  double mu=meanmu;
  double slope=meanslope;
  double sig=T1;

  double NewMean;
  double NewSlope;
  int Which1;



  //Important quantities we will need in the MCMC
  int  J=MeanInts.n_rows;
  double alpha=0;
  double U=0;
  double signew=0;
  double slopenew=0;
  double Munew=0;

  //For loop integers
  int j=0;
  int k=0;


  //Check to see if Dosestart = 0 otherwise rearrange for c++
  DoseStart=DoseStart-1;


  int StoreInx=0;
  //Vector Of group intercepts
  arma::vec a(J);
  //Vector of group slopes
  arma::vec b=a;
  //Vectors for Group MCMC



  //Vectors for proposals
  arma::vec aprop = a;
  arma::vec bprop=b;





  arma::vec NumA(J);
  arma::vec NumB(J);
  arma::vec IntA(J);
  arma::vec IntB(J);

  double NumMu = 2;
  double NumSlope=2;
  double IntMu=1;
  double IntSlope=1;

  arma::vec avar(J);
  arma::vec bvar(J);


  double muvar=1;
  double slopevar=1;
  double trialtime=0;
  NumericVector z9(2);
  NumericVector zprop(5);

  double NumSig=2;
  double IntSig=1;
  double sigvar=1;

  arma::vec etaG(J+1);
  arma::mat DoseProb(B1,(J+1)*nDose);
  arma::mat MeanDose(J+1,nDose);
  arma::vec sigstore(B1);
  arma::vec mustore(B1);
  arma::vec slopestore(B1);
  arma::mat astore(B1,J);
  arma::mat bstore=astore;
  double eta1=0;
  double DEC=0;

  arma::vec GroupMem(J);
  arma::vec GroupMemProp=GroupMem;

  arma::vec INVEC(J);
  INVEC.zeros();

  int eta2=0;

  arma::mat MeanVec((J+1),nDose);
  arma::vec NTox(J+1);
  arma::vec OptDose(J+1);
  arma::vec StoppedGroups(J+1);
  arma::vec SuspendGroups(J+1);
  arma::vec GLast(J+1);
  arma::vec nTreated(J+1);
  double stopped=0;
  arma::vec TrialTimes(nSims);
  arma::mat OptimalDoses(nSims,(J+1));
  arma::mat NTOX = OptimalDoses;
  arma::vec GroupVec(nDose);
  arma::vec TriedGroups(J+1);



  arma::mat DoseStore(nSims,Nmax);
  arma::mat GroupStore(nSims,Nmax);
  arma::mat DoseTried((J+1),nDose);
  arma::vec nTreated1(J+1);
  nTreated1.zeros();
INVEC.zeros();

arma::vec INVECNEW=INVEC;
arma::vec Y2;
arma::vec I2;
arma::vec Groups2;
arma::vec Doses2;
arma::vec GetGroup=StoppedGroups;


if(Param1.n_rows==2){
//Two Subgroups

  //Wrap nreps around
  for(rep=0;rep<nSims;rep++){




    if(rep%100==0){

      Rf_PrintValue(wrap(rep));

      Rprintf("Simulations Finished");

    }

    //Wrap i around

    MeanVec.zeros();
    //Reset all trial vectors
    trialtime=0;
    stopped=0;
    OptDose.zeros();
    StoppedGroups.zeros();
    SuspendGroups.zeros();
    TriedGroups.zeros();
    Y.zeros();
    I.zeros();
    Times.zeros();
    ACC.zeros();
    Groups.zeros();
    Doses.zeros();
    GLast.zeros()+10000;
    nTreated.zeros();
    DoseTried.zeros();

    NTox.zeros();




    //Enroll the first patient

    Group=Sample2(groupprob);

    Groups[0]=Group;

    Doses[0]=Dose[DoseStart[Group]];

    nTreated[Group]=nTreated[Group]+1;

    ACC[0]=0;

    if(Family==0){
      Times(0)=R::rexp(Param1(Group,DoseStart[Group]));
    }



    if(Family==1){
      Times[0]=R::rgamma(Param1(Group,DoseStart[Group]),1/Param2(Group,DoseStart[Group]));
    }

    if(Family==2){
      Times[0]=R::rlnorm(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
    }


    if(Family==3){
      if(as_scalar(arma::randu(1))<Param1(Group,DoseStart[Group])){
        Times[0]=as_scalar(arma::randu(1))*T1;
      }else{
        Times[0]=T1+1;
      }
    }

    if(Family==4){
      Times[0]=R::rweibull(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
    }

    DoseTried(Group,DoseStart(Group))=1;


    if(DoseStart[Group]==0){
      nTreated1[Group]=nTreated1[Group]+1;
    }

    for(i=1;i<Nmax;i++){
      //Start the Trial, Let's do this

      trialtime=trialtime+R::rexp(1/Accrue);
      Group=Sample2(groupprob);

      DEC=0;


      //Now its not in a stopped group or a subgroup so let's do our MCMC





      //Enrolling in a new subgroup
      if(TriedGroups[Group]==0){

        nTreated[Group]=nTreated[Group]+1;

        if(DoseStart[Group]==0){
          nTreated1[Group]=nTreated1[Group]+1;
        }

        Groups[i]=Group;

        Doses[i]=Dose[DoseStart[Group]];

        ACC[i]=trialtime;

        if(Family==0){
          Times(i)=R::rexp(Param1(Group,DoseStart[Group]));
        }



        if(Family==1){
          Times[i]=R::rgamma(Param1(Group,DoseStart[Group]),1/Param2(Group,DoseStart[Group]));
        }

        if(Family==2){
          Times[i]=R::rlnorm(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
        }


        if(Family==3){
          if(as_scalar(arma::randu(1))<Param1(Group,DoseStart[Group])){
            Times[i]=as_scalar(arma::randu(1))*T1;
          }else{
            Times[i]=T1+1;
          }
        }

        if(Family==4){
          Times[i]=R::rweibull(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
        }


        DoseTried(Group,DoseStart(Group))=1;
        TriedGroups[Group]=1;

      }else{
        //Not in a new subgroup, let's do our MCMC

        DEC=0;

        while(DEC==0){
          //If this patient group is

          I.zeros();
          Y.zeros();

          //Setup Y
          for(k=0;k<i;k++){
            Y[k]=min1(min1(Times[k],trialtime-ACC[k]),T1);
            if(Times[k]==Y[k]){
              I[k]=1;
            }
          }


          //Inner MCMC loop, only compute neccessary quantities. Reset counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB = IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          avar=avar.zeros()+1;
          bvar=bvar.zeros()+1;
          muvar=1;
          slopevar=1;
          NumSig=2;
          IntSig=1;
          sigvar=1;
          a.zeros();
          b.zeros();
          MeanVec.zeros();
          GroupVec.zeros();
          mu=meanmu;
          slope=meanslope;
   INVEC.zeros();
   INVECNEW.zeros();
   INVEC = INVEC+1;
   a=MeanInts;
   b=MeanSlopes;

          for(m=0;m<B;m++){

            if(m<(B/2 + 2)){
              if(m%100==0){


                for(k=0;k<(J+1);k++){
                  if((IntA[k]/NumA[k])>.5){
                    avar[k]=avar[k]*2;
                  }

                  if((IntA[k]/NumA[k])<.2){
                    avar[k]=avar[k]/2;
                  }








                  if((IntB[k]/NumB[k])>.5){
                    bvar[k]=bvar[k]*2;
                  }

                  if((IntB[k]/NumB[k])<.2){
                    bvar[k]=bvar[k]/2;
                  }



                }


                if((IntMu/NumMu)>.5){
                  muvar=muvar*2;
                }

                if((IntMu/NumMu)<.2){
                  muvar=muvar/2;
                }


                if((IntSlope/NumSlope)>.5){
                  slopevar=slopevar*2;
                }

                if((IntSlope/NumSlope)<.2){
                  slopevar=slopevar/2;
                }





                NumA=NumA.zeros()+2;
                NumB=NumB.zeros()+2;
                IntA=IntA.zeros()+1;
                IntB=IntB.zeros()+1;

                NumMu = 2;
                NumSlope=2;
                IntMu=1;
                IntSlope=1;
                IntSig=1;
                NumSig=2;




              }
            }



            for(k=0;k<aprop.n_rows;k++){

              if(INVEC[k]>0){

              aprop=a;
              aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);






              alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


              U=log(as_scalar(arma::randu(1)));


              if(U<alpha){
                a(k)=aprop(k);
                IntA[k]=IntA[k]+1;
              }

              NumA=NumA+1;


              }
            }
            //Now do a MH step for the slope

            for(k=0;k<aprop.size();k++){

              if(INVEC[k]>0){
              bprop=b;
              bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];



              alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1  +  Like1(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


              U=log(as_scalar(arma::randu(1)));

              z9[0]=alpha;

              //  Rf_PrintValue(wrap(z9));

              if(U<alpha){
                b[k]=bprop[k];
                IntB[k]=IntB[k]+1;
              }

              }

            }


            NumB=NumB+1;



            //Spike And Slab

            for(k=0;k<J;k++){
              if(INVEC[k]==0){
                //ADD Move
          aprop[k]=as_scalar(arma::randn(1))/10;
          bprop[k]=as_scalar(arma::randn(1))/10;


                //Density when 0 is automatically 1 so no proposal needed
                alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);

                alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);

                U=log(as_scalar(arma::randu(1)));

                z9[0]=alpha;

                //  Rf_PrintValue(wrap(z9));

                if(U<alpha){
                  b[k]=bprop[k];
                  a[k]=aprop[k];
                  INVEC[k]=1;
                }


              }else{
                //Delete Move

                //Delete Move
                aprop[k]=0;
                bprop[k]=0;


                //Density when 0 is automatically 1 so no proposal needed
                alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);

                alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

                U=log(as_scalar(arma::randu(1)));



                if(U<alpha){
                  b[k]=0;
                  a[k]=0;
                  INVEC[k]=0;
                }


              }







            }







            Munew = mu + as_scalar(arma::randn(1))*muvar;



            alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint+   Like1(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);




            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              mu=Munew;
              IntMu=IntMu+1;
            }

            NumMu=NumMu+1;




            //Sample the new slope and the new \beta_g coefficients jointly




            slopenew = slope + as_scalar(arma::randn(1))*slopevar;



            //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


            //Proposal ratio for beta

            alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


            alpha =   alpha+  Like1(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


            U=log(as_scalar(arma::randu(1)));



            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              slope=slopenew;
              IntSlope=IntSlope+1;
            }
            NumSlope=NumSlope+1;





            if(m>(B-B1-1)){
              //Make OptDose and StoppedGroups
              StoreInx=m-B1;

              for(j=0;j<J;j++){;
                astore(StoreInx,j)=a(j);
                bstore(StoreInx,j)=b[j];
              }

              mustore[StoreInx]=mu;
              slopestore[StoreInx]=slope;
              sigstore[StoreInx]=sig;


              //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
              //First nDose


              for(k=0;k<(J+1);k++){

                for(j=0;j<nDose;j++){
                  eta1=0;

                  if(k>0){
                    eta1 = a[k-1]+exp(b[k-1]+slope)*Dose[j]+mu;
                  }

                  if(k==0){
                    eta1=mu+exp(slope)*Dose[j];

                  }

                  //For all groups add the intercept and slope*Groups

                  eta1=exp(eta1);

                  if(eta1>100000){
                    DoseProb(StoreInx, k*nDose+j) = 1;

                  }else{
                    DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                  }



                }


              }
              //Now we have a matrix containing our sampled dose probabilities









            }

            //End MCMC
          }



          StoppedGroups.zeros();

          //Do we have a stopped group?
          for(k=0;k<(J+1);k++){
            eta2=0;

            for(j=0;j<DoseProb.n_rows;j++){
              //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
              if(DoseProb(j,k*nDose)>Target[k]){
                eta2++;
              }


            }



            if(eta2>(Upper[k]*DoseProb.n_rows)){
              StoppedGroups[k]=1;
            }

          }




          //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group
          for(k=0;k<(J+1);k++){
            if(nTreated1[k]<3){
              StoppedGroups[k]=0;
            }

            if(nTreated1[k]>2){
              GetGroup.zeros();
              GetGroup[k]=1;

              Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
              I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
              Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);

              if(GetFullyFollowed(Y2,I2,Doses2,Dose,T1)<3){
                StoppedGroups[k]=0;
              }

            }


          }






          //    Rf_PrintValue(wrap(StoppedGroups));
          if(sum(StoppedGroups)==(J+1)){
            //Let's try the groups separately to make sure we aren't stopping unneccessarily
          StoppedGroups.zeros();

          for(g=0;g<(J+1);g++){
            GetGroup.zeros();
            GetGroup[g]=1;

            Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
            I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
            Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);



        Which1 = g;


            if(Which1==0){
              NewMean = meanmu;
              NewSlope=meanslope;
            }else{
              NewMean = meanmu+MeanInts[Which1-1];
              NewSlope=meanslope+MeanSlopes[Which1-1];
            }





            //Inner MCMC loop, only compute neccessary quantities. Reset counters
            NumA=NumA.zeros()+2;
            NumB=NumB.zeros()+2;
            IntA=IntA.zeros()+1;
            IntB = IntB.zeros()+1;
            NumMu = 2;
            NumSlope=2;
            IntMu=1;
            IntSlope=1;
            avar=avar.zeros()+1;
            bvar=bvar.zeros()+1;
            muvar=1;
            slopevar=1;
            NumSig=2;
            IntSig=1;
            sigvar=1;
            a.zeros();
            b.zeros();
            MeanVec.zeros();
            GroupVec.zeros();
            mu=NewMean;
            slope=NewSlope;


            for(m=0;m<B;m++){

              if(m<(B/2 + 2)){
                if(m%100==0){




                  if((IntMu/NumMu)>.5){
                    muvar=muvar*2;
                  }

                  if((IntMu/NumMu)<.2){
                    muvar=muvar/2;
                  }


                  if((IntSlope/NumSlope)>.5){
                    slopevar=slopevar*2;
                  }

                  if((IntSlope/NumSlope)<.2){
                    slopevar=slopevar/2;
                  }







                  NumMu = 2;
                  NumSlope=2;
                  IntMu=1;
                  IntSlope=1;
                  IntSig=1;
                  NumSig=2;




                }
              }









              Munew = mu + as_scalar(arma::randn(1))*muvar;



              alpha= .5*pow((mu-NewMean),2)/varint -.5*pow((Munew-NewMean),2)/varint+   Like2(Y2, I2,  Doses2,  Munew,slope, sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope,  sig,Y2.n_rows);




              U=log(as_scalar(arma::randu(1)));

              if(U<alpha){
                mu=Munew;
                IntMu=IntMu+1;
              }

              NumMu=NumMu+1;




              //Sample the new slope and the new \beta_g coefficients jointly




              slopenew = slope + as_scalar(arma::randn(1))*slopevar;



              //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


              //Proposal ratio for beta

              alpha =.5*pow((slope-NewSlope),2)/varbeta -.5*pow((slopenew-NewSlope),2)/varbeta;


              alpha =   alpha+  Like2(Y2, I2,  Doses2,   mu,slopenew,sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope, sig,Y2.n_rows);


              U=log(as_scalar(arma::randu(1)));



              U=log(as_scalar(arma::randu(1)));

              if(U<alpha){
                slope=slopenew;
                IntSlope=IntSlope+1;
              }
              NumSlope=NumSlope+1;





              if(m>(B-B1-1)){
                //Make OptDose and StoppedGroups
                StoreInx=m-B1;



                //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
                //First nDose


                k=Which1;

                for(j=0;j<nDose;j++){
                  eta1=0;




                  eta1=mu+exp(slope)*Dose[j];



                  //For all groups add the intercept and slope*Groups

                  eta1=exp(eta1);

                  if(eta1>100000){
                    DoseProb(StoreInx, k*nDose+j) = 1;

                  }else{
                    DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                  }



                }



                //Now we have a matrix containing our sampled dose probabilities









              }

              //End MCMC
            }
            //Is our very last group stopped?

            k=Which1;
            eta2=0;

            for(j=0;j<DoseProb.n_rows;j++){
              //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
              if(DoseProb(j,k*nDose)>Target[k]){
                eta2++;
              }


            }


            if(eta2>(Upper[Which1]*DoseProb.n_rows)){
              StoppedGroups[Which1]=1;
            }






          }





          }




          //Should we stop the trial?
          if(sum(StoppedGroups)==(J+1)){
            stopped=1;
            break;
          }









          if(StoppedGroups[Group]==0){
            //Enroll this patient
            DEC=1;


            if(sum(StoppedGroups)>0){
              //If we have a stopped group, REMOVE IT

              Y2=ReturnStoppedY(Y,I,Groups,StoppedGroups,i);
              I2=ReturnStoppedI(Y,I,Groups,StoppedGroups,i);
              Groups2=ReturnStoppedGroups(Y,I,Groups,StoppedGroups,i);
              Doses2=ReturnStoppedY(Doses,I,Groups,StoppedGroups,i);








          //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group
          for(k=0;k<(J+1);k++){
            if(nTreated1[k]<3){
              StoppedGroups[k]=0;
            }

            if(nTreated1[k]>2){
              GetGroup.zeros();
              GetGroup[k]=1;

              Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
              I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
              Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);

              if(GetFullyFollowed(Y2,I2,Doses2,Dose,T1)<3){
                StoppedGroups[k]=0;
              }

            }


          }





          if(sum(StoppedGroups)==(J+1)){
            stopped=1;

            Rf_PrintValue(wrap(StoppedGroups));
            Rprintf("It Stopped");

            Rf_PrintValue(wrap(Y2));
            Rf_PrintValue(wrap(I2));
            Rf_PrintValue(wrap(Doses2));


            break;
          }




            }





          }else{

            trialtime=trialtime+R::rexp(1/Accrue);
            Group=Sample2(groupprob);

          }

//Do We have ANOTHER STOPPED GROUP?
//Should we stop the trial?







if(sum(StoppedGroups)==(J+1)){
  stopped=1;

  Rf_PrintValue(wrap(StoppedGroups));
  Rprintf("It Stopped");

  Rf_PrintValue(wrap(Y2));
  Rf_PrintValue(wrap(I2));
  Rf_PrintValue(wrap(Doses2));


  break;
}


        }

        //Ok now let's enroll this patient, Whats the Optimal Dose

        //If Stopped==1, Get out of the trial
        if(stopped==1){
          break;
        }



        if(MEANS==0){
          for(k=0;k<(J+1);k++){
            for(j=0;j<nDose;j++){
              MeanVec(k,j)=sum(DoseProb.col(k*nDose+j))/DoseProb.n_rows;

            }

          }



        }


        //Determine Optimal Dose
          k=Group;
          eta1=MinVec(abs(MeanVec.row(k).t()-Target(k)));
          j=0;




          GroupVec = abs(MeanVec.row(k).t()-Target(k));

          while(GroupVec[j]!=eta1 ){
            j++;
            //Loop Proceeds until the minimum difference is reached
          }

          OptDose[k]=j;




        //Now we have the subgroup specific optimal doses in the vector OptDose



        //Now let's see if the dose is bigger than the largest dose.
        if(OptDose[Group]>=nDose){
          OptDose[Group]=nDose-1;
        }



        //Is this dose bigger than the dose higher than what's been tried?
        if(DoseTried(Group,OptDose[Group])==0){
          j=0;

          while(DoseTried(Group,j)==1){
            j++;
          }

          OptDose[Group]=j;

        }






        //Assign Patient to subgroup, optdose, add accrual time.

        Groups[i]=Group;
        if(OptDose[Group]==0){
          //If we assign one at the lowest dose add this
          nTreated1[Group]=nTreated1[Group]+1;
        }
        nTreated[Group]=nTreated[Group]+1;

        //If this is the 3rd patient treated in this subgroup with the lowest dose, suspend accrual for TGroup
        if(nTreated[Group]==3){
          SuspendGroups[Group]=1;
          GLast[Group]=trialtime;
        }



        Doses[i]=Dose[OptDose[Group]];


        ACC[i]=trialtime;
        //Based on the Group and Optimal dose, randomly generate the next patient Data


        if(Family==0){
          Times(i)=R::rexp(Param1(Group,OptDose[Group]));
        }



        if(Family==1){
          Times[i]=R::rgamma(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
        }

        if(Family==2){
          Times[i]=R::rlnorm(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
        }


        if(Family==3){
          if(as_scalar(arma::randu(1))<Param1(Group,OptDose[Group])){
            Times[i]=as_scalar(arma::randu(1))*T1;
          }else{
            Times[i]=T1+1;
          }
        }

        if(Family==4){
          Times[i]=R::rweibull(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
        }


        DoseTried(Group,OptDose(Group))=1;












      }




      //End of the i loop, now we have ALL the patients
    }



    //Now let's make the final decision.
    if(stopped==1){
      //The trial Ended Early, i-1 is the index of the last patient enrolled.


      for(k=0;k<i;k++){
        for(j=0;j<(J+1);j++){
          if(Groups[k]==j){
            NTox[j]=NTox[j]+I[k];
          }
        }
      }


      //All Optimal Doses are -1
      for(k=0;k<(J+1);k++){
        OptDose[k]=-1;
      }


    }else{
      //The trial did NOT end early. Follow patients out until end of trial and get the optimal doses
      trialtime=trialtime+T1;



      Y.zeros();
      I.zeros();

      //Setup Y
      for(k=0;k<Nmax;k++){

        Y[k]=min1(Times[k],T1);


        if(Times[k]<T1){
          I[k]=1;
        }
      }





      //Inner MCMC loop, only compute neccessary quantities. Reset counters
      NumA=NumA.zeros()+2;
      NumB=NumB.zeros()+2;
      IntA=IntA.zeros()+1;
      IntB = IntB.zeros()+1;
      NumMu = 2;
      NumSlope=2;
      IntMu=1;
      IntSlope=1;
      avar=avar.zeros()+1;
      bvar=bvar.zeros()+1;
      muvar=1;
      slopevar=1;
      NumSig=2;
      IntSig=1;
      sigvar=1;
      a.zeros();
      b.zeros();
      MeanVec.zeros();
      GroupVec.zeros();
      mu=meanmu;
      slope=meanslope;
INVEC.zeros();
INVEC = INVEC+1;
a=MeanInts;
b=MeanSlopes;

      for(m=0;m<B;m++){

        if(m<(B/2 + 2)){
          if(m%100==0){


            for(k=0;k<(J+1);k++){
              if((IntA[k]/NumA[k])>.5){
                avar[k]=avar[k]*2;
              }

              if((IntA[k]/NumA[k])<.2){
                avar[k]=avar[k]/2;
              }








              if((IntB[k]/NumB[k])>.5){
                bvar[k]=bvar[k]*2;
              }

              if((IntB[k]/NumB[k])<.2){
                bvar[k]=bvar[k]/2;
              }



            }


            if((IntMu/NumMu)>.5){
              muvar=muvar*2;
            }

            if((IntMu/NumMu)<.2){
              muvar=muvar/2;
            }


            if((IntSlope/NumSlope)>.5){
              slopevar=slopevar*2;
            }

            if((IntSlope/NumSlope)<.2){
              slopevar=slopevar/2;
            }





            NumA=NumA.zeros()+2;
            NumB=NumB.zeros()+2;
            IntA=IntA.zeros()+1;
            IntB=IntB.zeros()+1;

            NumMu = 2;
            NumSlope=2;
            IntMu=1;
            IntSlope=1;
            IntSig=1;
            NumSig=2;




          }
        }



        for(k=0;k<aprop.n_rows;k++){

          if(INVEC[k]>0){

            aprop=a;
            aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);






            alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


            U=log(as_scalar(arma::randu(1)));


            if(U<alpha){
              a(k)=aprop(k);
              IntA[k]=IntA[k]+1;
            }

            NumA=NumA+1;


          }
        }
        //Now do a MH step for the slope

        for(k=0;k<aprop.size();k++){

          if(INVEC[k]>0){
            bprop=b;
            bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];



            alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1  +  Like1(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


            U=log(as_scalar(arma::randu(1)));

            z9[0]=alpha;

            //  Rf_PrintValue(wrap(z9));

            if(U<alpha){
              b[k]=bprop[k];
              IntB[k]=IntB[k]+1;
            }

          }

        }


        NumB=NumB+1;



        //Spike And Slab

        for(k=0;k<J;k++){
          if(INVEC[k]==0){
            //ADD Move
            aprop[k]=as_scalar(arma::randn(1))/10;
            bprop[k]=as_scalar(arma::randn(1))/10;


            //Density when 0 is automatically 1 so no proposal needed
            alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);

            alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);

            U=log(as_scalar(arma::randu(1)));

            z9[0]=alpha;

            //  Rf_PrintValue(wrap(z9));

            if(U<alpha){
              b[k]=bprop[k];
              a[k]=aprop[k];
              INVEC[k]=1;
            }


          }else{
            //Delete Move

            //Delete Move
            aprop[k]=0;
            bprop[k]=0;


            //Density when 0 is automatically 1 so no proposal needed
            alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);

            alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

            U=log(as_scalar(arma::randu(1)));



            if(U<alpha){
              b[k]=0;
              a[k]=0;
            }


          }







        }







        Munew = mu + as_scalar(arma::randn(1))*muvar;



        alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint+   Like1(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);




        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          mu=Munew;
          IntMu=IntMu+1;
        }

        NumMu=NumMu+1;




        //Sample the new slope and the new \beta_g coefficients jointly




        slopenew = slope + as_scalar(arma::randn(1))*slopevar;



        //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


        //Proposal ratio for beta

        alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


        alpha =   alpha+  Like1(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


        U=log(as_scalar(arma::randu(1)));



        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          slope=slopenew;
          IntSlope=IntSlope+1;
        }
        NumSlope=NumSlope+1;





        if(m>(B-B1-1)){
          //Make OptDose and StoppedGroups
          StoreInx=m-B1;

          for(j=0;j<J;j++){;
            astore(StoreInx,j)=a(j);
            bstore(StoreInx,j)=b[j];
          }

          mustore[StoreInx]=mu;
          slopestore[StoreInx]=slope;
          sigstore[StoreInx]=sig;


          //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
          //First nDose


          for(k=0;k<(J+1);k++){

            for(j=0;j<nDose;j++){
              eta1=0;

              if(k>0){
                eta1 = a[k-1]+exp(b[k-1]+slope)*Dose[j]+mu;
              }

              if(k==0){
                eta1=mu+exp(slope)*Dose[j];

              }

              //For all groups add the intercept and slope*Groups

              eta1=exp(eta1);

              if(eta1>100000){
                DoseProb(StoreInx, k*nDose+j) = 1;

              }else{
                DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

              }



            }


          }
          //Now we have a matrix containing our sampled dose probabilities









        }

        //End MCMC
      }



      //Get Optimal Dose


      StoppedGroups.zeros();
      //Get MeanVec
      //reset to means




      //Now we have the mean vector



      //Do we have a stopped group?
      for(k=0;k<(J+1);k++){
        eta2=0;

        for(j=0;j<DoseProb.n_rows;j++){
          //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
          if(DoseProb(j,k*nDose)>Target[k]){
            eta2++;
          }


        }







        if(eta2>(Upper[k]*DoseProb.n_rows)){
          StoppedGroups[k]=1;
        }

      }






    if(sum(StoppedGroups)>0){
      //If we have a stopped group, REMOVE IT

      Y2=ReturnStoppedY(Y,I,Groups,StoppedGroups,i);
      I2=ReturnStoppedI(Y,I,Groups,StoppedGroups,i);
      Groups2=ReturnStoppedGroups(Y,I,Groups,StoppedGroups,i);
      Doses2=ReturnStoppedY(Doses,I,Groups,StoppedGroups,i);









      //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group

      //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group
      for(k=0;k<(J+1);k++){
        if(nTreated1[k]<3){
          StoppedGroups[k]=0;
        }

        if(nTreated1[k]>2){
        GetGroup.zeros();
        GetGroup[k]=1;

        Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
        I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
        Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);

        if(GetFullyFollowed(Y2,I2,Doses2,Dose,T1)<3){
          StoppedGroups[k]=0;
        }

        }


      }




      //If all groups are stopped, let's do separate trials JUST TO MAKE SURE ONE GROUP IS NOT DOMINATING THE DECISION TO STOP!!


      //    Rf_PrintValue(wrap(StoppedGroups));
      if(sum(StoppedGroups)==(J+1)){
        //Let's try the groups separately to make sure we aren't stopping unneccessarily
        StoppedGroups.zeros();

        for(g=0;g<(J+1);g++){
          GetGroup.zeros();
          GetGroup[g]=1;

          Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
          I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
          Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);



          Which1 = g;


          if(Which1==0){
            NewMean = meanmu;
            NewSlope=meanslope;
          }else{
            NewMean = meanmu+MeanInts[Which1-1];
            NewSlope=meanslope+MeanSlopes[Which1-1];
          }





          //Inner MCMC loop, only compute neccessary quantities. Reset counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB = IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          avar=avar.zeros()+1;
          bvar=bvar.zeros()+1;
          muvar=1;
          slopevar=1;
          NumSig=2;
          IntSig=1;
          sigvar=1;
          a.zeros();
          b.zeros();
          MeanVec.zeros();
          GroupVec.zeros();
          mu=NewMean;
          slope=NewSlope;


          for(m=0;m<B;m++){

            if(m<(B/2 + 2)){
              if(m%100==0){




                if((IntMu/NumMu)>.5){
                  muvar=muvar*2;
                }

                if((IntMu/NumMu)<.2){
                  muvar=muvar/2;
                }


                if((IntSlope/NumSlope)>.5){
                  slopevar=slopevar*2;
                }

                if((IntSlope/NumSlope)<.2){
                  slopevar=slopevar/2;
                }







                NumMu = 2;
                NumSlope=2;
                IntMu=1;
                IntSlope=1;
                IntSig=1;
                NumSig=2;




              }
            }









            Munew = mu + as_scalar(arma::randn(1))*muvar;



            alpha= .5*pow((mu-NewMean),2)/varint -.5*pow((Munew-NewMean),2)/varint+   Like2(Y2, I2,  Doses2,  Munew,slope, sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope,  sig,Y2.n_rows);




            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              mu=Munew;
              IntMu=IntMu+1;
            }

            NumMu=NumMu+1;




            //Sample the new slope and the new \beta_g coefficients jointly




            slopenew = slope + as_scalar(arma::randn(1))*slopevar;



            //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


            //Proposal ratio for beta

            alpha =.5*pow((slope-NewSlope),2)/varbeta -.5*pow((slopenew-NewSlope),2)/varbeta;


            alpha =   alpha+  Like2(Y2, I2,  Doses2,   mu,slopenew,sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope, sig,Y2.n_rows);


            U=log(as_scalar(arma::randu(1)));



            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              slope=slopenew;
              IntSlope=IntSlope+1;
            }
            NumSlope=NumSlope+1;





            if(m>(B-B1-1)){
              //Make OptDose and StoppedGroups
              StoreInx=m-B1;



              //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
              //First nDose


              k=Which1;

              for(j=0;j<nDose;j++){
                eta1=0;




                eta1=mu+exp(slope)*Dose[j];



                //For all groups add the intercept and slope*Groups

                eta1=exp(eta1);

                if(eta1>100000){
                  DoseProb(StoreInx, k*nDose+j) = 1;

                }else{
                  DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                }



              }



              //Now we have a matrix containing our sampled dose probabilities









            }

            //End MCMC
          }
          //Is our very last group stopped?

          k=Which1;
          eta2=0;

          for(j=0;j<DoseProb.n_rows;j++){
            //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
            if(DoseProb(j,k*nDose)>Target[k]){
              eta2++;
            }


          }


          if(eta2>(Upper[Which1]*DoseProb.n_rows)){
            StoppedGroups[Which1]=1;
          }






        }





      }







    }








      if(MEANS==0){
        for(k=0;k<(J+1);k++){
          for(j=0;j<nDose;j++){
            MeanVec(k,j)=sum(DoseProb.col(k*nDose+j))/DoseProb.n_rows;

          }

        }



      }




      //Determine Optimal Doses
      for(k=0;k<(J+1);k++){

        if(StoppedGroups(k)==1){
          OptDose(k)=-1;
        }else{
          eta1=MinVec(abs(MeanVec.row(k).t()-Target[k]));
          j=0;

          GroupVec = abs(MeanVec.row(k).t()-Target[k]);

          //  Rf_PrintValue(wrap(MeanVec));


          while(GroupVec[j]!=eta1 ){
            j++;
            //Loop Proceeds until the minimum difference is reached
          }

          OptDose[k]=j;
        }


      }




      //Is this dose bigger than the dose higher than what's been tried?
      for(k=0;k<(J+1);k++){
        if(StoppedGroups[k]==0){
          for(j=0;j<(J+1);j++){
            if(OptDose[j]>=nDose){
              OptDose[j]=nDose-1;
            }
          }




        }



      }

      //Is this dose bigger than the dose higher than what's been tried?
      for(k=0;k<(J+1);k++){
        if(StoppedGroups(k)==0){
          if(DoseTried(k,OptDose[k])==0){
            j=0;

            while(DoseTried(k,j)==1){
              j++;
            }

            OptDose[k]=j;

          }
        }
      }






      //Num Toxicities
      for(k=0;k<Nmax;k++){
        for(j=0;j<(J+1);j++){
          if(Groups[k]==j){
            NTox[j]=NTox[j]+I[k];
          }
        }
      }


    }

    OptDose=OptDose+1;


    TrialTimes(rep)=trialtime;
    for(j=0;j<(J+1);j++){
      OptimalDoses(rep,j)=OptDose(j);
      NTOX(rep,j)=NTox(j);
    }




    //Fill in Doses and Groups Treated
    for(j=0;j<Nmax;j++){
      DoseStore(rep,j)=Doses(j);
      GroupStore(rep,j)=Groups(j);
    }




    //End Trial
  }
}else{
  //More than Two Subgroups


  arma::vec INVEC(J);
  INVEC.zeros();

  int eta2=0;

  arma::mat MeanVec((J+1),nDose);
  arma::vec NTox(J+1);
  arma::vec OptDose(J+1);
  arma::vec StoppedGroups(J+1);
  arma::vec SuspendGroups(J+1);
  arma::vec GLast(J+1);
  arma::vec nTreated(J+1);
  double stopped=0;
  arma::vec TrialTimes(nSims);
  arma::mat OptimalDoses(nSims,(J+1));
  arma::mat NTOX = OptimalDoses;
  arma::vec GroupVec(nDose);
  arma::vec TriedGroups(J+1);



  arma::mat DoseStore(nSims,Nmax);
  arma::mat GroupStore(nSims,Nmax);
  arma::mat DoseTried((J+1),nDose);
  arma::vec nTreated1(J+1);
  nTreated1.zeros();
  INVEC.zeros();

  arma::vec INVECNEW=INVEC;
  arma::vec Y2;
  arma::vec I2;
  arma::vec Groups2;
  arma::vec Doses2;
  arma::vec GetGroup=StoppedGroups;


  //Wrap nreps around
  for(rep=0;rep<nSims;rep++){




    if(rep%100==0){

      Rf_PrintValue(wrap(rep));

      Rprintf("Simulations Finished");

    }

    //Wrap i around

    MeanVec.zeros();
    //Reset all trial vectors
    trialtime=0;
    stopped=0;
    OptDose.zeros();
    StoppedGroups.zeros();
    SuspendGroups.zeros();
    TriedGroups.zeros();
    Y.zeros();
    I.zeros();
    Times.zeros();
    ACC.zeros();
    Groups.zeros();
    Doses.zeros();
    GLast.zeros()+10000;
    nTreated.zeros();
    DoseTried.zeros();

    NTox.zeros();




    //Enroll the first patient

    Group=Sample2(groupprob);

    Groups[0]=Group;

    Doses[0]=Dose[DoseStart[Group]];

    nTreated[Group]=nTreated[Group]+1;

    ACC[0]=0;

    if(Family==0){
      Times(0)=R::rexp(Param1(Group,DoseStart[Group]));
    }



    if(Family==1){
      Times[0]=R::rgamma(Param1(Group,DoseStart[Group]),1/Param2(Group,DoseStart[Group]));
    }

    if(Family==2){
      Times[0]=R::rlnorm(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
    }


    if(Family==3){
      if(as_scalar(arma::randu(1))<Param1(Group,DoseStart[Group])){
        Times[0]=as_scalar(arma::randu(1))*T1;
      }else{
        Times[0]=T1+1;
      }
    }

    if(Family==4){
      Times[0]=R::rweibull(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
    }

    DoseTried(Group,DoseStart(Group))=1;


    if(DoseStart[Group]==0){
      nTreated1[Group]=nTreated1[Group]+1;
    }

    for(i=1;i<Nmax;i++){
      //Start the Trial, Let's do this

      DEC=0;

      trialtime=trialtime+R::rexp(1/Accrue);
      Group=Sample2(groupprob);
      //Now its not in a stopped group or a subgroup so let's do our MCMC





      //Enrolling in a new subgroup
      if(TriedGroups[Group]==0){

        nTreated[Group]=nTreated[Group]+1;

        if(DoseStart[Group]==0){
          nTreated1[Group]=nTreated1[Group]+1;
        }

        Groups[i]=Group;

        Doses[i]=Dose[DoseStart[Group]];

        ACC[i]=trialtime;

        if(Family==0){
          Times(i)=R::rexp(Param1(Group,DoseStart[Group]));
        }



        if(Family==1){
          Times[i]=R::rgamma(Param1(Group,DoseStart[Group]),1/Param2(Group,DoseStart[Group]));
        }

        if(Family==2){
          Times[i]=R::rlnorm(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
        }


        if(Family==3){
          if(as_scalar(arma::randu(1))<Param1(Group,DoseStart[Group])){
            Times[i]=as_scalar(arma::randu(1))*T1;
          }else{
            Times[i]=T1+1;
          }
        }

        if(Family==4){
          Times[i]=R::rweibull(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
        }


        DoseTried(Group,DoseStart(Group))=1;
        TriedGroups[Group]=1;

      }else{
        //Not in a new subgroup, let's do our MCMC

        DEC=0;

        while(DEC==0){
          //If this patient group is

          I.zeros();
          Y.zeros();

          //Setup Y
          for(k=0;k<i;k++){
            Y[k]=min1(min1(Times[k],trialtime-ACC[k]),T1);
            if(Times[k]==Y[k]){
              I[k]=1;
            }
          }


          //Inner MCMC loop, only compute neccessary quantities. Reset counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB = IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          avar=avar.zeros()+1;
          bvar=bvar.zeros()+1;
          muvar=1;
          slopevar=1;
          NumSig=2;
          IntSig=1;
          sigvar=1;
          a.zeros();
          b.zeros();
          MeanVec.zeros();
          GroupVec.zeros();
          mu=meanmu;
          slope=meanslope;
          INVEC.zeros();
          INVECNEW.zeros();
          INVEC = INVEC+1;
          a=MeanInts;
          b=MeanSlopes;


          for(m=0;m<GroupMem.n_rows;m++){
            GroupMem[m]=m+1;
          }

          GroupMemProp=GroupMem;


          for(m=0;m<B;m++){

            if(m<(B/2 + 2)){
              if(m%100==0){


                for(k=0;k<(J+1);k++){
                  if((IntA[k]/NumA[k])>.5){
                    avar[k]=avar[k]*2;
                  }

                  if((IntA[k]/NumA[k])<.2){
                    avar[k]=avar[k]/2;
                  }








                  if((IntB[k]/NumB[k])>.5){
                    bvar[k]=bvar[k]*2;
                  }

                  if((IntB[k]/NumB[k])<.2){
                    bvar[k]=bvar[k]/2;
                  }



                }


                if((IntMu/NumMu)>.5){
                  muvar=muvar*2;
                }

                if((IntMu/NumMu)<.2){
                  muvar=muvar/2;
                }


                if((IntSlope/NumSlope)>.5){
                  slopevar=slopevar*2;
                }

                if((IntSlope/NumSlope)<.2){
                  slopevar=slopevar/2;
                }





                NumA=NumA.zeros()+2;
                NumB=NumB.zeros()+2;
                IntA=IntA.zeros()+1;
                IntB=IntB.zeros()+1;

                NumMu = 2;
                NumSlope=2;
                IntMu=1;
                IntSlope=1;
                IntSig=1;
                NumSig=2;




              }
            }



            for(k=0;k<aprop.n_rows;k++){

              if(INVEC[k]>0){

                aprop=a;
                aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);






                alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


                U=log(as_scalar(arma::randu(1)));


                if(U<alpha){
                  a(k)=aprop(k);
                  IntA[k]=IntA[k]+1;
                }

                NumA=NumA+1;


              }
            }
            //Now do a MH step for the slope

            for(k=0;k<aprop.size();k++){

              if(INVEC[k]>0){
                bprop=b;
                bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];



                alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1  +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


                U=log(as_scalar(arma::randu(1)));

                z9[0]=alpha;

                //  Rf_PrintValue(wrap(z9));

                if(U<alpha){
                  b[k]=bprop[k];
                  IntB[k]=IntB[k]+1;
                }

              }

            }


            NumB=NumB+1;



            //Spike And Slab


            U=as_scalar(arma::randu(1));


            if((3*U)<1){
              //Add All or delete all

              if(sum(INVEC)==(J)){
                //Auto Delete Move




                aprop.zeros();
                bprop.zeros();

                GroupMemProp=GroupMem;
                //Now Randomly Draw Our New Group Assignment


                INVECNEW.zeros();




                GroupMemProp.zeros();




                alpha=0;
                for(k=0;k<bprop.n_rows;k++){
                  alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

                }





                //Density when 0 is automatically 1 so no proposal needed
                alpha =  alpha + LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



                U=log(as_scalar(arma::randu(1)));



                if(U<alpha){
                  b.zeros();
                  a.zeros();
                  GroupMem.zeros();
                  INVEC.zeros();

                  //If Any other GroupMEM are in this group, we need to delete them.



                }






              }else{


                if(sum(INVEC)==0){
                  //Autho Add Move

                  //ADD Move

                  for(k=0;k<aprop.n_rows;k++){
                    aprop[k]=as_scalar(arma::randn(1));
                    bprop[k]=as_scalar(arma::randn(1))/10;
                    GroupMemProp[k]=k+1;
                  }



                  //Density when 0 is automatically 1 so no proposal needed

                  alpha =  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



                  for(k=0;k<aprop.n_rows;k++){


                    alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]) -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 ;

                  }




                  alpha=alpha + log(sum(INVEC)+1);

                  U=log(as_scalar(arma::randu(1)));

                  z9[0]=alpha;

                  //  Rf_PrintValue(wrap(z9));

                  if(U<alpha){

                    for(k=0;k<aprop.n_rows;k++){
                      b[k]=bprop[k];
                      a[k]=aprop[k];
                      INVEC[k]=1;
                      GroupMem[k]=k+1;
                    }


                  }




                }else{
                  //Randomly Decide to add or Delete

                  U=as_scalar(arma::randu(1));


                  if(U<.5){
                    //Add


                    //ADD Move

                    for(k=0;k<aprop.n_rows;k++){
                      aprop[k]=as_scalar(arma::randn(1));
                      bprop[k]=as_scalar(arma::randn(1))/10;
                      GroupMemProp[k]=k+1;
                    }



                    //Density when 0 is automatically 1 so no proposal needed

                    alpha =  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



                    for(k=0;k<aprop.n_rows;k++){


                      alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]) -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 ;

                    }




                    alpha=alpha + log(sum(INVEC)+1);

                    U=log(as_scalar(arma::randu(1)));

                    z9[0]=alpha;

                    //  Rf_PrintValue(wrap(z9));

                    if(U<alpha){

                      for(k=0;k<aprop.n_rows;k++){
                        b[k]=bprop[k];
                        a[k]=aprop[k];
                        INVEC[k]=1;
                        GroupMem[k]=k+1;
                      }


                    }



                  }else{
                    //Delete All




                    aprop.zeros();
                    bprop.zeros();

                    GroupMemProp=GroupMem;
                    //Now Randomly Draw Our New Group Assignment


                    INVECNEW.zeros();




                    GroupMemProp.zeros();




                    alpha=0;
                    for(k=0;k<bprop.n_rows;k++){
                      alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

                    }





                    //Density when 0 is automatically 1 so no proposal needed
                    alpha =  alpha + LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



                    U=log(as_scalar(arma::randu(1)));



                    if(U<alpha){
                      b.zeros();
                      a.zeros();
                      GroupMem.zeros();
                      INVEC.zeros();

                      //If Any other GroupMEM are in this group, we need to delete them.



                    }




                  }





                }




              }









            }else{



              //Should we swap?
              if((sum(INVEC)==J) || (sum(INVEC)==0)){
                //Add/Delete



                k=Sample1(J);

                if(INVEC[k]==0){
                  //ADD Move
                  aprop=a;
                  bprop=b;
                  aprop[k]=as_scalar(arma::randn(1));
                  bprop[k]=as_scalar(arma::randn(1))/10;

                  GroupMemProp=GroupMem;
                  GroupMemProp[k]=k+1;


                  //Density when 0 is automatically 1 so no proposal needed
                  alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

                  alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);


                  alpha=alpha + log(sum(INVEC)+1);

                  U=log(as_scalar(arma::randu(1)));

                  z9[0]=alpha;

                  //  Rf_PrintValue(wrap(z9));

                  if(U<alpha){
                    b[k]=bprop[k];
                    a[k]=aprop[k];
                    INVEC[k]=1;
                    GroupMem[k]=k+1;
                  }


                }else{
                  //Delete Move
                  aprop=a;
                  bprop=b;
                  //Delete Move
                  aprop[k]=0;
                  bprop[k]=0;

                  GroupMemProp=GroupMem;
                  //Now Randomly Draw Our New Group Assignment


                  INVECNEW=INVEC;
                  INVECNEW[k]=0;



                  GroupMemProp[k]=GetNewGroup(INVECNEW);










                  //Density when 0 is automatically 1 so no proposal needed
                  alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

                  alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

                  alpha=alpha-log(sum(INVEC)+1);

                  U=log(as_scalar(arma::randu(1)));



                  if(U<alpha){
                    b[k]=0;
                    a[k]=0;
                    GroupMem=GroupMemProp;
                    INVEC[k]=0;

                    //If Any other GroupMEM are in this group, we need to delete them.

                    for(j=0;j<GroupMem.n_rows;j++){
                      if(GroupMem[j]==(k+1)){
                        GroupMem[j]=GroupMem[k];
                      }
                    }



                  }


                }


              }else{




                if((3*U)>2){
                  //Swap Move
                  //Randomly Pick One In and One Out

                  IntIN = GetIn(INVEC);
                  IntOUT= GetOut(INVEC);

                  //Now IntIN and IntOUT contain the entries that are currently in or out





                }else{
                  //Add/Delete



                  k=Sample1(J);

                  if(INVEC[k]==0){
                    //ADD Move
                    aprop=a;
                    bprop=b;
                    aprop[k]=as_scalar(arma::randn(1));
                    bprop[k]=as_scalar(arma::randn(1))/10;

                    GroupMemProp=GroupMem;
                    GroupMemProp[k]=k+1;


                    //Density when 0 is automatically 1 so no proposal needed
                    alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

                    alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);


                    alpha=alpha + log(sum(INVEC)+1);

                    U=log(as_scalar(arma::randu(1)));

                    z9[0]=alpha;

                    //  Rf_PrintValue(wrap(z9));

                    if(U<alpha){
                      b[k]=bprop[k];
                      a[k]=aprop[k];
                      INVEC[k]=1;
                      GroupMem[k]=k+1;
                    }


                  }else{
                    //Delete Move
                    aprop=a;
                    bprop=b;
                    //Delete Move
                    aprop[k]=0;
                    bprop[k]=0;

                    GroupMemProp=GroupMem;
                    //Now Randomly Draw Our New Group Assignment


                    INVECNEW=INVEC;
                    INVECNEW[k]=0;



                    GroupMemProp[k]=GetNewGroup(INVECNEW);










                    //Density when 0 is automatically 1 so no proposal needed
                    alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

                    alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

                    alpha=alpha-log(sum(INVEC)+1);

                    U=log(as_scalar(arma::randu(1)));



                    if(U<alpha){
                      b[k]=0;
                      a[k]=0;
                      GroupMem=GroupMemProp;
                      INVEC[k]=0;

                      //If Any other GroupMEM are in this group, we need to delete them.

                      for(j=0;j<GroupMem.n_rows;j++){
                        if(GroupMem[j]==(k+1)){
                          GroupMem[j]=GroupMem[k];
                        }
                      }



                    }


                  }


                }



              }

            }






            //Swap Group

            for(k=0;k<J;k++){
              if(INVEC[k]==0){

                GroupMemProp=GroupMem;
                //Now Randomly Draw Our New Group Assignment


                INVECNEW=INVEC;
                INVECNEW[k]=0;



                GroupMemProp[k]=GetNewGroup(INVECNEW);










                //Density when 0 is automatically 1 so no proposal needed
                alpha = LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



                U=log(as_scalar(arma::randu(1)));



                if(U<alpha){
                  GroupMem[k]=GroupMemProp[k];
                }
              }

            }






            Munew = mu + as_scalar(arma::randn(1))*muvar;



            alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint+   LikeMULTI(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);




            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              mu=Munew;
              IntMu=IntMu+1;
            }

            NumMu=NumMu+1;




            //Sample the new slope and the new \beta_g coefficients jointly




            slopenew = slope + as_scalar(arma::randn(1))*slopevar;



            //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


            //Proposal ratio for beta

            alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


            alpha =   alpha+  LikeMULTI(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


            U=log(as_scalar(arma::randu(1)));



            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              slope=slopenew;
              IntSlope=IntSlope+1;
            }
            NumSlope=NumSlope+1;


            if(m>(B-B1-1)){
              //Make OptDose and StoppedGroups
              StoreInx=m-B1;

              for(j=0;j<J;j++){;
                astore(StoreInx,j)=INVEC(j);
                bstore(StoreInx,j)=GroupMem[j];
              }

              mustore[StoreInx]=mu;
              slopestore[StoreInx]=slope;
              sigstore[StoreInx]=sig;


              //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
              //First nDose

              aprop.zeros();
              bprop.zeros();
              for(j=0;j<J;j++){
                for(k=0;k<J;k++){
                  if(GroupMem[j]==(k+1)){
                    aprop[j]=a[k];
                    bprop[j]=b[k];
                  }
                }
              }


              for(k=0;k<(J+1);k++){

                for(j=0;j<nDose;j++){
                  eta1=0;

                  if(k>0){
                    eta1 = aprop[k-1]+exp(bprop[k-1]+slope)*Dose[j]+mu;
                  }

                  if(k==0){
                    eta1=mu+exp(slope)*Dose[j];

                  }

                  //For all groups add the intercept and slope*Groups

                  eta1=exp(eta1);

                  if(eta1>100000){
                    DoseProb(StoreInx, k*nDose+j) = 1;

                  }else{
                    DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                  }



                }


              }
              //Now we have a matrix containing our sampled dose probabilities









            }

            //End MCMC






          }

          ///STOPHERE


          StoppedGroups.zeros();

          //Do we have a stopped group?
          for(k=0;k<(J+1);k++){
            eta2=0;

            for(j=0;j<DoseProb.n_rows;j++){
              //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
              if(DoseProb(j,k*nDose)>Target[k]){
                eta2++;
              }


            }



            if(eta2>(Upper[k]*DoseProb.n_rows)){
              StoppedGroups[k]=1;
            }

          }




          //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group
          for(k=0;k<(J+1);k++){
            if(nTreated1[k]<3){
              StoppedGroups[k]=0;
            }

            if(nTreated1[k]>2){
              GetGroup.zeros();
              GetGroup[k]=1;

              Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
              I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
              Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);

              if(GetFullyFollowed(Y2,I2,Doses2,Dose,T1)<3){
                StoppedGroups[k]=0;
              }

            }


          }






          //    Rf_PrintValue(wrap(StoppedGroups));
          if(sum(StoppedGroups)==(J+1)){
            //Let's try the groups separately to make sure we aren't stopping unneccessarily
            StoppedGroups.zeros();

            for(g=0;g<(J+1);g++){
              GetGroup.zeros();
              GetGroup[g]=1;

              Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
              I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
              Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);



              Which1 = g;


              if(Which1==0){
                NewMean = meanmu;
                NewSlope=meanslope;
              }else{
                NewMean = meanmu+MeanInts[Which1-1];
                NewSlope=meanslope+MeanSlopes[Which1-1];
              }





              //Inner MCMC loop, only compute neccessary quantities. Reset counters
              NumA=NumA.zeros()+2;
              NumB=NumB.zeros()+2;
              IntA=IntA.zeros()+1;
              IntB = IntB.zeros()+1;
              NumMu = 2;
              NumSlope=2;
              IntMu=1;
              IntSlope=1;
              avar=avar.zeros()+1;
              bvar=bvar.zeros()+1;
              muvar=1;
              slopevar=1;
              NumSig=2;
              IntSig=1;
              sigvar=1;
              a.zeros();
              b.zeros();
              MeanVec.zeros();
              GroupVec.zeros();
              mu=NewMean;
              slope=NewSlope;


              for(m=0;m<B;m++){

                if(m<(B/2 + 2)){
                  if(m%100==0){




                    if((IntMu/NumMu)>.5){
                      muvar=muvar*2;
                    }

                    if((IntMu/NumMu)<.2){
                      muvar=muvar/2;
                    }


                    if((IntSlope/NumSlope)>.5){
                      slopevar=slopevar*2;
                    }

                    if((IntSlope/NumSlope)<.2){
                      slopevar=slopevar/2;
                    }







                    NumMu = 2;
                    NumSlope=2;
                    IntMu=1;
                    IntSlope=1;
                    IntSig=1;
                    NumSig=2;




                  }
                }









                Munew = mu + as_scalar(arma::randn(1))*muvar;



                alpha= .5*pow((mu-NewMean),2)/varint -.5*pow((Munew-NewMean),2)/varint+   Like2(Y2, I2,  Doses2,  Munew,slope, sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope,  sig,Y2.n_rows);




                U=log(as_scalar(arma::randu(1)));

                if(U<alpha){
                  mu=Munew;
                  IntMu=IntMu+1;
                }

                NumMu=NumMu+1;




                //Sample the new slope and the new \beta_g coefficients jointly




                slopenew = slope + as_scalar(arma::randn(1))*slopevar;



                //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


                //Proposal ratio for beta

                alpha =.5*pow((slope-NewSlope),2)/varbeta -.5*pow((slopenew-NewSlope),2)/varbeta;


                alpha =   alpha+  Like2(Y2, I2,  Doses2,   mu,slopenew,sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope, sig,Y2.n_rows);


                U=log(as_scalar(arma::randu(1)));



                U=log(as_scalar(arma::randu(1)));

                if(U<alpha){
                  slope=slopenew;
                  IntSlope=IntSlope+1;
                }
                NumSlope=NumSlope+1;





                if(m>(B-B1-1)){
                  //Make OptDose and StoppedGroups
                  StoreInx=m-B1;



                  //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
                  //First nDose


                  k=Which1;

                  for(j=0;j<nDose;j++){
                    eta1=0;




                    eta1=mu+exp(slope)*Dose[j];



                    //For all groups add the intercept and slope*Groups

                    eta1=exp(eta1);

                    if(eta1>100000){
                      DoseProb(StoreInx, k*nDose+j) = 1;

                    }else{
                      DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                    }



                  }



                  //Now we have a matrix containing our sampled dose probabilities









                }

                //End MCMC
              }
              //Is our very last group stopped?

              k=Which1;
              eta2=0;

              for(j=0;j<DoseProb.n_rows;j++){
                //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
                if(DoseProb(j,k*nDose)>Target[k]){
                  eta2++;
                }


              }


              if(eta2>(Upper[Which1]*DoseProb.n_rows)){
                StoppedGroups[Which1]=1;
              }






            }





          }




          //Should we stop the trial?
          if(sum(StoppedGroups)==(J+1)){
            stopped=1;
            break;
          }









          if(StoppedGroups[Group]==0){
            //Enroll this patient
            DEC=1;


            if(sum(StoppedGroups)>0){
              //If we have a stopped group, REMOVE IT

              Y2=ReturnStoppedY(Y,I,Groups,StoppedGroups,i);
              I2=ReturnStoppedI(Y,I,Groups,StoppedGroups,i);
              Groups2=ReturnStoppedGroups(Y,I,Groups,StoppedGroups,i);
              Doses2=ReturnStoppedY(Doses,I,Groups,StoppedGroups,i);




              if(sum(StoppedGroups)==J){
                //We just have one group left

                for(k=0;k<StoppedGroups.n_rows;k++){
                  if(StoppedGroups[k]==0){
                    Which1=k;
                  }
                }


                if(Which1==0){
                  NewMean = meanmu;
                  NewSlope=meanslope;
                }else{
                  NewMean = meanmu+MeanInts[Which1-1];
                  NewSlope=meanslope+MeanSlopes[Which1-1];
                }


                //Inner MCMC loop, only compute neccessary quantities. Reset counters
                NumA=NumA.zeros()+2;
                NumB=NumB.zeros()+2;
                IntA=IntA.zeros()+1;
                IntB = IntB.zeros()+1;
                NumMu = 2;
                NumSlope=2;
                IntMu=1;
                IntSlope=1;
                avar=avar.zeros()+1;
                bvar=bvar.zeros()+1;
                muvar=1;
                slopevar=1;
                NumSig=2;
                IntSig=1;
                sigvar=1;
                a.zeros();
                b.zeros();
                MeanVec.zeros();
                GroupVec.zeros();
                mu=NewMean;
                slope=NewSlope;


                for(m=0;m<B;m++){

                  if(m<(B/2 + 2)){
                    if(m%100==0){




                      if((IntMu/NumMu)>.5){
                        muvar=muvar*2;
                      }

                      if((IntMu/NumMu)<.2){
                        muvar=muvar/2;
                      }


                      if((IntSlope/NumSlope)>.5){
                        slopevar=slopevar*2;
                      }

                      if((IntSlope/NumSlope)<.2){
                        slopevar=slopevar/2;
                      }







                      NumMu = 2;
                      NumSlope=2;
                      IntMu=1;
                      IntSlope=1;
                      IntSig=1;
                      NumSig=2;




                    }
                  }









                  Munew = mu + as_scalar(arma::randn(1))*muvar;



                  alpha= .5*pow((mu-NewMean),2)/varint -.5*pow((Munew-NewMean),2)/varint+   Like2(Y2, I2,  Doses2,  Munew,slope, sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope,  sig,Y2.n_rows);




                  U=log(as_scalar(arma::randu(1)));

                  if(U<alpha){
                    mu=Munew;
                    IntMu=IntMu+1;
                  }

                  NumMu=NumMu+1;




                  //Sample the new slope and the new \beta_g coefficients jointly




                  slopenew = slope + as_scalar(arma::randn(1))*slopevar;



                  //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


                  //Proposal ratio for beta

                  alpha =.5*pow((slope-NewSlope),2)/varbeta -.5*pow((slopenew-NewSlope),2)/varbeta;


                  alpha =   alpha+  Like2(Y2, I2,  Doses2,   mu,slopenew,sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope, sig,Y2.n_rows);


                  U=log(as_scalar(arma::randu(1)));



                  U=log(as_scalar(arma::randu(1)));

                  if(U<alpha){
                    slope=slopenew;
                    IntSlope=IntSlope+1;
                  }
                  NumSlope=NumSlope+1;





                  if(m>(B-B1-1)){
                    //Make OptDose and StoppedGroups
                    StoreInx=m-B1;



                    //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
                    //First nDose


                    k=Which1;

                    for(j=0;j<nDose;j++){
                      eta1=0;




                      eta1=mu+exp(slope)*Dose[j];



                      //For all groups add the intercept and slope*Groups

                      eta1=exp(eta1);

                      if(eta1>100000){
                        DoseProb(StoreInx, k*nDose+j) = 1;

                      }else{
                        DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                      }



                    }



                    //Now we have a matrix containing our sampled dose probabilities









                  }

                  //End MCMC
                }
                //Is our very last group stopped?

                k=Which1;
                eta2=0;

                for(j=0;j<DoseProb.n_rows;j++){
                  //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
                  if(DoseProb(j,k*nDose)>Target[k]){
                    eta2++;
                  }


                }


                if(eta2>(Upper[Which1]*DoseProb.n_rows)){
                  StoppedGroups[Which1]=1;
                }






              }




              //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group
              for(k=0;k<(J+1);k++){
                if(nTreated1[k]<3){
                  StoppedGroups[k]=0;
                }

                if(nTreated1[k]>2){
                  GetGroup.zeros();
                  GetGroup[k]=1;

                  Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
                  I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
                  Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);

                  if(GetFullyFollowed(Y2,I2,Doses2,Dose,T1)<3){
                    StoppedGroups[k]=0;
                  }

                }


              }





              if(sum(StoppedGroups)==(J+1)){
                stopped=1;

                Rf_PrintValue(wrap(StoppedGroups));
                Rprintf("It Stopped");

                Rf_PrintValue(wrap(Y2));
                Rf_PrintValue(wrap(I2));
                Rf_PrintValue(wrap(Doses2));


                break;
              }




            }





          }else{

            trialtime=trialtime+R::rexp(1/Accrue);
            Group=Sample2(groupprob);

          }

          //Do We have ANOTHER STOPPED GROUP?
          //Should we stop the trial?







          if(sum(StoppedGroups)==(J+1)){
            stopped=1;

            Rf_PrintValue(wrap(StoppedGroups));
            Rprintf("It Stopped");

            Rf_PrintValue(wrap(Y2));
            Rf_PrintValue(wrap(I2));
            Rf_PrintValue(wrap(Doses2));


            break;
          }


        }

        //Ok now let's enroll this patient, Whats the Optimal Dose

        //If Stopped==1, Get out of the trial
        if(stopped==1){
          break;
        }



        if(MEANS==0){
          for(k=0;k<(J+1);k++){
            for(j=0;j<nDose;j++){
              MeanVec(k,j)=sum(DoseProb.col(k*nDose+j))/DoseProb.n_rows;

            }

          }



        }


        //Determine Optimal Dose
        k=Group;
        eta1=MinVec(abs(MeanVec.row(k).t()-Target(k)));
        j=0;




        GroupVec = abs(MeanVec.row(k).t()-Target(k));

        while(GroupVec[j]!=eta1 ){
          j++;
          //Loop Proceeds until the minimum difference is reached
        }

        OptDose[k]=j;




        //Now we have the subgroup specific optimal doses in the vector OptDose



        //Now let's see if the dose is bigger than the largest dose.
        if(OptDose[Group]>=nDose){
          OptDose[Group]=nDose-1;
        }



        //Is this dose bigger than the dose higher than what's been tried?
        if(DoseTried(Group,OptDose[Group])==0){
          j=0;

          while(DoseTried(Group,j)==1){
            j++;
          }

          OptDose[Group]=j;

        }






        //Assign Patient to subgroup, optdose, add accrual time.

        Groups[i]=Group;
        if(OptDose[Group]==0){
          //If we assign one at the lowest dose add this
          nTreated1[Group]=nTreated1[Group]+1;
        }
        nTreated[Group]=nTreated[Group]+1;

        //If this is the 3rd patient treated in this subgroup with the lowest dose, suspend accrual for TGroup
        if(nTreated[Group]==3){
          SuspendGroups[Group]=1;
          GLast[Group]=trialtime;
        }



        Doses[i]=Dose[OptDose[Group]];


        ACC[i]=trialtime;
        //Based on the Group and Optimal dose, randomly generate the next patient Data


        if(Family==0){
          Times(i)=R::rexp(Param1(Group,OptDose[Group]));
        }



        if(Family==1){
          Times[i]=R::rgamma(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
        }

        if(Family==2){
          Times[i]=R::rlnorm(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
        }


        if(Family==3){
          if(as_scalar(arma::randu(1))<Param1(Group,OptDose[Group])){
            Times[i]=as_scalar(arma::randu(1))*T1;
          }else{
            Times[i]=T1+1;
          }
        }

        if(Family==4){
          Times[i]=R::rweibull(Param1(Group,OptDose[Group]),Param2(Group,OptDose[Group]));
        }


        DoseTried(Group,OptDose(Group))=1;












      }




      //End of the i loop, now we have ALL the patients
    }



    //Now let's make the final decision.
    if(stopped==1){
      //The trial Ended Early, i-1 is the index of the last patient enrolled.


      for(k=0;k<i;k++){
        for(j=0;j<(J+1);j++){
          if(Groups[k]==j){
            NTox[j]=NTox[j]+I[k];
          }
        }
      }


      //All Optimal Doses are -1
      for(k=0;k<(J+1);k++){
        OptDose[k]=-1;
      }


    }else{
      //The trial did NOT end early. Follow patients out until end of trial and get the optimal doses
      trialtime=trialtime+T1;



      Y.zeros();
      I.zeros();

      //Setup Y
      for(k=0;k<Nmax;k++){

        Y[k]=min1(Times[k],T1);


        if(Times[k]<T1){
          I[k]=1;
        }
      }





      //Inner MCMC loop, only compute neccessary quantities. Reset counters
      NumA=NumA.zeros()+2;
      NumB=NumB.zeros()+2;
      IntA=IntA.zeros()+1;
      IntB = IntB.zeros()+1;
      NumMu = 2;
      NumSlope=2;
      IntMu=1;
      IntSlope=1;
      avar=avar.zeros()+1;
      bvar=bvar.zeros()+1;
      muvar=1;
      slopevar=1;
      NumSig=2;
      IntSig=1;
      sigvar=1;
      a.zeros();
      b.zeros();
      MeanVec.zeros();
      GroupVec.zeros();
      mu=meanmu;
      slope=meanslope;
      INVEC.zeros();
      INVECNEW.zeros();
      INVEC = INVEC+1;
      a=MeanInts;
      b=MeanSlopes;


      for(m=0;m<GroupMem.n_rows;m++){
        GroupMem[m]=m+1;
      }

      GroupMemProp=GroupMem;

      Rf_PrintValue(wrap(INVEC));
      Rf_PrintValue(wrap(GroupMem));

      for(m=0;m<B;m++){

        if(m<(B/2 + 2)){
          if(m%100==0){


            for(k=0;k<(J+1);k++){
              if((IntA[k]/NumA[k])>.5){
                avar[k]=avar[k]*2;
              }

              if((IntA[k]/NumA[k])<.2){
                avar[k]=avar[k]/2;
              }








              if((IntB[k]/NumB[k])>.5){
                bvar[k]=bvar[k]*2;
              }

              if((IntB[k]/NumB[k])<.2){
                bvar[k]=bvar[k]/2;
              }



            }


            if((IntMu/NumMu)>.5){
              muvar=muvar*2;
            }

            if((IntMu/NumMu)<.2){
              muvar=muvar/2;
            }


            if((IntSlope/NumSlope)>.5){
              slopevar=slopevar*2;
            }

            if((IntSlope/NumSlope)<.2){
              slopevar=slopevar/2;
            }





            NumA=NumA.zeros()+2;
            NumB=NumB.zeros()+2;
            IntA=IntA.zeros()+1;
            IntB=IntB.zeros()+1;

            NumMu = 2;
            NumSlope=2;
            IntMu=1;
            IntSlope=1;
            IntSig=1;
            NumSig=2;




          }
        }



        for(k=0;k<aprop.n_rows;k++){

          if(INVEC[k]>0){

            aprop=a;
            aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);






            alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


            U=log(as_scalar(arma::randu(1)));


            if(U<alpha){
              a(k)=aprop(k);
              IntA[k]=IntA[k]+1;
            }

            NumA=NumA+1;


          }
        }
        //Now do a MH step for the slope

        for(k=0;k<aprop.size();k++){

          if(INVEC[k]>0){
            bprop=b;
            bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];



            alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1  +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


            U=log(as_scalar(arma::randu(1)));

            z9[0]=alpha;

            //  Rf_PrintValue(wrap(z9));

            if(U<alpha){
              b[k]=bprop[k];
              IntB[k]=IntB[k]+1;
            }

          }

        }


        NumB=NumB+1;



        //Spike And Slab


        U=as_scalar(arma::randu(1));


        if((3*U)<1){
          //Add All or delete all

          if(sum(INVEC)==(J)){
            //Auto Delete Move




            aprop.zeros();
            bprop.zeros();

            GroupMemProp=GroupMem;
            //Now Randomly Draw Our New Group Assignment


            INVECNEW.zeros();




            GroupMemProp.zeros();




            alpha=0;
            for(k=0;k<bprop.n_rows;k++){
              alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

            }





            //Density when 0 is automatically 1 so no proposal needed
            alpha =  alpha + LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



            U=log(as_scalar(arma::randu(1)));



            if(U<alpha){
              b.zeros();
              a.zeros();
              GroupMem.zeros();
              INVEC.zeros();

              //If Any other GroupMEM are in this group, we need to delete them.



            }






          }else{


            if(sum(INVEC)==0){
              //Autho Add Move

              //ADD Move

              for(k=0;k<aprop.n_rows;k++){
                aprop[k]=as_scalar(arma::randn(1));
                bprop[k]=as_scalar(arma::randn(1))/10;
                GroupMemProp[k]=k+1;
              }



              //Density when 0 is automatically 1 so no proposal needed

              alpha =  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



              for(k=0;k<aprop.n_rows;k++){


                alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]) -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 ;

              }




              alpha=alpha + log(sum(INVEC)+1);

              U=log(as_scalar(arma::randu(1)));

              z9[0]=alpha;

              //  Rf_PrintValue(wrap(z9));

              if(U<alpha){

                for(k=0;k<aprop.n_rows;k++){
                  b[k]=bprop[k];
                  a[k]=aprop[k];
                  INVEC[k]=1;
                  GroupMem[k]=k+1;
                }


              }




            }else{
              //Randomly Decide to add or Delete

              U=as_scalar(arma::randu(1));


              if(U<.5){
                //Add


                //ADD Move

                for(k=0;k<aprop.n_rows;k++){
                  aprop[k]=as_scalar(arma::randn(1));
                  bprop[k]=as_scalar(arma::randn(1))/10;
                  GroupMemProp[k]=k+1;
                }



                //Density when 0 is automatically 1 so no proposal needed

                alpha =  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



                for(k=0;k<aprop.n_rows;k++){


                  alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]) -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 ;

                }




                alpha=alpha + log(sum(INVEC)+1);

                U=log(as_scalar(arma::randu(1)));

                z9[0]=alpha;

                //  Rf_PrintValue(wrap(z9));

                if(U<alpha){

                  for(k=0;k<aprop.n_rows;k++){
                    b[k]=bprop[k];
                    a[k]=aprop[k];
                    INVEC[k]=1;
                    GroupMem[k]=k+1;
                  }


                }



              }else{
                //Delete All




                aprop.zeros();
                bprop.zeros();

                GroupMemProp=GroupMem;
                //Now Randomly Draw Our New Group Assignment


                INVECNEW.zeros();




                GroupMemProp.zeros();




                alpha=0;
                for(k=0;k<bprop.n_rows;k++){
                  alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

                }





                //Density when 0 is automatically 1 so no proposal needed
                alpha =  alpha + LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



                U=log(as_scalar(arma::randu(1)));



                if(U<alpha){
                  b.zeros();
                  a.zeros();
                  GroupMem.zeros();
                  INVEC.zeros();

                  //If Any other GroupMEM are in this group, we need to delete them.



                }




              }





            }




          }









        }else{



          //Should we swap?
          if((sum(INVEC)==J) || (sum(INVEC)==0)){
            //Add/Delete



            k=Sample1(J);

            if(INVEC[k]==0){
              //ADD Move
              aprop=a;
              bprop=b;
              aprop[k]=as_scalar(arma::randn(1));
              bprop[k]=as_scalar(arma::randn(1))/10;

              GroupMemProp=GroupMem;
              GroupMemProp[k]=k+1;


              //Density when 0 is automatically 1 so no proposal needed
              alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

              alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);


              alpha=alpha + log(sum(INVEC)+1);

              U=log(as_scalar(arma::randu(1)));

              z9[0]=alpha;

              //  Rf_PrintValue(wrap(z9));

              if(U<alpha){
                b[k]=bprop[k];
                a[k]=aprop[k];
                INVEC[k]=1;
                GroupMem[k]=k+1;
              }


            }else{
              //Delete Move
              aprop=a;
              bprop=b;
              //Delete Move
              aprop[k]=0;
              bprop[k]=0;

              GroupMemProp=GroupMem;
              //Now Randomly Draw Our New Group Assignment


              INVECNEW=INVEC;
              INVECNEW[k]=0;



              GroupMemProp[k]=GetNewGroup(INVECNEW);










              //Density when 0 is automatically 1 so no proposal needed
              alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

              alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

              alpha=alpha-log(sum(INVEC)+1);

              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b[k]=0;
                a[k]=0;
                GroupMem=GroupMemProp;
                INVEC[k]=0;

                //If Any other GroupMEM are in this group, we need to delete them.

                for(j=0;j<GroupMem.n_rows;j++){
                  if(GroupMem[j]==(k+1)){
                    GroupMem[j]=GroupMem[k];
                  }
                }



              }


            }


          }else{




            if((3*U)>2){
              //Swap Move
              //Randomly Pick One In and One Out

              IntIN = GetIn(INVEC);
              IntOUT= GetOut(INVEC);

              //Now IntIN and IntOUT contain the entries that are currently in or out





            }else{
              //Add/Delete



              k=Sample1(J);

              if(INVEC[k]==0){
                //ADD Move
                aprop=a;
                bprop=b;
                aprop[k]=as_scalar(arma::randn(1));
                bprop[k]=as_scalar(arma::randn(1))/10;

                GroupMemProp=GroupMem;
                GroupMemProp[k]=k+1;


                //Density when 0 is automatically 1 so no proposal needed
                alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

                alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);


                alpha=alpha + log(sum(INVEC)+1);

                U=log(as_scalar(arma::randu(1)));

                z9[0]=alpha;

                //  Rf_PrintValue(wrap(z9));

                if(U<alpha){
                  b[k]=bprop[k];
                  a[k]=aprop[k];
                  INVEC[k]=1;
                  GroupMem[k]=k+1;
                }


              }else{
                //Delete Move
                aprop=a;
                bprop=b;
                //Delete Move
                aprop[k]=0;
                bprop[k]=0;

                GroupMemProp=GroupMem;
                //Now Randomly Draw Our New Group Assignment


                INVECNEW=INVEC;
                INVECNEW[k]=0;



                GroupMemProp[k]=GetNewGroup(INVECNEW);










                //Density when 0 is automatically 1 so no proposal needed
                alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

                alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

                alpha=alpha-log(sum(INVEC)+1);

                U=log(as_scalar(arma::randu(1)));



                if(U<alpha){
                  b[k]=0;
                  a[k]=0;
                  GroupMem=GroupMemProp;
                  INVEC[k]=0;

                  //If Any other GroupMEM are in this group, we need to delete them.

                  for(j=0;j<GroupMem.n_rows;j++){
                    if(GroupMem[j]==(k+1)){
                      GroupMem[j]=GroupMem[k];
                    }
                  }



                }


              }


            }



          }

        }






        //Swap Group

        for(k=0;k<J;k++){
          if(INVEC[k]==0){

            GroupMemProp=GroupMem;
            //Now Randomly Draw Our New Group Assignment


            INVECNEW=INVEC;
            INVECNEW[k]=0;



            GroupMemProp[k]=GetNewGroup(INVECNEW);










            //Density when 0 is automatically 1 so no proposal needed
            alpha = LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



            U=log(as_scalar(arma::randu(1)));



            if(U<alpha){
              GroupMem[k]=GroupMemProp[k];
            }
          }

        }






        Munew = mu + as_scalar(arma::randn(1))*muvar;



        alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint+   LikeMULTI(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);




        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          mu=Munew;
          IntMu=IntMu+1;
        }

        NumMu=NumMu+1;




        //Sample the new slope and the new \beta_g coefficients jointly




        slopenew = slope + as_scalar(arma::randn(1))*slopevar;



        //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


        //Proposal ratio for beta

        alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


        alpha =   alpha+  LikeMULTI(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


        U=log(as_scalar(arma::randu(1)));



        U=log(as_scalar(arma::randu(1)));

        if(U<alpha){
          slope=slopenew;
          IntSlope=IntSlope+1;
        }
        NumSlope=NumSlope+1;


        if(m>(B-B1-1)){
          //Make OptDose and StoppedGroups
          StoreInx=m-B1;

          for(j=0;j<J;j++){;
            astore(StoreInx,j)=INVEC(j);
            bstore(StoreInx,j)=GroupMem[j];
          }

          mustore[StoreInx]=mu;
          slopestore[StoreInx]=slope;
          sigstore[StoreInx]=sig;


          //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
          //First nDose

          aprop.zeros();
          bprop.zeros();
          for(j=0;j<J;j++){
            for(k=0;k<J;k++){
              if(GroupMem[j]==(k+1)){
                aprop[j]=a[k];
                bprop[j]=b[k];
              }
            }
          }


          for(k=0;k<(J+1);k++){

            for(j=0;j<nDose;j++){
              eta1=0;

              if(k>0){
                eta1 = aprop[k-1]+exp(bprop[k-1]+slope)*Dose[j]+mu;
              }

              if(k==0){
                eta1=mu+exp(slope)*Dose[j];

              }

              //For all groups add the intercept and slope*Groups

              eta1=exp(eta1);

              if(eta1>100000){
                DoseProb(StoreInx, k*nDose+j) = 1;

              }else{
                DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

              }



            }


          }
          //Now we have a matrix containing our sampled dose probabilities









        }

        //End MCMC






      }

      Rf_PrintValue(wrap(INVEC));
      Rf_PrintValue(wrap(GroupMem));
      //Get Optimal Dose


      StoppedGroups.zeros();
      //Get MeanVec
      //reset to means




      //Now we have the mean vector



      //Do we have a stopped group?
      for(k=0;k<(J+1);k++){
        eta2=0;

        for(j=0;j<DoseProb.n_rows;j++){
          //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
          if(DoseProb(j,k*nDose)>Target[k]){
            eta2++;
          }


        }







        if(eta2>(Upper[k]*DoseProb.n_rows)){
          StoppedGroups[k]=1;
        }

      }






      if(sum(StoppedGroups)>0){
        //If we have a stopped group, REMOVE IT

        Y2=ReturnStoppedY(Y,I,Groups,StoppedGroups,i);
        I2=ReturnStoppedI(Y,I,Groups,StoppedGroups,i);
        Groups2=ReturnStoppedGroups(Y,I,Groups,StoppedGroups,i);
        Doses2=ReturnStoppedY(Doses,I,Groups,StoppedGroups,i);




        if(sum(StoppedGroups)==J){
          //We just have one group left

          for(k=0;k<StoppedGroups.n_rows;k++){
            if(StoppedGroups[k]==0){
              Which1=k;
            }
          }


          if(Which1==0){
            NewMean = meanmu;
            NewSlope=meanslope;
          }else{
            NewMean = meanmu+MeanInts[Which1-1];
            NewSlope=meanslope+MeanSlopes[Which1-1];
          }

          z9[0]=NewMean;
          z9[1]=NewSlope;
          Rf_PrintValue(z9);

          //Inner MCMC loop, only compute neccessary quantities. Reset counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB = IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          avar=avar.zeros()+1;
          bvar=bvar.zeros()+1;
          muvar=1;
          slopevar=1;
          NumSig=2;
          IntSig=1;
          sigvar=1;
          a.zeros();
          b.zeros();
          MeanVec.zeros();
          GroupVec.zeros();
          mu=NewMean;
          slope=NewSlope;


          for(m=0;m<B;m++){

            if(m<(B/2 + 2)){
              if(m%100==0){




                if((IntMu/NumMu)>.5){
                  muvar=muvar*2;
                }

                if((IntMu/NumMu)<.2){
                  muvar=muvar/2;
                }


                if((IntSlope/NumSlope)>.5){
                  slopevar=slopevar*2;
                }

                if((IntSlope/NumSlope)<.2){
                  slopevar=slopevar/2;
                }







                NumMu = 2;
                NumSlope=2;
                IntMu=1;
                IntSlope=1;
                IntSig=1;
                NumSig=2;




              }
            }









            Munew = mu + as_scalar(arma::randn(1))*muvar;



            alpha= .5*pow((mu-NewMean),2)/varint -.5*pow((Munew-NewMean),2)/varint+   Like2(Y2, I2,  Doses2,  Munew,slope, sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope,  sig,Y2.n_rows);




            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              mu=Munew;
              IntMu=IntMu+1;
            }

            NumMu=NumMu+1;




            //Sample the new slope and the new \beta_g coefficients jointly




            slopenew = slope + as_scalar(arma::randn(1))*slopevar;



            //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


            //Proposal ratio for beta

            alpha =.5*pow((slope-NewSlope),2)/varbeta -.5*pow((slopenew-NewSlope),2)/varbeta;


            alpha =   alpha+  Like2(Y2, I2,  Doses2,   mu,slopenew,sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope, sig,Y2.n_rows);


            U=log(as_scalar(arma::randu(1)));



            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              slope=slopenew;
              IntSlope=IntSlope+1;
            }
            NumSlope=NumSlope+1;





            if(m>(B-B1-1)){
              //Make OptDose and StoppedGroups
              StoreInx=m-B1;

              for(j=0;j<J;j++){;
                astore(StoreInx,j)=INVEC(j);
                bstore(StoreInx,j)=GroupMem[j];
              }

              mustore(StoreInx)=mu;
              slopestore(StoreInx)=slope;

              //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
              //First nDose


              k=Which1;

              for(j=0;j<nDose;j++){
                eta1=0;




                eta1=mu+exp(slope)*Dose[j];



                //For all groups add the intercept and slope*Groups

                eta1=exp(eta1);

                if(eta1>100000){
                  DoseProb(StoreInx, k*nDose+j) = 1;

                }else{
                  DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                }



              }



              //Now we have a matrix containing our sampled dose probabilities









            }

            //End MCMC
          }
          //Is our very last group stopped?

          k=Which1;
          eta2=0;

          for(j=0;j<DoseProb.n_rows;j++){
            //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
            if(DoseProb(j,k*nDose)>Target[k]){
              eta2++;
            }


          }


          if(eta2>(Upper[Which1]*DoseProb.n_rows)){
            StoppedGroups[Which1]=1;
          }





        }

      }




      //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group

      //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group
      for(k=0;k<(J+1);k++){
        if(nTreated1[k]<3){
          StoppedGroups[k]=0;
        }

        if(nTreated1[k]>2){
          GetGroup.zeros();
          GetGroup[k]=1;

          Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
          I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
          Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);

          if(GetFullyFollowed(Y2,I2,Doses2,Dose,T1)<3){
            StoppedGroups[k]=0;
          }

        }


      }




      //If all groups are stopped, let's do separate trials JUST TO MAKE SURE ONE GROUP IS NOT DOMINATING THE DECISION TO STOP!!


      //    Rf_PrintValue(wrap(StoppedGroups));
      if(sum(StoppedGroups)==(J+1)){
        //Let's try the groups separately to make sure we aren't stopping unneccessarily
        StoppedGroups.zeros();

        for(g=0;g<(J+1);g++){
          GetGroup.zeros();
          GetGroup[g]=1;

          Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
          I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
          Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);



          Which1 = g;


          if(Which1==0){
            NewMean = meanmu;
            NewSlope=meanslope;
          }else{
            NewMean = meanmu+MeanInts[Which1-1];
            NewSlope=meanslope+MeanSlopes[Which1-1];
          }





          //Inner MCMC loop, only compute neccessary quantities. Reset counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB = IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          avar=avar.zeros()+1;
          bvar=bvar.zeros()+1;
          muvar=1;
          slopevar=1;
          NumSig=2;
          IntSig=1;
          sigvar=1;
          a.zeros();
          b.zeros();
          MeanVec.zeros();
          GroupVec.zeros();
          mu=NewMean;
          slope=NewSlope;


          for(m=0;m<B;m++){

            if(m<(B/2 + 2)){
              if(m%100==0){




                if((IntMu/NumMu)>.5){
                  muvar=muvar*2;
                }

                if((IntMu/NumMu)<.2){
                  muvar=muvar/2;
                }


                if((IntSlope/NumSlope)>.5){
                  slopevar=slopevar*2;
                }

                if((IntSlope/NumSlope)<.2){
                  slopevar=slopevar/2;
                }







                NumMu = 2;
                NumSlope=2;
                IntMu=1;
                IntSlope=1;
                IntSig=1;
                NumSig=2;




              }
            }









            Munew = mu + as_scalar(arma::randn(1))*muvar;



            alpha= .5*pow((mu-NewMean),2)/varint -.5*pow((Munew-NewMean),2)/varint+   Like2(Y2, I2,  Doses2,  Munew,slope, sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope,  sig,Y2.n_rows);




            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              mu=Munew;
              IntMu=IntMu+1;
            }

            NumMu=NumMu+1;




            //Sample the new slope and the new \beta_g coefficients jointly




            slopenew = slope + as_scalar(arma::randn(1))*slopevar;



            //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


            //Proposal ratio for beta

            alpha =.5*pow((slope-NewSlope),2)/varbeta -.5*pow((slopenew-NewSlope),2)/varbeta;


            alpha =   alpha+  Like2(Y2, I2,  Doses2,   mu,slopenew,sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope, sig,Y2.n_rows);


            U=log(as_scalar(arma::randu(1)));



            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              slope=slopenew;
              IntSlope=IntSlope+1;
            }
            NumSlope=NumSlope+1;





            if(m>(B-B1-1)){
              //Make OptDose and StoppedGroups
              StoreInx=m-B1;



              //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
              //First nDose


              k=Which1;

              for(j=0;j<nDose;j++){
                eta1=0;




                eta1=mu+exp(slope)*Dose[j];



                //For all groups add the intercept and slope*Groups

                eta1=exp(eta1);

                if(eta1>100000){
                  DoseProb(StoreInx, k*nDose+j) = 1;

                }else{
                  DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                }



              }



              //Now we have a matrix containing our sampled dose probabilities









            }

            //End MCMC
          }
          //Is our very last group stopped?

          k=Which1;
          eta2=0;

          for(j=0;j<DoseProb.n_rows;j++){
            //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
            if(DoseProb(j,k*nDose)>Target[k]){
              eta2++;
            }


          }


          if(eta2>(Upper[Which1]*DoseProb.n_rows)){
            StoppedGroups[Which1]=1;
          }






        }





      }




      if(MEANS==0){
        for(k=0;k<(J+1);k++){
          for(j=0;j<nDose;j++){
            MeanVec(k,j)=sum(DoseProb.col(k*nDose+j))/DoseProb.n_rows;

          }

        }


      }

      Rf_PrintValue(wrap(MeanVec));




      //Determine Optimal Doses
      for(k=0;k<(J+1);k++){

        if(StoppedGroups(k)==1){
          OptDose(k)=-1;
        }else{
          eta1=MinVec(abs(MeanVec.row(k).t()-Target[k]));
          j=0;

          GroupVec = abs(MeanVec.row(k).t()-Target[k]);

          //  Rf_PrintValue(wrap(MeanVec));


          while(GroupVec[j]!=eta1 ){
            j++;
            //Loop Proceeds until the minimum difference is reached
          }

          OptDose[k]=j;
        }


      }




      //Is this dose bigger than the dose higher than what's been tried?
      for(k=0;k<(J+1);k++){
        if(StoppedGroups[k]==0){
          for(j=0;j<(J+1);j++){
            if(OptDose[j]>=nDose){
              OptDose[j]=nDose-1;
            }
          }




        }



      }

      //Is this dose bigger than the dose higher than what's been tried?
      for(k=0;k<(J+1);k++){
        if(StoppedGroups(k)==0){
          if(DoseTried(k,OptDose[k])==0){
            j=0;

            while(DoseTried(k,j)==1){
              j++;
            }

            OptDose[k]=j;

          }
        }
      }


    }




    //Num Toxicities
    for(k=0;k<Nmax;k++){
      for(j=0;j<(J+1);j++){
        if(Groups[k]==j){
          NTox[j]=NTox[j]+I[k];
        }
      }
    }



    OptDose=OptDose+1;


    TrialTimes(rep)=trialtime;
    for(j=0;j<(J+1);j++){
      OptimalDoses(rep,j)=OptDose(j);
      NTOX(rep,j)=NTox(j);
    }


















    //Fill in Doses and Groups Treated
    for(j=0;j<Nmax;j++){
      DoseStore(rep,j)=Doses(j);
      GroupStore(rep,j)=Groups(j);
    }



  }



}




  List z1 = List::create(OptimalDoses,NTOX,TrialTimes,DoseStore,GroupStore,Y2,Doses2,mustore,slopestore,astore,bstore,Y,I,Doses,Groups);





  return(z1);




}






//[[Rcpp::export]]
arma::vec GetDose1( arma::vec Y, //Vector of Times
                    arma::vec I, //Vector of Toxicity Indicators
                    arma::vec Doses, //Standardized Dose Values given to each patient
                    arma::vec Groups, //SubGroup Membership
                    arma::mat DoseTried, //Matrix of Tried doses for each simulation
                    double T1, //Reference Time
                    arma::vec Target, //Target Toxicity Probability
                    arma::vec Upper, //Vector of Upper Cutoffs
                    arma::vec Dose, //Standardized values of Doses considered
                    double meanmu, //This is the prior mean of the baseline int
                    double meanslope, //Prior mean of baseline slope
                    arma::vec MeanInts, //Prior Means of Intercept Parameters
                    arma::vec MeanSlopes, //Prior Means of Slope Parameters
                    double varint, //Prior Variance of Intercept Parameters
                    double varbeta //Prior Variance of Slope Parameters
){


  int g=0;

  int StoreInx=0;

  //Change this to be an input
  arma::vec PROBIN(MeanInts.n_rows);
  for(g=0;g<PROBIN.n_rows;g++){
    PROBIN(g)=.9;
  }

  int MEANS=0;

  int IntIN=0;
  int IntOUT=0;

  //Group Slope Prior Var
  double varbeta1=varbeta;
  //Group Int Prior Var
  double varint1=varint;



  int Group=0;
  //Important For loop integer
  int m =0; //For the inner MCMC looprep
  int i=0; //For the patient index
  int rep=0; //For the simulation repetition.


  //First Fill out Target Vector
  int nDose=Dose.n_rows;


  //Important Integer Quantities
  int B=2000; //Number of iterations for MCMC


  int B1=B/2;

  //Make List objects we'll use





  //Innitialize parameters for MCMC
  double mu=meanmu;
  double slope=meanslope;
  double sig=T1;

  double NewMean;
  double NewSlope;
  int Which1;



  //Important quantities we will need in the MCMC
  int  J=MeanInts.n_rows;
  double alpha=0;
  double U=0;
  double signew=0;
  double slopenew=0;
  double Munew=0;

  //For loop integers
  int j=0;
  int k=0;


  //Vector Of group intercepts
  arma::vec a(J);
  //Vector of group slopes
  arma::vec b=a;
  //Vectors for Group MCMC



  //Vectors for proposals
  arma::vec aprop = a;
  arma::vec bprop=b;





  arma::vec NumA(J);
  arma::vec NumB(J);
  arma::vec IntA(J);
  arma::vec IntB(J);

  double NumMu = 2;
  double NumSlope=2;
  double IntMu=1;
  double IntSlope=1;

  arma::vec avar(J);
  arma::vec bvar(J);


  double muvar=1;
  double slopevar=1;
  double trialtime=0;
  NumericVector z9(2);
  NumericVector zprop(5);

  double NumSig=2;
  double IntSig=1;
  double sigvar=1;

  arma::vec etaG(J+1);
  arma::mat DoseProb(B1,(J+1)*nDose);
  arma::mat MeanDose(J+1,nDose);
  arma::vec sigstore(B1);
  arma::vec mustore(B1);
  arma::vec slopestore(B1);
  arma::mat astore(B1,J);
  arma::mat bstore=astore;
  double eta1=0;
  double DEC=0;

  arma::vec GroupMem(J);
  arma::vec GroupMemProp=GroupMem;

  arma::vec INVEC(J);
  INVEC.zeros();

  int eta2=0;

  arma::mat MeanVec((J+1),nDose);
  arma::vec NTox(J+1);
  arma::vec OptDose(J+1);
  arma::vec StoppedGroups(J+1);
  arma::vec SuspendGroups(J+1);
  arma::vec GLast(J+1);
  arma::vec nTreated(J+1);
  double stopped=0;

  arma::vec GroupVec(nDose);
  arma::vec TriedGroups(J+1);




  arma::vec nTreated1(J+1);
  nTreated1.zeros();



  arma::vec INVECNEW=INVEC;
  arma::vec Y2;
  arma::vec I2;
  arma::vec Groups2;
  arma::vec Doses2;
  arma::vec GetGroup=StoppedGroups;


  if(MeanInts.n_rows==1){
    //Two Subgroups







    //Inner MCMC loop, only compute neccessary quantities. Reset counters
    NumA=NumA.zeros()+2;
    NumB=NumB.zeros()+2;
    IntA=IntA.zeros()+1;
    IntB = IntB.zeros()+1;
    NumMu = 2;
    NumSlope=2;
    IntMu=1;
    IntSlope=1;
    avar=avar.zeros()+1;
    bvar=bvar.zeros()+1;
    muvar=1;
    slopevar=1;
    NumSig=2;
    IntSig=1;
    sigvar=1;
    a.zeros();
    b.zeros();
    MeanVec.zeros();
    GroupVec.zeros();
    mu=meanmu;
    slope=meanslope;
    INVEC.zeros();
    INVEC = INVEC+1;
    a=MeanInts;
    b=MeanSlopes;

    for(m=0;m<B;m++){

      if(m<(B/2 + 2)){
        if(m%100==0){


          for(k=0;k<(J+1);k++){
            if((IntA[k]/NumA[k])>.5){
              avar[k]=avar[k]*2;
            }

            if((IntA[k]/NumA[k])<.2){
              avar[k]=avar[k]/2;
            }








            if((IntB[k]/NumB[k])>.5){
              bvar[k]=bvar[k]*2;
            }

            if((IntB[k]/NumB[k])<.2){
              bvar[k]=bvar[k]/2;
            }



          }


          if((IntMu/NumMu)>.5){
            muvar=muvar*2;
          }

          if((IntMu/NumMu)<.2){
            muvar=muvar/2;
          }


          if((IntSlope/NumSlope)>.5){
            slopevar=slopevar*2;
          }

          if((IntSlope/NumSlope)<.2){
            slopevar=slopevar/2;
          }





          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB=IntB.zeros()+1;

          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          IntSig=1;
          NumSig=2;




        }
      }



      for(k=0;k<aprop.n_rows;k++){

        if(INVEC[k]>0){

          aprop=a;
          aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);






          alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


          U=log(as_scalar(arma::randu(1)));


          if(U<alpha){
            a(k)=aprop(k);
            IntA[k]=IntA[k]+1;
          }

          NumA=NumA+1;


        }
      }
      //Now do a MH step for the slope

      for(k=0;k<aprop.size();k++){

        if(INVEC[k]>0){
          bprop=b;
          bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];



          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1  +  Like1(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


          U=log(as_scalar(arma::randu(1)));

          z9[0]=alpha;

          //  Rf_PrintValue(wrap(z9));

          if(U<alpha){
            b[k]=bprop[k];
            IntB[k]=IntB[k]+1;
          }

        }

      }


      NumB=NumB+1;



      //Spike And Slab

      for(k=0;k<J;k++){
        if(INVEC[k]==0){
          //ADD Move
          aprop[k]=as_scalar(arma::randn(1))/10;
          bprop[k]=as_scalar(arma::randn(1))/10;


          //Density when 0 is automatically 1 so no proposal needed
          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);

          alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);

          U=log(as_scalar(arma::randu(1)));

          z9[0]=alpha;

          //  Rf_PrintValue(wrap(z9));

          if(U<alpha){
            b[k]=bprop[k];
            a[k]=aprop[k];
            INVEC[k]=1;
          }


        }else{
          //Delete Move

          //Delete Move
          aprop[k]=0;
          bprop[k]=0;


          //Density when 0 is automatically 1 so no proposal needed
          alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);

          alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

          U=log(as_scalar(arma::randu(1)));



          if(U<alpha){
            b[k]=0;
            a[k]=0;
          }


        }







      }







      Munew = mu + as_scalar(arma::randn(1))*muvar;



      alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint+   Like1(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);




      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        mu=Munew;
        IntMu=IntMu+1;
      }

      NumMu=NumMu+1;




      //Sample the new slope and the new \beta_g coefficients jointly




      slopenew = slope + as_scalar(arma::randn(1))*slopevar;



      //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


      //Proposal ratio for beta

      alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


      alpha =   alpha+  Like1(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


      U=log(as_scalar(arma::randu(1)));



      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        slope=slopenew;
        IntSlope=IntSlope+1;
      }
      NumSlope=NumSlope+1;





      if(m>(B-B1-1)){
        //Make OptDose and StoppedGroups
        StoreInx=m-B1;

        for(j=0;j<J;j++){;
          astore(StoreInx,j)=a(j);
          bstore(StoreInx,j)=b[j];
        }

        mustore[StoreInx]=mu;
        slopestore[StoreInx]=slope;
        sigstore[StoreInx]=sig;


        //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
        //First nDose


        for(k=0;k<(J+1);k++){

          for(j=0;j<nDose;j++){
            eta1=0;

            if(k>0){
              eta1 = a[k-1]+exp(b[k-1]+slope)*Dose[j]+mu;
            }

            if(k==0){
              eta1=mu+exp(slope)*Dose[j];

            }

            //For all groups add the intercept and slope*Groups

            eta1=exp(eta1);

            if(eta1>100000){
              DoseProb(StoreInx, k*nDose+j) = 1;

            }else{
              DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

            }



          }


        }
        //Now we have a matrix containing our sampled dose probabilities









      }

      //End MCMC
    }



    //Get Optimal Dose


    StoppedGroups.zeros();
    //Get MeanVec
    //reset to means




    //Now we have the mean vector



    //Do we have a stopped group?
    for(k=0;k<(J+1);k++){
      eta2=0;

      for(j=0;j<DoseProb.n_rows;j++){
        //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
        if(DoseProb(j,k*nDose)>Target[k]){
          eta2++;
        }


      }







      if(eta2>(Upper[k]*DoseProb.n_rows)){
        StoppedGroups[k]=1;
      }

    }






    if(sum(StoppedGroups)>0){
      //If we have a Stopped Group, Double Check that 3 patients have been fully evaluated

      Y2=ReturnStoppedY(Y,I,Groups,StoppedGroups,i);
      I2=ReturnStoppedI(Y,I,Groups,StoppedGroups,i);
      Groups2=ReturnStoppedGroups(Y,I,Groups,StoppedGroups,i);
      Doses2=ReturnStoppedY(Doses,I,Groups,StoppedGroups,i);









      //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group

      for(k=0;k<(J+1);k++){
        if(nTreated1[k]<3){
          StoppedGroups[k]=0;
        }

        if(nTreated1[k]>2){
          GetGroup.zeros();
          GetGroup[k]=1;

          Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
          I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
          Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);

          if(GetFullyFollowed(Y2,I2,Doses2,Dose,T1)<3){
            StoppedGroups[k]=0;
          }

        }


      }




      //If all groups are stopped, let's do separate trials JUST TO MAKE SURE ONE GROUP IS NOT DOMINATING THE DECISION TO STOP!!


      //    Rf_PrintValue(wrap(StoppedGroups));
      if(sum(StoppedGroups)==(J+1)){
        //Let's try the groups separately to make sure we aren't stopping unneccessarily
        StoppedGroups.zeros();

        for(g=0;g<(J+1);g++){
          GetGroup.zeros();
          GetGroup[g]=1;

          Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
          I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
          Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);



          Which1 = g;


          if(Which1==0){
            NewMean = meanmu;
            NewSlope=meanslope;
          }else{
            NewMean = meanmu+MeanInts[Which1-1];
            NewSlope=meanslope+MeanSlopes[Which1-1];
          }





          //Inner MCMC loop, only compute neccessary quantities. Reset counters
          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB = IntB.zeros()+1;
          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          avar=avar.zeros()+1;
          bvar=bvar.zeros()+1;
          muvar=1;
          slopevar=1;
          NumSig=2;
          IntSig=1;
          sigvar=1;
          a.zeros();
          b.zeros();
          MeanVec.zeros();
          GroupVec.zeros();
          mu=NewMean;
          slope=NewSlope;


          for(m=0;m<B;m++){

            if(m<(B/2 + 2)){
              if(m%100==0){




                if((IntMu/NumMu)>.5){
                  muvar=muvar*2;
                }

                if((IntMu/NumMu)<.2){
                  muvar=muvar/2;
                }


                if((IntSlope/NumSlope)>.5){
                  slopevar=slopevar*2;
                }

                if((IntSlope/NumSlope)<.2){
                  slopevar=slopevar/2;
                }







                NumMu = 2;
                NumSlope=2;
                IntMu=1;
                IntSlope=1;
                IntSig=1;
                NumSig=2;




              }
            }









            Munew = mu + as_scalar(arma::randn(1))*muvar;



            alpha= .5*pow((mu-NewMean),2)/varint -.5*pow((Munew-NewMean),2)/varint+   Like2(Y2, I2,  Doses2,  Munew,slope, sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope,  sig,Y2.n_rows);




            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              mu=Munew;
              IntMu=IntMu+1;
            }

            NumMu=NumMu+1;




            //Sample the new slope and the new \beta_g coefficients jointly




            slopenew = slope + as_scalar(arma::randn(1))*slopevar;



            //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


            //Proposal ratio for beta

            alpha =.5*pow((slope-NewSlope),2)/varbeta -.5*pow((slopenew-NewSlope),2)/varbeta;


            alpha =   alpha+  Like2(Y2, I2,  Doses2,   mu,slopenew,sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope, sig,Y2.n_rows);


            U=log(as_scalar(arma::randu(1)));



            U=log(as_scalar(arma::randu(1)));

            if(U<alpha){
              slope=slopenew;
              IntSlope=IntSlope+1;
            }
            NumSlope=NumSlope+1;





            if(m>(B-B1-1)){
              //Make OptDose and StoppedGroups
              StoreInx=m-B1;



              //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
              //First nDose


              k=Which1;

              for(j=0;j<nDose;j++){
                eta1=0;




                eta1=mu+exp(slope)*Dose[j];



                //For all groups add the intercept and slope*Groups

                eta1=exp(eta1);

                if(eta1>100000){
                  DoseProb(StoreInx, k*nDose+j) = 1;

                }else{
                  DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

                }



              }



              //Now we have a matrix containing our sampled dose probabilities









            }

            //End MCMC
          }
          //Is our very last group stopped?

          k=Which1;
          eta2=0;

          for(j=0;j<DoseProb.n_rows;j++){
            //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
            if(DoseProb(j,k*nDose)>Target[k]){
              eta2++;
            }


          }


          if(eta2>(Upper[Which1]*DoseProb.n_rows)){
            StoppedGroups[Which1]=1;
          }






        }





      }







    }








    if(MEANS==0){
      for(k=0;k<(J+1);k++){
        for(j=0;j<nDose;j++){
          MeanVec(k,j)=sum(DoseProb.col(k*nDose+j))/DoseProb.n_rows;

        }

      }



    }




    //Determine Optimal Doses
    for(k=0;k<(J+1);k++){

      if(StoppedGroups(k)==1){
        OptDose(k)=-1;
      }else{
        eta1=MinVec(abs(MeanVec.row(k).t()-Target[k]));
        j=0;

        GroupVec = abs(MeanVec.row(k).t()-Target[k]);

        //  Rf_PrintValue(wrap(MeanVec));


        while(GroupVec[j]!=eta1 ){
          j++;
          //Loop Proceeds until the minimum difference is reached
        }

        OptDose[k]=j;
      }


    }




    //Is this dose bigger than the dose higher than what's been tried?
    for(k=0;k<(J+1);k++){
      if(StoppedGroups[k]==0){
        for(j=0;j<(J+1);j++){
          if(OptDose[j]>=nDose){
            OptDose[j]=nDose-1;
          }
        }




      }



    }

    //Is this dose bigger than the dose higher than what's been tried?
    for(k=0;k<(J+1);k++){
      if(StoppedGroups(k)==0){
        if(DoseTried(k,OptDose[k])==0){
          j=0;

          while(DoseTried(k,j)==1){
            j++;
          }

          OptDose[k]=j;

        }
      }
    }







    OptDose=OptDose+1;





  }else{
    //More than Two Subgroups



    int eta2=0;

    arma::mat MeanVec((J+1),nDose);
    arma::vec NTox(J+1);
    arma::vec OptDose(J+1);
    arma::vec StoppedGroups(J+1);
    arma::vec SuspendGroups(J+1);
    arma::vec GLast(J+1);
    arma::vec nTreated(J+1);
    double stopped=0;
    arma::vec GroupVec(nDose);
    arma::vec TriedGroups(J+1);




    arma::mat DoseTried((J+1),nDose);
    arma::vec nTreated1(J+1);
    nTreated1.zeros();
    INVEC.zeros();

    arma::vec INVECNEW=INVEC;
    arma::vec Y2;
    arma::vec I2;
    arma::vec Groups2;
    arma::vec Doses2;
    arma::vec GetGroup=StoppedGroups;












    //Inner MCMC loop, only compute neccessary quantities. Reset counters
    NumA=NumA.zeros()+2;
    NumB=NumB.zeros()+2;
    IntA=IntA.zeros()+1;
    IntB = IntB.zeros()+1;
    NumMu = 2;
    NumSlope=2;
    IntMu=1;
    IntSlope=1;
    avar=avar.zeros()+1;
    bvar=bvar.zeros()+1;
    muvar=1;
    slopevar=1;
    NumSig=2;
    IntSig=1;
    sigvar=1;
    a.zeros();
    b.zeros();
    MeanVec.zeros();
    GroupVec.zeros();
    mu=meanmu;
    slope=meanslope;
    INVEC.zeros();
    INVECNEW.zeros();
    INVEC = INVEC+1;
    a=MeanInts;
    b=MeanSlopes;


    for(m=0;m<GroupMem.n_rows;m++){
      GroupMem[m]=m+1;
    }

    GroupMemProp=GroupMem;


    for(m=0;m<B;m++){

      if(m<(B/2 + 2)){
        if(m%100==0){


          for(k=0;k<(J+1);k++){
            if((IntA[k]/NumA[k])>.5){
              avar[k]=avar[k]*2;
            }

            if((IntA[k]/NumA[k])<.2){
              avar[k]=avar[k]/2;
            }








            if((IntB[k]/NumB[k])>.5){
              bvar[k]=bvar[k]*2;
            }

            if((IntB[k]/NumB[k])<.2){
              bvar[k]=bvar[k]/2;
            }



          }


          if((IntMu/NumMu)>.5){
            muvar=muvar*2;
          }

          if((IntMu/NumMu)<.2){
            muvar=muvar/2;
          }


          if((IntSlope/NumSlope)>.5){
            slopevar=slopevar*2;
          }

          if((IntSlope/NumSlope)<.2){
            slopevar=slopevar/2;
          }





          NumA=NumA.zeros()+2;
          NumB=NumB.zeros()+2;
          IntA=IntA.zeros()+1;
          IntB=IntB.zeros()+1;

          NumMu = 2;
          NumSlope=2;
          IntMu=1;
          IntSlope=1;
          IntSig=1;
          NumSig=2;




        }
      }



      for(k=0;k<aprop.n_rows;k++){

        if(INVEC[k]>0){

          aprop=a;
          aprop(k)=as_scalar(arma::randn(1))*avar[k] + a(k);






          alpha =    .5*pow((a[k]-MeanInts[k]),2)/varint1 -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


          U=log(as_scalar(arma::randu(1)));


          if(U<alpha){
            a(k)=aprop(k);
            IntA[k]=IntA[k]+1;
          }

          NumA=NumA+1;


        }
      }
      //Now do a MH step for the slope

      for(k=0;k<aprop.size();k++){

        if(INVEC[k]>0){
          bprop=b;
          bprop[k] = b[k] +   as_scalar(arma::randn(1))*bvar[k];



          alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((b[k]-MeanSlopes[k]),2)/varbeta1  +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


          U=log(as_scalar(arma::randu(1)));

          z9[0]=alpha;

          //  Rf_PrintValue(wrap(z9));

          if(U<alpha){
            b[k]=bprop[k];
            IntB[k]=IntB[k]+1;
          }

        }

      }


      NumB=NumB+1;



      //Spike And Slab


      U=as_scalar(arma::randu(1));


      if((3*U)<1){
        //Add All or delete all

        if(sum(INVEC)==(J)){
          //Auto Delete Move




          aprop.zeros();
          bprop.zeros();

          GroupMemProp=GroupMem;
          //Now Randomly Draw Our New Group Assignment


          INVECNEW.zeros();




          GroupMemProp.zeros();




          alpha=0;
          for(k=0;k<bprop.n_rows;k++){
            alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

          }





          //Density when 0 is automatically 1 so no proposal needed
          alpha =  alpha + LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



          U=log(as_scalar(arma::randu(1)));



          if(U<alpha){
            b.zeros();
            a.zeros();
            GroupMem.zeros();
            INVEC.zeros();

            //If Any other GroupMEM are in this group, we need to delete them.



          }






        }else{


          if(sum(INVEC)==0){
            //Autho Add Move

            //ADD Move

            for(k=0;k<aprop.n_rows;k++){
              aprop[k]=as_scalar(arma::randn(1));
              bprop[k]=as_scalar(arma::randn(1))/10;
              GroupMemProp[k]=k+1;
            }



            //Density when 0 is automatically 1 so no proposal needed

            alpha =  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



            for(k=0;k<aprop.n_rows;k++){


              alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]) -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 ;

            }




            alpha=alpha + log(sum(INVEC)+1);

            U=log(as_scalar(arma::randu(1)));

            z9[0]=alpha;

            //  Rf_PrintValue(wrap(z9));

            if(U<alpha){

              for(k=0;k<aprop.n_rows;k++){
                b[k]=bprop[k];
                a[k]=aprop[k];
                INVEC[k]=1;
                GroupMem[k]=k+1;
              }


            }




          }else{
            //Randomly Decide to add or Delete

            U=as_scalar(arma::randu(1));


            if(U<.5){
              //Add


              //ADD Move

              for(k=0;k<aprop.n_rows;k++){
                aprop[k]=as_scalar(arma::randn(1));
                bprop[k]=as_scalar(arma::randn(1))/10;
                GroupMemProp[k]=k+1;
              }



              //Density when 0 is automatically 1 so no proposal needed

              alpha =  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



              for(k=0;k<aprop.n_rows;k++){


                alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]) -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 ;

              }




              alpha=alpha + log(sum(INVEC)+1);

              U=log(as_scalar(arma::randu(1)));

              z9[0]=alpha;

              //  Rf_PrintValue(wrap(z9));

              if(U<alpha){

                for(k=0;k<aprop.n_rows;k++){
                  b[k]=bprop[k];
                  a[k]=aprop[k];
                  INVEC[k]=1;
                  GroupMem[k]=k+1;
                }


              }



            }else{
              //Delete All




              aprop.zeros();
              bprop.zeros();

              GroupMemProp=GroupMem;
              //Now Randomly Draw Our New Group Assignment


              INVECNEW.zeros();




              GroupMemProp.zeros();




              alpha=0;
              for(k=0;k<bprop.n_rows;k++){
                alpha=alpha +  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

              }





              //Density when 0 is automatically 1 so no proposal needed
              alpha =  alpha + LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b.zeros();
                a.zeros();
                GroupMem.zeros();
                INVEC.zeros();

                //If Any other GroupMEM are in this group, we need to delete them.



              }




            }





          }




        }









      }else{



        //Should we swap?
        if((sum(INVEC)==J) || (sum(INVEC)==0)){
          //Add/Delete



          k=Sample1(J);

          if(INVEC[k]==0){
            //ADD Move
            aprop=a;
            bprop=b;
            aprop[k]=as_scalar(arma::randn(1));
            bprop[k]=as_scalar(arma::randn(1))/10;

            GroupMemProp=GroupMem;
            GroupMemProp[k]=k+1;


            //Density when 0 is automatically 1 so no proposal needed
            alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

            alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);


            alpha=alpha + log(sum(INVEC)+1);

            U=log(as_scalar(arma::randu(1)));

            z9[0]=alpha;

            //  Rf_PrintValue(wrap(z9));

            if(U<alpha){
              b[k]=bprop[k];
              a[k]=aprop[k];
              INVEC[k]=1;
              GroupMem[k]=k+1;
            }


          }else{
            //Delete Move
            aprop=a;
            bprop=b;
            //Delete Move
            aprop[k]=0;
            bprop[k]=0;

            GroupMemProp=GroupMem;
            //Now Randomly Draw Our New Group Assignment


            INVECNEW=INVEC;
            INVECNEW[k]=0;



            GroupMemProp[k]=GetNewGroup(INVECNEW);










            //Density when 0 is automatically 1 so no proposal needed
            alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

            alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

            alpha=alpha-log(sum(INVEC)+1);

            U=log(as_scalar(arma::randu(1)));



            if(U<alpha){
              b[k]=0;
              a[k]=0;
              GroupMem=GroupMemProp;
              INVEC[k]=0;

              //If Any other GroupMEM are in this group, we need to delete them.

              for(j=0;j<GroupMem.n_rows;j++){
                if(GroupMem[j]==(k+1)){
                  GroupMem[j]=GroupMem[k];
                }
              }



            }


          }


        }else{




          if((3*U)>2){
            //Swap Move
            //Randomly Pick One In and One Out

            IntIN = GetIn(INVEC);
            IntOUT= GetOut(INVEC);

            //Now IntIN and IntOUT contain the entries that are currently in or out





          }else{
            //Add/Delete



            k=Sample1(J);

            if(INVEC[k]==0){
              //ADD Move
              aprop=a;
              bprop=b;
              aprop[k]=as_scalar(arma::randn(1));
              bprop[k]=as_scalar(arma::randn(1))/10;

              GroupMemProp=GroupMem;
              GroupMemProp[k]=k+1;


              //Density when 0 is automatically 1 so no proposal needed
              alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

              alpha = alpha -.5*pow((aprop[k]-MeanInts[k]),2)/varint1 + log(PROBIN[k])-log(1-PROBIN[k]);


              alpha=alpha + log(sum(INVEC)+1);

              U=log(as_scalar(arma::randu(1)));

              z9[0]=alpha;

              //  Rf_PrintValue(wrap(z9));

              if(U<alpha){
                b[k]=bprop[k];
                a[k]=aprop[k];
                INVEC[k]=1;
                GroupMem[k]=k+1;
              }


            }else{
              //Delete Move
              aprop=a;
              bprop=b;
              //Delete Move
              aprop[k]=0;
              bprop[k]=0;

              GroupMemProp=GroupMem;
              //Now Randomly Draw Our New Group Assignment


              INVECNEW=INVEC;
              INVECNEW[k]=0;



              GroupMemProp[k]=GetNewGroup(INVECNEW);










              //Density when 0 is automatically 1 so no proposal needed
              alpha =  .5*pow((b[k]-MeanSlopes[k]),2)/varbeta1 +  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, aprop, bprop, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);

              alpha = alpha +.5*pow((a[k]-MeanInts[k]),2)/varint1 + log(1-PROBIN[k])-log(PROBIN[k]);

              alpha=alpha-log(sum(INVEC)+1);

              U=log(as_scalar(arma::randu(1)));



              if(U<alpha){
                b[k]=0;
                a[k]=0;
                GroupMem=GroupMemProp;
                INVEC[k]=0;

                //If Any other GroupMEM are in this group, we need to delete them.

                for(j=0;j<GroupMem.n_rows;j++){
                  if(GroupMem[j]==(k+1)){
                    GroupMem[j]=GroupMem[k];
                  }
                }



              }


            }


          }



        }

      }






      //Swap Group

      for(k=0;k<J;k++){
        if(INVEC[k]==0){

          GroupMemProp=GroupMem;
          //Now Randomly Draw Our New Group Assignment


          INVECNEW=INVEC;
          INVECNEW[k]=0;



          GroupMemProp[k]=GetNewGroup(INVECNEW);










          //Density when 0 is automatically 1 so no proposal needed
          alpha = LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMemProp) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);



          U=log(as_scalar(arma::randu(1)));



          if(U<alpha){
            GroupMem[k]=GroupMemProp[k];
          }
        }

      }






      Munew = mu + as_scalar(arma::randn(1))*muvar;



      alpha= .5*pow((mu-meanmu),2)/varint -.5*pow((Munew-meanmu),2)/varint+   LikeMULTI(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);




      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        mu=Munew;
        IntMu=IntMu+1;
      }

      NumMu=NumMu+1;




      //Sample the new slope and the new \beta_g coefficients jointly




      slopenew = slope + as_scalar(arma::randn(1))*slopevar;



      //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


      //Proposal ratio for beta

      alpha =.5*pow((slope-meanslope),2)/varbeta -.5*pow((slopenew-meanslope),2)/varbeta;


      alpha =   alpha+  LikeMULTI(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,i,GroupMem) -  LikeMULTI(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i,GroupMem);


      U=log(as_scalar(arma::randu(1)));



      U=log(as_scalar(arma::randu(1)));

      if(U<alpha){
        slope=slopenew;
        IntSlope=IntSlope+1;
      }
      NumSlope=NumSlope+1;


      if(m>(B-B1-1)){
        //Make OptDose and StoppedGroups
        StoreInx=m-B1;

        for(j=0;j<J;j++){;
          astore(StoreInx,j)=INVEC(j);
          bstore(StoreInx,j)=GroupMem[j];
        }

        mustore[StoreInx]=mu;
        slopestore[StoreInx]=slope;
        sigstore[StoreInx]=sig;


        //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
        //First nDose

        aprop.zeros();
        bprop.zeros();
        for(j=0;j<J;j++){
          for(k=0;k<J;k++){
            if(GroupMem[j]==(k+1)){
              aprop[j]=a[k];
              bprop[j]=b[k];
            }
          }
        }


        for(k=0;k<(J+1);k++){

          for(j=0;j<nDose;j++){
            eta1=0;

            if(k>0){
              eta1 = aprop[k-1]+exp(bprop[k-1]+slope)*Dose[j]+mu;
            }

            if(k==0){
              eta1=mu+exp(slope)*Dose[j];

            }

            //For all groups add the intercept and slope*Groups

            eta1=exp(eta1);

            if(eta1>100000){
              DoseProb(StoreInx, k*nDose+j) = 1;

            }else{
              DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

            }



          }


        }
        //Now we have a matrix containing our sampled dose probabilities









      }

      //End MCMC






    }


    StoppedGroups.zeros();
    //Get MeanVec
    //reset to means




    //Now we have the mean vector



    //Do we have a stopped group?
    for(k=0;k<(J+1);k++){
      eta2=0;

      for(j=0;j<DoseProb.n_rows;j++){
        //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
        if(DoseProb(j,k*nDose)>Target[k]){
          eta2++;
        }


      }







      if(eta2>(Upper[k]*DoseProb.n_rows)){
        StoppedGroups[k]=1;
      }

    }









    //If we don't have at least 3 patients enrolled at the lowest dose, let's drop the Stopped Group
    for(k=0;k<(J+1);k++){
      if(nTreated1[k]<3){
        StoppedGroups[k]=0;
      }

      if(nTreated1[k]>2){
        GetGroup.zeros();
        GetGroup[k]=1;

        Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
        I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
        Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);

        if(GetFullyFollowed(Y2,I2,Doses2,Dose,T1)<3){
          StoppedGroups[k]=0;
        }

      }


    }




    //If all groups are stopped, let's do separate trials JUST TO MAKE SURE ONE GROUP IS NOT DOMINATING THE DECISION TO STOP!!


    //    Rf_PrintValue(wrap(StoppedGroups));
    if(sum(StoppedGroups)==(J+1)){
      //Let's try the groups separately to make sure we aren't stopping unneccessarily
      StoppedGroups.zeros();

      for(g=0;g<(J+1);g++){
        GetGroup.zeros();
        GetGroup[g]=1;

        Y2=ReturnStoppedY(Y,I,Groups,GetGroup,i);
        I2=ReturnStoppedI(Y,I,Groups,GetGroup,i);
        Doses2=ReturnStoppedY(Doses,I,Groups,GetGroup,i);



        Which1 = g;


        if(Which1==0){
          NewMean = meanmu;
          NewSlope=meanslope;
        }else{
          NewMean = meanmu+MeanInts[Which1-1];
          NewSlope=meanslope+MeanSlopes[Which1-1];
        }





        //Inner MCMC loop, only compute neccessary quantities. Reset counters
        NumA=NumA.zeros()+2;
        NumB=NumB.zeros()+2;
        IntA=IntA.zeros()+1;
        IntB = IntB.zeros()+1;
        NumMu = 2;
        NumSlope=2;
        IntMu=1;
        IntSlope=1;
        avar=avar.zeros()+1;
        bvar=bvar.zeros()+1;
        muvar=1;
        slopevar=1;
        NumSig=2;
        IntSig=1;
        sigvar=1;
        a.zeros();
        b.zeros();
        MeanVec.zeros();
        GroupVec.zeros();
        mu=NewMean;
        slope=NewSlope;


        for(m=0;m<B;m++){

          if(m<(B/2 + 2)){
            if(m%100==0){




              if((IntMu/NumMu)>.5){
                muvar=muvar*2;
              }

              if((IntMu/NumMu)<.2){
                muvar=muvar/2;
              }


              if((IntSlope/NumSlope)>.5){
                slopevar=slopevar*2;
              }

              if((IntSlope/NumSlope)<.2){
                slopevar=slopevar/2;
              }







              NumMu = 2;
              NumSlope=2;
              IntMu=1;
              IntSlope=1;
              IntSig=1;
              NumSig=2;




            }
          }









          Munew = mu + as_scalar(arma::randn(1))*muvar;



          alpha= .5*pow((mu-NewMean),2)/varint -.5*pow((Munew-NewMean),2)/varint+   Like2(Y2, I2,  Doses2,  Munew,slope, sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope,  sig,Y2.n_rows);




          U=log(as_scalar(arma::randu(1)));

          if(U<alpha){
            mu=Munew;
            IntMu=IntMu+1;
          }

          NumMu=NumMu+1;




          //Sample the new slope and the new \beta_g coefficients jointly




          slopenew = slope + as_scalar(arma::randn(1))*slopevar;



          //        slopenew = slopevar +   as_scalar(arma::randn(1))*slopevar;


          //Proposal ratio for beta

          alpha =.5*pow((slope-NewSlope),2)/varbeta -.5*pow((slopenew-NewSlope),2)/varbeta;


          alpha =   alpha+  Like2(Y2, I2,  Doses2,   mu,slopenew,sig,Y2.n_rows) -  Like2(Y2, I2,  Doses2,   mu,slope, sig,Y2.n_rows);


          U=log(as_scalar(arma::randu(1)));



          U=log(as_scalar(arma::randu(1)));

          if(U<alpha){
            slope=slopenew;
            IntSlope=IntSlope+1;
          }
          NumSlope=NumSlope+1;





          if(m>(B-B1-1)){
            //Make OptDose and StoppedGroups
            StoreInx=m-B1;



            //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
            //First nDose


            k=Which1;

            for(j=0;j<nDose;j++){
              eta1=0;




              eta1=mu+exp(slope)*Dose[j];



              //For all groups add the intercept and slope*Groups

              eta1=exp(eta1);

              if(eta1>100000){
                DoseProb(StoreInx, k*nDose+j) = 1;

              }else{
                DoseProb(StoreInx, k*nDose+j) = eta1/(1+eta1);

              }



            }



            //Now we have a matrix containing our sampled dose probabilities









          }

          //End MCMC
        }
        //Is our very last group stopped?

        k=Which1;
        eta2=0;

        for(j=0;j<DoseProb.n_rows;j++){
          //     eta2 = eta2 + DoseProb(j,k*nDose)>Target[k];
          if(DoseProb(j,k*nDose)>Target[k]){
            eta2++;
          }


        }


        if(eta2>(Upper[Which1]*DoseProb.n_rows)){
          StoppedGroups[Which1]=1;
        }






      }





    }




    if(MEANS==0){
      for(k=0;k<(J+1);k++){
        for(j=0;j<nDose;j++){
          MeanVec(k,j)=sum(DoseProb.col(k*nDose+j))/DoseProb.n_rows;

        }

      }


    }




    //Determine Optimal Doses
    for(k=0;k<(J+1);k++){

      if(StoppedGroups(k)==1){
        OptDose(k)=-1;
      }else{
        eta1=MinVec(abs(MeanVec.row(k).t()-Target[k]));
        j=0;

        GroupVec = abs(MeanVec.row(k).t()-Target[k]);

        //  Rf_PrintValue(wrap(MeanVec));


        while(GroupVec[j]!=eta1 ){
          j++;
          //Loop Proceeds until the minimum difference is reached
        }

        OptDose[k]=j;
      }


    }




    //Is this dose bigger than the dose higher than what's been tried?
    for(k=0;k<(J+1);k++){
      if(StoppedGroups[k]==0){
        for(j=0;j<(J+1);j++){
          if(OptDose[j]>=nDose){
            OptDose[j]=nDose-1;
          }
        }




      }



    }

    //Is this dose bigger than the dose higher than what's been tried?
    for(k=0;k<(J+1);k++){
      if(StoppedGroups(k)==0){
        if(DoseTried(k,OptDose[k])==0){
          j=0;

          while(DoseTried(k,j)==1){
            j++;
          }

          OptDose[k]=j;

        }
      }
    }










    OptDose=OptDose+1;























  }





  return(OptDose);




}











