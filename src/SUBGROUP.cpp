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


double Like(arma::vec Y,  //Vector of Toxicity Times
            arma::vec I,  // Vector of Censoring Indi
            arma::vec Dose, //Vector of Doses given to Patients
            arma::vec Group, //Vector of Group Membership, 0 is baseline
            double mu, // Reference Time
            double slope, //Baseline Slope for Dose Toxicity Relationship
            arma::vec a, //Vector of intercepts for non-baseline groups
            arma::vec b, //vector of slopes for non-baseline groups
            double T1 // Variance Parameter for Lognormal Distribution
){

  double LogL = 0; // Will store likelihood to return


  //First we calculate the eta term in the GLM for each patient.

  //Baseline vector
  arma::vec eta(Group.n_rows);

  eta.zeros();

  //Needed for For loops
  int m;
  int k=1;




  for(m=0;m<eta.n_rows;m++){
    for(k=0;k<a.n_rows;k++){
      if(Group[m]==(k+2)){
        eta[m] = a[k]+b[k]*Dose[m];
      }
    }

    //For all groups add the intercept and slope*Groups
    eta[m]=eta[m]+mu+slope*Dose[m];


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
    }

  }








  return(LogL);


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

  npats=npats-1;
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
      if(Group[m]==(k+2)){
        eta[m] = a[k]+b[k]*Dose[m];
      }
    }

    //For all groups add the intercept and slope*Groups
    eta[m]=eta[m]+mu+slope*Dose[m];


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
    }

  }








  return(LogL);


}




















//[[Rcpp::export]]
List GroupMCMCLN(arma::vec Y,  //Vector of Toxicity Times
                 arma::vec I,  // Vector of Censoring Indi
                 arma::vec Doses, //Vector of Doses given to Patients
                 arma::vec Groups, //Vector of Group Membership, 0 is baseline
                 double meanmu, //Reference Time
                 double meanslope, //Prior mean of the slope
                 arma::vec MeanInts , // Vector of prior means on group intercepts
                 arma::vec MeanSlopes, // Vector of prior means on group slopes
                 int B, // Number of iterations to perform
                 double T1
){


  //Innitialize parameters
  double mu=0;
  double slope=0;
  double sig=T1;

  arma::vec Dose=Doses;
  arma::vec Group=Groups;


  //Important quantities we will need
  int  J=MeanInts.n_rows;
  double alpha=0;
  double U=0;
  double signew=0;
  double slopenew=0;
  double Munew=0;

  //For loop integers
  int j=0;
  int k=0;
  int m=0;


  int StoreInx=0;
  //Vector Of group intercepts
  arma::vec a(J);
  //Vector of group slopes
  arma::vec b=a;
  //Vectors for Group MCMC

  for(m=0;m<a.n_rows;m++){
    a[m]=0;
    b[m]=0;
  }


  //Vectors for proposals
  arma::vec aprop = a;
  arma::vec bprop=b;


  //SetUp Storage Matrices
  int B1=B;

  arma::mat astore(B1,a.n_rows);
  arma::mat bstore =astore;
  arma::vec mustore(B1);
  arma::vec slopestore(B1);
  arma::vec sigstore(B1);


  double NumA=2;
  double NumB=2;
  double IntA=1;
  double IntB = 1;
  double NumMu = 2;
  double NumSlope=2;
  double IntMu=1;
  double IntSlope=1;

  double avar=1;
  double bvar=1;
  double muvar=1;
  double slopevar=1;
  NumericVector z9(2);
  NumericVector zprop(5);

  double NumSig=2;
  double IntSig=1;
  double sigvar=1;


  for(m=0;m<B;m++){

    if(m<(B/2 + 2)){
      if(m%250==0){
        if((IntA/NumA)>.8){
          avar=avar*2;
        }

        if((IntA/NumA)<.2){
          avar=avar/2;
        }


        if((IntSig/NumSig)>.8){
          sigvar=sigvar*2;
        }

        if((IntSig/NumSig)<.2){
          sigvar=sigvar/2;
        }




        if((IntB/NumB)>.8){
          bvar=bvar*2;
        }

        if((IntB/NumB)<.2){
          bvar=bvar/2;
        }


        if((IntMu/NumMu)>.8){
          muvar=muvar*2;
        }

        if((IntMu/NumMu)<.2){
          muvar=muvar/2;
        }


        if((IntSlope/NumSlope)>.8){
          slopevar=slopevar*2;
        }

        if((IntSlope/NumSlope)<.2){
          slopevar=slopevar/2;
        }



        zprop[0]=IntA/NumA;
        zprop[1]=IntB/NumB;
        zprop[2]=IntMu/NumMu;
        zprop[3]=IntSlope/NumSlope;
        zprop[4]=IntSig/NumSig;
        //   Rf_PrintValue(zprop);


        NumA=2;
        NumB=2;
        IntA=1;
        IntB = 1;
        NumMu = 2;
        NumSlope=2;
        IntMu=1;
        IntSlope=1;
        IntSig=1;
        NumSig=2;

        zprop[0]=muvar;
        zprop[1]=slopevar;
        zprop[2]=avar;
        zprop[3]=bvar;
        zprop[4]=sigvar;
        //  Rf_PrintValue(zprop);




      }
    }



    for(k=0;k<aprop.n_rows;k++){

      aprop=a;
      aprop(k)=as_scalar(arma::randn(1))*avar + a(k);






      alpha =    .5*pow((a[k]-MeanInts[k]),2) -.5*pow((aprop[k]-MeanInts[k]),2) +  Like(Y, I,  Dose, Group,  mu,slope, aprop, b, sig) -  Like(Y, I,  Dose, Group,  mu,slope, a, b, sig);


      U=log(as_scalar(arma::randu(1)));


      if(U<alpha){
        a(k)=aprop(k);
        IntA=IntA+1;
      }

      NumA=NumA+1;



    }
    //Now do a MH step for the slope

    for(k=0;k<aprop.size();k++){

      bprop=b;
      bprop[k]=max1(as_scalar(arma::randn(1))*bvar + b[k],-slope);

      //        alpha = log(1+pow((b[k]-MeanSlopes[k])/2.5,2)) - log(1+pow((bprop[k]-MeanSlopes[k])/2.5,2))  +  Like(Y, I,  Dose, Group,  mu,slope, a, bprop, sig) -  Like(Y, I,  Dose, Group,  mu,slope, a, b, sig);




      alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2) +.5*pow((b[k]-MeanSlopes[k]),2)  +  Like(Y, I,  Dose, Group,  mu,slope, a, bprop, sig) -  Like(Y, I,  Dose, Group,  mu,slope, a, b, sig);


      U=log(as_scalar(arma::randu(1)));

      z9[0]=alpha;

      //  Rf_PrintValue(wrap(z9));

      if(U<alpha){
        b[k]=bprop[k];
        IntB=IntB+1;
      }
      NumB=NumB+1;

    }





    Munew = mu + as_scalar(arma::randn(1))*muvar;



    alpha= .5*pow((mu-meanmu),2) -.5*pow((Munew-meanmu),2)+   Like(Y, I,  Dose, Group,  Munew,slope, a, b, sig) -  Like(Y, I,  Dose, Group,  mu,slope, a, b, sig);




    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      mu=Munew;
      IntMu=IntMu+1;
    }

    NumMu=NumMu+1;





    //Now do a MH step for the slope

    slopenew = max1(slope + as_scalar(arma::randn(1))*slopevar,.01);




    alpha = .5*pow((slope-meanslope),2) -.5*pow((slopenew-meanslope),2)     +    Like(Y, I,  Dose, Group,  mu,slopenew, a, b, sig) -  Like(Y, I,  Dose, Group,  mu,slope, a, b, sig);



    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      slope=slopenew;
      IntSlope=IntSlope+1;
    }
    NumSlope=NumSlope+1;




    if(m>(B-B1-1)){
      //Store Values in Matrix
      StoreInx = m-B+B1;

      for(j=0;j<J;j++){;
        astore(StoreInx,j)=a(j);
        bstore(StoreInx,j)=b[j];
      }

      mustore[StoreInx]=mu;
      slopestore[StoreInx]=slope;
      sigstore[StoreInx]=sig;

    }

  }









  List z1 = List::create(astore,bstore,mustore,slopestore,sigstore);


  return(z1);


}








int Sample1(int Num){
  double U=as_scalar(arma::randu(1))*Num;

  U=floor(U);

  return(U);

}





//[[Rcpp::export]]
List SimTrial1(int nSims, //Number of Simulations to Run
               int Nmax, //Max patients to enroll in trial
               double  T1, //Reference Time
               arma::vec Target, //Reference Probability Target (or vector)
               arma::vec Dose, //Standardized Vector of Doses
               arma::vec DoseStart, //Vector (or just one) Starting Dose Level
               double TGroup, //Time to suspend group accrual after 3 patients enrolled
               arma::vec Upper, //Thresholds for subgroup suspension
               double Accrue, //Expected Monthly accrual rate
               arma::vec groupprob, //Subgroup Enrollment Probabilities
               int Family, //Indicator of Distribution family to simulate from. Alphabetical order
               arma::mat Param1, //Group Specific Parameter matrix for each dose and subgroup to simulate from
               arma::mat Param2,
               double meanmu, //Prior mean of intercept
               double meanslope, //Prior mean of slope
               arma::vec MeanInts, //Prior means of group intercepts
               arma::vec MeanSlopes ){







int Group=0;
  //Important For loop integer
  int m =0; //For the inner MCMC looprep
  int i=0; //For the patient index
  int rep=0; //For the simulation repetition.


  //First Fill out Target Vector
int nDose=Dose.n_rows;


  //Important Integer Quantities
  int B=2000; //Number of iterations for MCMC

  //Make List objects we'll use

//Trial Vectors
arma::vec Y(Nmax);
arma::vec I(Nmax);
arma::vec Groups(Nmax);
arma::vec Doses(Nmax);
arma::vec Times(Nmax);
arma::vec ACC(Nmax);


  //Innitialize parameters for MCMC
  double mu=0;
  double slope=0;
  double sig=T1;



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





  double NumA=2;
  double NumB=2;
  double IntA=1;
  double IntB = 1;
  double NumMu = 2;
  double NumSlope=2;
  double IntMu=1;
  double IntSlope=1;

  double avar=1;
  double bvar=1;
  double muvar=1;
  double slopevar=1;
  double trialtime=0;
  NumericVector z9(2);
  NumericVector zprop(5);

  double NumSig=2;
  double IntSig=1;
  double sigvar=1;

  arma::vec etaG(J+1);
  arma::mat DoseProb(1000,(J+1)*nDose);
  arma::mat MeanDose(J+1,nDose);
  int B1=999;
  double eta1=0;
  double DEC=0;

  arma::mat MeanVec((J+1),nDose);
  arma::vec NTox(J+1);
  arma::vec OptDose(J+1);
  arma::vec StoppedGroups(J+1);
arma::vec SuspendGroups(J+1);
arma::vec GLast(J+1);
arma::vec nTreated(J+1);
double stopped;
arma::vec TrialTimes(nSims);
arma::mat OptimalDoses(nSims,(J+1));
arma::mat NTOX = OptimalDoses;
arma::vec GroupVec(nDose);



arma::mat DoseStore(nSims,Nmax);
arma::mat GroupStore(nSims,Nmax);
arma::mat DoseTried((J+1),nDose);




  //Wrap nreps around
for(rep=0;rep<nSims;rep++){


  //Wrap i around

MeanVec.zeros();
DoseProb.zeros();
  //Reset all trial vectors
  trialtime=0;
  stopped=0;
  OptDose.zeros();
  StoppedGroups.zeros();
  SuspendGroups.zeros();
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

  Group=Sample1((J+1));

  Groups[0]=Group;

  Doses[0]=Dose[DoseStart[Group]];

  nTreated[Group]=nTreated[Group]+1;

  ACC[0]=0;

  if(Family==0){
    Times(0)=R::rexp(Param1(Group,DoseStart[Group]));
  }



if(Family==1){
  Times[0]=R::rgamma(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
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


  for(i=1;i<Nmax;i++){
    //Start the Trial, Let's do this



trialtime=trialtime+R::rexp(Accrue);
    Group=Sample1((J+1));

    //Is this a new group? If so we automatically enroll them at starting dose
    eta1=0;
    for(j=0;j<i;j++){
      if(Group==Groups[j]){
        eta1=1;
      }
    }



    //Should we unsuspend any groups?
   for(k=0;k<(J+1);k++){
     if((trialtime-GLast[k])>TGroup){
       //Followed the group long enough so let's unsuspend it.
       SuspendGroups[k]=0;
     }
   }
//Check to see if this is a new subgroup
    if(eta1==0){
      //Automatically enroll this patient


      nTreated[Group]=nTreated[Group]+1;

      Groups[i]=Group;

      Doses[i]=Dose[DoseStart[Group]];

      ACC[i]=trialtime;

      if(Family==0){
        Times(i)=R::rexp(Param1(Group,DoseStart[Group]));
      }



      if(Family==1){
        Times[i]=R::rgamma(Param1(Group,DoseStart[Group]),Param2(Group,DoseStart[Group]));
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


    }else{
     //Check if group is in a suspended group
      if(SuspendGroups[Group]==0){
        DEC=0;
        //while loop to skip certain patients.

while(DEC==0){
//If this patient group is

I.zeros();

      //Setup Y
      for(k=0;k<i;k++){
Y[k]=min1(min1(Times[k],trialtime-ACC[k]),T1);
        if(Times[k]<T1){
          I[k]=1;
        }
      }

  //Inner MCMC loop, only compute neccessary quantities. Reset counters
   NumA=2;
     NumB=2;
     IntA=1;
     IntB = 1;
     NumMu = 2;
     NumSlope=2;
     IntMu=1;
     IntSlope=1;
     avar=1;
     bvar=1;
     muvar=1;
     slopevar=1;
     NumSig=2;
     IntSig=1;
     sigvar=1;
     a.zeros();
     b.zeros();


  for(m=0;m<B;m++){

    if(m<(B/2 + 2)){
      if(m%250==0){
        if((IntA/NumA)>.8){
          avar=avar*2;
        }

        if((IntA/NumA)<.2){
          avar=avar/2;
        }


        if((IntSig/NumSig)>.8){
          sigvar=sigvar*2;
        }

        if((IntSig/NumSig)<.2){
          sigvar=sigvar/2;
        }




        if((IntB/NumB)>.8){
          bvar=bvar*2;
        }

        if((IntB/NumB)<.2){
          bvar=bvar/2;
        }


        if((IntMu/NumMu)>.8){
          muvar=muvar*2;
        }

        if((IntMu/NumMu)<.2){
          muvar=muvar/2;
        }


        if((IntSlope/NumSlope)>.8){
          slopevar=slopevar*2;
        }

        if((IntSlope/NumSlope)<.2){
          slopevar=slopevar/2;
        }



        zprop[0]=IntA/NumA;
        zprop[1]=IntB/NumB;
        zprop[2]=IntMu/NumMu;
        zprop[3]=IntSlope/NumSlope;
        zprop[4]=IntSig/NumSig;
        //   Rf_PrintValue(zprop);


        NumA=2;
        NumB=2;
        IntA=1;
        IntB = 1;
        NumMu = 2;
        NumSlope=2;
        IntMu=1;
        IntSlope=1;
        IntSig=1;
        NumSig=2;

        zprop[0]=muvar;
        zprop[1]=slopevar;
        zprop[2]=avar;
        zprop[3]=bvar;
        zprop[4]=sigvar;
        //  Rf_PrintValue(zprop);




      }
    }



    for(k=0;k<aprop.n_rows;k++){

      aprop=a;
      aprop(k)=as_scalar(arma::randn(1))*avar + a(k);






      alpha =    .5*pow((a[k]-MeanInts[k]),2) -.5*pow((aprop[k]-MeanInts[k]),2) +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


      U=log(as_scalar(arma::randu(1)));


      if(U<alpha){
        a(k)=aprop(k);
        IntA=IntA+1;
      }

      NumA=NumA+1;



    }
    //Now do a MH step for the slope

    for(k=0;k<aprop.size();k++){

      bprop=b;
      bprop[k]=max1(as_scalar(arma::randn(1))*bvar + b[k],-slope);





      alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2) +.5*pow((b[k]-MeanSlopes[k]),2)  +  Like1(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


      U=log(as_scalar(arma::randu(1)));

      z9[0]=alpha;

      //  Rf_PrintValue(wrap(z9));

      if(U<alpha){
        b[k]=bprop[k];
        IntB=IntB+1;
      }
      NumB=NumB+1;

    }





    Munew = mu + as_scalar(arma::randn(1))*muvar;



    alpha= .5*pow((mu-meanmu),2) -.5*pow((Munew-meanmu),2)+   Like1(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);




    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      mu=Munew;
      IntMu=IntMu+1;
    }

    NumMu=NumMu+1;





    //Now do a MH step for the slope

    slopenew = max1(slope + as_scalar(arma::randn(1))*slopevar,.01);




    alpha = .5*pow((slope-meanslope),2) -.5*pow((slopenew-meanslope),2)     +    Like1(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);



    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      slope=slopenew;
      IntSlope=IntSlope+1;
    }
    NumSlope=NumSlope+1;





    if(m>B1){
    //Make OptDose and StoppedGroups
StoreInx=m-B1-1;

    //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
    //First nDose


for(k=0;k<(J+1);k++){

  for(j=0;j<nDose;j++){
  eta1=0;

      if(k>0){
        eta1 = a[k-1]+b[k-1]*Dose[j];
      }


    //For all groups add the intercept and slope*Groups
    eta1=eta1+mu+slope*Dose[j];

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




  //Do we have a stopped group?
  for(k=0;k<(J+1);k++){
    eta1=0;

    for(j=0;j<1000;j++){
      eta1 = eta1 + DoseProb(k*nDose)>Target[k];
    }

    eta1=eta1/1000;
    if(eta1>Upper[k]){
      StoppedGroups[k]=1;
    }

  }



  //If we don't have at least 3 patients enrolled, let's drop the Stopped Group
for(k=0;k<(J+1);k++){
  if(nTreated[k]<3){
    StoppedGroups[k]=0;
  }
}

//Should we stop the trial?
if(sum(StoppedGroups)==(J+1)){
  break;
  stopped=1;
}




  if(StoppedGroups[Group]==0){
    //Enroll this patient
    DEC=1;
  }else{
    //Contiunue to Accrue patients in the trial
    trialtime=trialtime+R::rexp(Accrue);
    Group=Sample1((J+1));
  }



    }



      }else{
        //Try Next Patient
        //Contiunue to Accrue patients in the trial
        trialtime=trialtime+R::rexp(Accrue);
        Group=Sample1((J+1));



        //Should we unsuspend any groups?
        for(k=0;k<(J+1);k++){
          if((trialtime-GLast[k])>TGroup){
            //Followed the group long enough so let's unsuspend it.
            SuspendGroups[k]=0;
          }
        }



      }



//Assign Next Dose
//Whats the optimal dose?


//Get MeanVec
for(k=0;k<(J+1);k++){
  for(j=0;j<nDose;j++){
    MeanVec(k,j) = sum(DoseProb.col(k*nDose+j))/1000;
  }

}
//Now we have the mean vector



//Determine Optimal Doses
for(k=0;k<(J+1);k++){

  eta1=MinVec(abs(MeanVec.row(k).t()-Target(k)));
  j=0;




GroupVec = abs(MeanVec.row(k).t()-Target(k));

  while(GroupVec[j]!=eta1 ){
    j++;
    //z9[0]=j;
    //z9[1]=abs(MeanVec.row(k).t()-Target(k))[j];
    //Rf_PrintValue(z9);
    //Loop Proceeds until the minimum difference is reached
  }

  OptDose[k]=j;



}
//Now we have the subgroup specific optimal doses in the vector OptDose



//Now let's see if the dose is bigger than the largest dose.
if(OptDose[Group]>=nDose){
 OptDose[Group]=nDose-1;
}



//Is this dose bigger than the dose higher than what's been tried?
if(DoseTried(Group,OptDose[Group])==0){
  j=OptDose[Group];

  while(DoseTried(Group,j)==0){
    j=j-1;
  }

  OptDose[Group]=j+1;

}




//Assign Patient to subgroup, optdose, add accrual time.


Groups[i]=Group;
nTreated[Group]=nTreated[Group]+1;

//If this is the 3rd patient treated in this subgroup, suspend accrual for TGroup
if(nTreated[Group]==3){
  SuspendGroups[Group]=1;
  GLast[Group]=trialtime;
}



Doses[i]=Dose[OptDose[Group]];

nTreated[Group]=nTreated[Group]+1;

ACC[i]=trialtime;
//Based on the Group and Optimal dose, randomly generate the next patient Data


if(Family==0){
  Times(i)=R::rexp(Param1(Group,DoseStart[Group]));
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



  I.zeros();

  //Setup Y
  for(k=0;k<Nmax;k++){
    Y[k]=min1(Times[k],T1);
    if(Times[k]<T1){
      I[k]=1;
    }
  }





  //Inner MCMC loop, only compute neccessary quantities. Reset counters
  NumA=2;
  NumB=2;
  IntA=1;
  IntB = 1;
  NumMu = 2;
  NumSlope=2;
  IntMu=1;
  IntSlope=1;
  avar=1;
  bvar=1;
  muvar=1;
  slopevar=1;
  NumSig=2;
  IntSig=1;
  sigvar=1;
  a.zeros();
  b.zeros();
  MeanVec.zeros();
  GroupVec.zeros();


  for(m=0;m<B;m++){

    if(m<(B/2 + 2)){
      if(m%250==0){
        if((IntA/NumA)>.8){
          avar=avar*2;
        }

        if((IntA/NumA)<.2){
          avar=avar/2;
        }


        if((IntSig/NumSig)>.8){
          sigvar=sigvar*2;
        }

        if((IntSig/NumSig)<.2){
          sigvar=sigvar/2;
        }




        if((IntB/NumB)>.8){
          bvar=bvar*2;
        }

        if((IntB/NumB)<.2){
          bvar=bvar/2;
        }


        if((IntMu/NumMu)>.8){
          muvar=muvar*2;
        }

        if((IntMu/NumMu)<.2){
          muvar=muvar/2;
        }


        if((IntSlope/NumSlope)>.8){
          slopevar=slopevar*2;
        }

        if((IntSlope/NumSlope)<.2){
          slopevar=slopevar/2;
        }



        zprop[0]=IntA/NumA;
        zprop[1]=IntB/NumB;
        zprop[2]=IntMu/NumMu;
        zprop[3]=IntSlope/NumSlope;
        zprop[4]=IntSig/NumSig;
        //   Rf_PrintValue(zprop);


        NumA=2;
        NumB=2;
        IntA=1;
        IntB = 1;
        NumMu = 2;
        NumSlope=2;
        IntMu=1;
        IntSlope=1;
        IntSig=1;
        NumSig=2;

        zprop[0]=muvar;
        zprop[1]=slopevar;
        zprop[2]=avar;
        zprop[3]=bvar;
        zprop[4]=sigvar;
        //  Rf_PrintValue(zprop);




      }
    }



    for(k=0;k<aprop.n_rows;k++){

      aprop=a;
      aprop(k)=as_scalar(arma::randn(1))*avar + a(k);






      alpha =    .5*pow((a[k]-MeanInts[k]),2) -.5*pow((aprop[k]-MeanInts[k]),2) +  Like1(Y, I,  Doses, Groups,  mu,slope, aprop, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


      U=log(as_scalar(arma::randu(1)));


      if(U<alpha){
        a(k)=aprop(k);
        IntA=IntA+1;
      }

      NumA=NumA+1;



    }
    //Now do a MH step for the slope

    for(k=0;k<aprop.size();k++){

      bprop=b;
      bprop[k]=max1(as_scalar(arma::randn(1))*bvar + b[k],-slope);





      alpha =  -.5*pow((bprop[k]-MeanSlopes[k]),2) +.5*pow((b[k]-MeanSlopes[k]),2)  +  Like1(Y, I,  Doses, Groups,  mu,slope, a, bprop, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);


      U=log(as_scalar(arma::randu(1)));

      z9[0]=alpha;

      //  Rf_PrintValue(wrap(z9));

      if(U<alpha){
        b[k]=bprop[k];
        IntB=IntB+1;
      }
      NumB=NumB+1;

    }





    Munew = mu + as_scalar(arma::randn(1))*muvar;



    alpha= .5*pow((mu-meanmu),2) -.5*pow((Munew-meanmu),2)+   Like1(Y, I,  Doses, Groups,  Munew,slope, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);




    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      mu=Munew;
      IntMu=IntMu+1;
    }

    NumMu=NumMu+1;





    //Now do a MH step for the slope

    slopenew = max1(slope + as_scalar(arma::randn(1))*slopevar,.01);




    alpha = .5*pow((slope-meanslope),2) -.5*pow((slopenew-meanslope),2)     +    Like1(Y, I,  Doses, Groups,  mu,slopenew, a, b, sig,i) -  Like1(Y, I,  Doses, Groups,  mu,slope, a, b, sig,i);



    U=log(as_scalar(arma::randu(1)));

    if(U<alpha){
      slope=slopenew;
      IntSlope=IntSlope+1;
    }
    NumSlope=NumSlope+1;





    if(m>B1){
      //Make OptDose and StoppedGroups
      StoreInx=m-B1-1;

      //Fill in Big matrix of dose toxicity probabilities for each subgroup and dose
      //First nDose


      for(k=0;k<(J+1);k++){

        for(j=0;j<nDose;j++){
          eta1=0;

          if(k>0){
            eta1 = a[k-1]+b[k-1]*Dose[j];
          }


          //For all groups add the intercept and slope*Groups
          eta1=eta1+mu+slope*Dose[j];

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



  //Get MeanVec
  for(k=0;k<(J+1);k++){
    for(j=0;j<nDose;j++){
      MeanVec(k,j) = sum(DoseProb.col(k*nDose+j))/1000;
    }

  }
  //Now we have the mean vector

  //Determine Optimal Doses
  for(k=0;k<(J+1);k++){

    eta1=MinVec(abs(MeanVec.row(k).t()-Target[k]));
    j=0;

    GroupVec = abs(MeanVec.row(k).t()-Target[k]);

    while(GroupVec[j]!=eta1 ){
      j++;
      //Loop Proceeds until the minimum difference is reached
    }

    OptDose[k]=j;



  }


  //If any optimal doses bigger than biggest, let's drop it to nDose


  for(k=0;k<(J+1);k++){
    if(OptDose[k]>=nDose){
      OptDose[k]=nDose-1;
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





  List z1 = List::create(OptimalDoses,NTOX,TrialTimes,DoseStore,GroupStore);


  return(z1);




}







