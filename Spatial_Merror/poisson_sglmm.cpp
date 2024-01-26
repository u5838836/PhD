//-----------------------------
#include <math.h> 
#include <TMB.hpp> 
#include <boost/accumulators/statistics/density.hpp>

using namespace density;
template<class Type>

Type objective_function<Type>::operator() () {  
     //declare data
     DATA_MATRIX(w); 
     DATA_MATRIX(sp_basis1); //sp_basis1 does not have an intercept; other basis functions do
     DATA_MATRIX(sp_basis2);
     DATA_MATRIX(y);
     int N = y.rows();
     int dimSigma1 = sp_basis1.cols();
     int dimSigma2 = sp_basis2.cols();
     int num_basis1 = sp_basis1.cols();
     int num_basis2 = sp_basis2.cols();

     Type nll = Type(0.0); //Negative log-likelihood

	 //declare parameters
     PARAMETER_VECTOR(beta); //Matrix of fixed effects
     PARAMETER_VECTOR(cholSigma1); //Vector of Choleksy decomposition for spatial random effects covariance matrix
     PARAMETER_VECTOR(cholSigma2); //Vector of Choleksy decomposition for X random effects covariance matrix
     PARAMETER_VECTOR(theta1); //spatial random effects
     PARAMETER_VECTOR(X_thetas); //X random effects
     PARAMETER(logsigma_e); //standard deviation of measurement error

     //To create cholesky decompostion as lower triangular matrix
     matrix<Type> cholSigmamat1(dimSigma1,dimSigma1); 
     int Sigma_ind0 = 0;
     for(int l2=0; l2<(dimSigma1) ; l2++) { for(int l1=0; l1<(dimSigma1); l1++) { 
          if(l2 > l1) { cholSigmamat1(l1,l2)=0; }
          if(l2 == l1) { cholSigmamat1(l1,l2) = exp(cholSigma1(Sigma_ind0)); Sigma_ind0 += 1; }
          if(l2 < l1) { cholSigmamat1(l1,l2) = cholSigma1(Sigma_ind0); Sigma_ind0 += 1; } 
          } }
     matrix<Type> Sigma1 = cholSigmamat1*cholSigmamat1.transpose();
     ADREPORT(Sigma1); 
     
     //To create cholesky decompostion as lower triangular matrix
     matrix<Type> cholSigmamat2(dimSigma2,dimSigma2); 
     int Sigma_ind = 0;
     for(int l2=0; l2<dimSigma2 ; l2++) { for(int l1=0; l1<dimSigma2; l1++) { 
          if(l2 > l1) { cholSigmamat2(l1,l2)=0; } 
          if(l2 == l1) { cholSigmamat2(l1,l2) = exp(cholSigma2(Sigma_ind)); Sigma_ind += 1; } 
          if(l2 < l1) { cholSigmamat2(l1,l2) = cholSigma2(Sigma_ind); Sigma_ind += 1; } 
          } }
     matrix<Type> Sigma2 = cholSigmamat2*cholSigmamat2.transpose();
     ADREPORT(Sigma2); 
     
     

	 vector<Type> eta1(N);
     for(int i=0; i<N; i++) {
          eta1(i)=beta(0); 
          for(int l1=0; l1<num_basis2; l1++) { eta1(i) += beta(1)*sp_basis2(i,l1)*X_thetas(l1); }
          for(int l1=0; l1<num_basis1; l1++) { eta1(i) += sp_basis1(i,l1)*theta1(l1); }
                }
     
     
     
	 vector<Type> eta2(N);
     for(int i=0; i<N; i++) {
          eta2(i)=0; 
          for(int l1=0; l1<num_basis2; l1++) { eta2(i) += sp_basis2(i,l1)*X_thetas(l1); }
                }

	 //Random effects from multivariate normal distribution
     nll += MVNORM(Sigma1)(theta1); //Returns negative logL by default
     nll += MVNORM(Sigma2)(X_thetas);

	 //data likelihood
	 for(int i=0; i<N; i++) {
		 nll -= dpois(y(i,0),exp(eta1(i)),true);
		 nll -= dnorm(w(i,0), eta2(i), exp(logsigma_e) , true );
		  }

	 return nll; 
	 }




