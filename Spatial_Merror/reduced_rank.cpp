//-----------------------------
#include <math.h> 
#include <TMB.hpp> 
//#include <boost/accumulators/statistics/density.hpp>

using namespace density;
template<class Type>

Type objective_function<Type>::operator() () {  
     //declare data
     DATA_VECTOR(predX); 
     DATA_MATRIX(sp_basis); //sp_basis does not have an intercept; other basis functions do
     DATA_MATRIX(y);
     DATA_MATRIX(identity_matrix);
     DATA_INTEGER(true_rank);
     int N = y.rows();
     int dimSigma = sp_basis.cols();
     int num_basis = sp_basis.cols();

     Type nll = Type(0.0); //Negative log-likelihood

	 //declare parameters
     PARAMETER_VECTOR(beta); //Matrix of fixed effects
     PARAMETER_VECTOR(cholSigma); //Vector of Choleksy decomposition for spatial random effects covariance matrix
     PARAMETER_VECTOR(alpha_rho); //spatial random effects
     //PARAMETER(sigma_nugget_sq); //nugget effect variance
     
     
     //To create cholesky decompostion as lower triangular matrix
     matrix<Type> cholSigmamat(dimSigma,true_rank); 
     int Sigma_ind = 0;
     for(int l2=0; l2<(true_rank) ; l2++) { for(int l1=0; l1<(dimSigma); l1++) { 
          if(l2 > l1) { cholSigmamat(l1,l2)=0; }
          if(l2 == l1) { cholSigmamat(l1,l2) = exp(cholSigma(Sigma_ind)); Sigma_ind += 1; }
          if(l2 < l1) { cholSigmamat(l1,l2) = cholSigma(Sigma_ind); Sigma_ind += 1; } 
          } }
     matrix<Type> Sigma = cholSigmamat*cholSigmamat.transpose() + 0.0001*identity_matrix;
     
     ADREPORT(Sigma); 
     
     
	 vector<Type> eta(N);
     for(int i=0; i<N; i++) {
          eta(i)=beta(0); 
          eta(i) += beta(1)*predX(i);
          for(int l1=0; l1<num_basis; l1++) { eta(i) += sp_basis(i,l1)*alpha_rho(l1); }
                }
     
     

	 //Random effects from multivariate normal distribution
     nll += MVNORM(Sigma)(alpha_rho); //Returns negative logL by default

	 //data likelihood
	 for(int i=0; i<N; i++) {
		 nll -= dpois(y(i,0),exp(eta(i)),true);
		  }

	 return nll; 
	 }




