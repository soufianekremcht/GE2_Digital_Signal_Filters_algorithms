

// Realisé Par : Kremcht Soufiane et Youssef Elmrabet
// GE2 Signal 2  2021/2022


#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;


double KALMAN(double U);

// Kalman Filter function definition

// We Have the systeme represented as : x= A*x +B*u & y= C*x_p 

double KALMAN_Predict(double* x_hat,double* P,double A,double B,double u,double T){
	
	// U : Input
	// x_hat : Previous estimated value
	// P : Error Covariance (increased in Prediction phase)
	// A,B : System matrices
	// T : Sampling Time
	cout << "x_hat : "<< *x_hat << " P : " << *P << endl;
	
	double Q = 1;
	
	*x_hat = A*(*x_hat) +B*u;
	
	*P = *P+ T*(A*(*P)*A + Q);
	
	cout << "Predicted x_hat : "<< *x_hat << " Predicted  P : " << *P << endl;
}

double KALMAN_Update(double y,double* x_hat,double *P){
	/*
	R : Noise Covariance
	Q : initial estimated covariance
	P : initial error covariance (decreased in Update Phase)
	x_hat : predicted estimate state (not the true estimated state)
	 
	K : initial Kalman gain
	C : measurement map scalar
	y : Measurement
	
	
	*/
	double C = 1.00;
	double R = 4;
	double Q = 1;
	double K = 0;
	
	// The higher R gets , K get reduced which means that the filter
	// will remove more noise but it can slow down the filter
	
	
	// Getting Kalman Gain
	K = (*P)*C/(C*(*P)*C+R);
	
	// Update the estimated Values (true estimate state)
	
	*x_hat = (*x_hat) + K*(y-C*(*x_hat));
	
	
	// Update Error Covariance
	*P = (1-K*C)*(*P) +Q;
	
	
	return *x_hat;
}



int main() {
	
	
	double P = 10;
	// estimated mesurement
	double x_hat = 10;
	
    string A4,file;
    
    vector<double>U_N ;
    // number of lines
    int i=0;
    
    // Read Signal Value from a file
    
    ifstream ceoff("noisy.txt");
    // if the file is open
    
    if (ceoff.is_open()){

       
       while(!ceoff.eof())
       {
    
           getline(ceoff,A4,'\n');
           U_N.push_back(stod(A4));
           i+=1;
       }
       
       ceoff.close();
       cout<<"number of entries :"<<i<<endl;
       
   }
   

   double NOISY_SIGNAL_DATA[i];
   
   copy(U_N.begin(),U_N.end(),NOISY_SIGNAL_DATA);
   // let's filter out these NOISY measurements ONE AT A TIME
   // initialize empty array of filtered values
   double FILTRED_SIGNAL[i];  
   // Apply Kalman Filter
   
   for (int c=0;c<i;c++){
        
        // Prediction 
        KALMAN_Predict(&x_hat,&P,1,1,0,0.1);
        // Update
        FILTRED_SIGNAL[c] = KALMAN_Update(NOISY_SIGNAL_DATA[c],&x_hat,&P);
        
        //cout << "Valeur : " << endl;
        //cout<<FILTRED_SIGNAL[c]<<endl;
    }
    
    ofstream MyFile("noisy_resultat.txt");

   // Save Results to noisy_resultat.txt
    int size =  sizeof(FILTRED_SIGNAL)/sizeof(FILTRED_SIGNAL[0]);
    
  	for (int p = 0;p<size;p++){
  		MyFile << FILTRED_SIGNAL[p]<<endl;
	}
	MyFile.close();
     
   //system("pause");
   return 0;
}
