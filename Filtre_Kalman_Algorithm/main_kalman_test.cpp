
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>




using namespace std;



// Kalman Filter funciton header
double KALMAN(double U);

// Kalman Filter function definition

// U : 

double KALMAN(double U){
	/*
	R : Noise Covariance
	H : measurement map scalar
	Q : initial estimated covariance
	P : initial error covariance 
	U_hat : initial estimated state
	K : initial Kalman gain
	
	
	*/
	static const double R = 10;
	static const double H = 1.00;
	
	static double Q = 10;
	static double P = 0;
	static double U_hat = 0;
	static double K = 0;
	
	// The higher R gets , K get reduced which means that the filter
	// will remove more noise but it can slow down the filter
	
	// constants

	
	// Updating Kalman Gain
	K = P*H/(H*P*H+R);
	
	cout << "Updating U_hat " <<U_hat<<endl;
	// Update the estimated Values
	U_hat = U_hat + K*(U-H*U_hat);
	
	
	
	// Update Error Covariance
	P = (1-K*H)*P +Q;
	
	return U_hat;
}



int main() {
	
	
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
        
        FILTRED_SIGNAL[c]=KALMAN(NOISY_SIGNAL_DATA[c]);
        //cout << "Valeur : " << endl;
        
        cout<<FILTRED_SIGNAL[c]<<endl;
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
