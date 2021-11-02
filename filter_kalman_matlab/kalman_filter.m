function [U_hat,K,P] = kalman_filter(U,K, P)

    %   R : Noise Covariance
    % 	H : measurement map scalar
    % 	Q : initial estimated covariance
    % 	P : initial error covariance 
    % 	U_hat : initial estimated state
    % 	K : initial Kalman gain
    
    % We iterate throught the input_data & we set Ki & Pi & U_hat
    % 
    
  
    R = 40;
    H = 1.00;
    Q = 10;
    U_hat = 0;
    
    K = P.*H./(H.*P.*H+R);
    
    % Estimated Value
    U_hat = U_hat + K.*(U-H.*U_hat);
    
    P = (1-K.*H).*P + Q;
    
end