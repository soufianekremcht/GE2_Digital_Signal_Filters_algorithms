function [x_hat,K,P] = kalman_filter(y,x_hat,P,A,B,u,T)

    %   R : Noise Covariance
    % 	H : measurement map scalar
    % 	Q : initial estimated covariance
    % 	P : initial error covariance 
    % 	U_hat : initial estimated state
    % 	K : initial Kalman gain
    %   T : Sampling Time
    
    % We Predict the Estimate Value
    % We Update The Kalman Gain & Estimate Value based on the current mesurements
    
    % We consider that R & Q are time varing variables Or Constants
    R = rand()*10+5;
    H = 1.00;
    Q = 5+rand()*50;
    
    
    % Prediction
	x_hat = x_hat +B*u;
	P = P + T*(A*P+P*A + Q);
    K = P.*H./(H.*P.*H+R);
    
    % Update
    % Estimated Value
    x_hat = x_hat + K.*(y-H.*x_hat);
    P = (1-K.*H).*P + Q; 
    
end