clc;
clear all;
x = 1:1:34;

noisy_signal_file = fopen('noisy_temp_sensor_data.txt','r');
filtered_signal_file= fopen('filtered_temp_sensor_data.txt','w');

formatSpec = '%f\n';
sizeNOISY = Inf;
Message = "Noisy Temperature Sensor Data"
% read signal Data
[NOISY,count] = fscanf(noisy_signal_file,formatSpec,sizeNOISY)
FILTRED = zeros(count,1);
c = 1;
% initial states
% 	P : initial error covariance 
% 	K : initial Kalman gain
Ki = 0;
Pi = 10;
initial_estimate = 0;

for i = 1:count-1
    y = NOISY(i);
    % FILTRED : Estimated Value after Using Kalman Filter
    
    [FILTRED(c),Ki,Pi] = kalman_filter(NOISY(i),initial_estimate,Pi,[1],[1],[10],1);
    initial_estimate = FILTRED(c)
   
    % add FILTERED Value to Results File
    fprintf(filtered_signal_file,'%f\n',FILTRED(c));
    c = c+1;
end

% save Filtered Sensor Data into .txt File
fclose(noisy_signal_file);
fclose(filtered_signal_file);

% Plotting The Results
% X : Time (s)
x = [0:1/(count-1):1]
plot(x,FILTRED,x,NOISY);
title('kalman Filter Example')
legend('Filtered Signal ','Noisy Signal')

xlabel('Time (s)')
ylabel('Temperature Sensor Data (°C)')

