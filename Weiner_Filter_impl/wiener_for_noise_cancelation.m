function wiener_for_noise_cancelation


%Author: Alex Dytso
%Description: This on of the project that shows how to implement Wiener
%filter as noise cancellations. 
% Our desired response is d(n) ouptut is x(n). We have noise v(n), v1(n)
% and v2(n)
% that have the following relationship v1(n)=0.8*v1(n-1)+v(n) and v2(n)=-0.6v2(n-1)+v(n)
% where v(n) is AWGN. 
% The ouptut is x(n)=d(n)+v1(n)
clc, clear all, close all

%%%%%Part2

theta=2*pi*rand-pi; % generating random varible thetha between -pi and pi
w0=0.05*pi;
N=10000;
n=1:N;

dn=sin(n*w0+theta); %generating our siganl

vn=randn(1,N); %generating 10000 points of White Gaussian Noise


% generating v1n=0.8*v1(n-1)+v(n) and v2(n)=-0.6v2(n-1)+v(n)

v1(1)=vn(1);
v2(1)=vn(1);
for i=2:N
    v1(i)=0.1*v1(i-1)+vn(i);
    v2(i)=-0.05*v2(i-1)+vn(i);
    
end


%generating x(n)=d(n)+v1(n)

xn=dn+v1;


%%%%%Part 3

%ploting xn and dn 
figure(1)
plot(1:200,xn(1:200))
hold on
plot(1:200,dn(1:200),'r')
xlabel('n')
legend('xn','dn')
title('Plot of x(n) and d(n)')
%ploting v2n
figure(2)
plot(1:200,v2(1:200))

xlabel('n')

title('Plot of v_2(n)')
%%%%% part 4 Computing auto-correlation of v2n

%comparing computatation with my function to computation with matlab
rv2=AUTOCORR(v2(1:40)) %computation with my function
rv2m=xcorr(v2(1:40),'biased'); %computation with MatLab

figure(3)
stem(rv2)
hold on
stem(rv2m(20:39),'xr')
xlabel('k')
title('Autocorrlation of v_2(n)')
legend('My function', 'MatLab function')


%%%% part5 Computing cross-correlation between xn and v2n

%comparing computation via my function with computation via MatLab
rxv2=CROSSCORR(xn(1:50),v2(1:50))  %my function 
rxv2m=xcorr(xn(1:50),v2(1:50),'biased');  %MatLab function

figure(4)
stem(rxv2)
hold on
stem(rxv2m(20:39),'xr')
xlabel('k')
title('Cross-correlation between x(n) and v_2(n)')
legend('My function', 'MatLab function')

%%%% part 6  Finding Coefficient of Wiener filter of order 1

%Coefficient are compute by w=Rv2^(-1)*rxv2


Rv2=[rv2(1),rv2(2);rv2(2),rv2(1)] %Autocorrelation matrix of v2


w=(Rv2)^(-1)*rxv2(1:2)  %Coefficients of Wiener fileter of oreder 1


v1nEst=conv(w,v2); %filtering 

dest=xn-v1nEst(1:N); %subtracting estimate of v1 from xn
figure(5)

plot(dest(1:200))
hold on
plot(dn(1:200),'r')
xlabel('n')
title('Estimated d(n) vs actual d(n). Filter order 1')
legend('Estimate d(n)', 'Actual d(n)')
%%%%part 7  order 6
 % autocorrelation matrix 
Rv2ord6 = [rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7);
         rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6);
         rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5);
         rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4);
         rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3);
         rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2);
         rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1)]
     
     
     word6=inv(Rv2ord6)*rxv2(1:7)%computing optimal coefficients
     
     v1nEst6=conv(word6,v2);%filtering with Wiener filter
     
     
     
     dest6=xn-v1nEst6(1:10000);%computing the error
figure(6)
plot(dest6(1:200))
hold on
plot(dn(1:200),'r')
xlabel('n')
title('Estimated d(n) vs actual d(n). Filter order 6')
legend('Estimate d(n)', 'Actual d(n)')
%%%part 8

     % autocorrelation matrix 
     Rv2ord12=[rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13);
              rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12);
              rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11);
              rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10);
              rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9);
              rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8);
              rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7);
              rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6);
              rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5);
              rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4);
              rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3);
              rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2);
              rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1)]
          
     word12 = inv(Rv2ord12)*rxv2(1:13) %computing optimal coefficients
     
     v1nEst12f=filter(word12,1,v2); %filtering with Wiener filter 
    
     dest12=xn-v1nEst12f(1:10000); %computing the error
     
figure(7)
 plot(dest12(1:200))
 hold on

 plot(dn(1:200),'r')         
 xlabel('n')
title('Estimated d(n) vs actual d(n). Filter order 12')
legend('Estimate d(n) ', 'Actual d(n)')

%%%%part 9 (EXTRA FILTER ORDER 20)


Rv2ord20=    [rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13),rv2(14),rv2(15),rv2(16),rv2(17),rv2(18),rv2(19),rv2(20);
              rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13),rv2(14),rv2(15),rv2(16),rv2(17),rv2(18),rv2(19);
              rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13),rv2(14),rv2(15),rv2(16),rv2(17),rv2(18);
              rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13),rv2(14),rv2(15),rv2(16),rv2(17);
              rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13),rv2(14),rv2(15),rv2(16);
              rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13),rv2(14),rv2(15);
              rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13),rv2(14);
              rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12),rv2(13);
              rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11),rv2(12);
              rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10),rv2(11);
              rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9),rv2(10);
              rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8),rv2(9);
              rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7),rv2(8);
              rv2(14),rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6),rv2(7);
              rv2(15),rv2(14),rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5),rv2(6);
              rv2(16),rv2(15),rv2(14),rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4),rv2(5);
              rv2(17),rv2(16),rv2(15),rv2(14),rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3),rv2(4);
              rv2(18),rv2(17),rv2(16),rv2(15),rv2(14),rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2),rv2(3);
              rv2(19),rv2(18),rv2(17),rv2(16),rv2(15),rv2(14),rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1),rv2(2);
              rv2(20),rv2(19),rv2(18),rv2(17),rv2(16),rv2(15),rv2(14),rv2(13),rv2(12),rv2(11),rv2(10),rv2(9),rv2(8),rv2(7),rv2(6),rv2(5),rv2(4),rv2(3),rv2(2),rv2(1)];


           word19=inv(Rv2ord20)*rxv2(1:20)
     
     v1nEst19=conv(word19,v2);
     dest19=xn-v1nEst19(1:10000);
figure(8)
plot(dest19(1:200))
hold on
 plot(dn(1:200),'r')         
 xlabel('n')
title('Estimated d(n) vs actual d(n). Filter order 19')
legend('Estimate d(n)', 'Actual d(n)')
          
        

%ETRA minimum error

en1=(dn-dest).^2;

en6=(dn-dest6).^2;

en12=(dn-dest12).^2;

en19=(dn-dest19).^2;

J1=mean(en1)
J6=mean(en6)
J12=mean(en12)
J19=mean(en19)

%  figure(9)
%  plot(en1(1:100))
%  hold on
% plot(en6(1:100),'r')
%  plot(en12(1:100),'g')
% plot(en19(1:100),'y')

%%%% EXTRA fft

figure(10)

subplot(4,1,1)
plot(1:200,fft(dest(1:200)))
hold on

subplot(4,1,2)
plot(1:200,fft(dest6(1:200)))

subplot(4,1,3)
plot(1:200,fft(dest12(1:200)))

subplot(4,1,4)
plot(1:200,fft(dest19(1:200)))

end


%%%% Autocorrelation and Cross-correlation functions



function y=AUTOCORR(x) %this function computes autocorrelation
K=length(x);
for index=1:K
    
    y(index,:)=sum([zeros(1,index-1),x].*[x,zeros(1,index-1)])/(K);
end



end


function y=CROSSCORR(X,Z) %this function computes autocorrelation

x=[X,zeros(1,length(Z)-length(X))]; %Making length of x and z to be the same
z=[Z,zeros(1,length(X)-length(Z))];

N=length(z);

for index=1:N
    
    y(index,:)=sum([x,zeros(1,index-1)].*[zeros(1,index-1),z])/(N);
end



end
