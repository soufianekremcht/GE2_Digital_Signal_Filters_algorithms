clc;
clear all;
clear workspace;
close all;

% INPUTS
% x = noisy signal
% y = reference signalsinit
% N = filter order


% OUTPUTS
% xest = estimated signal
% b = Wiener filter coefficents
% MSE = mean squared error

[y Fs] = audioread('sound_sample_simon_sinek.mp3');
audio_plyer = audioplayer(y,Fs);
% Creating reduced noise and adding it with the original signal
rand_noise = y + (randn(size(y)))/64;

% se = [1/6 1/6 1/6 1/6 1/6 1/6];
% E = imfilter(rand_noise,se);

startLog = " Testing Wiener Algorithm is Starting"
K = wiener2(rand_noise);
endLog = "Finish Testing Wiener Algorithm"

subplot(2,3,1)
imshow(K)
title("GrayScale Image")

