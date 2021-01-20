clear all; clc; close all; format long e; tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We estimate the parameters of a SEIR-type epidemiological model by
%%%% using a maximum a posteriori estimator. All the estimation procedures
%%%% are carried out by LSQNONLIN, although the likelihood function (or
%%%% data misfit) is log-Poisson. The model parameters are estimated from
%%%% the daily records of infections and deaths.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup for the ODE solver and the least-square function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data

%%%% Daily Cases, Hospitalizations and Deaths
%%%% CHICAGO, USA:

DATA = importdata('data-by-day.csv');
data = DATA.data;
ndays = 5; %%% Keeping the last 10 days for forecast.
data = data(1:end-ndays,:);

% DATA2 = importdata('Chicago_COVID-19_Daily_Testing_-_By_Person_20201003.csv');
DATA2 = importdata('tests.csv');
dataB = DATA2.data;
dataB = dataB(1:end-ndays,:);


%%% We shall delete the last 10 days.

t_actual = 0:size(data,1);

% %%%% Smoothing the data - averaging every 7e consecutive days:
% data1 = data;
% for jj=7:size(data,1)
% for ii = 1:size(data,2) 
% data1(jj,ii) = mean(data(jj-6:jj,ii));
% end
% end
% 
% data2 = dataB;
% for jj=7:size(dataB,1)
% for ii = 1:size(dataB,2) 
% data2(jj,ii) = mean(dataB(jj-6:jj,ii));
% end
% end

data1 = data;
for jj=4:size(data,1)-3
for ii = 1:size(data,2) 
data1(jj,ii) = mean(data(jj-3:jj+3,ii));
end
end

%%%% Daily Performed Tests
data2 = dataB;
for jj=4:size(dataB,1)-3
for ii = 1:size(dataB,2) 
data2(jj,ii) = mean(dataB(jj-3:jj+3,ii));
end
end
l=0;
t_span = datetime(2020,02,29) + caldays(0:length(t_actual)-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AutoCorr = zeros; % plot auto-correlations: yes = 1, no = 0;
SaveFig = 1;%zeros;  % plot auto-correlations: yes = 1, no = 0;
H = 300; % height of the plots

day = 1;

%%%%%% The time series beginning and end.

%%%% Shorter Period - closer to the mean value

% %%%% Larger Period - not so close to the mean, may be nonstationary
% inicio = 136; % 14-Jul-2020
% fim = 260; % 15-Nov-2020

END = length(t_actual)-1;%inicio-1;

AUX = [79,185;186,241;242,size(data1,1)-1-12;155,192;193,252;253,size(data1,1)-1-12;155,size(data1,1)-1-12];
TC = [];
TCD = [];
for zz = 1:size(AUX,1)
inicio = AUX(zz,1);%154; % 01-Aug-2020
fim = AUX(zz,2);%219; % 01-Oct-2020
DataNYCWithoutAge_20210101;
TC = [TC;TotalCorrectedHosp];
TCD = [TCD;TotalCorrectedDeath];
end