clear all; clc; close all; format long e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% We estimate the parameters of a SEIR-type epidemiological model by
%%%% using a maximum a posteriori estimator. All the estimation procedures
%%%% are carried out by LSQNONLIN, although the likelihood function (or
%%%% data misfit) is log-Poisson. The model parameters are estimated from
%%%% the daily records of infections and deaths.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Setup for the ODE solver and the least-square function
tol = 1.e-6;  % ode solver tolerance
% set 'Stats','on' to get more info
% options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats','on');

% note: set Refine switch to avoid interpolation
options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
                                                         'off','Refine',1);
options2 = [];%optimset('MaxFunEvals',10000,'MaxIter',7000,'TolFun',...
                                                       1e-30,'TolX',1e-30);
options3 = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Experimental Data

%%%% Daily Cases, Hospitalizations and Deaths
%%%% CHICAGO, USA:

DATA = importdata('COVID-19_Daily_Cases__Deaths__and_Hospitalizations_2020201204.csv');
data = DATA.data;
ndays = 0; %%% Keeping the last 10 days for forecast.
data = data(1:end-ndays,:);

DATA2 = importdata('COVID-19_Hospital_Capacity_Metrics_20201204.csv');
dataB = DATA2.data;
dataB = dataB(1:end-ndays,:);


%%% We shall delete the last 10 days.

t_actual = 0:size(data,1);

%%%% Smoothing the data - averaging every 7e consecutive days:

%%%% Daily Infections, Hospitalizations and Deaths

data1 = data;
for jj=4:size(data,1)-3
for ii = 1:size(data,2) 
data1(jj,ii) = mean(data(jj-3:jj+3,ii));
end
end

%%%% Total Population:
N = 2705988;      % CHICAGO, USA

%%%% Population proportion on each age range:
Proportion = ones;
PropInfections = sum(data1(:,1))/N;
PropHosp = sum(data1(:,3))/sum(data1(:,1));
PropICU = 0.4;
PropDeath = sum(data1(:,2))/(0.4*sum(data1(:,3)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial parameters

I_M0 = 7;      % Potential initial infected and mild at t=0
E_0 = 0;       % initial exposed
I_H0 = 0;      % initial infected and hospitalized at t=0
I_I0 = 0;      % initial infected and in ICU at t=0
R_0 = 0;       % initial recovered 
D_0 = 0;       % initial dead

%--------------------------------------------------------------------------
%%%% Asymptomatic
p = ones-0.17;
params.p = p; %%% proportio of symptomatic individuals
q = 0.58;
params.q = q;  %%% reduction in transmissibility
%--------------------------------------------------------------------------

%  params is a structure used to pass parameters to the
%   ODE solver

S_0 = N-(I_M0+I_I0+I_H0+R_0+D_0+E_0);    % Suceptible pop.,  excluding initial infected 
params.N = N;  % N = total population

NumberOfAgeClasses = 1;  % total age ranges

yinit(1:NumberOfAgeClasses) = S_0*Proportion;
yinit(NumberOfAgeClasses+1:2*NumberOfAgeClasses) = E_0*Proportion;
yinit(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = (1-p)*I_M0*Proportion;
yinit(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = p*I_M0*Proportion;
yinit(4*NumberOfAgeClasses+1:5*NumberOfAgeClasses) = I_H0*Proportion;
yinit(5*NumberOfAgeClasses+1:6*NumberOfAgeClasses) = I_I0*Proportion;
yinit(6*NumberOfAgeClasses+1:7*NumberOfAgeClasses) = R_0*Proportion;
yinit(7*NumberOfAgeClasses+1:8*NumberOfAgeClasses) = D_0*Proportion;
yinit = yinit/N;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Model Parameters

params.sigma = 1./5.1;   % inverse of mean incubation time
params.NumberOfAgeClasses = NumberOfAgeClasses;

% Mean time until recovery
nu_M = 1./14;
nu_H = 1./12;
nu_I = 1./9;

%--------------------------------------------------------------------------
% Mean time until death
mu_M = zeros;
mu_H = zeros;
mu_I = 1./7;

%--------------------------------------------------------------------------
% Mean time until passing to another infectious Class
gamma_M = 1./4;
gamma_H = 1./3.5;

%--------------------------------------------------------------------------
% Proportion of individuals that will recover
p_M  = (1-PropHosp);
params.p_M = p_M; % In Mild conditions

p_H = 1-PropICU;%[1,0.8531,0.7217,0.7024,0.6224,0.5148,0.4665];
params.p_H = p_H; % Hospitalized individuals

p_I = 1-PropDeath;%[1,0.9524,0.9153,0.7704,0.7523,0.7095,0.3772];
params.p_I = p_I; % In ICU

%--------------------------------------------------------------------------
% Proportion of individuals that will die
q_M = zeros(1,NumberOfAgeClasses);
params.q_M = q_M;  % In Mild conditions

q_H = zeros(1,NumberOfAgeClasses);
params.q_H = q_H; % Hospitalized individuals

%--------------------------------------------------------------------------
%%% RATES
%%% Recovery Rate

Recovery_M = (nu_M+mu_M+gamma_M)*p_M;
Recovery_H = zeros;%(nu_H+mu_H+gamma_H)*p_H;
Recovery_I = (nu_I+mu_I)*p_I;

params.Recovery_A = 1/14;
params.Recovery_M = Recovery_M; % in Mild conditions
params.Recovery_H = Recovery_H; % Hospitalized individuals
params.Recovery_I = Recovery_I; % in ICU individuals

%%% Getting Worse Rate

GetWorse_M = (nu_M+mu_M+gamma_M)*(1-p_M-q_M);
GetWorse_H = ones;%(nu_H+mu_H+gamma_H)*(1-p_H-q_H);

params.GetWorse_M = GetWorse_M; % Mild conditions
params.GetWorse_H = GetWorse_H; % Hospitalized individuals

%%% Death Rate

Death_M = (nu_M+mu_M+gamma_M)*q_M;
Death_H = (nu_H+mu_H+gamma_H)*q_H;
Death_I = (nu_I+mu_I)*(1-p_I);

params.Death_M = Death_M; % in Mild conditions
params.Death_H = Death_H; % Hospitalized individuals
params.Death_I = Death_I; % in ICU individuals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Correcting Hospitalization and Death Rates

n = size(data1,1)-size(dataB,1);
Hosp = ones(length(t_actual),1);
Hosp(3:end) = min(1,data1(2:end,3)./data1(1:end-1,1));
Hosp = min(20,Hosp/GetWorse_M);
% ICU = ones(length(t_actual),1);
% ICU(n+2:end) = min(1,dataB(:,15)./data1(n:end-1,3));
% ICU = min(20,ICU/GetWorse_H);
Death = ones(size(t_actual));
Death(3:end) = min(1,data1(2:end,2)./data1(1:end-1,3))/Death_I;
params.factorWorse = @(t)interp1(t_actual,Hosp,t);
% params.factorICU = @(t)interp1(t_actual,ICU,t);
params.factorDeath = @(t)interp1(t_actual,Death,t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
% A priori value for the transmission parameter
R_zero = 1.4*2.8/0.1782;
gamma = 1/18;

%--------------------------------------------------------------------------
% Transmission parameters
beta = 2.2911*R_zero*gamma;

beta_M = ones;
params.beta_M = beta_M;      % In Mild conditions
params.beta_H = beta_M;  % Hospitalized individuals
params.beta_I = beta_M; % In ICU

params.b = 0.1;%0.5;
params.c = 0.01;%0.05;

paramsOld = params;
yinitOld = yinit;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 1 until 20:

day = 20;

LB = [1/N,1E-3];
UB = [50/N,10];
unknowns0=[1/N,beta];
priors = unknowns0;

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:

OF = @(unknowns)ObjFun_InitialPopBetaWithoutAge(t_actual(1:day),params,...
data1(1:day-1,:),options,priors,yinit,Proportion,PropInfections,unknowns);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);

I_M = unknowns(1);
yinit(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = (1-p)*I_M*PropInfections;
yinit(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = p*I_M*PropInfections;
yinit(1:NumberOfAgeClasses) = Proportion-I_M*PropInfections;
params.a = unknowns(2);

PARAMS1 = unknowns;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 1
%%%% until 19:

BETA = zeros(length(t_actual),1);
BETA(1) = ones; 
unknowns0 = ones;
priors = unknowns0;
yinit2 = yinit;
yb = zeros(length(t_actual),8*NumberOfAgeClasses);
yb(1,:) = yinit;
R0 = zeros(1,length(t_actual)-1);
for jj =1:day-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta4(t_actual2,params,unknowns2,...
                    data1(jj,:),options,priors,BETA(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-2,10,options2);
unknowns0 = unknowns;
priors = unknowns;
BETA(jj+1) = unknowns;
beta2 = @(t)interp1(t_actual2,BETA(jj:jj+1),t);
[~,y2] = ode45(@(t,y)seir_death_age_beta2(t,y,params,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))
R0(jj) = basic_reproduction_rate_beta(Proportion,params,BETA(jj+1),t_actual(jj+1))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating basic parameters from day 21 until 80:

dt = 20;
nmax = floor(length(t_actual)/dt);
PARAMS2 = zeros(nmax,1);
for ss = 1:nmax-1

LB = 1E-3;
UB = 10;
unknowns0= beta;
priors = unknowns0;
yinit2 = y2(end,:);

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:
if ss == nmax-1
aux = (day-1)+(ss-1)*dt+1:length(t_actual)-1;
aux2 = (day-1)+(ss-1)*dt+1:length(t_actual);    
else
aux = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt;
aux2 = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt+1;
end
OF = @(unknowns)ObjFun_BetaWithoutAge(t_actual(aux2),params,...
           data1(aux,:),options,priors,yinit2,unknowns);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
params.a = unknowns(1);

PARAMS2(ss,:) = unknowns;

unknowns0 = ones;
priors = unknowns0;
for jj = aux
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns)ObjFun_beta4(t_actual2,params,unknowns,...
                  data1(jj,:),options,priors,BETA(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-2,10,options2);
unknowns0 = unknowns;
priors = unknowns;
BETA(jj+1,:) = unknowns;
beta2 = @(t)interp1(t_actual2,BETA(jj:jj+1),t);
[~,y2] = ode45(@(t,y)seir_death_age_beta2(t,y,params,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2(end,:);
yb(jj+1,:) = yinit2;
disp(num2str(unknowns))
R0(jj) = basic_reproduction_rate_beta(Proportion,params,BETA(jj+1),t_actual(jj+1));
end
end

factor = zeros(length(t_actual),1);
factorD = factor;
factorW = params.factorWorse;
factorDeath = params.factorDeath;
for jj = 1:length(t_actual)
factor(jj) = factorW(t_actual(jj));
factorD(jj) = factorDeath(t_actual(jj));
end  

%% Final Number of Cases for each Age Range
sigma = params.sigma;
NewCases = p*sigma*yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses)'*N;

%% Total Number of Deaths for each day
NewDeaths = zeros(size(yb(:,1)));
NewHospB = zeros(size(yb(:,1),1),NumberOfAgeClasses);
NewHosp = zeros(size(yb(:,1)));
Deaths = zeros(length(yb(:,1)),NumberOfAgeClasses);
for jj=1:NumberOfAgeClasses
NewDeaths = NewDeaths + Death_M(jj)*yb(:,3*NumberOfAgeClasses+jj)...
+ Death_H(jj)*yb(:,4*NumberOfAgeClasses+jj)...
+ Death_I(jj)*factorD.*yb(:,5*NumberOfAgeClasses+jj);
NewHosp = NewHosp + GetWorse_M(jj)*factor.*yb(:,3*NumberOfAgeClasses+jj);
NewHospB(:,jj) = GetWorse_M(jj)*factor.*yb(:,3*NumberOfAgeClasses+jj);
Deaths(:,jj) = Death_M(jj)*yb(:,3*NumberOfAgeClasses+jj)...
+Death_H(jj)*yb(:,4*NumberOfAgeClasses+jj)...
+ Death_I(jj)*factorD.*yb(:,5*NumberOfAgeClasses+jj);
end

%% Total Hospitalizations for each Age Range
TotalHosp = sum(NewHospB)'*N; 
TotalInfections = ...
      sum(p*sigma*yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses))'*N;
TotalDeaths = sum(Deaths)'*N;

disp('  Infections   Hosp.    Deaths  ')
disp(num2str(round([TotalInfections,TotalHosp,TotalDeaths])))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bootstrapping_20201204;
save dataChicago_20201204b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plotting Results:
%%% CI Evaluation
%%% R0

aux = sort(R0Boot);
MedianR0 = median(R0Boot);
aux2 = round(0.0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95R0 = [min(aux);max(aux)];

%%% Beta
aux = sort(BETABoot);
aux2 = round(0.0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95BETA = [min(aux);max(aux)];

%%% New Cases:
aux = sort(NewCasesBoot');
aux2 = round(0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewCases = [min(aux);max(aux)]*N;

%%% New Deaths:
aux = sort(NewDeathsBoot');
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewDeaths = [min(aux);max(aux)]*N;

%%% New Hospitalizations:
aux = sort(NewHospBoot');
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewHosp = [min(aux);max(aux)]*N;

aux = sort(PARAMS1Boot);
aux2 = round(0.0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95PARAMS1 = [min(aux);max(aux)];

aux = sort(PARAMS2Boot);
aux2 = round(0.0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95PARAMS2 = [min(aux);max(aux)];

CI95R0t = CI95R0;
CI95BETAt = CI95BETA;
R0t = R0;
MR0 = median(R0Boot);
MR0t = median(R0Boot);
len = 3;
for jj = 1+len:length(t_actual)-(1+len) 
CI95R0t(1,jj) = mean(CI95R0(1,jj-len:jj+len));
CI95R0t(2,jj) = mean(CI95R0(2,jj-len:jj+len));
R0t(jj) = mean(R0(jj-len:jj+len));
MR0t(jj) = mean(MR0(jj-len:jj+len));
CI95BETAt(1,jj) = mean(CI95BETA(1,jj-len:jj+len));
CI95BETAt(2,jj) = mean(CI95BETA(2,jj-len:jj+len));
end

t_span = datetime(2020,3,01) + caldays(0:length(t_actual)-1);

H = [100 100 600 400];
%%% plotting Results:
figure
hold on
% grid on
box on
title('Daily Infections - Chicago')
h1=area(t_span,CI95NewCases(2,:),'linestyle',':','FaceColor',[255,160,122]/255);%[51,236,255]/255);
h2=area(t_span,CI95NewCases(1,:),'linestyle',':','FaceColor',[1,1,1]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
bar(t_span(2:end),data1(:,1),'FaceColor',[0 0.75 0.75],'EdgeColor','none')%[255,119,51]/255,'EdgeColor','none')
plot(t_span,NewCases,'r','LineWidth',2)
legend('Reported','Estimated','Location','NorthWest')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end-3)])
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'InfectionsChicagoAS.fig');
print('-dpng','InfectionsChicagoAS');



figure
hold on
% grid on
box on
title('Daily Deaths - Chicago')
h1=area(t_span,CI95NewDeaths(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewDeaths(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data1(:,2),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
plot(t_span,N*NewDeaths,'r','LineWidth',2)
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Reported','Estimated')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end-3)])
% ylim([0,1.3*max(data(:,3))])
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'DeathsChicagoAS.fig');
print('-dpng','DeathsChicagoAS');


%%%%%%%%%
figure
hold on
% grid on
box on
title('Daily Hospitalizations - Chicago')
h1=area(t_span,CI95NewHosp(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewHosp(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data1(:,3),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t_span,N*NewHosp,'r','LineWidth',2)
% plot([t_span(day),t_span(day)],[0,1.3*max(data(:,2))],'--k','LineWidth',2)
legend('Reported','Estimated')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end-3)])
% ylim([0,1.3*max(data(:,2))])
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'HospChicagoAS.fig');
print('-dpng','HospChicagoAS');

%%% R0

figure
hold on
grid on
box on
title('Chicago')

area(t_span(2:end),CI95R0(2,:),'linestyle',':','FaceColor',[51,236,255]/255)
area(t_span(2:end),CI95R0(1,:),'linestyle',':','FaceColor',[1,1,1])
plot(t_span(2:end),R0,'b','LineWidth',2)
plot(t_span(2:end),median(R0Boot),'r','LineWidth',2)
% plot([t_span(day),t_span(day)],[0,8],'--k','LineWidth',2)
plot(t_span(2:end),ones(size(t_span(2:end))),'k')
ylabel('Basic Reproduction Rate')
xlim([t_span(1),t_span(end-3)])
ylim([0,10])
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'R0ChicagoAS.fig');
print('-dpng','R0ChicagoAS');

figure
hold on
grid on
box on
title('Chicago')

h1=area(t_span(2:end),CI95R0t(2,:),'linestyle',':','FaceColor',[255,160,122]/255);%[51,236,255]/255);
h2=area(t_span(2:end),CI95R0t(1,:),'linestyle',':','FaceColor',[1,1,1]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t_span(2:end),R0t,'r','LineWidth',2)
plot(t_span(2:end),ones(size(t_span(2:end))),'k')
ylabel('Effective Reproduction Number')
xlim([t_span(1),t_span(end-3)])
ylim([0,10])
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'R0Chicago_smooth.fig');
print('-dpng','R0Chicago_smooth');



