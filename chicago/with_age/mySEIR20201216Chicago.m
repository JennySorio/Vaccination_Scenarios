% clear all; clc; close all; format long e; tic;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% We estimate the parameters of a SEIR-type epidemiological model by
% %%%% using a maximum a posteriori estimator. All the estimation procedures
% %%%% are carried out by LSQNONLIN, although the likelihood function (or
% %%%% data misfit) is log-Poisson. The model parameters are estimated from
% %%%% the daily records of infections and deaths.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Setup for the ODE solver and the least-square function
% tol = 1.e-6;  % ode solver tolerance
% options = odeset('AbsTol', tol,'RelTol',tol,'MaxOrder',5,'Stats',...
%                                                          'off','Refine',1);
% options2 = optimset('MaxFunEvals',10000,'MaxIter',7000,'TolFun',...
%                                                        1e-30,'TolX',1e-30);
% options3 = [];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Experimental Data
% 
% %%%% Daily Cases, Hospitalizations and Deaths
% %%%% NYC, USA:
% 
% DATA = importdata('COVID-19_Daily_Cases__Deaths__and_Hospitalizations_2020201204.csv');
% data = DATA.data;
% ndays = 0; %%% Keeping the last 10 days for forecast.
% data = data(1:end-ndays,:);
% 
% %%% We shall delete the last 10 days.
% 
% t_actual = 0:size(data,1);
% 
% %%%% Smoothing the data - averaging every 7e consecutive days:
% data2 = data;
% for jj=4:size(data,1)-3
% for ii = 1:size(data,2) 
% data2(jj,ii) = mean(data(jj-3:jj+3,ii));
% end
% end
% t_span = datetime(2020,2,29) + caldays(0:length(t_actual)-1);
% 
% %%%% Total Population:
% N = 2705988;      % CHICAGO, USA
% 
% 
% %%%% Population proportion on each age range:
% Proportion = [23,17.2,17.4,12.6,11,9.7,5.9,3.3]/100;
% Population = [308979+310269,462964,469300,340219,294864,261009,157790,88565];
% PropInfections = 1E5*sum(data(:,4:11))./Population;
% PropInfections = PropInfections/max(PropInfections);
% factor = [1.185916802,0.890100477,0.819893104,1.047450903,1.201461658,1.112851255,1.096864419,1.503612944];
% PropInfections = PropInfections.*factor;
% PropHosp = sum(data(:,40:47))./sum(data(:,4:11));
% factor = [8.75,9.647540984,9.525974026,9.639269406,9.426380368,9.054830287,9.288888889,12.39267016];
% PropHosp = factor.*PropHosp;
% PropICU = ones(size(PropHosp));
% PropDeath = sum(data(:,22:29))./sum(data(:,40:47));
% factor = [4.79616307,4.8,4.736842105,4.609756098,4.638554217,4.662650602,4.718918919,4.790123457];
% PropDeath = factor.*PropDeath;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Initial parameters
% 
% I_M0 = 7;      % Potential initial infected and mild at t=0
% I_A0 = 0;
% E_0 = 0;       % initial exposed
% I_H0 = 0;      % initial infected and hospitalized at t=0
% I_I0 = 0;      % initial infected and in ICU at t=0
% R_0 = 0;       % initial recovered 
% D_0 = 0;       % initial dead
% 
% %  params is a structure used to pass parameters to the
% %   ODE solver
% 
% S_0 = N-(I_A0+I_M0+I_I0+I_H0+R_0+D_0+E_0);    % Suceptible pop.,  excluding initial infected 
% params.N =   N;  % N = total population
% 
% NumberOfAgeClasses = length(Proportion);  % total age ranges
% 
% yinit(1:NumberOfAgeClasses) = S_0*Proportion;
% yinit(NumberOfAgeClasses+1:2*NumberOfAgeClasses) = E_0*Proportion;
% yinit(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = I_A0*Proportion;
% yinit(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = I_M0*Proportion;
% yinit(4*NumberOfAgeClasses+1:5*NumberOfAgeClasses) = I_H0*Proportion;
% yinit(5*NumberOfAgeClasses+1:6*NumberOfAgeClasses) = I_I0*Proportion;
% yinit(6*NumberOfAgeClasses+1:7*NumberOfAgeClasses) = R_0*Proportion;
% yinit(7*NumberOfAgeClasses+1:8*NumberOfAgeClasses) = D_0*Proportion;
% yinit = yinit/N;
% 
% %--------------------------------------------------------------------------
% %%%% Asymptomatic
% p = ones-0.17;
% params.p = p; %%% proportio of symptomatic individuals
% q = 0.58;
% params.q = q;  %%% reduction in transmissibility
% 
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Model Parameters
% 
% params.sigma = 1./5.1;   % inverse of mean incubation time
% params.NumberOfAgeClasses = NumberOfAgeClasses;
% 
% %--------------------------------------------------------------------------
% % Mean time until recovery
% nu_M = 1./14;
% nu_H = 1./12;
% nu_I = 1./9;
% 
% %--------------------------------------------------------------------------
% % Mean time until death
% mu_M = zeros;
% mu_H = zeros;
% mu_I = 1./7;
% 
% %--------------------------------------------------------------------------
% % Mean time until passing to another infectious Class
% gamma_M = 1./1.2;
% gamma_H = 1./3.5;
% 
% %--------------------------------------------------------------------------
% % Proportion of individuals that will recover
% p_M  = (1-PropHosp);
% params.p_M = p_M; % In Mild conditions
% 
% p_H = 1-PropICU;
% params.p_H = p_H; % Hospitalized individuals
% 
% p_I = 1-PropDeath;
% params.p_I = p_I; % In ICU
% 
% %--------------------------------------------------------------------------
% % Proportion of individuals that will die
% q_M = zeros(1,NumberOfAgeClasses);
% params.q_M = q_M;  % In Mild conditions
% 
% q_H = zeros(1,NumberOfAgeClasses);
% params.q_H = q_H; % Hospitalized individuals
% 
% %--------------------------------------------------------------------------
% %%% RATES
% %%% Recovery Rate
% params.Recovery_A = 1/14;
% Recovery_M = ones-PropHosp;%(nu_M+mu_M+gamma_M)*p_M;
% Recovery_H = (nu_H+mu_H+gamma_H)*p_H;
% Recovery_I = ones-PropDeath;%(nu_I+mu_I)*p_I;
% 
% params.Recovery_M = Recovery_M; % in Mild conditions
% params.Recovery_H = ones-PropICU;%Recovery_H; % Hospitalized individuals
% params.Recovery_I = Recovery_I; % in ICU individuals
% 
% %%% Getting Worse Rate
% 
% GetWorse_M = PropHosp;%(nu_M+mu_M+gamma_M)*(1-p_M-q_M);
% GetWorse_H = PropICU;%(nu_H+mu_H+gamma_H)*(1-p_H-q_H);
% 
% params.GetWorse_M = GetWorse_M; % Mild conditions
% params.GetWorse_H = GetWorse_H; % Hospitalized individuals
% 
% %%% Death Rate
% 
% Death_M = (nu_M+mu_M+gamma_M)*q_M;
% Death_H = (nu_H+mu_H+gamma_H)*q_H;
% Death_I = PropDeath;%(nu_I+mu_I)*(1-p_I);
% 
% params.Death_M = Death_M; % in Mild conditions
% params.Death_H = Death_H; % Hospitalized individuals
% params.Death_I = Death_I; % in ICU individuals
% 
% %%%%%%%%
% Hosp = ones(length(t_actual),1);
% Hosp(3:end) = min(1,data2(2:end,3)./data2(1:end-1,1));
% % Hosp = min(2,Hosp/sum(GetWorse_M.*Proportion));
% % Hosp = min(3,Hosp/mean(GetWorse_M));
% Death = ones(size(t_actual));
% Death(3:end) = min(2,data2(2:end,2)./data2(1:end-1,3));
% % Death = min(3,Death/mean(Death_I));
% params.factorWorse = @(t)factorWorse(t,t_actual,Hosp);
% params.factorDeath = @(t)factorDeath(t,t_actual,Death);
% 
% %--------------------------------------------------------------------------
% % A priori value for the transmission parameter
% R_zero = 1.4*2.8/0.1782;
% gamma = 1/18;
% 
% %--------------------------------------------------------------------------
% % Transmission parameters
% beta = 2.2911*R_zero*gamma;
% bb = [0.5,0.25,0.20,0.15,0.1,0.05,0.025];
% beta_M = 0.5*diag(PropInfections);
% for jj = 1:NumberOfAgeClasses-1
% beta_M(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
% end
% beta_M = beta_M + beta_M';
% params.beta_M = beta_M;      % In Mild conditions
% params.beta_H = beta_M;  % Hospitalized individuals
% params.beta_I = beta_M; % In ICU
% 
% params.a = ones;
% params.b = 0.1;
% params.c = 0.01;
% 
% paramsOld = params;
% yinitOld = yinit;
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Estimating basic parameters from day 1 until 20:
% 
% day = 20;
% 
% LB1 = [1/N,1E-3];
% UB1 = [50/N,10];
% unknowns10=[1/N,beta];
% priors1 = unknowns10;
% 
% unknowns20 = [0.5,0.25,0.20,0.15,0.1,0.05,0.025];
% priors2 = unknowns20;
% LB2 = zeros(size(unknowns20));
% UB2 = ones(size(unknowns20));
% 
% unknowns0 = [unknowns10,unknowns20]; % initial parameters
% LB = [LB1,LB2]; % lower bounds
% UB = [UB1,UB2]; % upper bounds
% 
% priors = [priors1,priors2];
% 
% %%% Estimating the transmission constant parameters (M,H,I), the initial
% %%% infecve population (I_M0) and the transmission matrix:
% 
% OF = @(unknowns)ObjFun_InitialPopBetaM5(t_actual(1:day),params,...
% data2(1:day-1,:),options,priors,yinit,Proportion,PropInfections,unknowns);
% unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
% 
% I_M = unknowns(1);
% yinit(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = (1-p)*I_M*PropInfections;
% yinit(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = p*I_M*PropInfections;
% yinit(1:NumberOfAgeClasses) = Proportion-I_M*PropInfections;
% params.a = unknowns(2);
% 
% bb = unknowns(3:end);
% beta_M = 0.5*diag(PropInfections);
% for jj = 1:NumberOfAgeClasses-1
% beta_M(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
% end
% beta_M = beta_M + beta_M';
% params.beta_M = params.a.*beta_M;
% params.beta_H = params.a.*beta_M;
% params.beta_I = params.a.*beta_M;
% 
% PARAMS1 = unknowns;
% %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Estimating a time-dependent transmission parameter (beta) from day 1
% %%%% until 19:
% 
% BETA = zeros(length(t_actual),1);
% BETA(1) = unknowns(2); 
% unknowns0 = unknowns(2);
% priors = unknowns0;
% yinit2 = yinit;
% yb = zeros(length(t_actual),8*NumberOfAgeClasses);
% yb(1,:) = yinit;
% for jj =1:day-1
% t_actual2 = t_actual(jj:jj+1);
% OF2 = @(unknowns2)ObjFun_beta5(t_actual2,params,unknowns2,...
%                     data2(jj,:),options,priors,BETA(jj),yinit2);
% unknowns = lsqnonlin(OF2,unknowns0,1E-3,100,options3);
% unknowns0 = unknowns;
% priors = unknowns;
% BETA(jj+1) = unknowns;
% beta2 = @(t)interp1(t_actual2,BETA(jj:jj+1),t);
% [~,y2] = ode45(@(t,y)seir_death_age_beta3(t,y,params,beta2),...
%                                                  t_actual2,yinit2,options);
% yinit2 = y2(end,:);
% yb(jj+1,:) = yinit2;
% disp(num2str(unknowns))
% end
% 
% % R0 = zeros(1,length(t_actual)-1);
% % for jj=1:day-1
% % % R0(jj) = basic_reproduction_rate_beta2(Proportion,params,BETA(jj+1),t_actual(jj+1))';
% % R0(jj) = basic_reproduction_rate_beta2(ones(size(Proportion)),params,BETA(jj+1),t_actual(jj+1))';
% % end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Estimating basic parameters from day 21 until 80:
% 
% 
% dt = 20;
% nmax = floor(length(t_actual)/dt);
% PARAMS2 = zeros(nmax,8);
% for ss = 1:nmax-1
% 
% LB1 = 1E-3;
% UB1 = 10;
% unknowns10 = beta;
% priors1 = unknowns10;
% 
% unknowns20 = [0.5,0.25,0.20,0.15,0.1,0.05,0.025];
% priors2 = unknowns20;
% LB2 = zeros(size(unknowns20));
% UB2 = ones(size(unknowns20));
% 
% unknowns0 = [unknowns10,unknowns20]; % initial parameters
% LB = [LB1,LB2]; % lower bounds
% UB = [UB1,UB2]; % upper bounds
% 
% priors = [priors1,priors2];
% yinit2 = y2(end,:);
% 
% %%% Estimating the transmission constant parameters (M,H,I), the initial
% %%% infecve population (I_M0) and the transmission matrix:
% if ss == nmax-1
% aux = (day-1)+(ss-1)*dt+1:length(t_actual)-1;
% aux2 = (day-1)+(ss-1)*dt+1:length(t_actual);    
% else
% aux = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt;
% aux2 = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt+1;
% end
% 
% OF = @(unknowns)ObjFun_BetaM5(t_actual(aux2),params,...
%            data2(aux,:),options,priors,yinit2,unknowns,PropInfections);
% unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
% params.a = unknowns(1);
% 
% 
% PARAMS2(ss,:) = unknowns;
% 
% bb = unknowns(2:end);
% beta_M = 0.5*diag(PropInfections);
% for jj = 1:NumberOfAgeClasses-1
% beta_M(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
% end
% beta_M = beta_M + beta_M';
% params.beta_M = params.a.*beta_M;
% params.beta_H = params.a.*beta_M;
% params.beta_I = params.a.*beta_M;
% 
% unknowns0 = ones;
% priors = unknowns0;
% for jj = aux
% t_actual2 = t_actual(jj:jj+1);
% OF2 = @(unknowns)ObjFun_beta5(t_actual2,params,unknowns,...
%                   data2(jj,:),options,priors,BETA(jj),yinit2);
% unknowns = lsqnonlin(OF2,unknowns0,1E-2,10,options2);
% unknowns0 = unknowns;
% priors = unknowns;
% BETA(jj+1,:) = unknowns;
% beta2 = @(t)interp1(t_actual2,BETA(jj:jj+1),t);
% [~,y2] = ode45(@(t,y)seir_death_age_beta3(t,y,params,beta2),...
%                                                  t_actual2,yinit2,options);
% yinit2 = y2(end,:);
% yb(jj+1,:) = yinit2;
% disp(num2str(unknowns))
% % R0(jj) = basic_reproduction_rate_beta(Proportion,params,BETA(jj+1),t_actual(jj+1));
% end
% end
% 
% factor = zeros(length(t_actual),1);
% factorD = factor;
% factorW = @(t)factorWorse(t,t_actual,Hosp);
% factorDea = @(t)factorDeath(t,t_actual,Death);
% for jj = 1:length(t_actual)
% factor(jj) = factorW(t_actual(jj));
% factorD(jj) = factorDea(t_actual(jj));
% end
% 
% %%%% Final Number of Cases for each Age Range
% sigma = params.sigma;
% NewCases = p*sigma*sum(yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses),2)*N;
% 
% %%%% Total Number of Deaths for each day
% 
% NewDeaths = zeros(size(yb(:,1)));
% NewHospB = zeros(size(yb(:,1),1),NumberOfAgeClasses);
% NewHosp = zeros(size(yb(:,1)));
% Deaths = zeros(length(yb(:,1)),NumberOfAgeClasses);
% for jj=1:NumberOfAgeClasses
% NewDeaths = NewDeaths + Death_M(jj)*yb(:,3*NumberOfAgeClasses+jj)...
% + Death_H(jj)*yb(:,4*NumberOfAgeClasses+jj)...
% + Death_I(jj)*factorD.*yb(:,5*NumberOfAgeClasses+jj);
% NewHosp = NewHosp + GetWorse_M(jj)*factor.*yb(:,3*NumberOfAgeClasses+jj);
% NewHospB(:,jj) = GetWorse_M(jj)*factor.*yb(:,3*NumberOfAgeClasses+jj);
% Deaths(:,jj) = Death_M(jj)*yb(:,3*NumberOfAgeClasses+jj)...
% +Death_H(jj)*yb(:,4*NumberOfAgeClasses+jj)...
% +Death_I(jj)*factorD.*yb(:,5*NumberOfAgeClasses+jj);
% end
% 
% 
% %%%% Total Hospitalizations for each Age Range
% TotalHosp = sum(NewHospB)'*N; 
% TotalInfections = ...
%       sum(p*sigma*yb(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses))'*N;
% TotalDeaths = sum(Deaths)'*N;
% 
% disp('  Infections   Hosp.    Deaths  ')
% disp(num2str(round([TotalInfections,TotalHosp,TotalDeaths])))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Bootstraping_20201215;
% % save dataChicago_20201216;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%% CI Evaluation
% % %%% R0
% % aux = sort(R0Boot);
% % aux2 = round(0.05*NSamples);
% % aux = aux(aux2+1:end-aux2,:);
% % CI95R0 = [min(aux);max(aux)];
% % 
% %%% Beta
% aux = sort(BETABoot);
% aux2 = round(0.05*NSamples);
% aux = aux(aux2+1:end-aux2,:);
% CI95BETA = [min(aux);max(aux)];
% 
% %%% New Cases:
% aux = sort(NewCasesBoot');
% aux2 = round(0.05*NSamples);
% aux = aux(aux2+1:end-aux2,:);
% CI95NewCases = [min(aux);max(aux)]*N;
% 
% %%% New Deaths:
% aux = sort(NewDeathsBoot');
% aux2 = round(0.05*NSamples);
% aux = aux(aux2+1:end-aux2,:);
% CI95NewDeaths = [min(aux);max(aux)]*N;
% 
% %%% New Hospitalizations:
% aux = sort(NewHospBoot');
% aux2 = round(0.05*NSamples);
% aux = aux(aux2+1:end-aux2,:);
% CI95NewHosp = [min(aux);max(aux)]*N;
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % aux = sort(PARAMS1Boot);
% % aux2 = round(0.05*NSamples);
% % aux = aux(aux2+1:end-aux2,:);
% % CI95PARAMS1 = [min(aux);max(aux)];
% % 
% % aux = sort(PARAMS2Boot);
% % aux2 = round(0.05*NSamples);
% % aux = aux(aux2+1:end-aux2,:);
% % CI95PARAMS2 = [min(aux);max(aux)];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% disp(['Elapsed Time: ',num2str(toc),' seconds.'])
% % R02 = R0;
% % for jj=4:size(data,1)-3
% % R02(jj) = mean(R02(jj-3:jj+3));
% % end
% % figure
% % hold on
% % grid on
% % box on
% % title('Basic Reproduction Parameter')
% % % area(t_span,[0,CI95R0(2,:)],'linestyle',':','FaceColor',[255,160,122]/255)
% % % area(t_span,[0,CI95R0(1,:)],'linestyle',':','FaceColor',[1,1,1])
% % plot(t_span,ones(size(t_span)),'k')
% % plot(t_span,[0,R0],'r','LineWidth',2)
% % ylabel('R(t)')
% % xlim([t_span(1),t_span(end)])
% % ylim([0,1.3*max(R0)])
% % plot([t_span(day),t_span(day)],[0,1.3*max(R0)],'--k','LineWidth',2)
% % xtickformat('dd-MMM')
% % set(gca,'FontSize',16,'FontName','Arial')
% % hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting Results:

figure
hold on
box on
title('Daily New Infections')
% h1=area(t_span,CI95NewCases(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
% h2=area(t_span,CI95NewCases(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data(:,1),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
% plot(t_span(2:end),data2(:,1),'b','LineWidth',2)
plot(t_span,NewCases,'r','LineWidth',2)
plot([t_span(day),t_span(day)],[0,1.3*max(data(:,1))],'--k','LineWidth',2)
% h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Reported','Estimated')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end)])
ylim([0,1.3*max(data(:,1))])
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
hold off


figure
hold on
box on
title('Daily New Deaths')
% h1=area(t_span,CI95NewDeaths(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
% h2=area(t_span,CI95NewDeaths(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data(:,2),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
plot(t_span,N*NewDeaths,'r','LineWidth',2)
% plot(t_span(2:end),data2(:,4),'b','LineWidth',2)
% h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
legend('Reported','Estimated')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end)])
ylim([0,1.3*max(data(:,4))])
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
hold off

%%%%%%%
figure
hold on
box on
title('Daily New Hospitalizations')
% h1=area(t_span,CI95NewHosp(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
% h2=area(t_span,CI95NewHosp(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data(:,3),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
% plot(t_span(2:end),data2(:,3),'b','LineWidth',2)
% h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t_span,N*NewHosp,'r','LineWidth',2)
plot([t_span(day),t_span(day)],[0,1.3*max(data(:,3))],'--k','LineWidth',2)
legend('Reported','Estimated')
ylabel('Number of Individuals')
xlim([t_span(1),t_span(end)])
ylim([0,1.3*max(data(:,3))])
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
hold off

figure
hold on
box on
title('Transmission Parameter')
% area(t_span,CI95BETA(2,:),'linestyle',':','FaceColor',[255,160,122]/255)
% area(t_span,CI95BETA(1,:),'linestyle',':','FaceColor',[1,1,1])
plot(t_span,BETA,'r','LineWidth',2)
ylabel('\beta(t)')
xlim([t_span(1),t_span(end)])
ylim([0,1.3*max(BETA)])
plot([t_span(day),t_span(day)],[0,1.3*max(BETA)],'--k','LineWidth',2)
xtickformat('dd-MMM')
set(gca,'FontSize',16,'FontName','Arial')
hold off
