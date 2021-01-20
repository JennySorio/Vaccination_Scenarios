%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bootstrap Sampling
% load('data');
NSamples = 200;
trajectories = zeros(length(data2(:,1)),NSamples);
for jj = 1:NSamples
trajectories(:,jj) = poissrnd(NewCases(2:end));
end


PARAMS1Boot = zeros(NSamples,9);
BETABoot = zeros(NSamples,length(BETA));
yinit2B = zeros(NSamples,7*NumberOfAgeClasses);

PARAMS2Boot = zeros(nmax,8,NSamples);
t_actualB = t_actual;

R0Boot = zeros(NSamples,length(t_actual)-1);
TotalHospBoot = zeros(NumberOfAgeClasses,NSamples); 
TotalInfectionsBoot = zeros(NumberOfAgeClasses,NSamples);
TotalDeathsBoot = zeros(NumberOfAgeClasses,NSamples);

NewCasesBoot = zeros(length(t_actual),NSamples);
NewDeathsBoot = zeros(length(t_actual),NSamples);
NewHospBoot = zeros(length(t_actual),NSamples);

for ll = 1:1%NSamples
% %%% Estimating the transmission constant parameters (M,H,I), the initial
% %%% infecve population (I_M0) and the transmission matrix:
params2=paramsOld;
yinitB = yinitOld;
t_actual2 = t_actual(1:day);
LB1 = [1/N,1E-3];
UB1 = [50/N,10];
unknowns10=[1/N,beta];
priors1 = unknowns10;

unknowns20 = [0.5,0.25,0.20,0.15,0.1,0.05,0.025];
priors2 = unknowns20;
LB2 = zeros(size(unknowns20));
UB2 = ones(size(unknowns20));

unknowns0 = [unknowns10,unknowns20]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds

priors = [priors1,priors2];

OF = @(unknowns)ObjFun_InitialPopBetaM5(t_actual2,params2,...
trajectories(1:day-1,ll),options,priors,yinitB,Proportion,PropInfections,unknowns);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);
I_M = unknowns(1);
yinitB(2*NumberOfAgeClasses+1:3*NumberOfAgeClasses) = (1-p)*I_M*PropInfections;
yinitB(3*NumberOfAgeClasses+1:4*NumberOfAgeClasses) = p*I_M*PropInfections;
yinitB(1:NumberOfAgeClasses) = Proportion-I_M*PropInfections;
params2.a = unknowns(2);
bb = unknowns(3:end);
beta_M = 0.5*diag(PropInfections);
for jj = 1:NumberOfAgeClasses-1
beta_M(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
end
beta_M = beta_M+beta_M';
params2.beta_M = params2.a.*beta_M;
params2.beta_H = params2.a.*beta_M;
params2.beta_I = params2.a.*beta_M;

PARAMS1Boot(ll,:) = unknowns;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Estimating a time-dependent transmission parameter (beta) from day 1
%%%% until 19:

Beta = zeros(length(t_actual),1);
Beta(1) = unknowns(2); 
unknowns0 = unknowns(2);
priors = unknowns0;
yinit2 = yinitB;
yb2 = zeros(length(t_actual),8*NumberOfAgeClasses);
yb2(1,:) = yinit;
for jj =1:day-1
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta5(t_actual2,params2,unknowns2,...
                    trajectories(jj,ll),options,priors,Beta(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-3,50,options3);
unknowns0 = unknowns;
priors = unknowns;
Beta(jj+1) = unknowns;
beta2 = @(t)interp1(t_actual2,Beta(jj:jj+1),t);
[~,y2B] = ode45(@(t,y)seir_death_age_beta3(t,y,params2,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2B(end,:);
yb2(jj+1,:) = yinit2;
disp(num2str(unknowns))
end

% AUX = zeros(size(t_actual2));
% for jj=1:day-1
% % AUX(jj) = basic_reproduction_rate_beta2(Proportion,params2,...
% %                                                Beta(jj+1),t_actual(jj+1))';
% AUX(jj) = basic_reproduction_rate_beta2(ones(size(Proportion)),params2,...
%                                                Beta(jj+1),t_actual(jj+1))';
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:

for ss = 1:nmax-1

LB1 = 1E-3;
UB1 = 10;
unknowns10 = beta;
priors1 = unknowns10;

unknowns20 = [0.5,0.25,0.20,0.15,0.1,0.05,0.025];
priors2 = unknowns20;
LB2 = zeros(size(unknowns20));
UB2 = ones(size(unknowns20));

unknowns0 = [unknowns10,unknowns20]; % initial parameters
LB = [LB1,LB2]; % lower bounds
UB = [UB1,UB2]; % upper bounds

priors = [priors1,priors2];
yinit2 = y2B(end,:);

%%% Estimating the transmission constant parameters (M,H,I), the initial
%%% infecve population (I_M0) and the transmission matrix:
if ss == nmax-1
aux = (day-1)+(ss-1)*dt+1:length(t_actual)-1;
aux2 = (day-1)+(ss-1)*dt+1:length(t_actual);    
else
aux = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt;
aux2 = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt+1;
end

OF = @(unknowns)ObjFun_BetaM5(t_actualB(aux2),params2,...
trajectories(aux,ll),options,priors,yinit2,unknowns,PropInfections);
unknowns = lsqnonlin(OF,unknowns0,LB,UB,options2);

params2.a = unknowns(1);
PARAMS2Boot(ss,:,ll) = unknowns;

bb = unknowns(2:end);
beta_M = 0.5*diag(PropInfections);
for jj = 1:NumberOfAgeClasses-1
beta_M(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
end
beta_M = beta_M + beta_M';
params2.beta_M = params2.a.*beta_M;
params2.beta_H = params2.a.*beta_M;
params2.beta_I = params2.a.*beta_M;

unknowns0 = ones;
priors = unknowns0;
for jj = aux
t_actual2 = t_actual(jj:jj+1);
OF2 = @(unknowns2)ObjFun_beta5(t_actual2,params2,unknowns2,...
                    trajectories(jj,ll),options,priors,Beta(jj),yinit2);
unknowns = lsqnonlin(OF2,unknowns0,1E-2,10,options2);
unknowns0 = unknowns;
priors = unknowns;
Beta(jj+1,:) = unknowns;
beta2 = @(t)interp1(t_actual2,Beta(jj:jj+1),t);
[~,y2B] = ode45(@(t,y)seir_death_age_beta3(t,y,params2,beta2),...
                                                 t_actual2,yinit2,options);
yinit2 = y2B(end,:);
yb2(jj+1,:) = yinit2;
disp(num2str(unknowns))
% R0(jj) = basic_reproduction_rate_beta(Proportion,params,BETA(jj+1),t_actual(jj+1));
end
end
BETABoot(ll,:) = Beta;
% for jj=day:length(t_actual)-1
% % AUX(jj) = basic_reproduction_rate_beta2(Proportion,params2,...
% %                                                 Beta(jj+1),t_actual(jj+1));
% AUX(jj) = basic_reproduction_rate_beta2(ones(size(Proportion)),params2,...
%                                                 Beta(jj+1),t_actual(jj+1));
% end
% R0Boot(ll,:) = AUX;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Evaluating Curves of New Infections, Hospitalizations etc.
NewCasesBoot(:,ll) = ...
             p*sigma*sum(yb2(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses),2);

factor = zeros(length(t_actual),1);
factorD = factor;
factorW = params2.factorWorse;
factorDeath = params2.factorDeath;
for jj = 1:length(t_actual)
factor(jj) = factorW(t_actual(jj));
factorD(jj) = factorDeath(t_actual(jj));
end          

%%% Total Number of Deaths for each day
NewHospBb = zeros(size(yb2(:,1),1),NumberOfAgeClasses);
Deathsb = zeros(length(yb2(:,1)),NumberOfAgeClasses);
for jj=1:NumberOfAgeClasses
NewDeathsBoot(:,ll) = NewDeathsBoot(:,ll)...
+ Death_M(jj)*yb2(:,3*NumberOfAgeClasses+jj)...
+ Death_H(jj)*yb2(:,4*NumberOfAgeClasses+jj)...
+ Death_I(jj)*factorD.*yb2(:,5*NumberOfAgeClasses+jj);
NewHospBoot(:,ll) = NewHospBoot(:,ll)...
+ GetWorse_M(jj)*factor.*yb2(:,3*NumberOfAgeClasses+jj);
NewHospBb(:,jj) = GetWorse_M(jj)*factor.*yb2(:,3*NumberOfAgeClasses+jj);
Deathsb(:,jj) = Death_M(jj)*yb2(:,3*NumberOfAgeClasses+jj)...
+Death_H(jj)*yb2(:,4*NumberOfAgeClasses+jj)...
+Death_I(jj)*factorD.*yb2(:,5*NumberOfAgeClasses+jj);
end

%%% Total Hospitalizations for each Age Range
TotalHospBoot(:,ll) = sum(NewHospBb)'*N; 
TotalInfectionsBoot(:,ll) = ...
      sum(p*sigma*yb2(:,NumberOfAgeClasses+1:2*NumberOfAgeClasses))'*N;
TotalDeathsBoot(:,ll) = sum(Deathsb)'*N;
end
