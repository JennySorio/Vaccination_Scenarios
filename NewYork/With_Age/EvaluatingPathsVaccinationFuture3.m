
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluating BootVacstrap predictions

Number = NumberOfAgeClasses;
% NSamples = 1;

TotalHospBootVac = zeros(Number,NSamples); 
TotalInfectionsBootVac = zeros(Number,NSamples);
TotalDeathsBootVac = zeros(Number,NSamples);

Vaccinated = zeros(length(t_actual2),NSamples); 
NewCasesBootVac = zeros(length(t_actual2),NSamples);
NewDeathsBootVac = zeros(length(t_actual2),NSamples);
NewHospBootVac = zeros(length(t_actual2),NSamples);


parfor ll = 1:NSamples
params2 = params;
p = params2.p;

factor = zeros(length(t_actual2),1);
factorD = factor;
factorW = params2.factorWorse;
factorDeath = params2.factorDeath;
for jj = 1:length(t_actual2)
factor(jj) = factorW(t_actual2(jj));
factorD(jj) = factorDeath(t_actual2(jj));
end  

% %%% Estimating the transmission constant parameters (M,H,I), the initial
% %%% infecve population (I_M0) and the transmission matrix:
yinitB = yinitOld;
I_M = PARAMS1Boot(ll,1);
yinitB = zeros(1,9*Number);
for jj = 3:9
yinitB((jj-1)*Number+1:jj*Number) = yinit((jj-2)*Number+1:(jj-1)*Number);
end
yinitB(Number+1:2*Number) = zeros;
yinitB(3*Number+1:4*Number) = (1-p)*I_M*PropInfections;
yinitB(4*Number+1:5*Number) = p*I_M*PropInfections;
yinitB(1:Number) = Proportion-I_M*PropInfections;
params2.a = PARAMS1Boot(ll,2);

bb = PARAMS1Boot(ll,3:size(PARAMS1Boot,2));
beta_M = 0.5*diag(PropInfections);
for jj = 1:NumberOfAgeClasses-1
beta_M(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
end
beta_M = beta_M+beta_M';
params2.beta_M = params2.a.*beta_M;
params2.beta_H = params2.a.*beta_M;
params2.beta_I = params2.a.*beta_M;

Beta = BETABoot(ll,:);
Beta2 = [Beta(1:end-l),mean(Beta(end-l-9:end-l))*ones(1,ndays+1)];
beta2 = @(t)interp1(t_actual2,Beta2,t);
yinit2 = yinitB;
yb2 = zeros(length(t_actual2),9*Number);
yb2(1,:) = yinitB;
auxA = zeros(1,length(Beta2)-1);
for jj =1:day-1
t_actualB = t_actual2(jj:jj+1);
% beta2 = @(t)interp1(t_actualB,Beta2(jj:jj+1),t);
[~,y2B] = ode45(@(t,y)seir_Vaccination(t,y,params2,beta2),...
                                                 t_actualB,yinit2,options);
yinit2 = y2B(end,:);
yb2(jj+1,:) = yinit2;
end

for ss = 1:nmax-1

if ss == nmax-1
aux = (day-1)+(ss-1)*dt+1:length(t_actual2)-1;
aux2 = (day-1)+(ss-1)*dt+1:length(t_actual2);    
else
aux = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt;
aux2 = (day-1)+(ss-1)*dt+1:(day-1)+ss*dt+1;
end
params2.a = PARAMS2Boot(ss,1,ll);
bb = PARAMS2Boot(ss,2:size(PARAMS2Boot,2),ll);
beta_M = 0.5*diag(PropInfections);
for jj = 1:NumberOfAgeClasses-1
beta_M(jj,jj+1:end) = PropInfections(jj)*bb(1:end-jj+1);
end
beta_M = beta_M+beta_M';
params2.beta_M = params2.a.*beta_M;
params2.beta_H = params2.a.*beta_M;
params2.beta_I = params2.a.*beta_M;

BETAB = Beta2;
value1 = mean(Beta(end-l-9:end-l));
value2 = mean(Beta(200-9:200));

for jj = aux
t_actualB = t_actual2(jj:jj+1);
[~,y2] = ode45(@(t,y)seir_Vaccination(t,y,params2,beta2),...
                                                 t_actualB,yinit2,options);
yinit2 = y2(end,:);
yb2(jj+1,:) = yinit2;
end
end


NewCasesBootVac(:,ll) = ...
    sigma*sum(yb2(:,2*Number+1:3*Number),2);

% %%%% Vaccinated Population:
Vaccinationb = Vaccination;
for jj = 1:length(t_actual2)
    if sum(yb2(jj,1:Number)) < 0.4
    Vaccinationb(:,jj) = zeros;
    end
end

for jj=1:Number
Vaccinated(:,ll) = Vaccinated(:,ll) + ...
    Vaccinationb(jj,:)'.*yb2(:,jj);
NewDeathsBootVac(:,ll) = NewDeathsBootVac(:,ll)...
+ Death_M(jj)*yb2(:,4*Number+jj)...
+ Death_H(jj)*yb2(:,5*Number+jj)...
+ Death_I(jj)*factorD.*yb2(:,6*Number+jj);
NewHospBootVac(:,ll) = NewHospBootVac(:,ll)...
+ GetWorse_M(jj)*factor.*yb2(:,4*Number+jj);
end
end
