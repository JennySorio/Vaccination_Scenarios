clear all; clc; close all; format long e;
load dataNYC_20201211;

ndays = 183;
l=3;
z=0;
t_actual2 = 0:length(t_actual)-l+ndays;
t_span = datetime(2020,3,01) + caldays(0:length(t_actual2)-1);


Hosp = ones(length(t_actual2),1);
Hosp(3:length(t_actual)) = min(1,data1(2:end,3)./data1(1:end-1,1));
Hosp(length(t_actual)-l:end) = mean(Hosp(length(t_actual)-l-z:length(t_actual)-l))*ones;
Hosp = min(20,Hosp/GetWorse_M);
Death = ones(size(t_actual2));
Death(3:length(t_actual)) = min(1,data1(2:end,4)./data1(1:end-1,3))/Death_I;
Death(length(t_actual)-l:end) = mean(Death(length(t_actual)-l-z:length(t_actual)-l))*ones;
params.factorDeath = @(t)interp1(t_actual2,Death,t);
params.factorWorse = @(t)interp1(t_actual2,Hosp,t);


tt = [215,246,276,307,338,366,397];
ZZ = zeros(length(tt),12);
disp('Cases    Hospitalizations         Deaths    Vaccinated')
for zz = 1:length(tt)
Vaccination = zeros(size(t_actual2));
Vaccination(t_actual2>=tt(zz)) = 0.95*0.01*ones;
params.VaccinationRate = @(t)interp1(t_actual2,Vaccination,t);
EvaluatingPathsVaccinationFuture3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux = sort(Vaccinated');
aux2 = round(0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95Vac = [min(aux);max(aux)]*N;


%%% New Cases:
aux = sort(NewCasesBootVac');
aux2 = round(0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewCasesVac = [min(aux);max(aux)]*N;

%%% New Deaths:
aux = sort(NewDeathsBootVac');
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewDeathsVac = [min(aux);max(aux)]*N;

%%% New Hospitalizations:
aux = sort(NewHospBootVac');
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95NewHospVac = [min(aux);max(aux)]*N;



disp(['Total Vaccinated: ',num2str(round([sum(median(N*Vaccinated(tt(zz):end,:),2)),sum(CI95Vac(1,tt(zz):end)),sum(CI95Vac(2,tt(zz):end))]))])
disp(num2str(round([sum(median(N*NewCasesBootVac(tt(1):end,:),2)),sum(CI95NewCasesVac(1,tt(1):end)),sum(CI95NewCasesVac(2,tt(1):end)),...
    sum(median(N*NewHospBootVac(tt(1):end,:),2)),sum(CI95NewHospVac(1,tt(1):end)),sum(CI95NewHospVac(2,tt(1):end)),...
    sum(median(N*NewDeathsBootVac(tt(1):end,:),2)),sum(CI95NewDeathsVac(1,tt(1):end)),sum(CI95NewDeathsVac(2,tt(1):end)),...
    sum(median(N*Vaccinated(tt(zz):end,:),2)),sum(CI95Vac(1,tt(zz):end)),sum(CI95Vac(2,tt(zz):end))])))
ZZ(zz,:) = round([sum(median(N*NewCasesBootVac(tt(1):end,:),2)),sum(CI95NewCasesVac(1,tt(1):end)),sum(CI95NewCasesVac(2,tt(1):end)),...
    sum(median(N*NewHospBootVac(tt(1):end,:),2)),sum(CI95NewHospVac(1,tt(1):end)),sum(CI95NewHospVac(2,tt(1):end)),...
    sum(median(N*NewDeathsBootVac(tt(1):end,:),2)),sum(CI95NewDeathsVac(1,tt(1):end)),sum(CI95NewDeathsVac(2,tt(1):end)),...
    sum(median(N*Vaccinated(tt(zz):end,:),2)),sum(CI95Vac(1,tt(zz):end)),sum(CI95Vac(2,tt(zz):end))]);
end


AA = round(ZZ(:,7:end));

H = [100 100 600 300];
figure
hold on
box on
title('New York City')
area(t_span(tt),AA(:,3),'linestyle',':','FaceColor',[255,160,122]/255);
area(t_span(tt),AA(:,2),'linestyle',':','FaceColor',[1,1,1]);
plot(t_span(tt),AA(:,1),'-sr','LineWidth',2)
ylabel('Number of Individuals')
xlabel('Starting Date')
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'AccumDeathVacNYCB.fig');
print('-dpng','AccumDeathVacNYCB');

figure
hold on
box on
title('New York City')
area(t_span(tt(2:end)),diff(AA(:,3)),'linestyle',':','FaceColor',[255,160,122]/255);
area(t_span(tt(2:end)),diff(AA(:,2)),'linestyle',':','FaceColor',[1,1,1]);
plot(t_span(tt(2:end)),diff(AA(:,1)),'-sr','LineWidth',2)
ylabel('Number of Individuals')
xlabel('Starting Date')
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'NewDeathVacNYCB.fig');
print('-dpng','NewDeathVacNYCB');


disp(num2str([sum(data(tt(1)-1:end,1)),sum(data(tt(1)-1:end,3)),sum(data(tt(1)-1:end,4))]));
disp(num2str([sum(data(tt(1)-1:tt(2)-1,1)),sum(data(tt(1)-1:tt(2)-1,3)),sum(data(tt(1)-1:tt(2)-1,4))]));
