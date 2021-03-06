clear all; clc; close all; format long e;
load dataChicago_20201204b;

ndays = 183;%203;
l=3;
z=0;
t_actual2 = 0:length(t_actual)-l+ndays;
t_span = datetime(2020,3,01) + caldays(0:length(t_actual2)-1);


Hosp = ones(length(t_actual2),1);
Hosp(3:length(t_actual)) = min(1,data1(2:end,3)./data1(1:end-1,1));
Hosp(length(t_actual)-l:end) = mean(Hosp(length(t_actual)-l-z:length(t_actual)-l))*ones;
Hosp = min(20,Hosp/GetWorse_M);
Death = ones(size(t_actual2));
Death(3:length(t_actual)) = min(1,data1(2:end,2)./data1(1:end-1,3))/Death_I;
Death(length(t_actual)-l:end) = mean(Death(length(t_actual)-l-z:length(t_actual)-l))*ones;
params.factorDeath = @(t)interp1(t_actual2,Death,t);
params.factorWorse = @(t)interp1(t_actual2,Hosp,t);


EvaluatingPathsFuture;

tt = [215,246,276,307,338,366,397];
ZZ = zeros(length(tt),9);
for zz = 1:length(tt)
Vaccination = zeros(size(t_actual2));
Vaccination(t_actual2>=tt(zz)) = 0.95*0.01*ones;
params.VaccinationRate = @(t)interp1(t_actual2,Vaccination,t);
EvaluatingPathsVaccinationFuture2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux = sort(Vaccinated');
aux2 = round(0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95Vac = [min(aux);max(aux)]*N;

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

H = [100 100 600 300];
%%% plotting Results:
figure
hold on
% grid on
box on
title('Daily Infections - Chicago')
h1=area(t_span,CI95NewCases(2,:),'linestyle',':','FaceColor',[255,160,122]/255,'FaceAlpha',1);%[51,236,255]/255);
h2=area(t_span,CI95NewCases(1,:),'linestyle',':','FaceColor',[1,1,1]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
h1=area(t_span,CI95NewCasesVac(2,:),'linestyle',':','FaceColor',[51,236,255]/255,'FaceAlpha',0.3);
h2=area(t_span,CI95NewCasesVac(1,:),'linestyle',':','FaceColor',[1,1,1]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t_span,median(NewCasesBoot*N,2),'r','LineWidth',2)
% plot(t_span(2:end),data1(:,1),'r','LineWidth',2)
plot(t_span,median(NewCasesBootVac*N,2),'b','LineWidth',2)
legend('Without Vaccination','With Vaccination','Location','NorthEast')
ylabel('Number of Individuals')
xlim([t_span(tt(zz)-10),t_span(end)])
h1=plot([t_span(tt(zz)),t_span(tt(zz))],[0,1.3*max(median(NewCasesBoot(tt(zz)-10:end,:)*N,2))],'k','LineWidth',1);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylim([0,1.3*max(median(NewCasesBoot(tt(zz)-10:end,:)*N,2))])
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['InfectionsChicagoVac',num2str(zz),'.fig']);
print('-dpng',['InfectionsChicagoVac2',num2str(zz)]);

figure
hold on
% grid on
box on
title('Daily Deaths - Chicago')
h1=area(t_span,CI95NewDeaths(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewDeaths(1,:),'linestyle',':','FaceColor',[1,1,1]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
h1=area(t_span,CI95NewDeathsVac(2,:),'linestyle',':','FaceColor',[51,236,255]/255);
h2=area(t_span,CI95NewDeathsVac(1,:),'linestyle',':','FaceColor',[1,1,1]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t_span,median(NewDeathsBoot*N,2),'r','LineWidth',2)
% plot(t_span(2:end),data1(:,2),'r','LineWidth',2)
plot(t_span,median(NewDeathsBootVac*N,2),'b','LineWidth',2)
legend('Without Vaccination','With Vaccination','Location','NorthEast')
ylabel('Number of Individuals')
xlim([t_span(tt(zz)-10),t_span(end)])
h1=plot([t_span(tt(zz)),t_span(tt(zz))],[0,1.3*max(median(NewDeathsBoot(tt(zz)-10:end,:)*N,2))],'k','LineWidth',1);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylim([0,1.3*max(median(NewDeathsBoot(tt(zz)-10:end,:)*N,2))])
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['DeathsChicagoVac2',num2str(zz),'.fig']);
print('-dpng',['DeathsChicagoVac2',num2str(zz)]);

%%%%%%%%%
figure
hold on
% grid on
box on
title('Daily Hospitalizations - Chicago')
h1=area(t_span,CI95NewHosp(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewHosp(1,:),'linestyle',':','FaceColor',[1,1,1]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
h1=area(t_span,CI95NewHospVac(2,:),'linestyle',':','FaceColor',[51,236,255]/255);
h2=area(t_span,CI95NewHospVac(1,:),'linestyle',':','FaceColor',[1,1,1]);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
plot(t_span,median(NewHospBoot*N,2),'r','LineWidth',2)
% plot(t_span(2:end),data1(:,3),'r','LineWidth',2)
plot(t_span,median(NewHospBootVac*N,2),'b','LineWidth',2)
legend('Without Vaccination','With Vaccination','Location','NorthEast')
ylabel('Number of Individuals')
xlim([t_span(tt(zz)-10),t_span(end)])
h1=plot([t_span(tt(zz)),t_span(tt(zz))],[0,1.3*max(median(NewHospBoot(tt(zz)-10:end,:)*N,2))],'k','LineWidth',1);
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylim([0,1.3*max(median(NewHospBoot(tt(zz)-10:end,:)*N,2))])
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['HospChicagoVac2',num2str(zz),'.fig']);
print('-dpng',['HospChicagoVac2',num2str(zz)]);


disp(['Total Vaccinated: ',num2str(round([sum(median(N*Vaccinated(tt(zz):end,:),2)),sum(CI95Vac(1,tt(zz):end)),sum(CI95Vac(2,tt(zz):end))]))])
disp('Cases    Hospitalizations         Deaths')
disp(num2str(round([sum(median(N*NewCasesBootVac(tt(1):end,:),2)),sum(CI95NewCasesVac(1,tt(1):end)),sum(CI95NewCasesVac(2,tt(1):end)),...
    sum(median(N*NewHospBootVac(tt(1):end,:),2)),sum(CI95NewHospVac(1,tt(1):end)),sum(CI95NewHospVac(2,tt(1):end)),...
    sum(median(N*NewDeathsBootVac(tt(1):end,:),2)),sum(CI95NewDeathsVac(1,tt(1):end)),sum(CI95NewDeathsVac(2,tt(1):end))])))
ZZ(zz,:) = round([sum(median(N*NewCasesBootVac(tt(1):end,:),2)),sum(CI95NewCasesVac(1,tt(1):end)),sum(CI95NewCasesVac(2,tt(1):end)),...
    sum(median(N*NewHospBootVac(tt(1):end,:),2)),sum(CI95NewHospVac(1,tt(1):end)),sum(CI95NewHospVac(2,tt(1):end)),...
    sum(median(N*NewDeathsBootVac(tt(1):end,:),2)),sum(CI95NewDeathsVac(1,tt(1):end)),sum(CI95NewDeathsVac(2,tt(1):end))]);
end


AA = round(ZZ(:,7:end));

%%%% Accumulated deaths
figure
hold on
box on
title('Chicago')
area(t_span(tt),AA(:,3),'linestyle',':','FaceColor',[255,160,122]/255);
area(t_span(tt),AA(:,2),'linestyle',':','FaceColor',[1,1,1]);
plot(t_span(tt),AA(:,1),'-sr','LineWidth',2)
ylabel('Number of Individuals')
xlabel('Starting Date')
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['AccumDeathVac',num2str(zz),'.fig']);
print('-dpng',['AccumDeathVac',num2str(zz)]);

%%%% Increment in accumulated deaths
figure
hold on
box on
title('Chicago')
area(t_span(tt(2:end)),diff(AA(:,3)),'linestyle',':','FaceColor',[255,160,122]/255);
area(t_span(tt(2:end)),diff(AA(:,2)),'linestyle',':','FaceColor',[1,1,1]);
plot(t_span(tt(2:end)),diff(AA(:,1)),'-sr','LineWidth',2)
ylabel('Number of Individuals')
xlabel('Starting Date')
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,'NewDeathVacChicago.fig');
print('-dpng','NewDeathVacChicago');
