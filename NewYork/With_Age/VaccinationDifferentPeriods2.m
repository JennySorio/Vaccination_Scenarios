clear all; clc; close all; format long e;
load data_20201218;


ndays = 183;%203;
l=13;
z=0;
t_actual2 = 0:length(t_actual)-l+ndays;
t_span = datetime(2020,3,01) + caldays(0:length(t_actual2)-1);


Hosp = ones(length(t_actual2),1);
Hosp(3:length(t_actual)) = min(1,data2(2:end,3)./data2(1:end-1,1));
Hosp(length(t_actual)-l:end) = mean(Hosp(length(t_actual)-l-z:length(t_actual)-l))*ones;
Death = ones(size(t_actual2));
Death(3:length(t_actual)) = min(1,data2(2:end,4)./data2(1:end-1,3));
Death(length(t_actual)-l:end) = mean(Death(length(t_actual)-l-z:length(t_actual)-l))*ones;
params.factorDeath = @(t)interp1(t_actual2,Death,t);
params.factorWorse = @(t)interp1(t_actual2,Hosp,t);


% EvaluatingPathsFuture;

tt = [215,246,276,307,338,366,397];
ZZ = zeros(length(tt),12);
disp('Cases    Hospitalizations         Deaths    Vaccinated')
for zz = 1:length(tt)-1
Vaccination = zeros(NumberOfAgeClasses,length(t_actual2));
Vaccination(end,t_actual2>=tt(zz)) = 0.95*0.01*ones;
Vaccination(2:end-1,t_actual2>=tt(zz+1)) = 0.95*0.01*ones;
params.VaccinationRate = @(t)vaccination(Vaccination,t_actual2,t,NumberOfAgeClasses);
EvaluatingPathsVaccinationFuture3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aux = sort(Vaccinated');
aux2 = round(0*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95Vac = [min(aux);max(aux)]*N;

% %%% New Cases:
% aux = sort(NewCasesBoot');
% aux2 = round(0*NSamples);
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

% H = [100 100 600 300];
% %%% plotting Results:
% figure
% hold on
% % grid on
% box on
% title('Daily Infections - New York City')
% % h1=area(t_span,CI95NewCasesVac(2,:),'linestyle',':','FaceColor',[51,236,255]/255,'FaceAlpha',0.3);
% % h2=area(t_span,CI95NewCasesVac(1,:),'linestyle',':','FaceColor',[1,1,1]);
% % h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % plot(t_span,median(NewCasesBoot*N,2),'r','LineWidth',2)
% plot(t_span,median(NewCasesBootVac*N,2),'b','LineWidth',2)
% % legend('Without Vaccination','With Vaccination','Location','NorthEast')
% ylabel('Number of Individuals')
% xlim([t_span(tt(zz)-10),t_span(end)])
% % h1=plot([t_span(tt(zz)),t_span(tt(zz))],[0,1.3*max(median(NewCasesBoot(tt(zz)-10:end,:)*N,2))],'k','LineWidth',1);
% % h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % ylim([0,1.3*max(median(NewCasesBoot(tt(zz)-10:end,:)*N,2))])
% xtickformat('dd-MMM')
% set(gcf,'Position',H)
% set(gca,'FontSize',16,'FontName','Arial')
% hold off
% % saveas(gcf,['InfectionsNYCVac',num2str(zz),'.fig']);
% % print('-dpng',['InfectionsNYCVac2',num2str(zz)]);
% % 
% % figure
% % hold on
% % % grid on
% % box on
% % title('Daily Deaths - New York City')
% % h1=area(t_span,CI95NewDeaths(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
% % h2=area(t_span,CI95NewDeaths(1,:),'linestyle',':','FaceColor',[1,1,1]);
% % h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % h1=area(t_span,CI95NewDeathsVac(2,:),'linestyle',':','FaceColor',[51,236,255]/255);
% % h2=area(t_span,CI95NewDeathsVac(1,:),'linestyle',':','FaceColor',[1,1,1]);
% % h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % plot(t_span,median(NewDeathsBoot*N,2),'r','LineWidth',2)
% % % plot(t_span(2:end),data1(:,2),'r','LineWidth',2)
% % plot(t_span,median(NewDeathsBootVac*N,2),'b','LineWidth',2)
% % legend('Without Vaccination','With Vaccination','Location','NorthEast')
% % ylabel('Number of Individuals')
% % xlim([t_span(tt(zz)-10),t_span(end)])
% % h1=plot([t_span(tt(zz)),t_span(tt(zz))],[0,1.3*max(median(NewDeathsBoot(tt(zz)-10:end,:)*N,2))],'k','LineWidth',1);
% % h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % ylim([0,1.3*max(median(NewDeathsBoot(tt(zz)-10:end,:)*N,2))])
% % xtickformat('dd-MMM')
% % set(gcf,'Position',H)
% % set(gca,'FontSize',16,'FontName','Arial')
% % hold off
% % saveas(gcf,['DeathsNYCVac2',num2str(zz),'.fig']);
% % print('-dpng',['DeathsNYCVac2',num2str(zz)]);
% % 
% %%%%%%%%%
% figure
% hold on
% % grid on
% box on
% title('Daily Hospitalizations - New York City')
% % h1=area(t_span,CI95NewHospVac(2,:),'linestyle',':','FaceColor',[51,236,255]/255);
% % h2=area(t_span,CI95NewHospVac(1,:),'linestyle',':','FaceColor',[1,1,1]);
% % h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % plot(t_span,median(NewHospBoot*N,2),'r','LineWidth',2)
% % plot(t_span(2:end),data1(:,3),'r','LineWidth',2)
% plot(t_span,median(NewHospBootVac*N,2),'b','LineWidth',2)
% plot(t_span,300*ones(length(t_actual2)),'k')
% plot(t_span,100*ones(length(t_actual2)),'k')
% legend('Without Vaccination','With Vaccination','Location','NorthEast')
% ylabel('Number of Individuals')
% xlim([t_span(tt(zz)-10),t_span(end)])
% % h1=plot([t_span(tt(zz)),t_span(tt(zz))],[0,1.3*max(median(NewHospBoot(tt(zz)-10:end,:)*N,2))],'k','LineWidth',1);
% % h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
% % ylim([0,1.3*max(median(NewHospBoot(tt(zz)-10:end,:)*N,2))])
% xtickformat('dd-MMM')
% set(gcf,'Position',H)
% set(gca,'FontSize',16,'FontName','Arial')
% hold off
% % saveas(gcf,['HospNYCVac2',num2str(zz),'.fig']);
% % print('-dpng',['HospNYCVac2',num2str(zz)]);


% disp(['Total Vaccinated: ',num2str(round([sum(median(N*Vaccinated(tt(zz):end,:),2)),sum(CI95Vac(1,tt(zz):end)),sum(CI95Vac(2,tt(zz):end))]))])
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
area(t_span(tt(1:end-1)),AA(1:end-1,3),'linestyle',':','FaceColor',[255,160,122]/255);
area(t_span(tt(1:end-1)),AA(1:end-1,2),'linestyle',':','FaceColor',[1,1,1]);
plot(t_span(tt(1:end-1)),AA(1:end-1,1),'-sr','LineWidth',2)
ylabel('Number of Individuals')
xlabel('Starting Date')
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['AccumDeathVacNYCB3.fig']);
print('-dpng',['AccumDeathVacNYCB3']);

figure
hold on
box on
title('New York City')
area(t_span(tt(3:end-2)),AA(3:end-2,3),'linestyle',':','FaceColor',[255,160,122]/255);
area(t_span(tt(3:end-2)),AA(3:end-2,2),'linestyle',':','FaceColor',[1,1,1]);
plot(t_span(tt(3:end-2)),AA(3:end-2,1),'-sr','LineWidth',2)
ylabel('Number of Individuals')
xlabel('Starting Date')
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['AccumDeathVac2NYCB3.fig']);
print('-dpng',['AccumDeathVac2NYCB3']);

figure
hold on
box on
title('New York City')
area(t_span(tt(2:end-1)),diff(AA(1:end-1,3)),'linestyle',':','FaceColor',[255,160,122]/255);
area(t_span(tt(2:end-1)),diff(AA(1:end-1,2)),'linestyle',':','FaceColor',[1,1,1]);
plot(t_span(tt(2:end-1)),diff(AA(1:end-1,1)),'-sr','LineWidth',2)
ylabel('Number of Individuals')
xlabel('Starting Date')
xtickformat('dd-MMM')
set(gcf,'Position',H)
set(gca,'FontSize',16,'FontName','Arial')
hold off
saveas(gcf,['NewDeathVacNYCB3.fig']);
print('-dpng',['NewDeathVacNYCB3']);

