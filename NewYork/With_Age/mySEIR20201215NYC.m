
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bootstraping_20201215;
save dataNYC_20201215;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% CI Evaluation
%%% Beta
aux = sort(BETABoot);
aux2 = round(0.05*NSamples);
aux = aux(aux2+1:end-aux2,:);
CI95BETA = [min(aux);max(aux)];

%%% New Cases:
aux = sort(NewCasesBoot');
aux2 = round(0.05*NSamples);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Elapsed Time: ',num2str(toc),' seconds.'])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting Results:

figure
hold on
box on
title('Daily New Infections')
h1=area(t_span,CI95NewCases(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewCases(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data(:,1),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
% plot(t_span(2:end),data2(:,1),'b','LineWidth',2)
plot(t_span,NewCases,'r','LineWidth',2)
plot([t_span(day),t_span(day)],[0,1.3*max(data(:,1))],'--k','LineWidth',2)
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
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
h1=area(t_span,CI95NewDeaths(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewDeaths(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data(:,4),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
plot(t_span,N*NewDeaths,'r','LineWidth',2)
% plot(t_span(2:end),data2(:,4),'b','LineWidth',2)
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
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
h1=area(t_span,CI95NewHosp(2,:),'linestyle',':','FaceColor',[255,160,122]/255);
h2=area(t_span,CI95NewHosp(1,:),'linestyle',':','FaceColor',[1,1,1]);
bar(t_span(2:end),data(:,3),'FaceColor',[0 0.75 0.75],'EdgeColor','none')
% plot(t_span(2:end),data2(:,3),'b','LineWidth',2)
h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
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
