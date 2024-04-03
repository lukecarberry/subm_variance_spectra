% Luke Carberry
% 2024 04 01 
% Processing code to load the variance spectra from the Landsat - MODIS
% matchups, calculate spectral slope and statistics, and plot them
% Plots Figures 3, 4, Table 3

% Pxx = SST spectra
% Pyy = Chl spectra 

% LS = Landsat
% MD = MODIS 
% LM = Landsat @ MODIS resolution

% W = wavenumbers for spectra
% ECS = energy-containing scale, used to find low-wavenumber bound for spectral slope calculation
% med = mean spectrum for each region 
% line = straight line for plotting spectral slope

% this code will have to be adapted if you are working with the json files
% on github
clear
cd("/your/data/directory/")
region = ["AK","MX","PC","SB"];

r = struct;numreg = 4;

for reg = 1:numreg
files = dir('*_matchup_spectra.mat');load(files(reg).name);

dts = NaT(length(oc.date),1);
for i = 1:length(oc.date)
    try
        dts(i,1) = datetime((oc.date{1,i}),'InputFormat','yyyy-MM-dd');
    catch
        continue
    end
end
dts = day(dts,'dayofyear');
oc.dts = dts;
r.oc(reg) = oc;
clear oc
end

%% stats and making Manuscript Table 3

LCHL_slope = cell2mat(horzcat(r.oc(:).Pyy_LS_ECS_slope));
LSST_slope = cell2mat(horzcat(r.oc(:).Pxx_LS_ECS_slope));
MCHL_slope = cell2mat(horzcat(r.oc(:).Pyy_MD_ECS_slope));
MSST_slope = cell2mat(horzcat(r.oc(:).Pxx_MD_ECS_slope));
LMCHL_slope = cell2mat(horzcat(r.oc(:).Pyy_LM_ECS_slope));
LMSST_slope = cell2mat(horzcat(r.oc(:).Pxx_LM_ECS_slope));

clear J
J = isoutlier(LCHL_slope,'grubbs') + isnan(LCHL_slope);
J = J + isoutlier(LSST_slope,'grubbs') + isnan(LSST_slope);
J = J + isoutlier(MCHL_slope,'grubbs') + isnan(MCHL_slope);
J = J + isoutlier(MSST_slope,'grubbs') + isnan(MSST_slope);
J = J + isoutlier(LMCHL_slope,'grubbs') + isnan(LMCHL_slope);
J = J + isoutlier(LMSST_slope,'grubbs') + isnan(LMSST_slope);

LCHL_slope(J>0) = NaN;LSST_slope(J>0) = NaN;
MCHL_slope(J>0) = NaN;MSST_slope(J>0) = NaN;
LMCHL_slope(J>0) = NaN;LMSST_slope(J>0) = NaN;

variable = ["Landsat Chl";"Landsat SST";"MODIS Chl";"MODIS SST";"L@M Chl";"L@M SST"];

Slope = [mean(LCHL_slope,"omitnan"); mean(LSST_slope,"omitnan"); mean(MCHL_slope,"omitnan"); mean(MSST_slope,"omitnan"); mean(LMCHL_slope,"omitnan"); mean(LMSST_slope,"omitnan")];
Slope_std = [std(LCHL_slope,[],"omitnan"); std(LSST_slope,[],"omitnan"); std(MCHL_slope,[],"omitnan"); std(MSST_slope,[],"omitnan"); std(LMCHL_slope,[],"omitnan"); std(LMSST_slope,[],"omitnan")];
Slope_min = [min(LCHL_slope); min(LSST_slope); min(MCHL_slope); min(MSST_slope); min(LMCHL_slope); min(LMSST_slope)];
Slope_max = [max(LCHL_slope); max(LSST_slope); max(MCHL_slope); max(MSST_slope); max(LMCHL_slope); max(LMSST_slope)];

manusc = table(variable,Slope_min,Slope_max,Slope,Slope_std);

% comparing average spectra slope across instrument and measurement
[h,p] = ttest2(LCHL_slope,LSST_slope)
[h,p] = ttest2(MCHL_slope,MSST_slope)
[h,p] = ttest2(LMCHL_slope,LMSST_slope)
[h,p] = ttest2(LCHL_slope,MCHL_slope)
[h,p] = ttest2(LSST_slope,MSST_slope)
[h,p] = ttest2(LMCHL_slope,MCHL_slope)
[h,p] = ttest2(LMSST_slope,MSST_slope)
[h,p] = ttest2(LCHL_slope,LMCHL_slope)
[h,p] = ttest2(LSST_slope,LMSST_slope)

%% Manuscript Figure 4
fg1 = figure(1),clf

subplot(221)
clear a h
for reg = 1:4
I = find(~isnan(r.oc(reg).dts));
a(1)=scatter(cell2mat(r.oc(reg).Pyy_LS_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_MD_ECS_slope(I)),50,'b','filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',.4);hold on
a(2)=scatter(cell2mat(r.oc(reg).Pyy_LS_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_LM_ECS_slope(I)),50,'r','filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',.5);hold on
end
plot(-5/3,-5/3,'ksq','markersize',9);
plot([-10 1],[-10 1],'--k');
xlim([-4.9 -0.3]);ylim([-4.9 -0.3]);axis square;grid;box on
ax = gca;ax.YTickMode= 'auto';ax.YTick = unique(round(ax.YTick));
ylabel('1.5 km Res. Chl (m^2)'),xlabel('Landsat Chl (m^2)')
set(gca, 'FontSize',16);%title('Landsat Chl')
h(1) = scatter(nan, nan, 50,'b','filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',.4 , 'DisplayName', 'MODIS');
h(2) = scatter(nan, nan, 50,'r','filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',.4 , 'DisplayName', 'L@M');
legend(h,'location','northwest');legend boxoff
% [obj] = legend(a,'MODIS','L@M','location','northwest','fontsize',20);
% objhl = findobj(objh, 'type', 'hggroup');
% set(objhl, 'Markersize', 20);

subplot(222)
for reg = 1:4
I = find(~isnan(r.oc(reg).dts));
a(1)=scatter(cell2mat(r.oc(reg).Pxx_LS_ECS_slope(I)),cell2mat(r.oc(reg).Pxx_MD_ECS_slope(I)),50,'b','filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',.4);hold on
a(2)=scatter(cell2mat(r.oc(reg).Pxx_LS_ECS_slope(I)),cell2mat(r.oc(reg).Pxx_LM_ECS_slope(I)),50,'r','filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',.5);hold on
end
plot(-5/3,-5/3,'ksq','markersize',9);
plot([-10 1],[-10 1],'--k');
xlim([-4.9 -0.3]);ylim([-4.9 -0.3]);axis square;grid;box on
ax = gca;ax.YTickMode= 'auto';ax.YTick = unique(round(ax.YTick));
xlabel('Landsat SST (m^2)'),ylabel('1.5 km Res. SST (m^2)')
set(gca, 'FontSize',16);%title('Landsat Chl')
% legend(a,'MODIS Point Concep.','L@M Point Concep.','MODIS Other','L@M Other','location','southeast'),legend boxoff

%
subplot(223)
for reg = 1:4
I = find(~isnan(r.oc(reg).dts));
scatter(cell2mat(r.oc(reg).Pyy_MD_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_LM_ECS_slope(I)),15,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k'),hold on
end
% for reg = 4
% I = find(~isnan(r.oc(reg).dts));
% scatter(cell2mat(r.oc(reg).Pyy_MD_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_LM_ECS_slope(I)),50,cmap(reg,:),'filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',1),hold on
% end
plot(-5/3,-5/3,'ksq','markersize',9);
plot([-10 1],[-10 1],'--k');
xlim([-4.9 -0.3]);ylim([-4.9 -0.3]);axis square;grid;box on
ax = gca;ax.YTickMode= 'auto';ax.YTick = unique(round(ax.YTick));
xlabel('MODIS Chl (m^2)'),ylabel('L@M Chl (m^2)')
set(gca, 'FontSize',16);
% legend('AK','MX','PC','SBC','-5/3','location','northwest','fontsize',18);legend('boxoff')

subplot(224)
for reg = 1:4
I = find(~isnan(r.oc(reg).dts));
scatter(cell2mat(r.oc(reg).Pxx_MD_ECS_slope(I)),cell2mat(r.oc(reg).Pxx_LM_ECS_slope(I)),15,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k'),hold on
end
% for reg = 4
% I = find(~isnan(r.oc(reg).dts));
% scatter(cell2mat(r.oc(reg).Pxx_MD_ECS_slope(I)),cell2mat(r.oc(reg).Pxx_LM_ECS_slope(I)),50,cmap(reg,:),'filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',1),hold on
% end
plot(-5/3,-5/3,'ksq','markersize',9);
plot([-10 1],[-10 1],'--k');
xlim([-4.9 -0.3]);ylim([-4.9 -0.3]);grid;axis square;box on
ax = gca;ax.YTickMode= 'auto';ax.YTick = unique(round(ax.YTick));
xlabel('MODIS SST (m^2)'),ylabel('L@M SST (m^2)')
set(gca, 'FontSize',16);
% legend('','','Other Matchups','Point Conception','-5/3','location','northwest');legend boxoff

charlbl =  compose("(%s)",('a':'z').'); 
AddLetters2Plots(fg2,charlbl,'FontSize',25,'VShift',-.05,'HShift',-.02,'FontWeight','normal')

% legend('AK','MX','PC','SBC','-5/3','location','northwest','fontsize',18);legend('boxoff')
% sgtitle('Slopes','fontsize',21','fontweight','bold')
set(gcf,'Position',[50 1200 900 800])

cd ("/your/save/directory/")
fstr = "figure_4.jpg";
    % print(gcf, '-djpeg100',fstr, '-r600');

%% Manuscript Figure 3

N = size(r.oc(numreg).Pyy_LS_ECS_slope,2);
k = linspace(3e-5,2e-4,10);
days = dts(~isnan(dts));

fg3 = figure(4),clf
for reg = 4

I = ~isnan(r.oc(reg).dts);

N = size(padcat(2,r.oc(reg).W_LS{:}),2);
N2 = size(padcat(2,r.oc(reg).W_LS{:}),1);
N3 = 1:N;
N3 = N3(~cellfun('isempty',r.oc(reg).W_LS));

subplot(3,3,1)
patch([3e-6 3e-5 3e-5 3e-6],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
patch([2e-4 5e-3 5e-3 2e-4],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
set(gca,'Xscale','log','Yscale','log');
for n = 1:sum(~cellfun('isempty',r.oc(reg).W_LS))
    try
        if J(n)==0
            l=loglog(r.oc(reg).W_LS{1,N3(n)},r.oc(reg).Pyy_LS_med{1,N3(n)},'color','k','LineWidth',1.5);l.Color(4) = 0.25;hold on
        else
            continue
        end
        % scatter(r.oc(reg).W_LS{1,N3(n)},r.oc(reg).Pyy_LS_med{1,N3(n)},NaN,ones(N2,1)*days(n),'filled');hold on
    catch
        continue
    end
end
grid on
xlabel('log Wavenumber (m^{-1})'),ylabel('log Power Spectral Density (m)')
set(gca, 'FontSize',16);title('Landsat Chl')
xlim([3e-6 5e-3]);ylim([2e0 6e4])
h=loglog(k,2e-8*k.^Slope(1),'k','linewidth',1.7);
label(h,num2str(round(Slope(1),2)),'location','left','horizontalalignment','right','fontsize',14)

subplot(3,3,2)
patch([3e-6 3e-5 3e-5 3e-6],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
patch([2e-4 5e-3 5e-3 2e-4],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
set(gca,'Xscale','log','Yscale','log');
for n = 1:sum(~cellfun('isempty',r.oc(reg).W_LS))
    try
        if J(n)==0
            l=loglog(r.oc(reg).W_MD{1,N3(n)},r.oc(reg).Pyy_MD_med{1,N3(n)},'color','k','LineWidth',1.5);l.Color(4) = 0.25;hold on
        else
            continue
        end
    catch
        continue
    end
end
grid on
xlabel('log Wavenumber (m^{-1})'),ylabel('log Power Spectral Density (m)')
set(gca, 'FontSize',16);title('MODIS Chl')
xlim([3e-6 5e-3]);ylim([2e0 6e4])
h=loglog(k,4e-7*k.^Slope(3),'k','linewidth',1.7);
label(h,num2str(round(Slope(3),2)),'location','left','horizontalalignment','right','fontsize',14)

subplot(3,3,3)
patch([3e-6 3e-5 3e-5 3e-6],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
patch([2e-4 5e-3 5e-3 2e-4],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
set(gca,'Xscale','log','Yscale','log');
for n = 1:sum(~cellfun('isempty',r.oc(reg).W_LS))
    try
        if J(n)==0
            l=loglog(r.oc(reg).W_MD{1,N3(n)},r.oc(reg).Pyy_LM_med{1,N3(n)},'color','k','LineWidth',1.5);l.Color(4) = 0.25;hold on
        else
            continue
        end
    catch
        continue
    end
end
grid on
xlabel('log Wavenumber (m^{-1})'),ylabel('log Power Spectral Density (m)')
set(gca, 'FontSize',16);title('L@M Chl')
xlim([3e-6 5e-3]);ylim([2e0 6e4])
h=loglog(k,2e-8*k.^Slope(5),'k','linewidth',1.7);
label(h,num2str(round(Slope(5),2)),'location','left','horizontalalignment','right','fontsize',14)
% set(gca,'Yscale','linear'),set(gca,'Xscale','linear')
% clim([0 365]);colormap(gca,cmap);cb=colorbar;
% cb.Position = cb.Position + [1e-10];cb.Position(1) = cb.Position(1)+.04;
% cb.Label.String = 'Day of year';

subplot(3,3,4)
patch([3e-6 3e-5 3e-5 3e-6],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
patch([2e-4 5e-3 5e-3 2e-4],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
set(gca,'Xscale','log','Yscale','log');
for n = 1:sum(~cellfun('isempty',r.oc(reg).W_LS))
    try
        if J(n) == 0
            l=loglog(r.oc(reg).W_LS{1,N3(n)},r.oc(reg).Pxx_LS_med{1,N3(n)},'color','k','LineWidth',1.5);l.Color(4) = 0.25;hold on
        else
            continue
        end
    catch
        continue
    end
end
grid on
xlabel('log Wavenumber (m^{-1})'),ylabel('log Power Spectral Density (m)')
set(gca, 'FontSize',16);title('Landsat SST')
xlim([3e-6 5e-3]);ylim([2e0 6e4])
h=loglog(k,7e-11*k.^Slope(2),'k','linewidth',1.7);
label(h,num2str(round(Slope(2),2)),'location','left','horizontalalignment','right','fontsize',14)
% set(gca,'Yscale','linear'),set(gca,'Xscale','linear')

subplot(3,3,5)
patch([3e-6 3e-5 3e-5 3e-6],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
patch([2e-4 5e-3 5e-3 2e-4],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
set(gca,'Xscale','log','Yscale','log');
for n = 1:sum(~cellfun('isempty',r.oc(reg).W_LS))
    try
        if J(n)==0
            l=loglog(r.oc(reg).W_MD{1,N3(n)},r.oc(reg).Pxx_MD_med{1,N3(n)},'color','k','LineWidth',1.5);l.Color(4) = .25;hold on
        else
            continue
        end
    catch
        continue
    end
end
grid on
xlabel('log Wavenumber (m^{-1})'),ylabel('log Power Spectral Density (m)')
set(gca, 'FontSize',16);title('MODIS SST')
xlim([3e-6 5e-3]);ylim([2e0 6e4])
h=loglog(k,4e-7*k.^Slope(4),'k','linewidth',1.7);
label(h,num2str(round(Slope(4),2)),'location','left','horizontalalignment','right','fontsize',14)
% set(gca,'Yscale','linear'),set(gca,'Xscale','linear')

subplot(3,3,6)
patch([3e-6 3e-5 3e-5 3e-6],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
patch([2e-4 5e-3 5e-3 2e-4],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
set(gca,'Xscale','log','Yscale','log');
for n = 1:sum(~cellfun('isempty',r.oc(reg).W_LS))
    try
        if J(n) == 0
            l=loglog(r.oc(reg).W_MD{1,N3(n)},r.oc(reg).Pxx_LM_med{1,N3(n)},'color','k','LineWidth',1.5);l.Color(4) = .25;hold on
        else
            continue
        end
    catch
        continue
    end
end
grid on
xlabel('log Wavenumber (m^{-1})'),ylabel('log Power Spectral Density (m)')
set(gca, 'FontSize',16);title('L@M SST')
xlim([3e-6 5e-3]);ylim([2e0 6e4])
h=loglog(k,1e-10*k.^Slope(6),'k','linewidth',1.7);
label(h,num2str(round(Slope(6),2)),'location','left','horizontalalignment','right','fontsize',14)
% set(gca,'Yscale','linear'),set(gca,'Xscale','linear')
% clim([0 365]);colormap(gca,cmap);cb=colorbar;
% cb.Position = cb.Position + [1e-10];cb.Position(1) = cb.Position(1)+.04;
% cb.Label.String = 'Day of year';

end

N = size(padcat(2,r.oc(reg).W_LS{:}),2);
N2 = size(padcat(2,r.oc(reg).W_LS{:}),1);
cm = copper(4);

cmap = [cm(1,:);cm(1,:);cm(1,:);cm(3,:)];

subplot(337)
for reg = 1:4
I = find(~isnan(r.oc(reg).dts));
scatter(cell2mat(r.oc(reg).Pxx_LS_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_LS_ECS_slope(I)),15,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k'),hold on
end
% for reg = 4
% I = find(~isnan(r.oc(reg).dts));
% scatter(cell2mat(r.oc(reg).Pxx_LS_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_LS_ECS_slope(I)),50,cmap(reg,:),'filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',1),hold on
% end
plot([-10 1],[-10 1],'--k');
plot(-5/3,-5/3,'ksq','markersize',9);
xlim([-4.9 -0.3]);ylim([-4.9 -0.3]);axis square;grid;box on
ax = gca;ax.YTickMode= 'auto';ax.YTick = unique(round(ax.YTick));
xlabel('Landsat SST (m^2)'),ylabel('Landsat Chl (m^2)')
set(gca, 'FontSize',16);title('Landsat Spectral Slopes')

subplot(338)
for reg = 1:4
I = find(~isnan(r.oc(reg).dts));
scatter(cell2mat(r.oc(reg).Pxx_MD_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_MD_ECS_slope(I)),15,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k'),hold on
end
% for reg = 4
% I = find(~isnan(r.oc(reg).dts));
% scatter(cell2mat(r.oc(reg).Pxx_MD_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_MD_ECS_slope(I)),50,cmap(reg,:),'filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',1),hold on
% end
plot([-10 1],[-10 1],'--k');
plot(-5/3,-5/3,'ksq','markersize',9);
xlim([-4.9 -0.3]);ylim([-4.9 -0.3]);axis square;grid;box on
ax = gca;ax.YTickMode= 'auto';ax.YTick = unique(round(ax.YTick));
xlabel('MODIS SST (m^2)'),ylabel('MODIS Chl (m^2)')
set(gca, 'FontSize',16);title('MODIS Spectral Slopes')

subplot(339)
for reg = 1:4
I = find(~isnan(r.oc(reg).dts));
scatter(cell2mat(r.oc(reg).Pxx_LM_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_LM_ECS_slope(I)),15,'k','filled','MarkerFaceAlpha',1,'MarkerEdgeColor','k'),hold on
end
% for reg = 4
% I = find(~isnan(r.oc(reg).dts));
% scatter(cell2mat(r.oc(reg).Pxx_LM_ECS_slope(I)),cell2mat(r.oc(reg).Pyy_LM_ECS_slope(I)),50,cmap(reg,:),'filled','MarkerEdgeColor','k','LineWidth',.5,'MarkerFaceAlpha',1),hold on
% end
plot(-5/3,-5/3,'ksq','markersize',9);
plot([-10 1],[-10 1],'--k');
xlim([-4.9 -0.3]);ylim([-4.9 -0.3]);axis square;grid;box on
ax = gca;ax.YTickMode= 'auto';ax.YTick = unique(round(ax.YTick));
xlabel('L@M SST (m^2)'),ylabel('L@M Chl (m^2)')
set(gca, 'FontSize',16);title('L@M Spectral Slopes')
% legend('','','Other Matchups','Point Conception','-5/3','location','southeast');legend boxoff
charlbl =  compose("(%s)",('a':'z').'); 
% AddLetters2Plots(fg2,charlbl,'FontSize',25,'VShift',-.02,'HShift',-.035,'FontWeight','normal')


set(gcf,'Position',[50 1200 1200 1100])
charlbl =  compose("(%s)",('a':'z').'); 
AddLetters2Plots(fg3,charlbl,'FontSize',20,'VShift',-.03,'HShift',-.02,'FontWeight','normal')


cd ("/your/save/directory/")
fstr = "figure_3.jpg";
    % print(gcf, '-djpeg100',fstr, '-r600');
