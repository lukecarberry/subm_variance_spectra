% Luke Carberry
% 2024 04 03
% Processing code to plot a single matchup of Landsat & MODIS imagery, and 
% calculate spectral slope

% Plots Figures 1, 2 from manuscript

% uses functions stored in "functions" folder on github, including these
% third-party functions:

% LargestRectangle.m - Peter Seibold (2020). Largest inscribed rectangle square or circle, MATLAB Central File Exchange.

% brewermap.m - Stephen23 (2024). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap/releases/tag/3.2.5), GitHub. Retrieved April 3, 2024.

% crit_interp_g - Travis Wiens (2024). Peak Interpolation (https://www.mathworks.com/matlabcentral/fileexchange/24465-peak-interpolation), MATLAB Central File Exchange. Retrieved April 3, 2024.

% inpaint_nans.m - John D'Errico (2024). inpaint_nans (https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans), MATLAB Central File Exchange. Retrieved April 3, 2024.

clear
region = {"PC"};r = 1;
addpath(genpath('~/your/folders/subm_variance_spectra/functions'));
cd("~/your/folders/subm_variance_spectra/single_matchup")


addpath(genpath('~/Library/Mobile Documents/com~apple~CloudDocs/Documents/6_Nick_Project/Tools'));
cd("/Users/lukecarberry/Documents/6_Nick_Project/Projects/LS_timeseries/manuscript/subm_variance_spectra/single_matchup")

files = dir('*.nc');
fname = files(1).name(1:13);
rsublat(1,:) = [33.8 35.1]; %Point Conception, CA region
rsublon(1,:) = [-122.4 -120.6];

% load nc files
chl_ls = ncread(files(1).name,'chl_ls');
lat_ls = ncread(files(2).name,'lat_ls');
lon_ls = ncread(files(3).name,'lon_ls');
sst_ls = ncread(files(4).name,'sst_ls');
cloudmask_ls = ncread(files(5).name,'cloudmask_ls'); cloudmask_ls(cloudmask_ls>0) = 1;

chl_md = ncread(files(6).name,'chl_md');
lat_md = ncread(files(6).name,'lat_md');
lon_md = ncread(files(6).name,'lon_md');
sst_md = ncread(files(6).name,'sst_md');
cloudmask_md = ncread(files(7).name,'cloudmask_md'); cloudmask_md(cloudmask_md>0) = 1;

% apply the cloudmasks to sst and chl
% for reference, these are computed using my code coherent_cloudmasks.m
% the basic formula is combine the Landsat flagged regions and the MODIS 
% flagged regions, then interpolating the Landsat flags down to MODIS 
% resolution (using 50% of MODIS pixel as threshold), adding these to the 
% MODIS flag array for the MODIS flags, and then re-interpolating that back
% to the Landsat grid. Matchups with fewer than 30 usable Landsat transects 
% were excluded
sst_ls(cloudmask_ls>0) = NaN;
chl_ls(cloudmask_ls>0) = NaN;
sst_md(cloudmask_md>0) = NaN;
chl_md(cloudmask_md>0) = NaN;

% calculate spatial resolution of Landsat and MODIS images
[dx, dy] = dist_in_km([lat_ls(round(size(lat_ls,1)/2),1), lon_ls(round(size(lat_ls,1)/2),1)],[lat_ls(round(size(lat_ls,1)/2),2), lon_ls(round(size(lat_ls,1)/2),2)]);
dXL = sqrt(dx^2+dy^2)*1000;
[dx, dy] = dist_in_km([lat_md(round(size(lat_md,1)/2),1), lon_md(round(size(lat_md,1)/2),1)],[lat_md(round(size(lat_md,1)/2),2), lon_md(round(size(lat_md,1)/2),2)]);
dXM = sqrt(dx^2+dy^2)*1000;

% subset the larger image
[~,~,sst_md] = subregion(lat_md,lon_md,sst_md,rsublat(r,:),rsublon(r,:));
[lat_md,lon_md,chl_md] = subregion(lat_md,lon_md,chl_md,rsublat(r,:),rsublon(r,:));
[~,~,sst_ls] = subregion(lat_ls,lon_ls,sst_ls,rsublat(r,:),rsublon(r,:));
[lat_ls,lon_ls,chl_ls] = subregion(lat_ls,lon_ls,chl_ls,rsublat(r,:),rsublon(r,:));

% fill small holes within Landsat image to account for speckling
chl_ls(chl_ls>30) = NaN; % [chl] above 30 mg/m3 is likely inaccurate
chl_ls = fill_holes(chl_ls,14,5); 
sst_ls = fill_holes(sst_ls,14,5);

%find the largest rectangle within the region and further subset
LRout = LargestRectangle(~isnan(chl_ls),0,0,0,0,0) 
chl_ls_small = chl_ls(LRout(3,2):LRout(4,2),LRout(2,1):LRout(3,1));
sst_ls_small = sst_ls(LRout(3,2):LRout(4,2),LRout(2,1):LRout(3,1));
lat_ls_small = lat_ls(LRout(3,2):LRout(4,2),LRout(2,1):LRout(3,1));
lon_ls_small = lon_ls(LRout(3,2):LRout(4,2),LRout(2,1):LRout(3,1));

lon_box = [min(lon_ls_small,[],'all') min(lon_ls_small,[],'all') max(lon_ls_small,[],'all') max(lon_ls_small,[],'all') min(lon_ls_small,[],'all')];
lat_box = [min(lat_ls_small,[],'all') max(lat_ls_small,[],'all') max(lat_ls_small,[],'all') min(lat_ls_small,[],'all') min(lat_ls_small,[],'all')];

%get rid of nans within the rectangle (sometimes they pop up on the edge)
chl_ls_small = inpaint_nans(chl_ls_small,5);
sst_ls_small = inpaint_nans(sst_ls_small,5);

ntp = 7; % number of tapers

% calculating spectra for Landsat
[Pyy_LS,W_LS] = pmtm(detrend(chl_ls_small')./std(chl_ls_small,[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(chl_ls_small'),1/dXL,'Tapers','sine');
Pyy_LS_med = mean(Pyy_LS,2,'omitnan');
[Pxx_LS,W_LS] = pmtm(detrend(sst_ls_small')./std(sst_ls_small,[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(chl_ls_small'),1/dXL,'Tapers','sine');
Pxx_LS_med = mean(Pxx_LS,2,'omitnan');

% rescaling the rectangle from Landsat to MODIS dimensions
rowresc = size(chl_ls,1)./size(chl_md,1);
colresc = size(chl_ls,2)./size(chl_md,2);

% filling small MODIS holes
chl_md = fill_holes(chl_md,3,5);
sst_md = fill_holes(sst_md,3,5);

% selecting MODIS region, using rescaled box from Landsat
chl_md_small = chl_md(ceil(LRout(3,2)/rowresc):round(LRout(4,2)/rowresc),ceil(LRout(2,1)/colresc):round(LRout(3,1)/colresc));
sst_md_small = sst_md(ceil(LRout(3,2)/rowresc):round(LRout(4,2)/rowresc),ceil(LRout(2,1)/colresc):round(LRout(3,1)/colresc));
lat_md_small = lat_md(ceil(LRout(3,2)/rowresc):round(LRout(4,2)/rowresc),ceil(LRout(2,1)/colresc):round(LRout(3,1)/colresc));
lon_md_small = lon_md(ceil(LRout(3,2)/rowresc):round(LRout(4,2)/rowresc),ceil(LRout(2,1)/colresc):round(LRout(3,1)/colresc));

lon_box_md = [min(lon_md_small,[],'all') min(lon_md_small,[],'all') max(lon_md_small,[],'all') max(lon_md_small,[],'all') min(lon_md_small,[],'all')];
lat_box_md = [min(lat_md_small,[],'all') max(lat_md_small,[],'all') max(lat_md_small,[],'all') min(lat_md_small,[],'all') min(lat_md_small,[],'all')];

% inpainting just if there's nans on the edge
chl_md_small = inpaint_nans(chl_md_small,5);
sst_md_small = inpaint_nans(sst_md_small,5);

% calculating MODIS spectra
ntp = 7;
[Pyy_MD,W_MD] = pmtm(detrend(chl_md_small')./std(chl_md_small,[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(chl_md_small'),1/dXM,'Tapers','sine');
Pyy_MD_med = mean(Pyy_MD,2,'omitnan');
[Pxx_MD,W_MD] = pmtm(detrend(sst_md_small')./std(sst_md_small,[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(chl_md_small'),1/dXM,'Tapers','sine');
Pxx_MD_med = mean(Pxx_MD,2,'omitnan');

% testing whether the spectra were calculated
colorsst = brewermap(6,'Oranges');
figure(2),clf
pc = loglog(W_MD(1:end-1),Pxx_MD(1:end-1,:),'-','LineWidth',1.5,'color',[colorsst(end,:) .25]);hold on


%%

% calculating energy containing scale (peak in variance-conserving spectrum), for spectral slope 
[~,I] = max(Pxx_LS_med(W_LS<2e-4).*W_LS(W_LS<2e-4));
[Pxx_LS_ECS,Pxx_LS_E,~]=crit_interp_g(Pxx_LS_med(I-1:I+1).*W_LS(I-1:I+1),W_LS(I-1:I+1));

[~,I] = max(Pyy_LS_med(W_LS<2e-4).*W_LS(W_LS<2e-4));
[Pyy_LS_ECS,Pyy_LS_E,~]=crit_interp_g(Pyy_LS_med(I-1:I+1).*W_LS(I-1:I+1),W_LS(I-1:I+1));

[~,I] = max(Pxx_MD_med(W_MD<2e-4).*W_MD(W_MD<2e-4));
[Pxx_MD_ECS,Pxx_MD_E,~]=crit_interp_g(Pxx_MD_med(I-1:I+1).*W_MD(I-1:I+1),W_MD(I-1:I+1));

[~,I] = max(Pyy_MD_med(W_MD<2e-4).*W_MD(W_MD<2e-4));
[Pyy_MD_ECS,Pyy_MD_E,~]=crit_interp_g(Pyy_MD_med(I-1:I+1).*W_MD(I-1:I+1),W_MD(I-1:I+1));

% calculating slopes
mindX = 4e-5;maxdX = 2e-4;Jl = find(mindX < W_LS & maxdX > W_LS);Jm = find(mindX < W_MD & maxdX > W_MD);
[Pxx_LS_slope,Pxx_LS_line] = spectra_linefit2(W_LS,Pxx_LS_med,mindX,maxdX);
[Pyy_LS_slope,Pyy_LS_line] = spectra_linefit2(W_LS,Pyy_LS_med,mindX,maxdX);
[Pxx_MD_slope,Pxx_MD_line] = spectra_linefit2(W_MD,Pxx_MD_med,mindX,maxdX);
[Pyy_MD_slope,Pyy_MD_line] = spectra_linefit2(W_MD,Pyy_MD_med,mindX,maxdX);

[Pxx_LS_ECS_slope,Pxx_LS_ECS_line] = spectra_linefit2(W_LS,Pxx_LS_med,Pxx_LS_ECS,maxdX);
[Pyy_LS_ECS_slope,Pyy_LS_ECS_line] = spectra_linefit2(W_LS,Pyy_LS_med,Pxx_LS_ECS,maxdX);  
[Pxx_MD_ECS_slope,Pxx_MD_ECS_line] = spectra_linefit2(W_MD,Pxx_MD_med,Pxx_MD_ECS,maxdX);
[Pyy_MD_ECS_slope,Pyy_MD_ECS_line] = spectra_linefit2(W_MD,Pyy_MD_med,Pyy_MD_ECS,maxdX);

%% averaging Landsat @ MODIS resolution

%scaling factor
dX = round(dXM/dXL);
if isinteger(dX/2)==1
    dX = dX+1;
end

% spatial smoothing of Landsat
h = fspecial('average',[dX dX]);
chl_lm1 = imfilter(log(inpaint_nans(chl_ls,5)),h,'replicate');
sst_lm1 = imfilter(inpaint_nans(sst_ls,5),h,'replicate');

% create new landsat @ modis (L@M) lat lon grid
lon = linspace(min(lon_ls,[],'all'),max(lon_ls,[],'all'),size(lon_md,2));
lat = linspace(min(lat_ls,[],'all'),max(lat_ls,[],'all'),size(lat_md,1))';
[lon_lm,lat_lm] = meshgrid(lon,lat);

% find indices to subsample smoothed landsat data at coarse grid
for i = 1:length(lat)
    [~,lat_loc(i,1)] = min(abs(lat_ls(:,1) - lat_lm(i,1)));
end
for j = 1:length(lon)
    [~,lon_loc(1,j)] = min(abs(lon_ls(1,:) - lon_lm(1,j)));
end

% subsample
chl_lm = exp(chl_lm1(lat_loc,lon_loc));
sst_lm = sst_lm1(lat_loc,lon_loc);
% make sure the nans match up 
chl_lm(isnan(chl_md)) = NaN;
sst_lm(isnan(chl_md)) = NaN;

% select the new rectangular L@M region
chl_lm_small = chl_lm(ceil(LRout(3,2)/rowresc):round(LRout(4,2)/rowresc),ceil(LRout(2,1)/colresc):round(LRout(3,1)/colresc));
sst_lm_small = sst_lm(ceil(LRout(3,2)/rowresc):round(LRout(4,2)/rowresc),ceil(LRout(2,1)/colresc):round(LRout(3,1)/colresc));
lat_lm_small = lat_lm(ceil(LRout(3,2)/rowresc):round(LRout(4,2)/rowresc),ceil(LRout(2,1)/colresc):round(LRout(3,1)/colresc));
lon_lm_small = lon_lm(ceil(LRout(3,2)/rowresc):round(LRout(4,2)/rowresc),ceil(LRout(2,1)/colresc):round(LRout(3,1)/colresc));


%% calculate spectra for L@M and plot
[dx, dy] = dist_in_km([lat_lm_small(round(size(lat_lm_small,1)/2),1), lon_lm_small(round(size(lat_lm_small,1)/2),1)],[lat_lm_small(round(size(lat_lm_small,1)/2),2), lon_lm_small(round(size(lat_lm_small,1)/2),2)]);
dXLM = sqrt(dx^2+dy^2)*1000;clear dx dy

[Pyy_LM,~] = pmtm(detrend(chl_lm_small')./std(chl_lm_small,[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(chl_lm_small'),1/dXLM,'Tapers','sine');
Pyy_LM_med = mean(Pyy_LM,2,'omitnan');
[Pxx_LM,W_LM] = pmtm(detrend(sst_lm_small')./std(sst_lm_small,[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(sst_lm_small'),1/dXLM,'Tapers','sine');
Pxx_LM_med = mean(Pxx_LM,2,'omitnan');

% calculating LM energy containing scale 
[~,I] = max(Pxx_LM_med(W_LM<2e-4).*W_LM(W_LM<2e-4));
[Pxx_LM_ECS,Pxx_LM_E,~]=crit_interp_g(Pxx_LM_med(I-1:I+1).*W_LM(I-1:I+1),W_LM(I-1:I+1));

[~,I] = max(Pyy_LM_med(W_LM<2e-4).*W_LM(W_LM<2e-4));
[Pyy_LM_ECS,Pyy_LM_E,~]=crit_interp_g(Pyy_LM_med(I-1:I+1).*W_LM(I-1:I+1),W_LM(I-1:I+1));

% calculating slopes
[Pxx_LM_ECS_slope,Pxx_LM_ECS_line] = spectra_linefit2(W_LM,Pxx_LM_med,Pxx_LM_ECS,maxdX);
[Pyy_LM_ECS_slope,Pyy_LM_ECS_line] = spectra_linefit2(W_LM,Pyy_LM_med,Pxx_LM_ECS,maxdX);  

colorchl = brewermap(6,'Blues');
colorsst = brewermap(6,'Oranges');

% plot a single transect
row = 1
figure(3),clf
subplot(2,3,1:2)
plot(lon_ls_small(1,:),chl_ls_small(row,:),'linewidth',2),hold on
plot(lon_lm_small(1,:),chl_lm_small(row,:),'linewidth',2)
plot(lon_md_small(1,:),chl_md_small(row,:),'linewidth',2)
legend('Landsat','L@M','MODIS','location','northwest'),title('Chl')
subplot(2,3,4:5)
plot(lon_ls_small(1,:),sst_ls_small(row,:),'linewidth',2),hold on
plot(lon_lm_small(1,:),sst_lm_small(row,:),'linewidth',2)
plot(lon_md_small(1,:),sst_md_small(row,:),'linewidth',2)
legend('Landsat','L@M','MODIS','location','northwest'),title('SST')

[pyy_ls,~] = pmtm(detrend(chl_ls_small(row,:)')./std(chl_ls_small(row,:),[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(chl_ls_small(row,:)'),1/dXL,'Tapers','sine');
[pyy_lm,~] = pmtm(detrend(chl_lm_small(row,:)')./std(chl_lm_small(row,:),[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(chl_lm_small(row,:)'),1/dXLM,'Tapers','sine');

[pxx_ls,w_ls] = pmtm(detrend(sst_ls_small(row,:)')./std(sst_ls_small(row,:),[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(sst_ls_small(row,:)'),1/dXL,'Tapers','sine');
[pxx_lm,w_lm] = pmtm(detrend(sst_lm_small(row,:)')./std(sst_lm_small(row,:),[],2,'omitnan')',(ntp:-1:1)/sum(1:ntp),length(sst_lm_small(row,:)'),1/dXLM,'Tapers','sine');

subplot(2,3,[3,6])
p(1) = loglog(w_ls,pxx_ls,'-','LineWidth',1.5,'color',colorsst(end,:));hold on
p(2) = loglog(w_lm,pxx_lm,'-','LineWidth',2.5,'color',colorsst(end-1,:));
p(3) = loglog(w_ls,pyy_ls,'-','LineWidth',1.5,'color',colorchl(end,:));
p(4) = loglog(w_lm,pyy_lm,'-','LineWidth',2.5,'color',colorchl(end-1,:));grid on
legend(p,"Landsat SST slope = "+num2str(round(Pxx_LS_ECS_slope,1)),"L@M SST slope = "+num2str(round(Pxx_LM_ECS_slope,1)),"Landsat Chl slope = "+num2str(round(Pyy_LS_ECS_slope,1)),"L@M Chl slope = "+num2str(round(Pyy_LM_ECS_slope,1)));legend boxoff

%% Figure 2 in manuscript

colorchl = brewermap(6,'Blues');
colorsst = brewermap(6,'Oranges');

figure(2),clf
patch([3e-6 Pxx_LS_ECS Pxx_LS_ECS 3e-6],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
patch([2e-4 5e-3 5e-3 2e-4],[1e0 1e0 6e4 6e4],'k','edgecolor','none');alpha(0.1),hold on
set(gca,'Xscale','log','Yscale','log');

p(1) = loglog(W_LS(1:end-1),Pyy_LS_med(1:end-1),'-','LineWidth',1.5,'color',colorchl(end,:));
p(2) = loglog(W_MD(1:end-1),Pyy_MD_med(1:end-1),'-','LineWidth',3.5,'color',colorchl(end-2,:));
p(3) = loglog(W_LM(1:end-1),Pyy_LM_med(1:end-1),'-','LineWidth',2.5,'color',colorchl(end-1,:));
p(4) = loglog(W_LS(1:end-1),Pxx_LS_med(1:end-1),'-','LineWidth',1.5,'color',colorsst(end,:));hold on
p(5) = loglog(W_MD(1:end-1),Pxx_MD_med(1:end-1),'-','LineWidth',3.5,'color',colorsst(end-2,:));
p(6) = loglog(W_LM(1:end-1),Pxx_LM_med(1:end-1),'-','LineWidth',2.5,'color',colorsst(end-1,:));
% loglog([3e-5;2e-4],4e-6*[3e-5;2e-4].^(-5/3),'-k')
% text([6e-5],[2.5e1],"k^{-5/3}",'HorizontalAlignment','right','fontsize',14)

xlabel('log Wavenumber (m^{-1})','FontSize',16);
ylabel('log Power Spectral Density (m)','FontSize',16);
xlim([3e-6 5e-3]);ylim([1e0 6e4]);grid minor
set(gca,'FontSize',16,'FontWeight','normal','XGrid','on','YGrid','on');%,'XDir','reverse');
% legend(p,"Landsat SST slope = "+num2str(round(Pxx_LS_ECS_slope,1)),"L@M SST slope = "+num2str(round(Pxx_LM_ECS_slope,1)),"Landsat Chl slope = "+num2str(round(Pyy_LS_ECS_slope,1)),"L@M Chl slope = "+num2str(round(Pyy_LM_ECS_slope,1)));legend boxoff
% xline([3e-5;2e-4],'--m')
legend(p(1:6),"Landsat Chl slope = "+num2str(round(Pyy_LS_ECS_slope,1)),"MODIS Chl slope = "+num2str(round(Pyy_MD_ECS_slope,1)),"L@M Chl slope = "+num2str(round(Pyy_LM_ECS_slope,1)),"Landsat SST slope = "+num2str(round(Pxx_LS_ECS_slope,1)),"MODIS SST slope = "+num2str(round(Pxx_MD_ECS_slope,1)),"L@M SST slope = "+num2str(round(Pxx_LM_ECS_slope,1)),'fontsize',18);%legend boxoff


set(gcf,'Position',[50 1200 600 550])
fstr = fname + "_6spec.jpg";
     % print(gcf, '-djpeg100',fstr, '-r600');

%% Figure 1 in manuscript

% ybound = [34.1 34.45];xbound = [-121.1 -120.7]; % small region to see detail
ybound = [rsublat(1) rsublat(2)];xbound = [rsublon(1) rsublon(2)];

fg1 = figure(1);clf
subplot(2,3,1)
h=pcolor(lon_ls,lat_ls,chl_ls);hold on
set(h,'linestyle','none');title('Landsat')
map=cmocean('-gray');colormap(gca,map);
set(gca, 'FontSize',16);caxis([.2 5]);
c=colorbar;c.FontSize=16;%c.Label.String = 'Chl (mg m^{-3})'; %labels the colorbar
xlabel('Longitude','FontSize',16); %labels x axis
ylabel('Latitude','FontSize',16); %labels y axiscolorbar;
daspect([1/cos(mean(lat_md,'all')*pi/180) 1 1]),set(gca,'ColorScale','log')
plot(lon_box,lat_box,'+r','MarkerSize',10)
ylim([ybound(1) ybound(2)]);xlim([xbound(1) xbound(2)])

subplot(2,3,2)
h=pcolor(lon_md,lat_md,chl_md);hold on
set(h,'linestyle','none');title('MODIS')
map=cmocean('-gray');colormap(gca,map);
set(gca, 'FontSize',16);caxis([.2 5]);
c=colorbar;c.FontSize=16;%c.Label.String = 'Chl (mg m^{-3})';c.Label.FontWeight = "bold";c.Label.Color = 'g' %labels the colorbar
xlabel('Longitude','FontSize',16); %labels x axis
ylabel('Latitude','FontSize',16); %labels y axiscolorbar;
daspect([1/cos(mean(lat_md,'all')*pi/180) 1 1]),set(gca,'ColorScale','log')
plot(lon_box_md,lat_box_md,'+r','MarkerSize',10)
ylim([ybound(1) ybound(2)]);xlim([xbound(1) xbound(2)])

subplot(2,3,3)
h=pcolor(lon_lm,lat_lm,chl_lm);hold on
set(h,'linestyle','none');title('Landsat @ MODIS Res.')
map=cmocean('-gray');colormap(gca,map);
set(gca, 'FontSize',16);caxis([.2 5]);
c=colorbar;c.Label.String = 'Chl (mg m^{-3})';c.FontSize=16;c.Label.FontWeight = "bold";c.Label.Color = 'k' %labels the colorbar
xlabel('Longitude','FontSize',16); %labels x axis
ylabel('Latitude','FontSize',16); %labels y axiscolorbar;
daspect([1/cos(mean(lat_md,'all')*pi/180) 1 1]),set(gca,'ColorScale','log')
plot(lon_box_md,lat_box_md,'+r','MarkerSize',10)
ylim([ybound(1) ybound(2)]);xlim([xbound(1) xbound(2)])

subplot(2,3,4)
h=pcolor(lon_ls,lat_ls,sst_ls);hold on
set(h,'linestyle','none')
map=cmocean('gray');colormap(gca,map);
set(gca, 'FontSize',16);caxis([13 17]);%caxis([12.5 18]);
c=colorbar;c.FontSize=16;%c.Label.String = 'SST (^{o}C)'; %labels the colorbar
xlabel('Longitude','FontSize',16); %labels x axis
ylabel('Latitude','FontSize',16); %labels y axiscolorbar;
daspect([1/cos(mean(lat_md,'all')*pi/180) 1 1])
plot(lon_box,lat_box,'+r','MarkerSize',10)
ylim([ybound(1) ybound(2)]);xlim([xbound(1) xbound(2)])

subplot(2,3,5)
h=pcolor(lon_md,lat_md,sst_md);hold on
set(h,'linestyle','none')
map=cmocean('gray');colormap(gca,map);
set(gca, 'FontSize',16);caxis([13 17]);%caxis([12.5 18]);
c=colorbar;c.FontSize=16;%c.Label.String = 'SST (^{o}C)';c.Label.FontWeight = "bold";c.Label.Color = 'r'; %labels the colorbar
xlabel('Longitude','FontSize',16); %labels x axis
ylabel('Latitude','FontSize',16); %labels y axiscolorbar;
daspect([1/cos(mean(lat_md,'all')*pi/180) 1 1])
plot(lon_box_md,lat_box_md,'+r','MarkerSize',10)
ylim([ybound(1) ybound(2)]);xlim([xbound(1) xbound(2)])

subplot(2,3,6)
h=pcolor(lon_lm,lat_lm,sst_lm);hold on
set(h,'linestyle','none');
map=cmocean('gray');colormap(gca,map);
set(gca, 'FontSize',16);caxis([13 17]);%caxis([12.5 18]);
c=colorbar;c.Label.String = 'SST (^{o}C)';c.FontSize=16;c.Label.FontWeight = "bold";c.Label.Color = 'k'; %labels the colorbar
xlabel('Longitude','FontSize',16); %labels x axis
ylabel('Latitude','FontSize',16); %labels y axiscolorbar;
daspect([1/cos(mean(lat_md,'all')*pi/180) 1 1])
plot(lon_box_md,lat_box_md,'+r','MarkerSize',10)
ylim([ybound(1) ybound(2)]);xlim([xbound(1) xbound(2)])

% charlbl =  compose("(%s)",('a':'z').'); 
% AddLetters2Plots(fg1,charlbl,'FontSize',25,'VShift',-.025,'HShift',-.03,'FontWeight','normal')

fstr = fname + "_same_color_L@M.jpg";
set(gcf,'Position',[50 1200 1400 800])

     % print(gcf, '-djpeg100',fstr, '-r600');

