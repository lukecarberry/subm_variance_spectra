% Luke Carberry
% 2024 04 04
% Processing code to create a coherent cloud mask across the four fields (Landsat + MODIS Chl
% & SST)

clear
addpath(genpath('~/your/folders/subm_variance_spectra/functions'));
cd("~/your/folders/landsat_imagery_mat_files/")
region = "SBC";

files = dir('*.mat');

% rsublat(1,:) = [58.0 59.1]; %AK region
% rsublon(1,:) = [-147 -144.6];
% rsublat(1,:) = [23.8 25.0]; %Baja region
% rsublon(1,:) = [-117.1 -115.5];
% rsublat(1,:) = [33.7 35.1]; %PC region
% rsublon(1,:) = [-122.6 -120.4];
rsublat(1,:) = [33.7 35.1]; %PC region
rsublon(1,:) = [-121 -119.1];

f=1;
for f = 1:length(files)
files(f).name

cd(files(f).folder)

dt = datetime(files(f).name(4:13) + "183400",'InputFormat','yyyy-MM-ddHHmmss','TimeZone','UTC');
str = files(f).name(4:13);

load(files(f).name,'chlOC','lat','lon','sst');

[~,~,sst_ls] = subregion(lat,lon,sst,rsublat(1,:),rsublon(1,:));
[lat_ls,lon_ls,chl_ls] = subregion(lat,lon,chlOC,rsublat(1,:),rsublon(1,:));

[dx, dy] = dist_in_km([lat_ls(round(size(lat_ls,1)/2),1), lon_ls(round(size(lat_ls,1)/2),1)],[lat_ls(round(size(lat_ls,1)/2),2), lon_ls(round(size(lat_ls,1)/2)+1,2)]);     
dXL = round(sqrt(dx^2+dy^2)*1000);

[~,~,sst_ls] = rotate_grid(lat_ls,lon_ls,sst_ls,[rsublat(1,1),rsublat(1,2)],[rsublon(1,1),rsublon(1,2)],dXL);
[lat_ls,lon_ls,chl_ls] = rotate_grid(lat_ls,lon_ls,chl_ls,[rsublat(1,1),rsublat(1,2)],[rsublon(1,1),rsublon(1,2)],dXL);

% if sum(~isnan(lrs(:,1))) < round(size(sstr,1)/5)
%     continue
% end

cd("~/your/folders/modis_imagery_OC_mat_files/")
dt2 = datetime(dt,'InputFormat','yyyyMMdd','Format','yyMMdd');

try
    [matchupfilename_OC] = load_modis_matchup(dt2);
    load(matchupfilename_OC,'chlr','latr','lonr')
    chl_md = chlr;lat_mdc = latr;lon_mdc = lonr;clear chlr latr lonr 
catch
    continue
end

cd("~/your/folders/modis_imagery_SST_mat_files/")
try
    matchupfilename_SST = strrep(matchupfilename_OC,'OC','SST');
    load(matchupfilename_SST,'sstr')
    sst_mdc = sstr;clear sstr latr lonr 
catch
    try
        [matchupfilename_SST] = load_modis_matchup(dt2);
        load(matchupfilename_SST,'sstr','latr','lonr')

        %interpolate sst to chl dimensions
        sst_md_obj = griddedInterpolant(lonr',latr',sstr');sst_md_obj.Method = 'nearest';
        sst_i = sst_md_obj(lon_mdc',lat_mdc');sst_mdc = sst_i';clear sstr_mdi2 sst_i sst_md_obj

        sst_mds = sstr;lat_mds = latr;lon_mds = lonr;clear sstr latr lonr 
    catch
        continue
    end
end

%removing outliers
sst_ls(sst_ls < mean(sst_ls,'all','omitnan')-2*std(sst_ls,[],'all','omitnan') | sst_ls > mean(sst_ls,'all','omitnan')+2*std(sst_ls,[],'all','omitnan')) = NaN;
% chl_ls(chl_ls < mean(chl_ls,'all','omitnan')-2*std(chl_ls,[],'all','omitnan') | chl_ls > mean(chl_ls,'all','omitnan')+2*std(chl_ls,[],'all','omitnan')) = NaN;
sst_mdc(sst_mdc < mean(sst_mdc,'all','omitnan')-2*std(sst_mdc,[],'all','omitnan') | sst_mdc > mean(sst_mdc,'all','omitnan')+2*std(sst_mdc,[],'all','omitnan')) = NaN;
% chl_md(chl_md < mean(chl_md,'all','omitnan')-2*std(chl_md,[],'all','omitnan') | chl_md > mean(chl_md,'all','omitnan')+2*std(chl_md,[],'all','omitnan')) = NaN;

chl_ls = fill_holes(chl_ls,16,5);sst_ls = fill_holes(sst_ls,16,5);
chl_md = fill_holes(chl_md,4,5);sst_mdc = fill_holes(sst_mdc,4,5);


%% adding the respective cloud masks

cloud_md = isnan(chl_md) + isnan(sst_mdc);
cloud_ls = isnan(chl_ls) + isnan(sst_ls);

% interpolating the landsat cloudmask

cloud_ls_obj = griddedInterpolant(lon_ls',lat_ls',cloud_ls');%sstr_md_obj.Method = 'nearest';
cloud_ls3 = cloud_ls_obj(lon_mdc',lat_mdc');cloud_ls2 = cloud_ls3'; clear cloud_ls3
cloud_ls3 = round(cloud_ls2);cloud_ls3(cloud_ls3 > 0) = 1;
cloud_ls4 = cloud_ls2;cloud_ls4(cloud_ls2 >= 0.4) = 1;cloud_ls4(cloud_ls2 < 0.4) = 0;

test = sum(cloud_ls4 - cloud_ls3,'all');

%% making the final cloud masks?

fullcloud4md = cloud_md + cloud_ls4;fullcloud4md(fullcloud4md > 0) = 1;
cloud_md_obj = griddedInterpolant(lon_mdc',lat_mdc',fullcloud4md');cloud_md_obj.Method = 'nearest';
cloud_lsi2 = cloud_md_obj(lon_ls',lat_ls');cloud_lsi = cloud_lsi2'; clear cloud_lsi2
cloud_lsi(cloud_lsi>0) = 1;

fullcloud4ls = cloud_ls + cloud_lsi;%fullcloud4ls(fullcloud4ls>0.5) = 1;

fullcloud4ls(fullcloud4ls>0) = 1;
fullcloud4md(fullcloud4md>0) = 1;

sst_ls(fullcloud4ls>0) = NaN;
chl_ls(fullcloud4ls>0) = NaN;
sst_mdc(fullcloud4md>0) = NaN;
chl_md(fullcloud4md>0) = NaN;

[lrs_ls,lrl_ls] = findgoodruns(chl_ls,30,37000);
[lrs_md,lrl_md] = findgoodruns(chl_md,1500,37000);

if sum(~isnan(lrs_ls(:,1))) < 30
    continue
end

%%
figure(1),clf
subplot(2,2,1)
h = pcolor(lon_ls,lat_ls,chl_ls);
set(h,'linestyle','none')
map=cmocean('algae');colormap(gca,map);caxis([median(chl_ls,'all','omitnan')-std(chl_ls,[],'all','omitnan') median(chl_ls,'all','omitnan')+2*std(chl_ls,[],'all','omitnan')]);
c=colorbar;%,title("Chl")
c.Label.String = 'Chl (mg m^{-3})';c.FontSize=14; %labels the colorbar
xlabel('Longitude'); %labels x axis
ylabel('Latitude'); %labels y axiscolorbar;
daspect([1/cos(mean(lat_ls,'all','omitnan')*pi/180) 1 1]);%set(gca,'ColorScale','log');
set(gca, 'FontSize',16);
set(gcf,'Position',[100 1000 1100 800])

subplot(2,2,2)
h = pcolor(lon_ls,lat_ls,sst_ls);
set(h,'linestyle','none')
map=cmocean('thermal');colormap(gca,map);caxis([median(sst_ls,'all','omitnan')-1.5 median(sst_ls,'all','omitnan')+1.5]);
c=colorbar;%,title("SST")
c.Label.String = 'SST (^{o}C)';c.FontSize=14; %labels the colorbar
xlabel('Longitude'); %labels x axis
ylabel('Latitude'); %labels y axiscolorbar;
daspect([1/cos(mean(lat_ls,'all','omitnan')*pi/180) 1 1]);%set(gca,'ColorScale','log');
set(gca, 'FontSize',16);

subplot(2,2,3)
h = pcolor(lon_mdc,lat_mdc,chl_md);
set(h,'linestyle','none')
map=cmocean('algae');colormap(gca,map);caxis([median(chl_md,'all','omitnan')-std(chl_md,[],'all','omitnan') median(chl_md,'all','omitnan')+2*std(chl_md,[],'all','omitnan')]);
c=colorbar;%,title("Chl")
c.Label.String = 'Chl (mg m^{-3})';c.FontSize=14; %labels the colorbar
xlabel('Longitude'); %labels x axis
ylabel('Latitude'); %labels y axiscolorbar;
daspect([1/cos(mean(lat_ls,'all','omitnan')*pi/180) 1 1]);%set(gca,'ColorScale','log');
set(gca, 'FontSize',16);

subplot(2,2,4)
h = pcolor(lon_mdc,lat_mdc,sst_mdc);
set(h,'linestyle','none');
map=cmocean('thermal');colormap(gca,map);caxis([median(sst_mdc,'all','omitnan')-1.5 median(sst_mdc,'all','omitnan')+1.5]);
c=colorbar;%,title(ltstr + " PT")
c.Label.String = 'SST (^{o}C)';c.FontSize=14; %labels the colorbar
xlabel('Longitude'); %labels x axis
ylabel('Latitude'); %labels y axiscolorbar;
daspect([1/cos(mean(lat_ls,'all','omitnan')*pi/180) 1 1]);%set(gca,'ColorScale','log');
set(gca, 'FontSize',16);

    cd("~/your/folders/matchup_scene_jpgs/")
    fstr = region + "_" + str + "_scene.jpg"
    print(gcf, '-djpeg100',fstr, '-r600');

%% make a fstr and save masks
% also need to optimize modis image selection re: landsat image
% load   run indices and save to the same file?

    cd("~/your/folders/matchup_mat_files/")
    fstr = region + "_" + str + "_flags.mat";
    save(fstr,"sst_mdc","sst_ls","fullcloud4ls","fullcloud4md","chl_md", ...
        "chl_ls","lat_ls","lon_ls","lat_mdc","lon_mdc","lrs_ls","lrl_ls","lrs_md","lrl_md","str")

end 









