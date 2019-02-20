% POEM output at all locations
% Hist fished prod only ran through 1999 (139 years)
% All values starting at mo 1669 are 9.969e+36 == NA

clear all
close all

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
harv = 'All_fish03';

fpath=['/Volumes/GFDL/NC/Matlab_new_size/' cfile '/'];

%% SP
ncid = netcdf.open([fpath 'Historic_' harv '_prod_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

[ni,nj] = size(prod);
nas = find(isnan(prod(1,:)));
nt = nas(1) - 1;

SP.prod = prod;
Sml_p.prod = prod(:,nt);
clear prod

%% SF
ncid = netcdf.open([fpath 'Historic_' harv '_prod_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

SF.prod = prod(:,1:nj);
Sml_f.prod = prod(:,nt);
clear prod

%% SD
ncid = netcdf.open([fpath 'Historic_' harv '_prod_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

SD.prod = prod;
Sml_d.prod = prod(:,nt);
clear prod

% MP
ncid = netcdf.open([fpath 'Historic_' harv '_prod_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MP.prod = prod;
MP.yield = yield;
Med_p.prod = prod(:,nt);
clear prod yield

% MF
ncid = netcdf.open([fpath 'Historic_' harv '_prod_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MF.prod = prod;
MF.yield = yield;
Med_f.prod = prod(:,nt);
clear prod yield

% MD
ncid = netcdf.open([fpath 'Historic_' harv '_prod_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

MD.prod = prod;
MD.yield = yield;
Med_d.prod = prod(:,nt);
clear prod yield

% LP
ncid = netcdf.open([fpath 'Historic_' harv '_prod_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

LP.prod = prod;
LP.yield = yield;
Lrg_p.prod = prod(:,nt);
clear prod yield

% LD
ncid = netcdf.open([fpath 'Historic_' harv '_prod_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

LD.prod = prod;
LD.yield = yield;
Lrg_d.prod = prod(:,nt);
clear prod yield

% Benthic material
ncid = netcdf.open([fpath 'Historic_' harv '_prod_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
BENT.bio = biomass(:,nt);
clear biomass

%% Take means and totals
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%Time
sp_tprod=mean(SP.prod,1);
sf_tprod=mean(SF.prod,1);
sd_tprod=mean(SD.prod,1);
mp_tprod=mean(MP.prod,1);
mf_tprod=mean(MF.prod,1);
md_tprod=mean(MD.prod,1);
lp_tprod=mean(LP.prod,1);
ld_tprod=mean(LD.prod,1);
b_tmean=mean(Bent.bio,1);

mf_tmy=mean(MF.yield,1);
mp_tmy=mean(MP.yield,1);
md_tmy=mean(MD.yield,1);
lp_tmy=mean(LP.yield,1);
ld_tmy=mean(LD.yield,1);

%% 50 yrs (1951-2000)
y = 1860+(1/12):(1/12):2005;
yr50=find(y>=1951 & y<2001);
sp_prod50=nanmean(SP.prod(:,yr50),2);
sf_prod50=nanmean(SF.prod(:,yr50),2);
sd_prod50=nanmean(SD.prod(:,yr50),2);
mp_prod50=nanmean(MP.prod(:,yr50),2);
mf_prod50=nanmean(MF.prod(:,yr50),2);
md_prod50=nanmean(MD.prod(:,yr50),2);
lp_prod50=nanmean(LP.prod(:,yr50),2);
ld_prod50=nanmean(LD.prod(:,yr50),2);
b_mean50=nanmean(Bent.bio(:,yr50),2);

mf_my50=nanmean(MF.yield(:,yr50),2);
mp_my50=nanmean(MP.yield(:,yr50),2);
md_my50=nanmean(MD.yield(:,yr50),2);
lp_my50=nanmean(LP.yield(:,yr50),2);
ld_my50=nanmean(LD.yield(:,yr50),2);

% 1990-1995 (Climatology)
lyr=find(y>=1990 & y<1995);
%Means
sp_prod5=mean(SP.prod(:,lyr),2);
sf_prod5=mean(SF.prod(:,lyr),2);
sd_prod5=mean(SD.prod(:,lyr),2);
mp_prod5=mean(MP.prod(:,lyr),2);
mf_prod5=mean(MF.prod(:,lyr),2);
md_prod5=mean(MD.prod(:,lyr),2);
lp_prod5=mean(LP.prod(:,lyr),2);
ld_prod5=mean(LD.prod(:,lyr),2);
b_mean5=mean(Bent.bio(:,lyr),2);

mf_my5=mean(MF.yield(:,lyr),2);
mp_my5=mean(MP.yield(:,lyr),2);
md_my5=mean(MD.yield(:,lyr),2);
lp_my5=mean(LP.yield(:,lyr),2);
ld_my5=mean(LD.yield(:,lyr),2);

% 1990 (Climatology)
yr1=find(y>=1990 & y<1991);
%Totals
sp_ptot=sum(SP.prod(:,yr1).*MNTH,2);
sf_ptot=sum(SF.prod(:,yr1).*MNTH,2);
sd_ptot=sum(SD.prod(:,yr1).*MNTH,2);
mp_ptot=sum(MP.prod(:,yr1).*MNTH,2);
mf_ptot=sum(MF.prod(:,yr1).*MNTH,2);
md_ptot=sum(MD.prod(:,yr1).*MNTH,2);
lp_ptot=sum(LP.prod(:,yr1).*MNTH,2);
ld_ptot=sum(LD.prod(:,yr1).*MNTH,2);
b_tot=sum(Bent.bio(:,yr1).*MNTH,2);

%% Every 5 years
st=1:60:length(time);
en=60:60:length(time);

for n=1:length(st)
    b_mean(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
    
    sp_prod(:,n)=nanmean(SP.prod(:,st(n):en(n)),2);
    sf_prod(:,n)=nanmean(SF.prod(:,st(n):en(n)),2);
    sd_prod(:,n)=nanmean(SD.prod(:,st(n):en(n)),2);
    mp_prod(:,n)=nanmean(MP.prod(:,st(n):en(n)),2);
    mf_prod(:,n)=nanmean(MF.prod(:,st(n):en(n)),2);
    md_prod(:,n)=nanmean(MD.prod(:,st(n):en(n)),2);
    lp_prod(:,n)=nanmean(LP.prod(:,st(n):en(n)),2);
    ld_prod(:,n)=nanmean(LD.prod(:,st(n):en(n)),2);
    
%     sp_rec(:,n)=nanmean(SP.rec(:,st(n):en(n)),2);
%     sf_rec(:,n)=nanmean(SF.rec(:,st(n):en(n)),2);
%     sd_rec(:,n)=nanmean(SD.rec(:,st(n):en(n)),2);
%     mp_rec(:,n)=nanmean(MP.rec(:,st(n):en(n)),2);
%     mf_rec(:,n)=nanmean(MF.rec(:,st(n):en(n)),2);
%     md_rec(:,n)=nanmean(MD.rec(:,st(n):en(n)),2);
%     lp_rec(:,n)=nanmean(LP.rec(:,st(n):en(n)),2);
%     ld_rec(:,n)=nanmean(LD.rec(:,st(n):en(n)),2);
    
%     mp_my(:,n)=nanmean(MP.yield(:,st(n):en(n)),2);
%     mf_my(:,n)=nanmean(MF.yield(:,st(n):en(n)),2);
%     md_my(:,n)=nanmean(MD.yield(:,st(n):en(n)),2);
%     lp_my(:,n)=nanmean(LP.yield(:,st(n):en(n)),2);
%     ld_my(:,n)=nanmean(LD.yield(:,st(n):en(n)),2);
end

%%
save([fpath 'Means_Historic_' harv '_prod_' cfile '.mat'],'time',...
    'y','yr50','yr1','lyr',...
    'sf_tprod','sp_tprod','sd_tprod',...
    'mf_tprod','mp_tprod','md_tprod',...
    'lp_tprod','ld_tprod','b_tmean',...
    'sf_prod50','sp_prod50','sd_prod50',...
    'mf_prod50','mp_prod50','md_prod50',...
    'lp_prod50','ld_prod50','b_mean50',...
    'sf_prod5','sp_prod5','sd_prod5',...
    'mf_prod5','mp_prod5','md_prod5',...
    'lp_prod5','ld_prod5','b_mean5',...
    'sf_ptot','sp_ptot','sd_ptot',...
    'mf_ptot','mp_ptot','md_ptot',...
    'lp_ptot','ld_ptot','b_tot',...
    'sf_prod','sp_prod','sd_prod',...
    'mf_prod','mp_prod','md_prod',...
    'lp_prod','ld_prod',...
    'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
    'mf_my5','mp_my5','md_my5','lp_my5','ld_my5',...
    'mf_my50','mp_my50','md_my50','lp_my50','ld_my50');

%save([fpath 'Means_Historic_' harv '_prod_' cfile '.mat']);

%%
load([fpath 'Means_Historic_' harv '_prod_' cfile '.mat']);
HistProdT(1,:)=sf_tprod;
HistProdT(2,:)=sp_tprod;
HistProdT(3,:)=sd_tprod;
HistProdT(4,:)=mf_tprod;
HistProdT(5,:)=mp_tprod;
HistProdT(6,:)=md_tprod;
HistProdT(7,:)=lp_tprod;
HistProdT(8,:)=ld_tprod;
HistProdT(9,:)=b_tmean;
save([fpath 'Means_Historic_',harv,'_prod_' cfile '.mat'],'HistProdT','-append');







