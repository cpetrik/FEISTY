function [TEeffM,TEeff_ATL,TEeff_LTLd,TEeff_HTLd] = hist_fished_effTEs_useDet_ensem(BE,mnpp,mdet,mmz_loss,mlz_loss,...
    mf_prod50,mp_prod50,md_prod50,lp_prod50,ld_prod50,fname,simname)

% FEISTY Hindcast globally, monthly means saved
% Transfer efficiency ("effective") 
% Use BE*det

Pmf=mf_prod50;
Pmp=mp_prod50;
Pmd=md_prod50;
Plp=lp_prod50;
Pld=ld_prod50;

Pmf(Pmf(:)<0) = 0;
Pmp(Pmp(:)<0) = 0;
Pmd(Pmd(:)<0) = 0;
Plp(Plp(:)<0) = 0;
Pld(Pld(:)<0) = 0;

AllM = Pmp+Pmf+Pmd;
AllL = Plp+Pld;

%% ltl as vectors
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
grid = csvread([cpath 'grid_csv.csv']);
ID = grid(:,1);

mnpp = mnpp(ID);
mdet = mdet(ID);
mmz_loss = mmz_loss(ID);
mlz_loss = mlz_loss(ID);

%% Effective TEs
% With BE*det instead of Bent
TEeffM = AllM ./ (BE*mdet + mmz_loss + mlz_loss); 

%TEeff_ATL = production_L/NPP
TEeff_ATL = AllL ./ mnpp;
TEeff_ATL(TEeff_ATL==-Inf) = NaN;
TEeff_ATL(TEeff_ATL==Inf) = NaN;
TEeff_ATL(TEeff_ATL<0) = NaN;

%TEeff_LTL = (production_benthic_invert+mesozoo_prod_to_fish)/NPP
TEeff_LTLd = (BE*mdet + mmz_loss + mlz_loss) ./ mnpp;
TEeff_LTLd(TEeff_LTLd==-Inf) = NaN;
TEeff_LTLd(TEeff_LTLd==Inf) = NaN;
TEeff_LTLd(TEeff_LTLd<0) = NaN;

%TEeff_HTL = production_L/(production_benthic_invert+mesozoo_prod_to_fish)
TEeff_HTLd = AllL ./ (BE*mdet + mmz_loss + mlz_loss); 
TEeff_HTLd(TEeff_HTLd<0) = NaN;

% TELTLd1 = real(TEeff_LTLd.^(1/1.25));
% TELTLd2 = real(TEeff_LTLd.^(1/1.5));
% 
% TEM = real(TEeffM.^(1/2));          %should this be 1/1?
% TEL = real(TEeff_ATL.^(1/4));       %should this be 1/3?
% TEHTLd = real(TEeff_HTLd.^(1/3));   %should this be 1/2?

%% save
save([fname '_TEeffDet_' simname '.mat'],'TEeffM',...
    'AllM','AllL','mmz_loss','mlz_loss','mnpp','mdet',...
    'TEeff_ATL','TEeff_LTLd','TEeff_HTLd');


end
