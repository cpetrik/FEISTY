%%%%!! RUN Climatol Single Locations
function Locs_clim_midw()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Make core parameters/constants (global)
param = [];
param = make_params_midw(param);

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/MIP/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Grid; choose where and when to run the model
load('/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/clim_complete/ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
load('/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/clim_grid_180x360_id_locs_area_dep.mat','ids');
ID = ids;
NX = length(ID);
param.NX = NX;
param.ID = ID;

%! Create a directory for output
Fdir = pwd;
addpath([Fdir '/naming/'])
[fname] = sub_fname_clim_mid(param);

%! Storage
Clim_Sml_f = NaN*ones(DAYS,26,NX);
Clim_Sml_m = NaN*ones(DAYS,26,NX);
Clim_Sml_p = NaN*ones(DAYS,26,NX);
Clim_Sml_l = NaN*ones(DAYS,26,NX);
Clim_Sml_d = NaN*ones(DAYS,26,NX);
Clim_Med_f = NaN*ones(DAYS,26,NX);
Clim_Med_m = NaN*ones(DAYS,26,NX);
Clim_Med_p = NaN*ones(DAYS,26,NX);
Clim_Med_l = NaN*ones(DAYS,26,NX);
Clim_Med_d = NaN*ones(DAYS,26,NX);
Clim_Lrg_p = NaN*ones(DAYS,26,NX);
Clim_Lrg_l = NaN*ones(DAYS,26,NX);
Clim_Lrg_d = NaN*ones(DAYS,26,NX);
Clim_Cobalt = NaN*ones(DAYS,5,NX);

S_Sml_f = NaN*ones(12*YEARS,26,NX);
S_Sml_m = NaN*ones(12*YEARS,26,NX);
S_Sml_p = NaN*ones(12*YEARS,26,NX);
S_Sml_l = NaN*ones(12*YEARS,26,NX);
S_Sml_d = NaN*ones(12*YEARS,26,NX);
S_Med_f = NaN*ones(12*YEARS,26,NX);
S_Med_m = NaN*ones(12*YEARS,26,NX);
S_Med_p = NaN*ones(12*YEARS,26,NX);
S_Med_l = NaN*ones(12*YEARS,26,NX);
S_Med_d = NaN*ones(12*YEARS,26,NX);
S_Lrg_p = NaN*ones(12*YEARS,26,NX);
S_Lrg_l = NaN*ones(12*YEARS,26,NX);
S_Lrg_d = NaN*ones(12*YEARS,26,NX);
S_Cobalt = NaN*ones(12*YEARS,5,NX);


%! Initialize
[Sml_f,Sml_m,Sml_p,Sml_l,Sml_d,Med_f,Med_m,Med_p,Med_l,Med_d,Lrg_p,Lrg_l,Lrg_d,BENT] = ...
        sub_init_fish_midw(ID);

%% %%%%%%%%%%%%%%%%%%%% Run the Model %%%%%%%%%%%%%%%%%%%%
%! Run model
MNT=0;
for YR = 1:YEARS % years

    num2str(YR)

    for DAY = 1:param.DT:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));

        [Sml_f,Sml_m,Sml_p,Sml_l,Sml_d,Med_f,Med_m,Med_p,Med_l,Med_d,...
        Lrg_p,Lrg_l,Lrg_d,BENT,ENVR] = ...
            sub_futbio_midw(ID,DY,COBALT,GRD,Sml_f,Sml_m,Sml_p,Sml_l,Sml_d,Med_f,...
            Med_m,Med_p,Med_l,Med_d,Lrg_p,Lrg_l,Lrg_d,BENT,param);

        %! Store last year
        %if (YR==YEARS)
        Clim_Sml_f(DY,:,:) = [Sml_f.bio Sml_f.enc_f Sml_f.enc_p Sml_f.enc_d Sml_f.enc_zm ...
            Sml_f.enc_zl Sml_f.enc_be Sml_f.con_f Sml_f.con_p Sml_f.con_d Sml_f.con_zm ...
            Sml_f.con_zl Sml_f.con_be Sml_f.I Sml_f.nu Sml_f.gamma Sml_f.die Sml_f.rep ...
            Sml_f.rec Sml_f.clev Sml_f.prod Sml_f.pred Sml_f.nmort Sml_f.met Sml_f.caught Sml_f.fmort]';
        Clim_Sml_m(DY,:,:) = [Sml_m.bio Sml_m.enc_f Sml_m.enc_p Sml_m.enc_d Sml_m.enc_zm Sml_m.enc_zl Sml_m.enc_be Sml_m.con_f Sml_m.con_p Sml_m.con_d Sml_m.con_zm Sml_m.con_zl Sml_m.con_be Sml_m.I Sml_m.nu Sml_m.gamma Sml_m.die Sml_m.rep Sml_m.rec Sml_m.clev Sml_m.prod Sml_m.pred Sml_m.nmort Sml_m.met Sml_m.caught Sml_m.fmort]';
        Clim_Sml_p(DY,:,:) = [Sml_p.bio Sml_p.enc_f Sml_p.enc_p Sml_p.enc_d Sml_p.enc_zm Sml_p.enc_zl Sml_p.enc_be Sml_p.con_f Sml_p.con_p Sml_p.con_d Sml_p.con_zm Sml_p.con_zl Sml_p.con_be Sml_p.I Sml_p.nu Sml_p.gamma Sml_p.die Sml_p.rep Sml_p.rec Sml_p.clev Sml_p.prod Sml_p.pred Sml_p.nmort Sml_p.met Sml_p.caught Sml_p.fmort]';
        Clim_Sml_l(DY,:,:) = [Sml_l.bio Sml_l.enc_f Sml_l.enc_p Sml_l.enc_d Sml_l.enc_zm Sml_l.enc_zl Sml_l.enc_be Sml_l.con_f Sml_l.con_p Sml_l.con_d Sml_l.con_zm Sml_l.con_zl Sml_l.con_be Sml_l.I Sml_l.nu Sml_l.gamma Sml_l.die Sml_l.rep Sml_l.rec Sml_l.clev Sml_l.prod Sml_l.pred Sml_l.nmort Sml_l.met Sml_l.caught Sml_l.fmort]';
        Clim_Sml_d(DY,:,:) = [Sml_d.bio Sml_d.enc_f Sml_d.enc_p Sml_d.enc_d Sml_d.enc_zm Sml_d.enc_zl Sml_d.enc_be Sml_d.con_f Sml_d.con_p Sml_d.con_d Sml_d.con_zm Sml_d.con_zl Sml_d.con_be Sml_d.I Sml_d.nu Sml_d.gamma Sml_d.die Sml_d.rep Sml_d.rec Sml_d.clev Sml_d.prod Sml_d.pred Sml_d.nmort Sml_d.met Sml_d.caught Sml_d.fmort]';
        Clim_Med_f(DY,:,:) = [Med_f.bio Med_f.enc_f Med_f.enc_p Med_f.enc_d Med_f.enc_zm Med_f.enc_zl Med_f.enc_be Med_f.con_f Med_f.con_p Med_f.con_d Med_f.con_zm Med_f.con_zl Med_f.con_be Med_f.I Med_f.nu Med_f.gamma Med_f.die Med_f.rep Med_f.rec Med_f.clev Med_f.prod Med_f.pred Med_f.nmort Med_f.met Med_f.caught Med_f.fmort]';
        Clim_Med_m(DY,:,:) = [Med_m.bio Med_m.enc_f Med_m.enc_p Med_m.enc_d Med_m.enc_zm Med_m.enc_zl Med_m.enc_be Med_m.con_f Med_m.con_p Med_m.con_d Med_m.con_zm Med_m.con_zl Med_m.con_be Med_m.I Med_m.nu Med_m.gamma Med_m.die Med_m.rep Med_m.rec Med_m.clev Med_m.prod Med_m.pred Med_m.nmort Med_m.met Med_m.caught Med_m.fmort]';
        Clim_Med_p(DY,:,:) = [Med_p.bio Med_p.enc_f Med_p.enc_p Med_p.enc_d Med_p.enc_zm Med_p.enc_zl Med_p.enc_be Med_p.con_f Med_p.con_p Med_p.con_d Med_p.con_zm Med_p.con_zl Med_p.con_be Med_p.I Med_p.nu Med_p.gamma Med_p.die Med_p.rep Med_p.rec Med_p.clev Med_p.prod Med_p.pred Med_p.nmort Med_p.met Med_p.caught Med_p.fmort]';
        Clim_Med_l(DY,:,:) = [Med_l.bio Med_l.enc_f Med_l.enc_p Med_l.enc_d Med_l.enc_zm Med_l.enc_zl Med_l.enc_be Med_l.con_f Med_l.con_p Med_l.con_d Med_l.con_zm Med_l.con_zl Med_l.con_be Med_l.I Med_l.nu Med_l.gamma Med_l.die Med_l.rep Med_l.rec Med_l.clev Med_l.prod Med_l.pred Med_l.nmort Med_l.met Med_l.caught Med_l.fmort]';
        Clim_Med_d(DY,:,:) = [Med_d.bio Med_d.enc_f Med_d.enc_p Med_d.enc_d Med_d.enc_zm Med_d.enc_zl Med_d.enc_be Med_d.con_f Med_d.con_p Med_d.con_d Med_d.con_zm Med_d.con_zl Med_d.con_be Med_d.I Med_d.nu Med_d.gamma Med_d.die Med_d.rep Med_d.rec Med_d.clev Med_d.prod Med_d.pred Med_d.nmort Med_d.met Med_d.caught Med_d.fmort]';
        Clim_Lrg_p(DY,:,:) = [Lrg_p.bio Lrg_p.enc_f Lrg_p.enc_p Lrg_p.enc_d Lrg_p.enc_zm Lrg_p.enc_zl Lrg_p.enc_be Lrg_p.con_f Lrg_p.con_p Lrg_p.con_d Lrg_p.con_zm Lrg_p.con_zl Lrg_p.con_be Lrg_p.I Lrg_p.nu Lrg_p.gamma Lrg_p.die Lrg_p.rep Lrg_p.rec Lrg_p.clev Lrg_p.prod Lrg_p.pred Lrg_p.nmort Lrg_p.met Lrg_p.caught Lrg_p.fmort]';
        Clim_Lrg_l(DY,:,:) = [Lrg_l.bio Lrg_l.enc_f Lrg_l.enc_p Lrg_l.enc_d Lrg_l.enc_zm Lrg_l.enc_zl Lrg_l.enc_be Lrg_l.con_f Lrg_l.con_p Lrg_l.con_d Lrg_l.con_zm Lrg_l.con_zl Lrg_l.con_be Lrg_l.I Lrg_l.nu Lrg_l.gamma Lrg_l.die Lrg_l.rep Lrg_l.rec Lrg_l.clev Lrg_l.prod Lrg_l.pred Lrg_l.nmort Lrg_l.met Lrg_l.caught Lrg_l.fmort]';
        Clim_Lrg_d(DY,:,:) = [Lrg_d.bio Lrg_d.enc_f Lrg_d.enc_p Lrg_d.enc_d Lrg_d.enc_zm Lrg_d.enc_zl Lrg_d.enc_be Lrg_d.con_f Lrg_d.con_p Lrg_d.con_d Lrg_d.con_zm Lrg_d.con_zl Lrg_d.con_be Lrg_d.I Lrg_d.nu Lrg_d.gamma Lrg_d.die Lrg_d.rep Lrg_d.rec Lrg_d.clev Lrg_d.prod Lrg_d.pred Lrg_d.nmort Lrg_d.met Lrg_d.caught Lrg_d.fmort]';
        Clim_Cobalt(DY,:,:) = [BENT.mass BENT.pred ENVR.fZm ENVR.fZl ENVR.fB]';
        %end

    end %Days

    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)];    % start of the month
    b = cumsum(MNTH);       % end of the month
    for i = 1:12
        MNT = MNT+1;        % Update monthly ticker
        S_Cobalt(MNT,:,:) = mean(Clim_Cobalt(a(i):b(i),:,:),1);
        S_Sml_f(MNT,:,:) = mean(Clim_Sml_f(a(i):b(i),:,:),1);
        S_Sml_m(MNT,:,:) = mean(Clim_Sml_m(a(i):b(i),:,:),1);
        S_Sml_p(MNT,:,:) = mean(Clim_Sml_p(a(i):b(i),:,:),1);
        S_Sml_l(MNT,:,:) = mean(Clim_Sml_l(a(i):b(i),:,:),1);
        S_Sml_d(MNT,:,:) = mean(Clim_Sml_d(a(i):b(i),:,:),1);
        S_Med_f(MNT,:,:) = mean(Clim_Med_f(a(i):b(i),:,:),1);
        S_Med_m(MNT,:,:) = mean(Clim_Med_m(a(i):b(i),:,:),1);
        S_Med_p(MNT,:,:) = mean(Clim_Med_p(a(i):b(i),:,:),1);
        S_Med_l(MNT,:,:) = mean(Clim_Med_l(a(i):b(i),:,:),1);
        S_Med_d(MNT,:,:) = mean(Clim_Med_d(a(i):b(i),:,:),1);
        S_Lrg_p(MNT,:,:) = mean(Clim_Lrg_p(a(i):b(i),:,:),1);
        S_Lrg_l(MNT,:,:) = mean(Clim_Lrg_l(a(i):b(i),:,:),1);
        S_Lrg_d(MNT,:,:) = mean(Clim_Lrg_d(a(i):b(i),:,:),1);
    end

end %Years

%%% Save
save([fname '_locs.mat'],...
'S_Sml_f','S_Sml_m','S_Sml_p','S_Sml_l','S_Sml_d',...
'S_Med_f','S_Med_m','S_Med_p','S_Med_l','S_Med_d',...
'S_Lrg_p','S_Lrg_l','S_Lrg_d','S_Cobalt')

end
