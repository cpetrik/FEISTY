%%%%!! RUN HISTORIC WITH FISHING FOR ALL LOCATIONS
function Historic_fished_gfdl_1meso_loop()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Grid
load('/Users/cpetrik/Dropbox/Princeton/POEM_2.0/CODE/Data/Data_grid_hindcast_NOTflipped.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;
param.NX = NX;
param.ID = ID;

%! How long to run the model
YEARS = 145;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Set fishing rate
param.frate = 0.3;
param.dfrate = param.frate/365.0;

lam = 0.55:0.025:0.675;
aM = 4.25:0.25:6;

for L = 1:length(lam)
    for m = 1:length(aM)
        
        param.Lambda = lam(L);
        param.amet = aM(m);
        
        %! Make core parameters/constants (global)
        param = make_parameters_1meso(param); % make core parameters/constants
        
        %! Create a directory for output
        [fname,simname] = sub_fname_hist_1meso(param);
        
        %! Storage variables
        S_Bent_bio = zeros(NX,DAYS);
        S_Mzoo_frac = zeros(NX,DAYS);
        
        S_Sml_f = zeros(NX,DAYS);
        S_Sml_p = zeros(NX,DAYS);
        S_Sml_d = zeros(NX,DAYS);
        S_Med_f = zeros(NX,DAYS);
        S_Med_p = zeros(NX,DAYS);
        S_Med_d = zeros(NX,DAYS);
        S_Lrg_p = zeros(NX,DAYS);
        S_Lrg_d = zeros(NX,DAYS);
        
        %! Dims of monthyl means
        nt = 12*YEARS;
        Mo_Sml_f.bio = NaN*ones(NX,nt);
        Mo_Sml_p.bio = NaN*ones(NX,nt);
        Mo_Sml_d.bio = NaN*ones(NX,nt);
        Mo_Med_f.bio = NaN*ones(NX,nt);
        Mo_Med_p.bio = NaN*ones(NX,nt);
        Mo_Med_d.bio = NaN*ones(NX,nt);
        Mo_Lrg_p.bio = NaN*ones(NX,nt);
        Mo_Lrg_d.bio = NaN*ones(NX,nt);
        Mo_Bent.bio = NaN*ones(NX,nt);
        Mo_Mzoo.frac = NaN*ones(NX,nt);
        
        
        %% ! Initialize
        %init_sim = simname;
        init_sim ='Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
        load(['/Volumes/FEISTY/NC/Matlab_new_size/',init_sim '/Last_mo_preindust_' init_sim '.mat']);
        BENT.mass = BENT.bio;
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_hist(ID,DAYS,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
        Med_d.td(1:NX) = 0.0;
        Lrg_d.td(1:NX) = 0.0;
        ENVR = sub_init_env_1meso(ID);
        
        %% %%%%%%%%%%%%%%%%%%%% Run the Model
        MNT = 0;
        %! Run model with no fishing
        for YR = 1:YEARS % years
            %! Load a year's COBALT data
            ti = num2str(YR+1860);
            load(['/Volumes/FEISTY/POEM_JLD/esm2m_hist/Data_ESM2Mhist_',ti,'.mat'],'COBALT');
            
            [num2str(param.Lambda),',',num2str(param.amet),',',num2str(YR),]
                
            for DAY = 1:param.DT:DAYS % days
                
                %%%! Future time step
                DY = int64(ceil(DAY));
                %[num2str(YR),' , ',num2str(mod(DY,365))]
                [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                    sub_futbio_1meso(ID,DY,COBALT,GRD,Sml_f,Sml_p,Sml_d,...
                    Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
                
                %! Store
                S_Bent_bio(:,DY) = BENT.mass;
                S_Mzoo_frac(:,DY) = ENVR.fZm;
                
                S_Sml_f(:,DY) = Sml_f.bio;
                S_Sml_p(:,DY) = Sml_p.bio;
                S_Sml_d(:,DY) = Sml_d.bio;
                S_Med_f(:,DY) = Med_f.bio;
                S_Med_p(:,DY) = Med_p.bio;
                S_Med_d(:,DY) = Med_d.bio;
                S_Lrg_p(:,DY) = Lrg_p.bio;
                S_Lrg_d(:,DY) = Lrg_d.bio;
                
            end %Days
            
            %! Calculate monthly means and save
            aa = (cumsum(MNTH)+1);
            a = [1,aa(1:end-1)]; % start of the month
            b = cumsum(MNTH); % end of the month
            for i = 1:12
                MNT = MNT+1; % Update monthly ticker
                
                Mo_Bent.bio(:,MNT) = mean(S_Bent_bio(:,a(i):b(i)),2);
                Mo_Mzoo.frac(:,MNT) = mean(S_Mzoo_frac(:,a(i):b(i)),2);
                
                Mo_Sml_f.bio(:,MNT) = mean(S_Sml_f(:,a(i):b(i)),2);
                Mo_Sml_p.bio(:,MNT) = mean(S_Sml_p(:,a(i):b(i)),2);
                Mo_Sml_d.bio(:,MNT) = mean(S_Sml_d(:,a(i):b(i)),2);
                Mo_Med_f.bio(:,MNT) = mean(S_Med_f(:,a(i):b(i)),2);
                Mo_Med_p.bio(:,MNT) = mean(S_Med_p(:,a(i):b(i)),2);
                Mo_Med_d.bio(:,MNT) = mean(S_Med_d(:,a(i):b(i)),2);
                Mo_Lrg_p.bio(:,MNT) = mean(S_Lrg_p(:,a(i):b(i)),2);
                Mo_Lrg_d.bio(:,MNT) = mean(S_Lrg_d(:,a(i):b(i)),2);
                
            end %Monthly mean
            
        end %Years
        
        %% Take mean of last 50 yrs and save
        %Time
        sp_tmean=nanmean(Mo_Sml_p.bio,1);
        sf_tmean=nanmean(Mo_Sml_f.bio,1);
        sd_tmean=nanmean(Mo_Sml_d.bio,1);
        mp_tmean=nanmean(Mo_Med_p.bio,1);
        mf_tmean=nanmean(Mo_Med_f.bio,1);
        md_tmean=nanmean(Mo_Med_d.bio,1);
        lp_tmean=nanmean(Mo_Lrg_p.bio,1);
        ld_tmean=nanmean(Mo_Lrg_d.bio,1);
        b_tmean=nanmean(Mo_Bent.bio,1);
        mz_tmfrac=nanmean(Mo_Mzoo.frac,1);
        
        y = 1860+(1/12):(1/12):2005;
        yr50=find(y>=1951 & y<2001);
        sp_mean50=nanmean(Mo_Sml_p.bio(:,yr50),2);
        sf_mean50=nanmean(Mo_Sml_f.bio(:,yr50),2);
        sd_mean50=nanmean(Mo_Sml_d.bio(:,yr50),2);
        mp_mean50=nanmean(Mo_Med_p.bio(:,yr50),2);
        mf_mean50=nanmean(Mo_Med_f.bio(:,yr50),2);
        md_mean50=nanmean(Mo_Med_d.bio(:,yr50),2);
        lp_mean50=nanmean(Mo_Lrg_p.bio(:,yr50),2);
        ld_mean50=nanmean(Mo_Lrg_d.bio(:,yr50),2);
        b_mean50=nanmean(Mo_Bent.bio(:,yr50),2);
        mz_mfrac50=nanmean(Mo_Mzoo.frac(:,yr50),2);
        
        MZ.over = nan*ones(size(Mo_Mzoo.frac));
        MZ.over(Mo_Mzoo.frac > 1) = ones;
        MZ.over(Mo_Mzoo.frac <= 1) = zeros;
        mz_ttf=nansum(MZ.over,1);
        mz_mtf50=nansum(MZ.over(:,yr50),2);
        
        save([fname,'_Means.mat'],...
            'sf_tmean','sp_tmean','sd_tmean',...
            'mf_tmean','mp_tmean','md_tmean',...
            'lp_tmean','ld_tmean','b_tmean',...
            'sf_mean50','sp_mean50','sd_mean50',...
            'mf_mean50','mp_mean50','md_mean50',...
            'lp_mean50','ld_mean50','b_mean50',...
            'mz_tmfrac','mz_mfrac50','mz_ttf','mz_mtf50')
        
    end
end


end
