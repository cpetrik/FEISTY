%%%%!! RUN Climatol FOR ALL LOCATIONS
function Climatol_fishing_diffF()

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/FEISTY/FEISTY_clim/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');

%! How long to run the model
YEARS = 150;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! choose where and when to run the model
load('ESM26_1deg_5yr_clim_191_195_grid.mat','GRD');
NX = length(GRD.Z);
ID = 1:NX;
param.NX = NX;
param.ID = ID;

load('/Volumes/FEISTY/NC/Matlab_new_size/bio_rates/LHS_diffF.mat','X');

% Search mult combinations of F on 3 functional types
for j = 19:length(X)
    
    %! Change individual parameters
    pset = X(j,:);
    param.MFsel = pset(1);
    param.LPsel = pset(2);
    param.LDsel = pset(3);
    
    %! Make core parameters/constants (global)
    param = make_parameters(param);
    
    %! Create a directory for output
    fname = sub_fname_Ftest(param);
    
    %! Storage variables
    % Dims
    nt = 12;
    
    S_Bent_bio = zeros(NX,DAYS);
    
    S_Sml_f = zeros(NX,DAYS);
    S_Sml_p = zeros(NX,DAYS);
    S_Sml_d = zeros(NX,DAYS);
    S_Med_f = zeros(NX,DAYS);
    S_Med_p = zeros(NX,DAYS);
    S_Med_d = zeros(NX,DAYS);
    S_Lrg_p = zeros(NX,DAYS);
    S_Lrg_d = zeros(NX,DAYS);
    
    S_Med_f_fish = zeros(NX,DAYS);
    S_Med_p_fish = zeros(NX,DAYS);
    S_Med_d_fish = zeros(NX,DAYS);
    S_Lrg_p_fish = zeros(NX,DAYS);
    S_Lrg_d_fish = zeros(NX,DAYS);
    
    Clim_Sml_f.bio = NaN*ones(NX,nt);
    Clim_Sml_p.bio = NaN*ones(NX,nt);
    Clim_Sml_d.bio = NaN*ones(NX,nt);
    Clim_Med_f.bio = NaN*ones(NX,nt);
    Clim_Med_p.bio = NaN*ones(NX,nt);
    Clim_Med_d.bio = NaN*ones(NX,nt);
    Clim_Lrg_p.bio = NaN*ones(NX,nt);
    Clim_Lrg_d.bio = NaN*ones(NX,nt);
    Clim_Bent.bio = NaN*ones(NX,nt);
    
    Clim_Med_f.yield = NaN*ones(NX,nt);
    Clim_Med_p.yield = NaN*ones(NX,nt);
    Clim_Med_d.yield = NaN*ones(NX,nt);
    Clim_Lrg_p.yield = NaN*ones(NX,nt);
    Clim_Lrg_d.yield = NaN*ones(NX,nt);
    
    %! Initialize
    [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,DAYS);
    Med_d.td(1:NX) = 0.0;
    Lrg_d.td(1:NX) = 0.0;
    
    %% %%%%%%%%%%%%%%%%%%%% Run the Model
    %! Run model with no fishing
    MNT=0;
    for YR = 1:YEARS % years
        [num2str(j),' , ', num2str(YR)]
        
        for DAY = 1:param.DT:DAYS % days
            
            %%%! Future time step
            DY = int64(ceil(DAY));
            %[num2str(j),' , ', num2str(YR),' , ', num2str(mod(DY,365))]
            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = ...
                sub_futbio(ID,DY,COBALT,GRD,Sml_f,Sml_p,Sml_d,...
                Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param.dfrate,param);
            
            if (YR==YEARS)
                S_Bent_bio(:,DY) = BENT.mass;
                
                S_Sml_f(:,DY) = Sml_f.bio;
                S_Sml_p(:,DY) = Sml_p.bio;
                S_Sml_d(:,DY) = Sml_d.bio;
                S_Med_f(:,DY) = Med_f.bio;
                S_Med_p(:,DY) = Med_p.bio;
                S_Med_d(:,DY) = Med_d.bio;
                S_Lrg_p(:,DY) = Lrg_p.bio;
                S_Lrg_d(:,DY) = Lrg_d.bio;
                
                S_Med_f_fish(:,DY) = Med_f.caught;
                S_Med_p_fish(:,DY) = Med_p.caught;
                S_Med_d_fish(:,DY) = Med_d.caught;
                S_Lrg_p_fish(:,DY) = Lrg_p.caught;
                S_Lrg_d_fish(:,DY) = Lrg_d.caught;
            end
            
        end %Days
        
    end %Years
    
    %! Calculate monthly means after final year and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        Clim_Bent.bio(:,i) = mean(S_Bent_bio(:,a(i):b(i)),2);
        Clim_Sml_f.bio(:,i) = mean(S_Sml_f(:,a(i):b(i)),2);
        Clim_Sml_p.bio(:,i) = mean(S_Sml_p(:,a(i):b(i)),2);
        Clim_Sml_d.bio(:,i) = mean(S_Sml_d(:,a(i):b(i)),2);
        Clim_Med_f.bio(:,i) = mean(S_Med_f(:,a(i):b(i)),2);
        Clim_Med_p.bio(:,i) = mean(S_Med_p(:,a(i):b(i)),2);
        Clim_Med_d.bio(:,i) = mean(S_Med_d(:,a(i):b(i)),2);
        Clim_Lrg_p.bio(:,i) = mean(S_Lrg_p(:,a(i):b(i)),2);
        Clim_Lrg_d.bio(:,i) = mean(S_Lrg_d(:,a(i):b(i)),2);
        
        Clim_Med_f.yield(:,i) = mean(S_Med_f_fish(:,a(i):b(i)),2);
        Clim_Med_p.yield(:,i) = mean(S_Med_p_fish(:,a(i):b(i)),2);
        Clim_Med_d.yield(:,i) = mean(S_Med_d_fish(:,a(i):b(i)),2);
        Clim_Lrg_p.yield(:,i) = mean(S_Lrg_p_fish(:,a(i):b(i)),2);
        Clim_Lrg_d.yield(:,i) = mean(S_Lrg_d_fish(:,a(i):b(i)),2);
        
    end %Monthly mean
    %%% Save
    save([fname,'.mat'],...
        'Clim_Sml_f','Clim_Sml_p','Clim_Sml_d','Clim_Med_f',...
        'Clim_Med_p','Clim_Med_d','Clim_Lrg_p','Clim_Lrg_d',...
        'Clim_Bent')
    
end

end
