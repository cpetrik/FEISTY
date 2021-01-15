%%%%!! RUN Climatol Single Locations
function Locs_clim_fishvec()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Make core parameters/constants (global)
param = [];
param = make_parameters_fishvec(param);

%! Setup Climatol (loop 5-year climatology of ESM2.6-COBALT)
load('/Volumes/FEISTY/POEM_JLD/esm26_hist/ESM26_1deg_5yr_clim_191_195_daily.mat','COBALT');

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
[fname] = sub_fname_fishvec(param);

%! Storage
Clim_bio = NaN*ones(DAYS,NX,11);
Clim_I = Clim_bio;
Clim_nu = Clim_bio;
Clim_gamma = Clim_bio;
Clim_die = Clim_bio;
Clim_rep = Clim_bio;
Clim_rec = Clim_bio;
Clim_clev = Clim_bio;
Clim_prod = Clim_bio;
Clim_pred = Clim_bio;
Clim_nmort = Clim_bio;
Clim_met = Clim_bio;
Clim_caught = Clim_bio;
Clim_fmort = Clim_bio;
Clim_Cobalt = NaN*ones(DAYS,NX,5);

S_bio = NaN*ones(12*YEARS,NX,11);
S_I = S_bio;
S_nu = S_bio;
S_gamma = S_bio;
S_die = S_bio;
S_rep = S_bio;
S_rec = S_bio;
S_clev = S_bio;
S_prod = S_bio;
S_pred = S_bio;
S_nmort = S_bio;
S_met = S_bio;
S_caught = S_bio;
S_fmort = S_bio;
S_Cobalt = NaN*ones(12*YEARS,NX,5);


%! Initialize
[fish,BENT] = sub_init_fishvec(ID,param);

%% %%%%%%%%%%%%%%%%%%%% Run the Model %%%%%%%%%%%%%%%%%%%%
%! Run model
MNT=0;
for YR = 1:YEARS % years
    
    num2str(YR)
    
    for DAY = 1:param.DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        %[num2str(YR),' , ', num2str(mod(DY,365))]
        
        [fish,BENT,ENVR] = sub_futbio_fishvec(ID,DY,COBALT,GRD,fish,BENT,param);
        
        %! Store last year
        %if (YR==YEARS)
        Clim_bio(DY,:,:) = fish.bio;
        Clim_I(DY,:,:) = fish.I;
        Clim_nu(DY,:,:) = fish.nu;
        Clim_gamma(DY,:,:) = fish.gamma;
        Clim_die(DY,:,:) = fish.die;
        Clim_rep(DY,:,:) = fish.rep;
        Clim_rec(DY,:,:) = fish.rec;
        Clim_clev(DY,:,:) = fish.clev;
        Clim_prod(DY,:,:) = fish.prod;
        Clim_pred(DY,:,:) = fish.pred;
        Clim_nmort(DY,:,:) = fish.nmort;
        Clim_met(DY,:,:) = fish.met;
        Clim_caught(DY,:,:) = fish.caught;
        Clim_fmort(DY,:,:) = fish.fmort;
        Clim_Cobalt(DY,:,:) = [BENT.mass BENT.pred ENVR.fZm ENVR.fZl ENVR.fB];
        %end
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        S_Cobalt(MNT,:,:) = mean(Clim_Cobalt(a(i):b(i),:,:),1);
        S_bio(MNT,:,:) = mean(Clim_bio(a(i):b(i),:,:),1);
        S_I(MNT,:,:) = mean(Clim_I(a(i):b(i),:,:),1);
        S_nu(MNT,:,:) = mean(Clim_nu(a(i):b(i),:,:),1);
        S_gamma(MNT,:,:) = mean(Clim_gamma(a(i):b(i),:,:),1);
        S_die(MNT,:,:) = mean(Clim_die(a(i):b(i),:,:),1);
        S_rep(MNT,:,:) = mean(Clim_rep(a(i):b(i),:,:),1);
        S_rec(MNT,:,:) = mean(Clim_rec(a(i):b(i),:,:),1);
        S_clev(MNT,:,:) = mean(Clim_clev(a(i):b(i),:,:),1);
        S_prod(MNT,:,:) = mean(Clim_prod(a(i):b(i),:,:),1);
        S_pred(MNT,:,:) = mean(Clim_pred(a(i):b(i),:,:),1);
        S_nmort(MNT,:,:) = mean(Clim_nmort(a(i):b(i),:,:),1);
        S_met(MNT,:,:) = mean(Clim_met(a(i):b(i),:,:),1);
        S_caught(MNT,:,:) = mean(Clim_caught(a(i):b(i),:,:),1);
        S_fmort(MNT,:,:) = mean(Clim_fmort(a(i):b(i),:,:),1);
    end
    
end %Years

%%% Save
save([fname '_locs.mat'],'S_Cobalt',...
    'S_bio','S_I','S_nu','S_gamma','S_die','S_rep','S_rec','S_clev',...
    'S_prod','S_pred','S_nmort','S_met','S_caught','S_fmort');

end
