%%%%!! RUN SPINUP FOR ONE LOCATION, ALL DEPTHS
function BATS_1D_spinup()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Set fishing rate
param.frate = 0.3;
param.dfrate = param.frate/365.0;
param.dfrateF = nan;
param.dfrateP = nan;
param.dfrateD = nan;

%! Make core parameters/constants
param = make_parameters(param);

%! Grids
vpath = '/Volumes/petrik-lab/Feisty/';
%vpath = '/project/Feisty/GCM_Data/';

%1-D
load([vpath 'Data_grid_ocean_cobalt_1D.mat'],'GRD');
GRD1 = GRD;

param.NX = GRD.NID;
param.ID = 1:GRD.NID;

%! How long to run the model
YEARS = 100;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
exper = 'BATS_spinup_COBALT2003';
opath = '/Volumes/petrik-lab/Feisty/NC/MOM6-1D/BATS_vert/offline_feisty/';
%opath = '/project/Feisty/NC/MOM6-1D/BATS_vert/offline_feisty/';
[fname,simname,sname] = sub_fname_spin(param,opath,exper);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ids,DAYS);

%! Storage
D_Sml_f = zeros(NX,DAYS);
D_Sml_p = zeros(NX,DAYS);
D_Sml_d = zeros(NX,DAYS);
D_Med_f = zeros(NX,DAYS);
D_Med_p = zeros(NX,DAYS);
D_Med_d = zeros(NX,DAYS);
D_Lrg_p = zeros(NX,DAYS);
D_Lrg_d = zeros(NX,DAYS);
D_Bent_bio = zeros(NX,DAYS);

S_Sml_f = zeros(NX,12*YEARS);
S_Sml_p = zeros(NX,12*YEARS);
S_Sml_d = zeros(NX,12*YEARS);
S_Med_f = zeros(NX,12*YEARS);
S_Med_p = zeros(NX,12*YEARS);
S_Med_d = zeros(NX,12*YEARS);
S_Lrg_p = zeros(NX,12*YEARS);
S_Lrg_d = zeros(NX,12*YEARS);
S_Bent_bio = zeros(NX,12*YEARS);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
%! Iterate forward in time 
load([vpath,'COBALT_2023_10_spinup_2003_subset.mat'],'COBALT');

MNT=0;
for YR = 1:YEARS % years
    num2str(YR)

    for DAY = 1:DAYS % days
        
        %%%! ticker
        DY = int64(ceil(DAY));
        
        %%%! Future time step
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio(DY,COBALT,GRD1,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
        
        %! Daily for one year
        D_Bent_bio(:,DY) = BENT.mass;
        D_Sml_f(:,DY) = Sml_f.bio;
        D_Sml_p(:,DY) = Sml_p.bio;
        D_Sml_d(:,DY) = Sml_d.bio;
        D_Med_f(:,DY) = Med_f.bio;
        D_Med_p(:,DY) = Med_p.bio;
        D_Med_d(:,DY) = Med_d.bio;
        D_Lrg_p(:,DY) = Lrg_p.bio;
        D_Lrg_d(:,DY) = Lrg_d.bio;
        
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        S_Sml_f(:,MNT) = mean(D_Sml_f(:,a(i):b(i)),1);
        S_Sml_p(:,MNT) = mean(D_Sml_p(:,a(i):b(i)),1);
        S_Sml_d(:,MNT) = mean(D_Sml_d(:,a(i):b(i)),1);
        S_Med_f(:,MNT) = mean(D_Med_f(:,a(i):b(i)),1);
        S_Med_p(:,MNT) = mean(D_Med_p(:,a(i):b(i)),1);
        S_Med_d(:,MNT) = mean(D_Med_d(:,a(i):b(i)),1);
        S_Lrg_p(:,MNT) = mean(D_Lrg_p(:,a(i):b(i)),1);
        S_Lrg_d(:,MNT) = mean(D_Lrg_d(:,a(i):b(i)),1);
    end
    
end %Years

%%% Save
save([fname '.mat'],...
'S_Sml_f','S_Sml_p','S_Sml_d','S_Med_f','S_Med_p','S_Med_d',...
'S_Lrg_p','S_Lrg_d')

end
