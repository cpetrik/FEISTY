% Read GFDL netcdfs
% Regridded by code merging Xiao's with 1 line from Matthias
% Need to be re-oriented to match other ISIMIP files

clear
close all


%% loop over files, each one year (12 mo)
yr = 1961:2010;
ttot = 12*length(yr); %total number of months across all files
M = 0;

for y = 1:length(yr)
    Y = yr(y);

    %% Loop over mos
    for t=1:nt
        M = M+1;

        testM = double(squeeze(mz100(:,:,t)));
        testL = double(squeeze(lz100(:,:,t)));

        Mtsplit1 = testM(1:720,:);
        Mtsplit2 = testM(721:end,:);
        Mtflip1 = fliplr(Mtsplit1);
        Mtflip2 = fliplr(Mtsplit2);
        Mtcomb = [Mtflip2;Mtflip1];

        Ltsplit1 = testL(1:720,:);
        Ltsplit2 = testL(721:end,:);
        Ltflip1 = fliplr(Ltsplit1);
        Ltflip2 = fliplr(Ltsplit2);
        Ltcomb = [Ltflip2;Ltflip1];

        nmdz_100(:,:,M) = Mtcomb;
        nlgz_100(:,:,M) = Ltcomb;

        clear Mtsplit1 Mtsplit2 Mtflip1 Mtflip2 Mtcomb
        clear Ltsplit1 Ltsplit2 Ltflip1 Ltflip2 Ltcomb

    end

    clear mz100 lz100 nmdz nlgz

end %yrs







