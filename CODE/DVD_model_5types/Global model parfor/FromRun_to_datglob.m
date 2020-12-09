
% double click RUN_TOT_global_date
% select datcell

datglobsub = zeros(37911,30);

for i = 1:37911
    datglobsub(i,:) = datcell{i,1};
end

datglob = zeros(37911,9);

datglob(:,1) = datglobsub(:,1);
datglob(:,2) = datglobsub(:,2);
datglob(:,3) = datglobsub(:,3);
datglob(:,4) = datglobsub(:,4);
datglob(:,5) = sum(datglobsub(:,5:8),2);
datglob(:,6) = sum(datglobsub(:,9:12),2);
datglob(:,7) = sum(datglobsub(:,13:18),2);
datglob(:,8) = sum(datglobsub(:,19:24),2);
datglob(:,9) = sum(datglobsub(:,25:30),2);

csvwrite('datglob.csv',datglob)

