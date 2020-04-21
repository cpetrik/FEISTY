function netcdf_read_fore_fished_cmax_ens(fname,simname)

  % FEISTY output at all locations

  %% SP
  ncid = netcdf.open([fname '_cmax_sml_p.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
  end
  netcdf.close(ncid);

  [ni,nt] = size(con);

  cmax = con./clev;
  SP.con  = con;
  SP.clev = clev;
  SP.cmax = cmax;
  clear cmax con clev

  %% SF
  ncid = netcdf.open([fname '_cmax_sml_f.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
  end
  netcdf.close(ncid);

  cmax = con./clev;
  SF.con  = con(:,1:nt);
  SF.clev = clev(:,1:nt);
  SF.cmax = cmax(:,1:nt);
  clear cmax con clev

  % SD
  ncid = netcdf.open([fname '_cmax_sml_d.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
  end
  netcdf.close(ncid);

  cmax = con./clev;
  SD.con  = con;
  SD.clev = clev;
  SD.cmax = cmax;
  clear cmax con clev

  %% MP
  ncid = netcdf.open([fname '_cmax_med_p.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
  end
  netcdf.close(ncid);

  cmax = con./clev;
  MP.con  = con;
  MP.clev = clev;
  MP.cmax = cmax;
  clear cmax con clev

  % MF
  ncid = netcdf.open([fname '_cmax_med_f.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
  end
  netcdf.close(ncid);

  cmax = con./clev;
  MF.con  = con;
  MF.clev = clev;
  MF.cmax = cmax;
  clear cmax con clev

  % MD
  ncid = netcdf.open([fname '_cmax_med_d.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
  end
  netcdf.close(ncid);

  cmax = con./clev;
  MD.con  = con;
  MD.clev = clev;
  MD.cmax = cmax;
  clear cmax con clev

  % LP
  ncid = netcdf.open([fname '_cmax_lrg_p.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
  end
  netcdf.close(ncid);

  cmax = con./clev;
  LP.con  = con;
  LP.clev = clev;
  LP.cmax = cmax;
  clear cmax con clev

  % LD
  ncid = netcdf.open([fname '_cmax_lrg_d.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' >= 9.9692e+36) = NaN;']);
  end
  netcdf.close(ncid);

  cmax = con./clev;
  LD.con  = con;
  LD.clev = clev;
  LD.cmax = cmax;
  clear cmax con clev

  % Benthic material (Time)
  ncid = netcdf.open([fname '_bent.nc'],'NC_NOWRITE');
  [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
  for i = 1:nvars
      varname = netcdf.inqVar(ncid, i-1);
      eval([ varname ' = netcdf.getVar(ncid,i-1);']);
      eval([ varname '(' varname ' == 99999) = NaN;']);
  end
  netcdf.close(ncid);
  clear biomass

  %% Take means and totals

  %Time
  sp_tcmax=mean(SP.cmax,1);
  sf_tcmax=mean(SF.cmax,1);
  sd_tcmax=mean(SD.cmax,1);
  mp_tcmax=mean(MP.cmax,1);
  mf_tcmax=mean(MF.cmax,1);
  md_tcmax=mean(MD.cmax,1);
  lp_tcmax=mean(LP.cmax,1);
  ld_tcmax=mean(LD.cmax,1);

  %% 50 yrs (2051-2100)
  yr50=time((end-(50*12)+1):end);
  sp_cmax50=nanmean(SP.cmax(:,yr50),2);
  sf_cmax50=nanmean(SF.cmax(:,yr50),2);
  sd_cmax50=nanmean(SD.cmax(:,yr50),2);
  mp_cmax50=nanmean(MP.cmax(:,yr50),2);
  mf_cmax50=nanmean(MF.cmax(:,yr50),2);
  md_cmax50=nanmean(MD.cmax(:,yr50),2);
  lp_cmax50=nanmean(LP.cmax(:,yr50),2);
  ld_cmax50=nanmean(LD.cmax(:,yr50),2);

  %% Every 5 years
  st=1:60:length(time);
  en=60:60:length(time);

  for n=1:length(st)
      sp_cmax(:,n)=nanmean(SP.cmax(:,st(n):en(n)),2);
      sf_cmax(:,n)=nanmean(SF.cmax(:,st(n):en(n)),2);
      sd_cmax(:,n)=nanmean(SD.cmax(:,st(n):en(n)),2);
      mp_cmax(:,n)=nanmean(MP.cmax(:,st(n):en(n)),2);
      mf_cmax(:,n)=nanmean(MF.cmax(:,st(n):en(n)),2);
      md_cmax(:,n)=nanmean(MD.cmax(:,st(n):en(n)),2);
      lp_cmax(:,n)=nanmean(LP.cmax(:,st(n):en(n)),2);
      ld_cmax(:,n)=nanmean(LD.cmax(:,st(n):en(n)),2);

  end

  %%
  save([fname '_Means_cmax_' simname '.mat'],'time',...
      'yr50',...
      'sf_tcmax','sp_tcmax','sd_tcmax',...
      'mf_tcmax','mp_tcmax','md_tcmax',...
      'lp_tcmax','ld_tcmax',...
      'sf_cmax50','sp_cmax50','sd_cmax50',...
      'mf_cmax50','mp_cmax50','md_cmax50',...
      'lp_cmax50','ld_cmax50',...
      'sf_cmax','sp_cmax','sd_cmax',...
      'mf_cmax','mp_cmax','md_cmax',...
      'lp_cmax','ld_cmax');

end
