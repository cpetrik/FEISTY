% diff e-ratio calcs

% Laws et al. 2000, https://doi.org/10.1029/1999GB001229
% empirical model
PE_Laws00=(-0.02*SST)+0.6186

% Laws et al. 2000, https://doi.org/10.1029/1999GB001229
% food web model
PE_LawsFW=efintrpC(CHL,SST,PP)

% Laws et al. 2011, https://doi.org/10.4319/lom.2011.9.593 
PE_Laws11=0.04756*(0.78-(0.43*SST/30))*(PP)^0.307

% Henson et al. 2011, http://dx.doi.org/10.1029/2011GL046735 
PE_Henson11=0.23*exp(-0.08*SST);

% Dunne et al., 2005, https://doi.org/10.1029/2004GB002390 
PE_Dunne05=(-0.0081*SST)+(0.0806*log(CHL))+0.426;
PE_Dunne05(find(PE_Dunne05<0.04))=0.04;
PE_Dunne05(find(PE_Dunne05>0.72))=0.72;
