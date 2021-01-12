function[theta, depthDay, depthNight, avlocDay, avlocNight] = calcpreference(param)

   % calculate size-preference matrix
   %
   for i = 1:length(param.ixR)
    EncR(i,:) = zeros(length(param.wc),1)';
   end

   for i = 1:length(param.ixFish)
    thetaP = sqrt(pi/2)*param.sigma.* (...
     erf( (log(param.wu)-log(param.wc(param.ixFish(i))./param.beta))./(sqrt(2)*param.sigma) ) ...
     - erf( (log(param.w)-log(param.wc(param.ixFish(i))./param.beta))./(sqrt(2)*param.sigma) ));

    Enc(i,:) = thetaP./(log(param.wu)-log(param.w));
   end
   prefer = vertcat(EncR,Enc); 

   % calculate overlap from depth distribution
   %
   sigma = 10; % width of initial distribution
   tau = 10;  % increase in width
   sigmap = sigma + tau*log10(param.wc/param.wc(1)); % width for each size class
   xrange = linspace(0,param.bottom,(param.bottom+1));
   dvm = param.photic + 500; % vertical migration depth
   dvm(param.bottom < (param.photic + 500)) = param.bottom; % migration to bottom in intermediate habitats
   dvm(param.bottom <= param.mesop) = 0; % no migration in shallow habitats
   
   % first stages as juvenile/adult for predators
   [~,ixjuv] = min(abs(param.sizes-param.smat));
   [~,ixadult] = min(abs(param.sizes-param.lmat));
   
   % zooplankton night
   xloc = 0; 
   zp_n = (1./(sqrt(2*3.14*sigmap(param.ixR(1:2)).^2))) .* ...
            exp(-((xrange-xloc).^2' ./(2*sigmap(param.ixR(1:2)).^2)));
   zp_n = zp_n*diag(1./sum(zp_n,1));
        
   % zooplankton day (half at surface, half at dvm depth
   xloc = dvm;
   ix = param.ixR(1:2);
   zp_d = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc).^2' ./(2*sigmap(ix).^2)));
   zp_d = zp_d*diag(1./sum(zp_d,1));
   zp_d = (zp_n + zp_d)/2;

   % benthos small and large (at bottom with width sigma)
   xloc = param.bottom; 
   bent = (1/(sqrt(2*3.14*sigma^2))) .* exp(-((xrange-xloc).^2/(2*sigma^2)));
   bent = bent /sum(bent); 
   bent_dn = [bent' bent']; 

   % small pelagic fish (day + night) 
   xloc = 0; 
   ix = param.ix1(1):param.ix2(1);
   spe = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc).^2' ./(2*sigmap(ix).^2)));
   spel_dn = spe*diag(1./sum(spe,1));
   
   % meso pelagic night
   mpel_n = spel_dn;
   
   % meso pelagic day (dvm)
   xloc = dvm; 
   ix = param.ix1(2):param.ix2(2);
   mpe = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc).^2' ./(2*sigmap(ix).^2)));
   mpel_d = mpe*diag(1./sum(mpe,1));
      
   % large pelagic fish night (surface)
   xloc = 0; 
   ix = param.ix1(3):param.ix2(3);
   lpe = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc).^2' ./(2*sigmap(ix).^2)));
   lpe_n = lpe*diag(1./sum(lpe,1));
   
   % large pelagic fish day (surface + dvm)
   xloc = zeros(param.nstage,1)';
   xloc(ixadult:end) = dvm;
   ix = param.ix1(3):param.ix2(3);
   lpe = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc').^2' ./(2*sigmap(ix).^2)));
   lpe_d = lpe*diag(1./sum(lpe,1));
   lpe_d = (lpe_d + lpe_n)/2;
   
   % bathypelagic night (adults in midwater, others at surface)
   xloc = zeros(param.nstage,1)';
   xloc(ixadult:end) = dvm;
   ix = param.ix1(4):param.ix2(4);
   bpe = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc').^2' ./(2*sigmap(ix).^2)));
   bpel_n = bpe*diag(1./sum(bpe,1)); 
       
   % bathypelagic day (dvm)
   xloc = repmat(dvm,param.nstage,1)';
   ix = param.ix1(4):param.ix2(4);
   bpe = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc').^2' ./(2*sigmap(ix).^2)));
   bpel_d = bpe*diag(1./sum(bpe,1)); 
   
   % demersal fish night
   xloc = zeros(param.nstage,1)';
   xloc(ixjuv:end) = param.bottom;
   ix = param.ix1(5):param.ix2(5);
   dem = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc').^2' ./(2*sigmap(ix).^2)));
   dem_n = dem*diag(1./sum(dem,1));
   
   % demersal fish day
   demmig = dvm;
   demmig((param.bottom - dvm) >= 1200) = dvm + (param.bottom-dvm-1200);
   demmig((param.bottom - dvm) >= 1500) = param.bottom;
   xloc(ixadult:end) = demmig;
   ix = param.ix1(5):param.ix2(5);
   dem = (1./(sqrt(2*3.14*sigmap(ix).^2))) .* ...
            exp(-((xrange-xloc').^2' ./(2*sigmap(ix).^2)));
   dem_d = dem*diag(1./sum(dem,1));
 
   % calculate overlap during day
   %
   depthDay = [zp_d bent_dn spel_dn mpel_d lpe_d bpel_d dem_d];
   dayout = magic(param.ixFish(end)).*0;
   for i = 1:param.ixFish(end)
    test =  min(depthDay(1:end,i) , depthDay(1:end,1:end));
    dayout(1:end,i) = sum(test,1)'; 
   end
 
   % calculate overlap during night
   %
   depthNight = [zp_n bent_dn spel_dn mpel_n lpe_n bpel_n dem_n];
   nightout = magic(param.ixFish(end)).*0;
   for i = 1:param.ixFish(end)
    test =  min(depthNight(1:end,i) , depthNight(1:end,1:end));
    nightout(1:end,i) = sum(test,1)'; 
   end
   
   % visual predation is good at light, bad in the dark)
   visualpred = [(param.ix1(1):param.ix2(1)) (param.ix1(3):param.ix2(3))];
   dayout(visualpred,:) = dayout(visualpred,:)*param.visual;
   nightout(visualpred,:) = nightout(visualpred,:)*(2-param.visual);  
   
   % pelagic predators have limited vision in twilight zone during day
   pelpred = param.ix1(3):param.ix2(3);
   pelpred = pelpred(ixadult:end);
   preytwi = [param.ix1(2):param.ix2(2) param.ix1(4):param.ix2(4)];
   dayout(pelpred,preytwi) =  dayout(pelpred,preytwi)/param.visual*(2-param.visual);
   
   % average overlap during the whole day
   location = (dayout+nightout).*0.5;

   % calculate combined preference matrix   
   theta = prefer.*location;

   % change specific interactions
   %
   % benthivory is a strategy (larvae + pelagics do not eat benthos)
   idx_be = 5:(param.ix1(5)+(ixjuv-2)); % all larvae + pelagics
   theta(idx_be,3:4) = 0;   
   
   % small demersals are less preyed on
   idx_smd = param.ix1(5)+(ixjuv-1):param.ix1(5)+(ixadult-2);
   theta(idx_be,idx_smd)= theta(idx_be,idx_smd)*.25;
   
   % demersals do not eat much zooplankton
   idx_dems = param.ix1(5)+(ixjuv-1):param.ix2(5);
   theta(idx_dems,1:2)=  theta(idx_dems,1:2)*0;
   
   % provide benefit to forage and mesopelagic fish (predator avoidance)
   pred1 = param.ix1(3)+(ixadult-1):param.ix2(3);
   pred2 = param.ix1(4)+(ixadult-1):param.ix2(4);
   pred3 = param.ix1(5)+(ixadult-1):param.ix2(5);
   prey1 = param.ix1(1)+(ixjuv-1):param.ix2(1);
   prey2 = param.ix1(2)+(ixjuv-1):param.ix2(2);
   idx_predat = [pred1 pred2 pred3];
   idx_prey= [prey1 prey2];
   theta(idx_predat,idx_prey) = theta(idx_predat,idx_prey)*0.80;  
        
   % remove feeding on its own stage
   %idown = param.ix2(5);
   %idb = linspace(1,idown*idown,idown);
   %theta(theta==idb)=0;
   
   % calculate center of vertical distribution during the day
   [~,idi] = max(depthDay);
   avlocDay = xrange(idi);
   %
   % calculate center of vertical distribution during the night
   [~,idi] = max(depthNight);
   avlocNight = xrange(idi);
   
