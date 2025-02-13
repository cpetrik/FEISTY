% Test q10 fns


M_s = 10^((log10(0.001)+log10(0.5))/2);
M_m = 10^((log10(0.5)+log10(250))/2);
M_l = 10^((log10(250)+log10(125000))/2);

wgt = M_l;
temp = -2:33;
cmax = (exp(0.063*(temp-10.0)) .* 20 .* wgt^(-0.25));
cmax1 = (exp(0.0405*(temp-10.0)) .* 20 .* wgt^(-0.25));
cmax2 = (exp(0.0693*(temp-10.0)) .* 20 .* wgt^(-0.25));
cmax3 = (exp(0.0916*(temp-10.0)) .* 20 .* wgt^(-0.25));

%R(q10) = R .* q10 .* (temp-10.0)./10;
c2 = 20 .* wgt^(-0.25) .* 2 .^ ((temp-10.0)./10);
c1 = 20 .* wgt^(-0.25) .* 1.5 .^ ((temp-10.0)./10);
c3 = 20 .* wgt^(-0.25) .* 2.5 .^ ((temp-10.0)./10);

figure(1)
plot(temp,c1,'-r'); hold on;
plot(temp,c2,'-b'); hold on;
plot(temp,c3,'-m'); hold on;
plot(temp,cmax,'ko'); hold on;
plot(temp,cmax1,'ro'); hold on;
plot(temp,cmax2,'bo'); hold on;
plot(temp,cmax3,'mo'); hold on;

%% Find q10
%ks = [0.0405:0.01:0.1216];
ks = [0.0405,0.0555,0.0705,0.0855,0.1005,0.1155,0.1305]';
q10s = exp(10*ks)

%% BE as temp-dep
%Want BE to decr with temp
eff = 0.025;
BE = exp(0.063*(temp-10.0)) .* eff;
mx = exp(0.063*(33.1-10.0));
BE2 = (mx - exp(0.063*(temp-10.0))) .* eff;

figure(2)
plot(temp,BE,'-b'); hold on;
plot(temp,BE2,'-r'); hold on;

%% Comp my temp scaling to J&C 15
temp2 = temp+273;
Tref = 283;
E=0.6;
k=8.62e-5;
oth = (exp(0.063*(temp-10.0)));
met = (exp(0.0855*(temp-10.0)));
jc = exp((-1*E/k)*((1./temp2)-(1./Tref)));

figure(3)
plot(temp,oth,'-k'); hold on;
plot(temp,met,'-b'); hold on;
plot(temp,jc,'-r'); hold on;

c=(E/k)/Tref;
kjc=c./temp2;
q10jc = exp(10*kjc);


