f0 = dlmread("flux0.dat",",");
f = dlmread("flux.dat",",");

wvls = 1./f(:,1);
R_meep = -f(:,2)./f(:,2);

eps_quartz = @(l) 1+(0.6961663*l.^2)./(l.^2-0.0684043^2)+(0.4079426*l.^2)./(l.^2-0.1162414^2)+(0.8974794*l.^2)./(l.^2-9.896161^2);
R_fresnel = @(l) abs((1-eps_quartz(l).^0.5)./(1+eps_quartz(l).^0.5)).^2;

plot(wvls,R_meep,'bo-',wvls,R_fresnel(wvls),'rs-');
xlabel("wavelength (μm)");
ylabel("reflectance");
legend("meep","analytic")
