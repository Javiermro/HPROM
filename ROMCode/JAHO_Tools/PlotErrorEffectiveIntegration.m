function PlotErrorEffectiveIntegration(igausEFF,error_INTEGRATION_nz)

if nargin == 0
    load('tmp.mat')
end

figure(1789)
hold on
%  error_INTEGRATION_nz = error_INTEGRATION_nz(error_INTEGRATION_nz>0) ;
ngausEFF =  max(igausEFF) ;
errorPLOT = zeros(ngausEFF,1) ;
for i=1:ngausEFF
   igg =  find(igausEFF==i)  ;
   if ~isempty(igg)
       errorPLOT(i) = min([error_INTEGRATION_nz(igg);errorPLOT(1:i-1)]); 
   else
       errorPLOT(i)  =errorPLOT(i-1) ;
   end
end

plot(errorPLOT)
xlabel('Number of points')
ylabel('norm(RESID/RESID_INI) (%)')

figure(1714)
hold on
plot(log10(errorPLOT/100))
xlabel('Number of points')
ylabel('log(RESID/RESID_INI)')

figure(3001)
hold on
plot(igausEFF)
xlabel('Number of input points')
ylabel('Number of points with w>0')