% Test goodness of fit of the Gumbel distribution for modeling flood return
% periods. 
%
% Meg D. Fowler, 2017-10-09
%

%% Calculate return period 

%Fit data to extreme value distribution (Gumbel)
periods = [2,3,5,7,10:5:125];
eulers = 0.57721;       %Euler's constant, used in computing Gumbel Distribution

%Array of K-values used to find discharge associated with T-year floods
K      = -(sqrt(6)/pi)*(eulers + log(log(periods./(periods-1))));

caseNames = {'control','full','radiation','physiology'};

%Calculate area-weighted differences between each case for return period
CaMaFile = '/Users/meganfowler/gp_fuse/Flooding-physiology/MatlabData/fldare1161.nc'; 
clon = ncread(CaMaFile,'lon');      %Longitudes in original data file
clat = ncread(CaMaFile, 'lat');     %Latitudes in original data file 

for icase=1:1
   %recName = ['/Users/meganfowler/Documents/MATLAB/FloodProject/CaMaFlood/Physiology/30yr/',caseNames{icase},'_outflw.mat'];
   %load(recName);
   load('/Users/meganfowler/gp_fuse/Flooding-physiology/MatlabData/control_outflw.mat'); 
   
   %Old way of computing Gumbel dist - Method of Moments
   avgOutflw(icase,:,:) = avgVar_full;
   qbar(icase,:,:) = nanmean(maxYearly_full,3);
   sig(icase,:,:)  = std(maxYearly_full,0,3);  %flag=0 means default normalization by N-1; Use =1 to normalize by N
   %Compute discharge of extreme flood events with return periods of T
    for p=1:numel(K)
        qReturn(icase,:,:,p) = qbar(icase,:,:) + sig(icase,:,:)*K(p);
    end
   
   %Use Maximum Likelihood Estimation for parameters
   for i=1:numel(clon)
      for j=1:numel(clat)
          %rec   = squeeze(maxYearly_full(i,j,:));    %Isolate lat/lon record
          rec   = squeeze(maxYearly_full(i,j,:))*(1e9)/86400;    %Isolate lat/lon record in m^3/s
          
%           %% Chi-squared GOF test
%           if numel(find(isnan(rec)==0))==30
%               EV  = evfit(rec);
%               pdEV = makedist('ExtremeValue','mu',EV(1),'sigma',EV(2));
%               hEV_chi(icase,i,j) = chi2gof(rec,'CDF',pdEV);              
%               
%               GEV = gevfit(rec);
%               pdGEV = makedist('GeneralizedExtremeValue','k',GEV(1),'sigma',GEV(2),'mu',GEV(3));            
%               hGEV_chi(icase,i,j) = chi2gof(rec,'CDF',pdGEV); 
%           else
%               hEV_chi(icase,i,j) = NaN;
%               hGEV_chi(icase,i,j) = NaN;
%           end
              
          
%           %% Anderson-Darling tests
%           if numel(find(isnan(rec)==0))>=4 
%               hEV(icase,i,j) = adtest(rec,'Distribution','ev','alpha',0.1); %Test at 90% CI
%               
%               %GEV PD defined above, if not commented out
%               %pd = makedist('gev');     %Should this be made without specifying k,sigma,mu? 
%               hGEV(icase,i,j) = adtest(rec,'Distribution',pdGEV,'alpha',0.1);
%           else
%               hEV(icase,i,j) = NaN;
%               hGEV(icase,i,j) = NaN;
%           end
          
          %% PPCC tests
          if sum(isnan(rec))==0
              gumb  = evfit(rec);    %Estimate parameters of Gumbel dist
              mu    = gumb(1);         %Location Parameter
              sigma = gumb(2);      %Scale Parameter
              qs    = sort(rec);     %Sort observations to compute PPCC
              n     =length(qs);         
              p     =(1:n)/(n+1);        %Evaluate at range of values from 0:1
              qhat  = evinv(p,mu,sigma); %Inverse CDF of Gumbel
              ppcc(icase,i,j) = corr(qs,qhat'); %PPCC (via hydrology.usu.edu HW solutions)
              
              % Using method of moments instead...
              sigmaOld = (sqrt(6)/pi)*sig(icase,i,j);
              muOld    = qbar(icase,i,j) - sigmaOld*eulers;
              qhatOld = evinv(p,muOld,sigmaOld);
              ppccOld(icase,i,j) = corr(qs,qhatOld');

              % How about using L-moments? 
              M100 = 0;
              M110 = 0;
              recSort = sort(rec);
              for ir=1:numel(recSort)
                  M100 = M100+recSort(ir);
                  M110 = M110 + ((ir-1)/(numel(recSort)-1) * recSort(ir));
              end
              M100 = (1/numel(recSort))*M100;
              M110 = (1/numel(recSort))*M110;
              L1 = M100;
              L2 = 2*M110 - M100;
              Lmom_sigma = L2/log10(2);
              Lmom_mu = L1 - Lmom_sigma*eulers;

              qhatLmom = evinv(p,Lmom_mu,Lmom_sigma);
              ppccLmom(icase,i,j) = corr(qs,qhatLmom');
              
              % Do things improve if using the generalized extreme value
              % distriution? 
              gev      = gevfit(rec);
              k        = gev(1);
              sigmaGEV = gev(2);
              muGEV    = gev(3);
              qhatGEV  = gevinv(p,k,sigmaGEV,muGEV);
              ppccGEV(icase,i,j) = corr(qs,qhatGEV'); 

          else
              ppcc(icase,i,j)     = NaN;
              ppccOld(icase,i,j)  = NaN;
              ppccLmom(icase,i,j) = NaN;
              ppccGEV(icase,i,j)  = NaN; 
          end
          
          
          
      end
   end

    %clearvars maxYearly_full
    
end
    
% i100 = find(periods==100);
% qMan = squeeze(qReturn(1,:,:,i100));
% qAuto = squeeze(qReturnAuto(1,:,:,i100));
% diffPct = (qAuto-qMan)./qMan;

%Test case for manual calculation of PPCC
i=500;j=500;
rec   = squeeze(maxYearly_full(i,j,:));    %Isolate lat/lon record

%Sort in order 
rec = sort(rec);
b = pi/(std(rec)*sqrt(6));
a = eulers - nanmean(rec)*b;
for i=1:numel(rec)
    fHat(i) = (i-0.44)/(numel(rec)+0.12);
end

%OR, compute estimates via MLE and try again...
gumb  = evfit(rec);    %Estimate parameters of Gumbel dist
mu    = gumb(1);         %Location Parameter
sigma = gumb(2);      %Scale Parameter

%M = (-a - log(-log(fHat)))./b;
M  = (-mu - log(-log(fHat)))./sigma;

topFull=0; bottom1Full=0; bottom2Full=0; 
for i=1:numel(rec)
    top = (rec(i)-nanmean(rec))*(M(i)-nanmean(M));
    topFull = top+topFull;
    
    bottom1 = (rec(i)-nanmean(rec))^2;
    bottom1Full = bottom1+bottom1Full;
    
    bottom2 = (M(i)-nanmean(M))^2; 
    bottom2Full = bottom2+bottom2Full;
end
rHat = topFull/(bottom1Full*bottom2Full);

%Test estimating parameters via L-moments
M100 = 0;
M110 = 0;
rec = sort(rec);
for i=1:numel(rec)
    M100 = M100+rec(i);
    M110 = M110 + ((i-1)/(numel(rec)-1) * rec(i));
end

M100 = (1/numel(rec)).*M100;
M110 = (1/numel(rec)).*M110;
L1 = M100;
L2 = 2*M110 - M100;

Lmom_b = L2/log10(2);
Lmom_a = L1 - Lmom_b*eulers;
%% Plot PPCC

PPCC_color = [1 0 0;...         %red
              1 0.5 0;...       %orange 
              0 1 0;...         %green
              0 1 1;...         %cyan
              (51/255) (153/255) 1;... %Med blue
              0 0 1;...         %Blue
              0 0 (153/255)];   %Dark blue
          
ppcc1 = squeeze(ppcc(1,:,:));
ppcc2 = squeeze(ppccOld(1,:,:));
ppcc3 = squeeze(ppccGEV(1,:,:));
          
colors = NaN(size(ppcc1));
colors(ppcc3 <=0.8)               = 1;
colors(ppcc3>0.8  & ppcc3<= 0.9)  = 2;
colors(ppcc3>0.9  & ppcc3<=0.95)  = 3;
colors(ppcc3>0.95 & ppcc3<=0.97)  = 4;
colors(ppcc3>0.97 & ppcc3<=0.98)  = 5;
colors(ppcc3>0.98 & ppcc3<=0.99)  = 6;
colors(ppcc3>0.99)                = 7;

coast = load('coast.mat');
figure;
contourf(clon,clat,colors','LineColor','none');   %Change in return period from non-sp to SP 
hold on;
plot(coast.long, coast.lat, 'k', coast.long+360, coast.lat, 'k');
hold off
axis('xy','equal',[-180 180 -90 90]);
title('PPCC of GEV Dist in Control');
my_color=PPCC_color;
colormap(my_color);
%caxis([min(c) max(c)]);
hcb = colorbar('location', 'southoutside');
%set(hcb,'XTick',round(c));
hcb.Label.String = 'PPCC Test';
set(hcb,'ytick',linspace(1,7,8))
set(hcb,'yticklabel',{'<0.80','0.80','0.90','0.95','0.97','0.98','0.99','>0.99'});
set(gca,'FontSize',24);
set(hcb,'FontSize',24);
% ----------------- End Full




%% Use maximum likelihood estimate in Matlab to estimate distribution parameters instead 

%TEST CASE 
% i=500; j=500;
% rec = squeeze(maxYearly_full(i,j,:));
% gumb = evfit(rec);
% qs = sort(rec);
% mu = gumb(1); sigma = gumb(2);
% n=length(qs);
% p=(1:n)/(n+1);
% qhat = evinv(p,mu,sigma);
% ppcc = corr(qs,qhat');
% T = 100;
% q100 = mu - (sigma*log(-log(1-(1/T))));


