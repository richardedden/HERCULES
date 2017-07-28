runme_iv;
%%%%%%%
%I guess the first step is to set up a 2D basis set which has dimensions
%(n_pointsx4)xN_met
%i.e. Hadamard spectra are concatenated.
metabolites={'Cho','Cr','NAA', 'Lac', 'GABA'};
%n_metabolites=28;
%weights=ones(1,n_metabolites);
weights=[1 1 7 2 3 3 10 3 3 3 1 1 1 1 4 1 17 17 5 5 3 1 1 1 4 1 0.1 ];

% 
 

load ../SIMULATIONS.nosync/HERMES5_TE80_EP20_20170525/sim_IndividualMetabs.mat;
%metabolites={'Cho','Cr','NAA', 'Lac', 'GABA'};
%%%
metabolites={'Asc','Asp', 'Cr', 'GABA', 'Gln1','Gln2','Glu','GPC1','GPC2','GPC3','GSH1','GSH2','GSH3','H20','Ins','Lac','NAA1','NAA2','NAAG1','NAAG2','NAAG3','PCh1','PCh2','PCr1','PCr2','Scyllo','bHG'};
n_metabolites=27;
%weights=ones(1,n_metabolites);`
clear ANames
ANames=who('A_*');
%%%
BasisSet=zeros(length(ppm),n_metabolites,4);

eval(['ppm=' ANames{1} '.ppm;']);
%Restrict fit range

ppm(1)
ppm(end)
size(ppm)
%stophere
fit_range = ppm >= 1 & ppm <= 4.5;
ppm=ppm(fit_range);
BasisSet=zeros(sum(fit_range),n_metabolites,4);

for ii=1:n_metabolites
    ii
    eval(['BasisSet(:,ii,1)=' ANames{ii} '.specs(fit_range);']);
    eval(['BasisSet(:,ii,2)=B' ANames{ii}(2:end) '.specs(fit_range);']);
    eval(['BasisSet(:,ii,3)=C' ANames{ii}(2:end) '.specs(fit_range);']);
    eval(['BasisSet(:,ii,4)=D' ANames{ii}(2:end) '.specs(fit_range);']);
    dummy=squeeze(BasisSet(:,ii,:))*hadamard(4);
    BasisSet(:,ii,:)=dummy;
end



BasisSet=real(permute(BasisSet,[1,3,2]));
BasisSet=reshape(BasisSet,[length(BasisSet(:))/n_metabolites n_metabolites]);



% %Set up a fake spectrum
% %simulate noise of same size as the spectrum
% Noise=0.4*randn(length(ppm)*4,1);
% weights_sim=weights.*(1+0.1*randn(n_metabolites,1).');
% 
% p1=0.05;
% p2=0.66;
% width=p1-p1*p2;
% widthG=p1*p2;
% 
% lineshape = conv(LorentzianModel_herc([1 width 0 0 0],0.5:-0.01:-0.5).',GaussModel_herc([1 widthG 0 0 0],0.5:-0.01:-0.5).','same');
% lineshape =lineshape/sum(lineshape);
% baseline_c=[0 0.01 0.04 -0.2];
% baseline_m=[0.01 0.0 0.004 -0.03];
% baseline=[(ppm-3)*baseline_m(1)+baseline_c(1) (ppm-3)*baseline_m(2)+baseline_c(2) (ppm-3)*baseline_m(3)+baseline_c(3) (ppm-3)*baseline_m(4)+baseline_c(4)].';
% %Add additional basline curvature
% baseline=baseline+4*([8*sin(ppm*pi/3+1).'; 4*sin(ppm*pi/2.3+2).'; 4*sin(ppm*pi/8+2).'; -4*sin(ppm*pi/1.58+2).']);
% baseline=baseline';
% size(baseline')
% size(Noise)
% size(conv(sum(BasisSet.*repmat(weights_sim,[size(BasisSet,1) 1]),2),lineshape,'same'))
% Sim_Spectrum = conv(sum(BasisSet.*repmat(weights_sim,[size(BasisSet,1) 1]),2),lineshape,'same')+Noise+baseline';
% %Sim_Spectrum = sum(BasisSet.*repmat(weights_sim,[size(BasisSet,1) 1]),2)+Noise+baseline;
% Sim_Spectrum =real(Sim_Spectrum);
%NEXT STEP IS TO INTEGRATE MORE PARAMETERS


%StackPlot(Sim_Spectrum(end:-1:1,:),40)

%pause(10)
% hold on
% plot(lineshape)
% plot(Sim_Spectrum);
% 
% pause(5)
%set(gca,'XDir','reverse');

Sim_Spectrum =real(M.spec.H4(end:-1:1,:));
%Sim_Spectrum =Sim_Spectrum(end:-1:1,:);
Sim_Spectrum=Sim_Spectrum(:);


%Work on a linear simulation
%Drop in code from GannetFit
size(Sim_Spectrum)

weights_init=[weights zeros(1,8) 0.1 0.5 0];


        lb = [zeros(1,n_metabolites) ones(1,8)*(-10000)  0.01 0]; %NP; our bounds are 0.03 less due to creatine shift
        ub = [ones(1,n_metabolites)*20 ones(1,8)*10000 1 1];
        options = optimset('lsqcurvefit');
        options = optimset(options,'Display','off','TolFun',1e-10,'Tolx',1e-10,'MaxIter',1e5);
        nlinopts = statset('nlinfit');
        nlinopts = statset(nlinopts, 'MaxIter', 1e5);
 
        Fit_Spectrum2=Hercules_model(weights_init,BasisSet);

        Sim_Spectrum=Sim_Spectrum/max(Sim_Spectrum(:))*1E6;

%  StackPlot(Sim_Spectrum,40)
%  hold on
%  StackPlot(Fit_Spectrum2,40,'r')
% % %StackPlot(Fit_Spectrum,1,'b')
%  hold off     
% 
plot(1:(sum(fit_range)*4),Sim_Spectrum,1:(sum(fit_range)*4),Fit_Spectrum2);
%stophere
% pause(10)
        %Fitting happens here
         [Fitted_weights,resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) Hercules_model(xdummy,ydummy), ...
            weights_init, BasisSet,Sim_Spectrum, ...
            lb,ub,options);

       % Sim_Spectrum = sum(BasisSet.*repmat(weights_sim,[size(BasisSet,1) 1]),2)+Noise;

%plot(1:2004,Sim_Spectrum,1:2004,sum(BasisSet.*repmat(Fitted_weights,[size(BasisSet,1) 1]),2));

% Sim_Spectrum=reshape(Sim_Spectrum,[,4]);
% %linear=linspace(-0.5,0.5,size(BasisSet,1)/4);
% %fit_baseline = [linear*Fitted_weights(10)+Fitted_weights(6) linear*Fitted_weights(11)+Fitted_weights(7) linear*Fitted_weights(12)+Fitted_weights(8) linear*Fitted_weights(13)+Fitted_weights(9)]
% 
 Fit_Spectrum=reshape(Hercules_model(Fitted_weights,BasisSet),[length(ppm),4]);;
% %hold on
% %StackPlot(Sim_Spectrum,1)
% %StackPlot(Fit_Spectrum,1,'r')
% %hold off
spline_vault=zeros([length(ppm),4,5]);
%Introduce iterative baseline spline
weights_init=Fitted_weights;
for ii=1:5
    Sim_Spectrum=reshape(Sim_Spectrum,[length(ppm),4]);
    
    
    Residuals=Sim_Spectrum-Fit_Spectrum;
    %focus on fitting the first bit

    
% 
    spline1=fit(ppm',Residuals(1:length(ppm))','smoothingspline','SmoothingParam',0.5);
    spline2=fit(ppm',Residuals((1:length(ppm))+length(ppm))','smoothingspline','SmoothingParam',0.7);
    spline3=fit(ppm',Residuals((1:length(ppm))+length(ppm)*2)','smoothingspline','SmoothingParam',0.5);
    spline4=fit(ppm',Residuals((1:length(ppm))+length(ppm)*3)','smoothingspline','SmoothingParam',0.5);
% 
%     
    spline_vault(:,1,ii)=feval(spline1,ppm);
    spline_vault(:,2,ii)=feval(spline2,ppm);
    spline_vault(:,3,ii)=feval(spline3,ppm);
    spline_vault(:,4,ii)=feval(spline4,ppm);
    
    %plot(ppm,spline_vault(:,2,ii),ppm,Residuals((1:length(ppm))+length(ppm)));
    %title(['Residual SD:' num2str(std(Residuals(:)))]);
    descent(ii)=std(Residuals(:));
    %pause(4)
%     %plot(ppm,Residuals(1:501));
%     %hold on
%     %plot(ppm,feval(spline,ppm))
%     %hold off
    Sim_Spectrum=Sim_Spectrum-squeeze(spline_vault(:,:,ii));
    Sim_Spectrum=Sim_Spectrum(:);
%     size(Sim_Spectrum)
    %loop through the fit again
    if(ii>1)
       weights_init=Fitted_weights(ii-1,:);
       %weights_init(end)=0.5;
       %weights_init(end-1)=0.1;
    end
    [Fitted_weights(ii,:),resnorm,residg] = lsqcurvefit(@(xdummy,ydummy) Hercules_model(xdummy,ydummy), ...
                weights_init, BasisSet,Sim_Spectrum, ...
                lb,ub,options);

%           
%      
% 
%     
%     %Fit_Spectrum2(1:501)=Fit_Spectrum2(1:501)+feval(spline,ppm)';
Fit_Spectrum=reshape(Hercules_model(Fitted_weights(ii,:),BasisSet),[length(ppm),4]);;
end
Sim_Spectrum=reshape(Sim_Spectrum,[length(ppm),4]);
Sim_Spectrum=Sim_Spectrum+squeeze(sum(spline_vault,3)); 
Fit_Spectrum2=reshape(Hercules_model(Fitted_weights,BasisSet),[length(ppm),4])+squeeze(sum(spline_vault,3));
%ii=1;
Fit_Spectrum2=reshape(Hercules_model(Fitted_weights(ii,:),BasisSet),[length(ppm),4])+squeeze(sum(spline_vault,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

figure(2)
Sim_Spectrum_display=Sim_Spectrum;
Fit_Spectrum_display=Fit_Spectrum2;
Sim_Spectrum_display(:,2:4)=Sim_Spectrum_display(:,2:4)*5;
Fit_Spectrum_display(:,2:4)=Fit_Spectrum_display(:,2:4)*5;
%Sim_Spectrum_display(Sim_Spectrum_display<-180)=-180;
%Fit_Spectrum_display(Fit_Spectrum_display<-180)=-180;
StackPlot(Sim_Spectrum_display,3000000)
hold on
StackPlot(Fit_Spectrum_display,3000000,'r')
spline_vault_display=squeeze(sum(spline_vault,3));
spline_vault_display(:,2:4)=spline_vault_display(:,2:4)*5;
StackPlot(spline_vault_display,3000000,'g')
hold off     
text(400,220,'* 0.2')

resnorm

% orig_width=-200;
% BasisSet=zeros(length(ppm),4,length(weights));
% %simpleCho
% BasisSet(:,1,1)=LorentzianModel([0.8 orig_width 3.2 0 0],ppm);
% %simpleCr
% BasisSet(:,1,2)=LorentzianModel([1 orig_width 3.0 0 0],ppm);
% %simpleNAA
% BasisSet(:,1,3)=LorentzianModel([1.2 orig_width 2.0 0 0],ppm);
% BasisSet(:,3,3)=LorentzianModel([-1.2 orig_width 2.0 0 0],ppm);
% %simpleLac
% BasisSet(:,2,4)=LorentzianModel([0.2 orig_width 1.33 0 0],ppm)+LorentzianModel([0.2 orig_width 1.28 0 0],ppm);
% %simpleGABA
% BasisSet(:,1,5)=LorentzianModel([0.2 orig_width 3.02 0 0],ppm);
% BasisSet(:,3,5)=LorentzianModel([0.1 orig_width 3.02+0.06 0 0],ppm)+LorentzianModel([0.1 orig_width 3.02-0.06 0 0],ppm)+LorentzianModel([0.05 -200 3.02 0 0],ppm);
% 
% 
% 
% BasisSet=reshape(BasisSet,[length(ppm)*4,length(weights)]);

%Model_Spectrum = sum(BasisSet.*repmat(weights,[size(BasisSet,1) 1]),2);
%plot(Model_Spectrum);