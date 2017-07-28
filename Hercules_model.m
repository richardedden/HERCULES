function F = Hercules_model(weights,basis)
%% Function for basis set fitting
%weights
%%%%for now weights(1:5) are metabolite 'concentrations'.
%%%%for now weights(6:13) are baseline.
%weights(14) is an additional line broadening parameter
%need to add in an interpolative x-shift

size(basis,2)

%make baseline
linear=linspace(-0.5,0.5,size(basis,1)/4);
baseline = [linear*weights(size(basis,2)+5)+weights(size(basis,2)+1) linear*weights(size(basis,2)+6)+weights(size(basis,2)+2) linear*weights(size(basis,2)+7)+weights(size(basis,2)+3) linear*weights(size(basis,2)+8)+weights(size(basis,2)+4)].';

%make lineshape
param1=weights(size(basis,2)+9);
param2=weights(size(basis,2)+10);

ppm=0.5:-0.01:-0.5;


sigma=(param2*param1)/2/sqrt(2*log(2));
Gauss3=exp(-(ppm.^2)/2/sigma^2)/(sigma*sqrt(2*pi));

Lorentz3=(((ppm.^2)+(0.5*(param1-(param2*param1)))^2).^-1)/pi*0.5*(param1-(param2*param1));

lineshape=conv(Lorentz3,Gauss3,'same');

%plot(ppm,Gauss3)
%pause(3)

%plot(ppm,Lorentz3)
%pause(3)
%plot(ppm,lineshape)
%pause(10)
%stop here
%lineshape = LorentzianModel([1 width 0 0 0],0.5:-0.01:-0.5).'.* GaussModel([1 widthG 0 0 0],0.5:-0.01:-0.5).';
%lineshape =lineshape/sum(lineshape);
%size(lineshape)
%size(sum(basis.*repmat(weights(1:5),[size(basis,1) 1]),2))
F = conv(sum(basis.*repmat(weights(1:size(basis,2)),[size(basis,1) 1]),2),lineshape,'same')+baseline;

F = fshift(F,weights(size(basis,2)+11));
%F = sum(basis.*repmat(weights(1:5),[size(basis,1) 1]),2)+baseline;

end
