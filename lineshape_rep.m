%Gaussian
x_bit=-1:0.0001:1;


G_width=0.2;
sigma=G_width/2/sqrt(2*log(2));
Gauss3=exp(-(x_bit.^2)/2/sigma^2)/(sigma*sqrt(2*pi));
plot(x_bit,Gauss3)


%Lorentzian
GAMMA_cap=0.2;
Lorentz3=(((x_bit.^2)+(0.5*GAMMA_cap)^2).^-1)/pi*0.5*GAMMA_cap;

GAMMA_Voigt=G_width/sqrt(1+8*pi^2*log(2)*(GAMMA_cap/G_width)^2)




% Since the equation in the paper is complicated, approximate combined
% width to be:

ii=0;
for gg=0.01:0.01:0.2
    for ll=0.01:0.01:0.2
  ii=ii+1;  
    G_width=gg;
    GAMMA_cap=ll;

    param1=G_width+GAMMA_cap;
    param2=G_width/(G_width+GAMMA_cap);

sigma=(param2*param1)/2/sqrt(2*log(2));
Gauss3=exp(-(x_bit.^2)/2/sigma^2)/(sigma*sqrt(2*pi));

Lorentz3=(((x_bit.^2)+(0.5*(param1-(param2*param1)))^2).^-1)/pi*0.5*(param1-(param2*param1));

Voigt=conv(Lorentz3,Gauss3,'same');
%max(Voigt)
%plot(Voigt)
%pause(10)
Voigt=Voigt/max(real(Voigt));
%plot(Voigt)
%pause(10)
%max(Voigt)
peak_top=Voigt>0.5;
width(ii)=sum(peak_top);
approx(ii)=gg+ll;
%plot(Voigt)
%pause(10)
    end
end
plot(approx,width*0.0001,'o')

