 function[H_LS]=LS_CE(Y,Xp,pilot_loc,Nfft,Nps,int_opt)
 %LSchannelestimationfunction
 % Inputs:
 %      Y           = Frequency-domainreceivedsignal
 %      Xp          = Pilotsignal
 %      pilot_loc   = Pilotlocation
 %      N           = FFTsize
 %      Nps         = Pilotspacing
 %      int_opt     = 'linear' or 'spline'
 % Output:
 %      H_LS =LSChannelestimate

 Np=Nfft/Nps; k=1:Np;

 LS_est(k)=Y(pilot_loc(k))./Xp(k); %LSchannelestimation

 if lower(int_opt(1))=='l'
    method = 'linear';
 else
    method = 'spline';
 end

 %Linear/Splineinterpolation
 H_LS=interpolate(LS_est,pilot_loc,Nfft,method);