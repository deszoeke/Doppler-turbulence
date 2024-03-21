function [Smoothed, Frq]=specsmooth(SS,DD)

%*****************************************************
% Smooths spectra.
% [Smoothed, Frq]=specsmooth(SS,DD)
%
% Inputs:  SS is the spectrum to be smoothed
%          DD is the sampling rate
% Outputs: Smoothed is the smoothed spectrum
%          Frq is the frequency of smoothed points
%******************************************************

jj=1;

%Larger number in Rdi results in less points in the smoothed spectrum

%Rdi=.2;      % 31 points output from 4097 input
%Rdi=.175;    % 35 points output from 4097 input
%Rdi=.15;     % 41 points output from 4097 input
%Rdi=.125;    % 47 points output from 4097 input
%Rdi=.123;    % 49 points output from 4097 input
Rdi=.118;    % 50 points output from 4097 input


Mj=1;
Di=0;
Slength=length(SS);
ii=0;
Dfx=DD/(Slength*2);        %DD==sampling rate
while ii<Slength;
   if ii==Slength;
      break;
   end;
   Di=round(Rdi*(Mj+Di));
   Smoothed(jj)=0;
   for ii=Mj:Mj+Di;
      if ii==Slength;
         break;
      end;
      Smoothed(jj)=Smoothed(jj)+SS(ii+1);
   end;
   Smoothed(jj)=Smoothed(jj)./(Di+1);
   Frq(jj)=(Mj+Di/2)*Dfx;
   Mj=Mj+Di+1;
   jj=jj+1;
end;
Smoothed=Smoothed(1:length(Frq)-1);
Frq=Frq(1:length(Frq)-1);
