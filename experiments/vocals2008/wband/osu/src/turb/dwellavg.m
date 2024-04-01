function x=dwellavg(tt,xx,t,dwell)
% x=dwellavg(tt,xx,t,dwell);
% Average xx(tt) over new time intervals t. tt and t are the beginning
% of time intervals. The average conserves integrals. Specify dwell, the
% time over which to average the signal.

% SPdeS 8 June 2007

ni=length(t);
x=NaN+zeros(ni,1);

nii=length(tt);

sti=max(1,find(t>tt(1),1,'first')-1);
%t(ni+1)=2*t(ni)-t(ni-1); don't make up end of interval at last point

for i=sti:ni-1
    ii=max(1,find(tt<t(i),1,'last')); % subscript of the time interval tt straddling the beginning of interval t
    if isempty(ii) % no data xx before t(i)
        ii=1;
    end
    s=0; dts=0;
    while (tt(ii)<min(t(i+1),t(i)+dwell) && ii<nii )
        if isfinite(xx(ii))
            dt=min([tt(ii+1),t(i+1),t(i)+dwell])-max(tt(ii),t(i));
            s=s+xx(ii)*dt;
            dts=dts+dt;
        end
        ii=ii+1;
    end
    if dts==0
        x(i)=NaN;
    else
        x(i)=s/dts;
    end
end
i=ni; % t(i+1) does not exist, but dwell does
    ii=find(tt<t(i),1,'last');
    s=0; dts=0;
    while (tt(ii)<t(i)+dwell && ii<nii && isfinite(xx(ii)))
        dt=min(tt(ii+1),t(i)+dwell)-max(tt(ii),t(i));
        s=s+xx(ii)*dt;
        dts=dts+dt;
        ii=ii+1;
    end
    if dts==0
        x(i)=NaN;
    else
        x(i)=s/dts;
    end
%
return
% The last xx(nii) is not averaged though tt(nii) is used as an end for the
% previous interval [tt(nii-1) tt(nii)].

if 0 % plot for visual inspection
    figure
    ttp=[tt(1:end-1) tt(2:end)]';
    tp=[t(1:end-1) t(2:end)]';
    xxp=[xx(1:end-1) xx(1:end-1)]';
    xp=[x(1:end-1) x(1:end-1)]';
    ttp=ttp(:); tp=tp(:); xxp=xxp(:); xp=xp(:);
    plot(ttp(:),xxp(:),'rx-',tp(:),xp(:),'b.-')
end