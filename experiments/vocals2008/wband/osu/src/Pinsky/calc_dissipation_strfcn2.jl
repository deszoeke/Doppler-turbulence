"""
Calculate cloud top TKE dissipation from Pinsky-adjusted Doppler velocity
structure functions.
Simon de Szoeke 2013.04.26
2018.07

yet to be implemented:
Vertically localize the structure function.
Assume noise is the same at every height in the 10 minute ensemble.
"""

using MAT, NetCDF
using PyPlot


## set environment and constants
#cd("/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband")
#cd("/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband")
if ismatch(r"squall",gethostname())
    stem="~/Data/cruises/VOCALS_2008/RHB"
    way_Pinsky_retrieval="/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband/Pinsky/retrieval/"
elseif ismatch(r"fog",gethostname())
    stem="/Volumes/Pegasus2/cruises/VOCALS_2008/RHB"
    way_Pinsky_retrieval="/Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/"
end

year=2008
yyyy=@sprintf("%04d",year)
stday=318
enday=335
fs=3.5; #Hz
crosscomp=3/4
kolmogorov=0.54 # Matches atmospheric boundary layer estimates and Sreenivasan 1995
# kolmogorov=0.5 # probably only 1 digit of precision, Sreenivasan 1995
C1prime=4/3*kolmogorov # as in Pope eqn 6.243
factr=crosscomp/kolmogorov # 1/C1prime, used for my dissipation calculation

# universal constant for 2nd order structure function()
# longitudinal structure function constant C2 from Pope 6.2 [p.193] after Saddoughi and Veeravalli [1994]
C2ll=2.0 # 2.1, up to 2.2? Wiles et al. 2006
factrz=1/C2ll
factrx=3/4/C2ll
# Kolmogorov [viscous] scale for atmospheric boundary layer is()
# eta=(nu^3/epsilon); epsilon~3e-4 m^2/s^2/s, kinematic viscosity nu=2e-5 m^2/s
# --> eta= 2.3 mm
#     25 m = 1000 eta;  1000 m = 440 000 eta

dten=10*60; # seconds per window [10 min()]
diten=floor(Int,fs*dten); # samples per window

## structure function parameters
nmin=100; # min number of Doppler velocities to calc structure function in a 10 min iterval

"statistics functions"
nanmean(x) = mean(filter(!isnan,x)) # easy, not as fast as Images.meanfinite
nanmean(x,d) = mapslices(nanmean,x,d)
nanstd(x) = std(filter(!isnan,x))
nanstd(x,d) = mapslices(nanstd,x,d)
anom(x)=broadcast(-,x,nanmean(x))
anom(x,d)=broadcast(-,x,nanmean(x,d))
nancov(x,y) = nanmean(anom(x).*anom(y))
nancov(x,y,d) = nanmean(anom(x,d).*anom(y,d),d)
nancorr(x,y) = nancov(x,y)./(nanstd(x).*nanstd(y))
nancorr(x,y,d) = nancov(x,y,d)./(nanstd(x,d).*nanstd(y,d))


h=matread("$(way_Pinsky_retrieval)Pinsky_wRetrieval2008_318_23.mat")["h"]; # read height from a data file

"""
Calculate 2nd order longitudinal stucture function D of W, update count and D at bigind.
  calc_2dord_strfcn_10min!(count, D, bigind, W[ii10,:])
"""
function calc_2dord_strfcn_10min!(count, D, bigind, W)

nw=sum(isfinite.(W),1)
if maximum(nw)>100
    nmax,knmax=(maximum(nw),indmax(nw))
    # choose the highest vertical level ik that has enough data coverage
    fink=find(nw.>100)
    ik=fink[end]  # highest
    hk=fink[1]; # lowest
        
    # compute dissipation by longitudinal 2nd order stucture function
    # if there is 10 min cloud Doppler velocity data (A cloud top is identified.)
    if ~isempty(ik) && ik<100 # A cloud top is identified.
        ktop[bigind]=min(ik+1,100)
        kk=(max(hk-1,4):ktop[bigind])
        w10=W[:,kk]; # subset w10 and filter out extreme values
        w10[abs.(w10-nanmean(w10)).>4*nanstd(w10)]=NaN
        #iiwt=isfinite.(w10)
        nj=length(kk)
        kmax[bigind]=nj
        for j=1:nj # could be 1:nj-1 # index relative to 10min window cloud top k=ktop+(j-1)
            # k1=nj-j+1; # indexes w10[ ,k1]
            iij0=isfinite.(w10[:,nj-j+1]); # j=1 -> cloud top
            for l=1:nj-j # structure fcn displacement index
                iijl=isfinite.(w10[:,nj-j+1-l])
                iijj=iij0 .& iijl
                count[bigind,j,l]=sum(iijj)
                dw=w10[iijj, nj-j+1]-w10[iijj, nj-j+1-l]
                D[bigind,j,l]=mean(dw.*dw)
            end
        end
    end # there is 10 min cloud Doppler velocity data                
end # maximum(nw)>100
end

## loop length numbers
nhr=24*(enday-stday+1)
n10min=nhr*6
## allocate output data variables
jmax=83; lmax=jmax-1; # 72,71 # actually jmax could be jmax-1 too
ydayw10=NaN+zeros(n10min,1)
ktop=zeros(Int8,n10min,1)
kmax=zeros(Int8,n10min,1)
count=zeros(Int16,n10min,jmax,lmax)
D=NaN+zeros(n10min,jmax,lmax)

## loop over files by day, hour
dirlist=readdir(way_Pinsky_retrieval)
bigind=0
for iday=stday:enday
    @printf("\n%3d ",iday)
    for idhr=0:23
        @printf(" %02d ",idhr)

        retrfile=@sprintf("Pinsky_wRetrieval%04i_%03i_%02i.mat",year,iday,idhr)
        if ~isempty(retrfile) && contains(==,dirlist,retrfile)
            matf=matread("$(way_Pinsky_retrieval)$(retrfile)")
            W=matf["W"]
            nt=size(W,1)
            # calculate cloud-top dissipation from inertial range of Pinsky-adjusted vel spectra
            # loop through, interp and calc spectra for 10-min chunks
            for ind=1:round(Int,nt/diten) # index of the 10 minute interval in the hour
                bigind=24*6*(iday-stday)+6*idhr+ind
                i10=(ind-1)*diten+1; # first index in the 10 minute interval
                ydayw10[bigind,1]=iday+(idhr+(ind-1)/6)/24; # current time of this 10-minute window
                #Utop=double(interp1(yday10,U,ydayw10[bigind]+5/(24*60))); # interp ship-rel wind speed to the middle of the current time window
                # find the max coverage of nonnoise returns --> nmax
                ii10=i10:min(i10+diten-1,nt)
                calc_2dord_strfcn_10min!(count, D, bigind, W[ii10,:]) #updates count, D at bigind                
            end # 10 min
        end # there is an hourly file
    end # hr
end # day

# truncate data
ydayw10=ydayw10[1:bigind]
ktop=ktop[1:bigind]
kmax=kmax[1:bigind]
jmax=maximum(kmax); lmax=jmax-1
count=count[1:bigind,1:jmax,1:lmax]
D=D[1:bigind,1:jmax,1:lmax]

# save as a matlab file
matwrite("D_strfcn.mat", Dict( "ydayw10" => ydayw10,
                               "ktop" => ktop,
	                       "kmax" => kmax,
                     	       "jmax" => jmax,
	                       "count" => count,
                               "D" => D              )
        )

# calculate dissipation from the longtudinal structure function
# preliminaries
dh1=h[50]-h[49]
m=(0:jmax-1)'.+(1:lmax) # m=j-1+l; midpoint index
mid=m*0.5*dh1;                # midpoint displacement below cloud top [m]
r23=(dh1*(1:lmax)').^(2/3); # displacement[71], top[72] 
r23_2=repmat(r23',1,jmax) # 2D version of r^(2/3)
# solve D = N + A*r^(2/3) for noise N and coefficient A
# x=lsqlin[Cee,Dee,Aye,bee] minimizes the residual Cee*x-Dee subject to A*x<=b
# Aye=-eye(2); # limiting conditions for solution
# bee=zeros(2,1)
Aeq=[]; beq=[]; # no equality constraint
lb=[0; 0]; # noise and dissipation must be positive
ub=[Inf; Inf]
x0=[9e-2; 6e-3]; # not used for interior-point, but must supply [] argument
#options = optimoptions["lsqlin','Algorithm','interior-point','display()','none"]

fii=find(sum(sum(isfinite.(D),2),3).>100); # times that have reasonable ensembles of structure functions()

itime=1042; # test time
    nj=kmax[fii[itime]]
    D2=D[fii[itime],:,:]'
    count2=count[fii[itime],:,:]'
    ix=mapslices(indmax,(count2.>=200).*D2,1)
    mx=[ ((count2.>=200).*D2)[ix[i],i] for i in 1:size(D2)[2] ]
    D2filt=D2; # D2[displacement, top position]
    for i=1:size(D2,1)
        D2filt[(ix[i]+1):end,i]=NaN
    end
    
clf()
subplot(2,1,1)
plot(r23',D2,marker=".",linewidth=0)
plot(r23',D2filt)
#plot(r23_2[ix,:][1,:,1],mx,marker="^")
plot(r23[ix],mx,marker="^")
xlim(0, 50)
ylabel("structure function D = N + Ar\$^{2/3}\$")
subplot(2,1,2)
semilogy(r23',count2,".-")
xlim(0, 50)
xlabel("structure function displacement r\$^{2/3}\$ [m\$^{2/3}\$]")
ylabel("count points in structure function")

# require at least 200 points in the structure function
# cut D for displacements greater than global maximum(D)

noise=NaN+zeros(n10min,lmax)
A=NaN+zeros(n10min,lmax)
for itime=1:length(fii) # pick a time
    #nj=kmax[fii[itime]]
    D2=D[fii[itime],:,:]'
    count2=count[fii[itime],:,:]'
    ix=mapslices(indmax,(count2.>=200).*D2,1)
    #mx=[ ((count2.>=200).*D2)[ix[i],i] for i in 1:size(D2)[2] ]
    D2filt=D2; # D2[displacement, top position]
    for i=1:size(D2,1)
        D2filt[(ix[i]+1):end,i]=NaN
    end
    jjll=isfinite.(D2filt) .& (r23_2.>10.0); # displacement, position; exclude first r bc adjacent range gates are correlated
    #jjll=count2>10
    # good structure functions are found at D2[jjll]
    nr=sum(jjll,1)'
    njjll=sum(nr)
    # solve for a uniform noise and all A_k simultaneously
    if njjll>=5
        Dee=D2filt[jjll] # D2[jjll]
        # build the matrix Cee that can be solved for a single value of noise and different A for each
	# midpoint height
	iinr=nr.>0; # some columns excluded by conditions above, so need to remove from Cee
        nj=sum(iinr)+1
	nrend=cumsum(nr[iinr])
	ic=collect(1:njjll)
        jc=zeros(size(ic))
        jc[1:nrend[1]]=1
        for i=2:length(nrend)
            jc[(nrend[i-1]+1):nrend[i]]=i # put ith r23_2 column jc the ith
        end
        #Cee=sparse(ic,jc+1,r23_2[jjll], njjll,nj, 2*njjll) # for matlab to allocate
	Cee=sparse(ic,jc+1,r23_2[jjll])
        Cee[:,1]=1.0
        Aye=-eye(sum(iinr)+1)
        bee=zeros(sum(iinr)+1,1)
#        x=lsqlin[Cee,Dee,Aye,bee,Aeq,beq,lb,ub,[],options]; # x is() [noise A_k]
#         x=lsqlin[Cee,Dee,Aye,bee]; # x is() [noise A_k]
	x = (Cee'*Cee) \ (Cee'*Dee); # left divides
        # if noise<0, set it to zero and recompute the slope
        if x[1]<0
            x[1]=0
            meanr23=sum(Cee[:,2:end],1)./sum(Cee[:,2:end].>0,1)
            meanD  =nanmean(D2filt[r23'.>10.0,:],1)
            x[2:end]=meanD[isfinite.(meanD)]./meanr23[:]
        end
        noise[fii[itime]]=x[1]
        for pickm=1:nj-1 # indexes displacement
            A[fii[itime],pickm]=x[pickm+1]
        end
    end
end

PyPlot.plt[:hist](noise[isfinite.(noise)],linspace(0,.5,50))

clf()
plot(r23',D2,marker=".",linewidth=0)
plot(r23',D2filt)
plot(r23_2[ix,:][1,:,1],mx,marker="^",linewidth=0)
xlim(0, 50)
ylabel("structure function D = N + Ar\$^{2/3}\$")
xlabel("structure function displacement r\$^{2/3}\$ [m\$^{2/3}\$]")

epsz=(A*factrz).^1.5
semilogy(ydayw10,epsz,marker=".")

ax=pcolor(ydayw10,-mid[:,1],log10(epsz)')
caxis(-5, -2.5)
ax[:set_facecolor]((0.7, 0.7, 0.7))
colorbar()
title("dissipation from longitudinal structure function log\$_{10}\$[m\$^2\$ s\$^{-3}\$]")
xlabel("yearday 2008")
ylabel("midpoint displacement from cloud top [m]")
savefig("epsz_cloudtoprel.eps")

# reconstruct dissipation on absolute height grid
dh2=dh1/2
habs=flipdim(h[maximum(ktop)]-dh2*(1:size(epsilonz,2))', 1); # increasing
epsilonz=NaN(length(ydayw10),maximum(ktop)+size(mid,1)-1)
fl=find(isfinite(epsz)); # going in the direction that sparse is better than full(), but not quite there yet
for fli=1:length(fl) # very fast
    l=fl[fli]
    [i,em]=ind2sub(size(epsz),l)
    epsilonz[i,double(ktop[i])-em] = epsz[l];   
end

ax=pcolor(ydayw10,habs,log10(epsilonz)')
caxis((-5 -2.5))
ylim((400, 1300))
colorbar()
ax[:set_facecolor]((0.7, 0.7 0.7))

# compare with spectral dissipation calculation
# at every other vertical level of the strfcn calc.

spectr=matread("radarturb3.mat")
find(abs(h-habs[end-1])<10)
#habs[125]==h[76]
# koffs=(length(habs)-1) - find(abs(h-habs[end-1])<10); # offset at top
koffs=find(abs(h-habs[1])<5); # offset at base, 1-based, =14
# habs[1] == h[koffs]
plot(koffs-1+(1:length(habs)/2),habs[1:2:end],'.')
nk=fix(length(habs)/2); # 63
# sum( abs( h[koffs+(0:nk-1)] - habs[1:2:end] )>5 ) # 0
# h[koffs+(0:nk-1)]==habs[1:2:end]
ii=(any(isfinite(epsilonz),2))
x=spectr["epsilonk"][ii,koffs+(0:nk-1)]
y=epsilonz[ii,1:2:end]
loglog(x[:],y[:],marker=".",markersize=1)
plot([1e-6 1],[1e-6 1],color="k")
plot([1e-6 1],2*[1e-6 1],color="k",linestyle="--")
ylabel("structure function dissipation")
xlabel("spectral dissipation")
axis([1e-6 1e-2 1e-6 1e-2])
#[nanmedian[x[:]] nanmean[x[:]]; nanmedian[y[:]] nanmean[y[:]]]
# not well correlated, but distributions agree pretty well in the
# mean(),median()
# strfcn dissipation is not well localized, and I haven't conditioned on
# regions of nearly constant dissipation.

# ...or could do by interpolation
# xi=find(any(isfinite(epsilonz),2))
# interp2(xi,habs


plot([0; 120],x[1]+x[2]*[0; 120],"r-")
hold on
plot(r23_2[jjllmm],Dee,'.')

clf()
pcolor(dh1*(1:60),dh1*(0:50-1),squeeze(D[fii[itime],1:50,1:60])); shading flat
cx=caxis()
contour(dh1*(1:60),dh1*(0:50-1),mid[1:50,1:60],0:100:1000,'k')
set(gca(),"tickdir','out','color',0.7+[0 0 0],'fontsize",14)
xlabel("displacement [m]")
ylabel("high point distance below cloud top [m]")
caxis(cx)
colorbar("East")

plot(ydayw10[ktop>0],h[ktop[ktop>0]],marker=".")
