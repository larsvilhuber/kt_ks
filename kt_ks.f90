!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kt_ks.f90
!
! Fortran code for the KS solution of the Khan and Thomas (2008) 
! model, in the baseline version.
!
! 'Alternative Methods for Solving Heterogeneous Firm Models'
! Stephen Terry (2015)
!
! This Version : 12/18/15
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module modparams
implicit none
    
!This module contains fixed model parameter values, available
!to the rest of the code. It also declares constants which are made available
!to subprograms when set elsewhere. It also declares some variables
!threadprivate, which makes them unique to OpenMP parallel threads.

integer, parameter :: knum = 10 !grid size for idio capital
integer, parameter :: znum = 5  !grid size for idio productivity
integer, parameter :: anum = 5  !grid size for agg productivity
integer, parameter :: kbarnum = 10 !grid size for agg capital
integer, parameter :: kdensenum = 50 !grid size for idio capital on dist
integer, parameter :: numper = 2500 !number of periods in unconditional simulation
integer, parameter :: numdiscard = 500 !number of periods discarded from unconditional simulation 
integer, parameter :: seedint=2503 !random number seed
integer, parameter :: ainit=3 !initial gridpoint for agg prod in simulations
integer, parameter :: numsimIRF=2000 !number of IRF economies to simulate
integer, parameter :: numperIRF=50 !number of periods per IRF economy
integer, parameter :: shockperIRF=25 !period at which to shock each IRF economy
integer, parameter :: shockanumIRF = anum !grid point that you impose in IRF period
integer, parameter :: maxvfit = 1000 !max number of vf iters
integer, parameter :: maxaccelit = 50 !number of howard accelerations
integer, parameter :: maxpit=100 !max number of price iterations when clearing market
integer, parameter :: maxfcstit=30 !max number of GE forecast rule iterations
integer, parameter :: maxergit = 5000 !max number of iters to compute ergodic dist of agg prod
integer, parameter :: pcutoff = 17 !iterations at which to restart price clearing algorithms
integer, parameter :: nummicro = 7 !number of micro moments to compute
integer, parameter :: doVFI = 1 !complete VF iteration or read from files?
integer, parameter :: doIRF = 1 !do IRF simulation?

double precision, parameter :: alpha = 0.256 !capital elasticity
double precision, parameter :: nu = 0.640 !labor elasticity
double precision, parameter :: phi = 2.4 !disutility of labor
double precision, parameter :: xibar = 0.0083 !upper bound of capital AC dist
double precision, parameter :: beta = 0.977 !discount rate
double precision, parameter :: delta = 0.069 !capital depreciation rate
double precision, parameter :: kmin= 0.1 !min of idio capital grid
double precision, parameter :: kmax = 8.0 !max of idio capital grid
double precision, parameter :: kbarmin =1.25 !min of agg capital grid
double precision, parameter :: kbarmax = 2.0 !max of agg capital grid
double precision, parameter :: nstdevz = 2.0 !number of st dev's to span for idio prod discretization
double precision, parameter :: rhoz = 0.859  !persistence of idio prod shock
double precision, parameter :: sigmaz = 0.022 !st dev of shock to idio prod
double precision, parameter :: nstdeva = 2.0 !number of st dev's to span for agg prod discretization
double precision, parameter :: rhoa = 0.859  !persistence of agg prod shock
double precision, parameter :: sigmaa=0.014 !st dev of shock to agg prod
double precision, parameter :: shocksizeIRF=sigmaa !shock size in IRF simulations
double precision, parameter :: ergdistatol = 1e-5 !tolerance on computing ergodic distriution
double precision, parameter :: vferrortol=1e-4 !tolerance of VF iteration
double precision, parameter :: kprimeerrortol = 1e-4 !tolerance on policy convergence
double precision, parameter :: xistarerrortol = 1e-4 !tolerance on adjustment threshold convergence
double precision, parameter :: perrortol = 1e-4 !tolerance on mkt clearing/price
double precision, parameter :: brenttol = 1e-6 !tolerance on Brent optimization of adjusting capital
double precision, parameter :: fcsterrortol = 1e-3 !tolerance on fcst rule convergence
double precision, parameter :: fcstgain=0.5 !dampening parameter for fcst rule update
double precision, parameter :: plb=2.0 !lb for price clearing
double precision, parameter :: pub=2.4 !ub for price clearing
double precision, parameter :: pwindow=0.01 !initial window around pfcst for price clearing

!some stuff to be available globally
double precision :: k0(knum),z0(znum),a0(anum),kbar0(kbarnum),kdense0(kdensenum),pr_mat_z(znum,znum),pr_mat_a(anum,anum),&
					pr_mat(anum*znum,anum*znum),kfcstmat(anum,2),pfcstmat(anum,2),kbarfcstwgts(anum,kbarnum),&
					kbarfcstvals(anum,kbarnum),pfcstvals(anum,kbarnum),V(znum,anum,knum,kbarnum),&
					Vold(znum,anum,knum,kbarnum),Vna(znum,anum,knum,kbarnum),Va(znum,anum,knum,kbarnum),&
					kprime(znum,anum,knum,kbarnum),xistar(znum,anum,knum,kbarnum),V2old(znum,anum,knum,kbarnum),&
					kprime2(znum,anum,knum,kbarnum),xistar2(znum,anum,knum,kbarnum),kprimeold(znum,anum,knum,kbarnum),&
					xistarold(znum,anum,knum,kbarnum),distkz(znum,kdensenum,numper),kprimep(znum,kdensenum),&
					xistarp(znum,kdensenum),padjustp(znum,kdensenum),perrorsim(numper),kbarfcstsim(numper),&
					kbarsim(numper),psim(numper),asimshock(numper),pfcstsim(numper),ysim(numper),isim(numper),&
					asimshockIRF(numperIRF,numsimIRF),arriveshockIRF(numsimIRF),kfcststore(anum,2,maxfcstit),&
					pfcststore(anum,2,maxfcstit),kfcstmatnew(anum,2),pfcstmatnew(anum,2),ergz0(znum),ergz0old(znum),&
					pr_mat_z_dum(znum,znum),pstore(anum,kbarnum),Kprimestore(anum,kbarnum),perrorstore(anum,kbarnum),&
					Kbarnoagg(anum),pnoagg(anum),Kbaraggrule(numper),paggrule(numper),ergdista(anum),&
					ergdistaold(anum),kbarsimIRF(numperIRF,numsimIRF,2),psimIRF(numperIRF,numsimIRF,2),&
					perrorsimIRF(numperIRF,numsimIRF,2),kbarfcstsimIRF(numperIRF,numsimIRF,2),&
					pfcstsimIRF(numperIRF,numsimIRF,2),ysimIRF(numperIRF,numsimIRF,2),isimIRF(numperIRF,numsimIRF,2),&
					distIRF(znum,kdensenum,numperIRF,numsimIRF,2),NsimIRF(numperIRF,numsimIRF,2),Nsim(numper),&
					MICROsim(nummicro,numper),kfcsterror,pfcsterror

integer :: asimpos(numper),asimposIRF(numperIRF,numsimIRF,2),kbarfcstinds(anum,kbarnum)

integer, allocatable :: seedarray(:)

!stuff to be available to subfunctions and which may need to be varied across OpenMP threads
double precision :: pval,wval,aval,weight
integer :: kbarfcstind,RHSact,RHSzct
!$omp threadprivate(pval,wval,aval,kbarfcstind,weight,RHSact,RHSzct)
    
end module modparams

program kt_ks
use omp_lib
use base_lib
use modparams
implicit none

integer :: zct,act,debugind,kct,zprimect,aprimect,kbarct,statect,stateprimect,&
    vfct,ct,zvalpos,avalpos,kvalpos,pfcstvalpos,kbarfcstvalpos,zctpos,&
    actpos,t,piter,kprimeind,kprimeindnoadj,fcstiter,perct,&
    adjustbias,accelit,pinitflag,numstates,seeddim,shockct,simct,ind,&
    kbarfcstindstore,momct
       
double precision :: start,finish,vferror,zval,kval,kbarval,kbarfcstval,kprimeval,&
    pfcstval,Vnaval,Vaval,wfcstval,Vnextval,xival,Vval,kprimeerror,xistarerror,&
    Vnextval1,Vnextval2,Yactualp,Iactualp,Cactualp,yval,perror,Kprimeactualp,&
    kprimevalnoadj,weightnoadj,xmean,ymean,x2mean,xymean,xval,&
    ergerror,kfcstbias,pfcstbias,p1val,p2val,p1error,p2error,shockprobIRF,&
    shockprobdenom,Nactualp,wgt,pvalstore,wvalstore,avalstore,weightstore,&
    fa,fb,fc,pvala,pvalb,pvalc,ival,padjust

start = omp_get_wtime()

open(13,file="kt_ks.txt")

!$omp parallel
write(*,*) "parallel hello to you."
!$omp end parallel

!!!!!!! INSERT PARAMETERS
numstates = znum*anum*knum*kbarnum

!write constants
open(8,file="constants.txt")
write(8,*) xibar,delta,numper,numdiscard,numperIRF,numsimIRF,shockperIRF,doIRF
close(8)


!capital grids
call linspace(k0,log(kmin),log(kmax),knum); k0 = exp(k0);
call linspace(kbar0,kbarmin,kbarmax,kbarnum);
call linspace(kdense0,log(kmin),log(kmax),kdensenum); kdense0=exp(kdense0);

!fcst rule initialization
kfcstmat(1,:) = (/   5.5703123796310788E-002, 0.81914773272391017 /)
kfcstmat(2,:) = (/   6.9665036089883670E-002, 0.81538518691927342 /)
kfcstmat(3,:) = (/   8.1809515956025064E-002, 0.81478097273136640 /)
kfcstmat(4,:) = (/   9.4412606856224934E-002, 0.81379804303600745 /)
kfcstmat(5,:) = (/  0.10996621006719427     , 0.80876644960499511 /)

kfcststore(:,:,:)=0.0; kfcststore(:,:,1)=kfcstmat;

pfcstmat(1,:) = (/  1.0156048668373361  ,   -0.40906585902130099 /)
pfcstmat(2,:) = (/ 0.99893966987725991  ,   -0.40622933253252513 /)
pfcstmat(3,:) = (/ 0.98065875147066461  ,   -0.40188307388436484 /)
pfcstmat(4,:) = (/ 0.96172566916304070  ,   -0.39616735766609229 /)
pfcstmat(5,:) = (/ 0.94552394727216549  ,   -0.39495042951610010 /)

pfcststore(:,:,:)=0.0; pfcststore(:,:,1)=pfcstmat;

!process coefficients into arrays used in solution below
call fcst_process()

!discretize exogenous processes A, z, and simulate A
call discretize_simulate()

!intialize the value function
Vold(:,:,:,:) = 0.0
V2old(:,:,:,:) = 0.0
kprimeold(:,:,:,:) = 0.0
xistarold(:,:,:,:) = 0.5

if (doVFI==0) then ;!read in the Vold and V2old matrices from data files
open(8,file="kprimeold.txt")
do zct=1,znum
do act=1,anum
do kct=1,knum
do kbarct=1,knum
read(8,*) kprimeold(zct,act,kct,kbarct)
end do
end do
end do
end do
close(8)
end if

!initialize kprime for Howard acceleration
if (doVFI==1) then
!$omp parallel private(zct,act,kct,kbarct)
!$omp do collapse(4)
do zct=1,znum
do act=1,anum
do kct=1,knum
do kbarct=1,kbarnum
 	   
    kprimeold(zct,act,kct,kbarct) = k0(kct)

end do !kbarct
end do !kct
end do !act
end do !zct
!$omp end do nowait
!$omp end parallel
end if


!!!INSERT FCST RULE LOOPING APPARATUS HERE
do fcstiter=1,maxfcstit


!!!!!!!!!!!!!!!!! GIVEN A FCST RULE, PERFORM VFI

do vfct=1,maxvfit
    
    !note that with howard iteration, it is kprime that's important, not Vold
    Vold(:,:,:,:) = 0.0; V2old(:,:,:,:) = 0.0;
    
    !here is howard acceleration step
    do accelit=1,maxaccelit
    
        !$omp parallel private(zct,act,kct,kbarct,zval,kval,kbarval,&
        !$omp& kbarfcstval,Vnaval,kprimeval,statect,&
        !$omp& zprimect,aprimect,stateprimect,Vnextval,Vaval,&
        !$omp& xival,Vval,Vnextval1,Vnextval2)
        !$omp do collapse(4)
        do zct=1,znum
        do act=1,anum
        do kct=1,knum
        do kbarct=1,kbarnum
            
            !determine states
            zval = z0(zct); aval = a0(act); kval = k0(kct); kbarval = kbar0(kbarct);
            
            !determine fcsts
            kbarfcstval = kbarfcstvals(act,kbarct)
            pval = pfcstvals(act,kbarct)
            wval = phi/pval;
            
            kbarfcstind = kbarfcstinds(act,kbarct)
            weight = kbarfcstwgts(act,kbarct)
            
            !! determine Vna(z,a,k,K)
            Vnaval = pval * freduced(zval,kval)
            kprimeval = max((1.0-delta)*kval,kmin)
            statect=(act-1)*znum+zct
            do zprimect=1,znum
            do aprimect=1,anum
                stateprimect = (aprimect-1)*znum+zprimect
                
                !now, need to evaluate Vold at (zprimect,aprimect,kprimeval,kbarfcstval)
                call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
                knum,kprimeval,Vnextval1)
                call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
                knum,kprimeval,Vnextval2)
                
                Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
                
                
                !now that next period's value function is evaluated, add to expectation
                Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
            end do
            end do
          
            
            
            
            !! determine Va(z,a,k,K)
            kprimeval = kprimeold(zct,act,kct,kbarct)
            Vaval = pval *  freduced(zval,kval)
            Vaval = Vaval - pval*(kprimeval - (1.0-delta)*kval)
            
            statect=(act-1)*znum+zct
            do zprimect=1,znum
            do aprimect=1,anum
                stateprimect = (aprimect-1)*znum+zprimect
                
                call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
                knum,kprimeval,Vnextval1)
                call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
                knum,kprimeval,Vnextval2)
                
                Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
                
                
                !now that next period's value function is evaluated, add to expectation
                Vaval = Vaval + beta * pr_mat(statect,stateprimect) * Vnextval
            end do
            end do
            
            xival = (Vaval - Vnaval)/phi
            Vval = -1.0 * phi * expecxi(xival)&
            + cdfxi(xival) * Vaval + (1.0 - cdfxi(xival)) * Vnaval;
            
            V(zct,act,kct,kbarct) = Vval
            
            !write(*,*) Vnaval,Vaval,Vval
            
        end do
        end do
        end do
        end do
        !$omp end do nowait
        !$omp end parallel
        
        Vold = V

        !given Vold from last iteration, spline the function along k dimension
        do zct=1,znum
        do act=1,anum
        do kbarct=1,kbarnum
            call spline(k0,Vold(zct,act,:,kbarct),knum,dble(1.0e30),dble(1.0e30),V2old(zct,act,:,kbarct))
        end do !kbarct
        end do !act
        end do !zct
    
        
    end do !accelct
    
    
    !$omp parallel private(zct,act,kct,kbarct,zval,kval,kbarval,&
    !$omp& kbarfcstval,Vnaval,kprimeval,statect,&
    !$omp& zprimect,aprimect,stateprimect,Vnextval,Vaval,&
    !$omp& xival,Vval,Vnextval1,Vnextval2)
    !$omp do collapse(4)
    do zct=1,znum
	do act=1,anum
	do kct=1,knum
	do kbarct=1,kbarnum
       
        !determine states
        zval = z0(zct); aval = a0(act); kval = k0(kct); kbarval = kbar0(kbarct);
        
        !determine fcsts
        kbarfcstval = kbarfcstvals(act,kbarct)
        pval = pfcstvals(act,kbarct)
        wval = phi/pval;
        
        kbarfcstind = kbarfcstinds(act,kbarct)
        weight = kbarfcstwgts(act,kbarct)
        
        !! determine Vna(z,a,k,K)
        Vnaval = pval * freduced(zval,kval)
        kprimeval = max((1.0-delta)*kval,kmin)
        statect=(act-1)*znum+zct
        do zprimect=1,znum
        do aprimect=1,anum
            stateprimect = (aprimect-1)*znum+zprimect
            
            !now, need to evaluate Vold at (zprimect,aprimect,kprimeval,kbarfcstval)
            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
            knum,kprimeval,Vnextval1)
            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
            knum,kprimeval,Vnextval2)
            
            Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
            
            
            !now that next period's value function is evaluated, add to expectation
            Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
        end do
        end do
        
        Vna(zct,act,kct,kbarct) = Vnaval
        !end of block computing Vnaval
      
        !now, this block determines Va(z,a,k,K)
        
        !first, you need to optimize the function wrt kprimeval
        RHSzct =zct
        RHSact = act
        Vaval = brent(k0(1),k0(knum/2),k0(knum),fnRHS,brenttol,kprimeval) !note this sets correct kprimeval, not correct vaval
        Vaval = -1.0 * Vaval
      
        !at this point, Vaval = -p k' + beta E V'
        !now, construct current period return remaining portions
        Vaval = Vaval + pval * (freduced(zval,kval)+(1.0-delta)*kval)
        
        Va(zct,act,kct,kbarct) = Vaval
        kprime(zct,act,kct,kbarct) = kprimeval
        
        xival = (Vaval - Vnaval)/phi
        xistar(zct,act,kct,kbarct) = xival
        !!end of block computuing Vaval
        
        !now, process the info to compute the value function
        Vval = -1.0 * phi* expecxi(xival)&
            + cdfxi(xival) * Vaval + (1.0 - cdfxi(xival)) * Vnaval;
  
        V(zct,act,kct,kbarct) = Vval
        

    end do !kbarct
    end do !kct
    end do !act
    end do !zct
    !$omp end do nowait
    !$omp end parallel


    vferror = maxval(abs( (log(V)-log(Vold) )) )    
    kprimeerror = maxval(abs(log(kprime)-log(kprimeold)))
    xistarerror = maxval(abs(log(xistar)-log(xistarold)))
    if (mod(vfct,1)==0) then
    write(13,*) "VF iter = ",vfct,"VF error = ",vferror
    write(13,*) "VF iter = ",vfct,"Kprime error = ",kprimeerror
    write(13,*) "VF iter = ",vfct,"Xistar error = ",xistarerror
    write(13,*) " "
    
    write(*,*) "VF iter = ",vfct,"VF error = ",vferror
    write(*,*) "VF iter = ",vfct,"Kprime error = ",kprimeerror
    write(*,*) "VF iter = ",vfct,"Xistar error = ",xistarerror
    write(*,*) " "
    
    
    
    end if
    if (doVFI==0.and.fcstiter==1) exit; !in this case, only need one run through to initialize the V2 and param matrices
    if (kprimeerror<kprimeerrortol .and. xistarerror<xistarerrortol) exit
    Vold=V
    kprimeold = kprime
    xistarold = xistar
end do !vfct

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished VFI in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished VFI in ",finish-start," seconds."


if (minval(kprime)<1.05*k0(1)) then
    write(13,*) "Lowest idio capital close to bottom of grid."
    write(*,*) "Lowest idio capital close to bottom of grid."
end if
if (maxval(kprime)>0.95*k0(knum)) then
    write(13,*) "Highest idio capital close to top of grid."
    write(*,*) "Highest idio capital close to top of grid."
end if


open(8,file="k0.txt")
do kct=1,knum; write(8,*) k0(kct); end do
close(8)

open(8,file="kbar0.txt")
do kbarct=1,kbarnum; write(8,*) kbar0(kbarct); end do
close(8)

open(8,file="V.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) V(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="Vna.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) Vna(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="Va.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) Va(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="xistar.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) xistar(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="kprime.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) kprime(zct,act,kct,kbarct)
end do; end do;end do; end do
close(8)

open(8,file="padjust.txt")
do zct=1,znum; do act=1,anum; do kct=1,knum; do kbarct=1,kbarnum
write(8,*) cdfxi(xistar(zct,act,kct,kbarct))
end do; end do;end do; end do
close(8)

!!!!END OF VF BLOCK

!NOW THAT THE VFI IS DONE, SIMULATE THE MODEL

distkz(:,:,:) = 0.0

!initialize
kbarsim(1) = 1.55
distkz(:,:,1) = 1.0
distkz(:,:,1) = distkz(:,:,1) / sum(distkz(:,:,1))
kbaraggrule(:) = 0.0
paggrule(:) = 0.0

!MICROsim(nummicro,numper)
MICROsim(:,:) = 0.0
! 1 = i/k
! 2 = stdev(i/k)
! 3 = P(inaction)
! 4 = P(i/k>=0.2)
! 5 = P(i/k<=-0.2)
! 6 = P(i/k > 0)
! 7 = P(i/k < 0)

do t=1,numper-1
!do t=1,1

    kprimep(:,:) = 0.0
    xistarp(:,:) = 0.0
    padjustp(:,:) = 0.0
    
    act=asimpos(t)
    aval = a0(act); avalstore = aval
    kbarval = kbarsim(t)
    kbarfcstval = kfcstmat(act,1) + kfcstmat(act,2)*log(kbarval); kbarfcstval = exp(kbarfcstval);
    pfcstval = pfcstmat(act,1) + pfcstmat(act,2)*log(kbarval); pfcstval = exp(pfcstval);


    kbarfcstind = kbarnum/2
    call hunt(kbar0,kbarnum,kbarfcstval,kbarfcstind)
    
    if (kbarfcstind<=0) then
        weight=0.0
        kbarfcstind=1
    else if (kbarfcstind>=1.and.kbarfcstind<=(kbarnum-1)) then
        weight = (kbarfcstval-kbar0(kbarfcstind))/(kbar0(kbarfcstind+1)-kbar0(kbarfcstind))
    else if (kbarfcstind>=kbarnum) then
        weight = 1.0
        kbarfcstind=kbarnum-1
    end if
    
    kbarfcstindstore = kbarfcstind
    weightstore = weight

    
    !note that we will in future iterate over p
    pvala = pfcstval - pwindow
    pvalb = pfcstval + pwindow
    pvalc = pvala + dble(0.67) * (pvalb-pvala) 
    
    do piter=1,maxpit
    
    !brent method for price optimization
    
    
    if (piter==1) pval = pvala
    if (piter==2) pval = pvalb
    if (piter==3) pval = pvalc
    
    if (piter==pcutoff) then
        !pvala = plb
        pvala = pfcstval - 4.0*pwindow
        pval = pvala
    else if (piter==pcutoff+1) then
        !pvalb = pub
        pvalb = pfcstval + 4.0*pwindow
        pval = pvalb
    else if (piter==pcutoff+2) then
        pvalc = pvala + dble(0.67) * (pvalb-pvala) 
            pval = pvalc
        end if
    
    
    if ((piter>3.and.piter<pcutoff).or.(piter>pcutoff+2)) then
            
                        !first, try inverse quadratic interpolation of the excess demand function
            pval = ( pvala * fb * fc ) / ( (fa - fb) * (fa - fc) ) &
                + ( pvalb * fa * fc ) / ( (fb - fa) * (fb - fc ) ) &
                + ( pvalc * fa * fb ) / ( (fc - fa) * (fc - fb ) )
       
            !if it lies within bounds, and isn't too close to the bounds, then done
            
            !o/w, take bisection step
            if ((minval( (/ abs(pvala - pval), abs(pvalb-pval) /) )<&
                    abs( (pvalb-pvala)/dble(9.0) ) ).or.(pval<pvala).or.(pval>pvalb))   then
                pval = (pvala + pvalb) / dble(2.0)
            
                    end if
        end if
    
    wval = phi/pval;

    pvalstore = pval
    wvalstore = wval

    !$omp parallel private(zct,kct,zval,kval,Vnaval,kprimeval,statect,stateprimect,&
    !$omp& aprimect,zprimect,Vnextval1,Vnextval2,Vnextval,Vaval,xival,Vval)
    !$omp do collapse(2)
    do zct=1,znum
	do kct=1,kdensenum
            
        aval = avalstore
        pval = pvalstore
        wval = wvalstore
        kbarfcstind = kbarfcstindstore
        weight = weightstore
            
        if (distkz(zct,kct,t)>0.0) then
        
        zval = z0(zct)
        kval = kdense0(kct)
        
        !block to determine Vnaval
        Vnaval = pval * freduced(zval,kval)
        kprimeval = max(kmin,(1.0-delta)*kval)
        
        statect=(act-1)*znum+zct
        do aprimect=1,anum
        do zprimect=1,znum
            stateprimect=(aprimect-1)*znum+zprimect
            

            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
            knum,kprimeval,Vnextval1)
            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
            knum,kprimeval,Vnextval2)
            
            Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
            
            Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
        
        end do !zprimect
        end do !aprimect
        
        !block to construct Vaval
        RHSzct =zct
        RHSact = act
        Vaval = brent(k0(1),k0(knum/2),k0(knum),fnRHS,brenttol,kprimeval) !note this sets correct kprimeval, not correct vaval
        Vaval = -1.0 * Vaval
        !at this point, Vaval = -p k' + beta E V'
        !now, construct current period return remaining portions
        Vaval = Vaval + pval * (freduced(zval,kval)+(1.0-delta)*kval)
        
        !now process this info
        xival = (Vaval - Vnaval)/phi
        Vval = -1.0 * phi * expecxi(xival) + cdfxi(xival)*Vaval &
            + (1.0-cdfxi(xival))*Vnaval
            
        kprimep(zct,kct)=kprimeval
        xistarp(zct,kct)=xival
        padjustp(zct,kct) = cdfxi(xival)
        

    end if !distkz tolerance
    
	end do !kct
	end do !zct
    !$omp end do nowait
    !$omp end parallel

    Yactualp = 0.0
    Iactualp = 0.0
    Kprimeactualp = 0.0
    Cactualp = 0.0
    Nactualp = 0.0
    
    !!$omp parallel private(zct,kct,zval,kval,kprimeval,xival,&
    !!$omp& yval) reduction(+:Yactualp,Iactualp,Kprimeactualp,Nactualp)
    !!$omp do collapse(2)
    do zct=1,znum
    do kct=1,kdensenum
         if (distkz(zct,kct,t)>0.0) then
        zval = z0(zct); 
        kval = kdense0(kct)
        kprimeval=kprimep(zct,kct)
        xival=xistarp(zct,kct)
        
        
        yval = yreduced(zval,kval)
        
        Yactualp = Yactualp  + distkz(zct,kct,t) * yval
        Iactualp = Iactualp + distkz(zct,kct,t) * cdfxi(xival) &
            * (kprimeval - (1.0-delta)*kval)
    
        Kprimeactualp = Kprimeactualp + distkz(zct,kct,t) * cdfxi(xival) * kprimeval &
            +  distkz(zct,kct,t) * (1.0 - cdfxi(xival)) * (1.0-delta)*kval
            
        Nactualp = Nactualp +  distkz(zct,kct,t) * (nreduced(zval,kval) + expecxi(xival))              
            
        end if
    end do !kct
    end do  !zct
    !!$omp end do nowait
    !!$omp end parallel

    Cactualp = Yactualp - Iactualp

        !are you initializing the brent?
        if (piter==1) fa = (1.0/pval) - Cactualp
        if (piter==2) fb = (1.0/pval) - Cactualp
        if (piter==3) fc = (1.0/pval) - Cactualp
        
        if (piter==pcutoff) fa = (1.0/pval) - Cactualp
        if (piter==pcutoff+1) fb = (1.0/pval) - Cactualp
        if (piter==pcutoff+2) fc = (1.0/pval) - Cactualp
        

    perror = 1.0/pval - Cactualp
    
    !if not restarting or initializing
       if ((piter>3.and.piter<pcutoff).or.(piter>pcutoff+2)) then 
            if (perror<0) then
                pvalc = pvalb; fc = fb;
                pvalb = pval; fb = perror;
                !pval a doesn't change
            else if (perror>=0) then
                pvalc = pvala; fc = fa;
                pvala = pval; fa = perror;
                !pval b doesn't change
            end if
        end if
    
    if (abs(perror)<perrortol) exit

    end do !piter
    if (mod(t,100)==1) then
        write(13,"(A,I5,A,I5,A,F7.5,A,F7.5,A,F7.5)") "t= ",t,", pit = ",piter,", err = ",perror,", p = ",pval,", aval = ",aval
        write(*,"(A,I5,A,I5,A,F7.5,A,F7.5,A,F7.5)") "t= ",t,", pit = ",piter,", err = ",perror,", p = ",pval,", aval = ",aval
    end if 
    
    !now, the value of p has been determined, in a market clearing fashion
    kbarsim(t+1)=Kprimeactualp
    psim(t) = pval
    perrorsim(t) = perror
    kbarfcstsim(t+1) = kbarfcstval
    pfcstsim(t) = pfcstval
    ysim(t) = Yactualp
    isim(t) = Iactualp
    Nsim(t) = Nactualp
    
    !DH stat series calculations
    if (t==(numdiscard)) then
        
        kbaraggrule(t) = kbarsim(t)
        paggrule(t) = pval
        
        kbaraggrule(t+1)=kbarfcstval
    
    else if (t>numdiscard) then
        kbaraggrule(t+1)=kfcstmat(act,1)+kfcstmat(act,2)*log(kbaraggrule(t))
        kbaraggrule(t+1) = exp(kbaraggrule(t+1))
        
        paggrule(t)=pfcstmat(act,1)+pfcstmat(act,2)*log(kbaraggrule(t))
        paggrule(t) = exp(paggrule(t))
        
    end if
    
    !now that market clearing has occured, insert distribution into next period's
    do zct=1,znum
    do kct=1,kdensenum
        if (distkz(zct,kct,t)>0.0) then
        zval = z0(zct);
        kval = kdense0(kct);
        kprimeval = kprimep(zct,kct)
        kprimevalnoadj = max(kmin,(1.0-delta)*kval)
        xival = xistarp(zct,kct)
        
        !bracket kprimeval (which is kprime given that you adjust)
        kprimeind=kct
        call hunt(kdense0,kdensenum,kprimeval,kprimeind)
        
        if (kprimeind<=0) then
            kprimeind = 1; weight = 0.0
        else if (kprimeind>=1.and.kprimeind<=(kdensenum-1)) then
            weight = (kprimeval - kdense0(kprimeind))/(kdense0(kprimeind+1) - kdense0(kprimeind));
        else if (kprimeind>=kdensenum) then
            kprimeind = kdensenum-1; weight = 1.0
        end if
        
        !bracket kprimevalnoadj (which is kprime given that you don't adjust)
        kprimeindnoadj=kct
        call hunt(kdense0,kdensenum,kprimevalnoadj,kprimeindnoadj)
        
        if (kprimeindnoadj<=0) then
            kprimeindnoadj = 1; weightnoadj = 0.0
        else if (kprimeindnoadj>=1.and.kprimeindnoadj<=(kdensenum-1)) then
            weightnoadj = (kprimevalnoadj - kdense0(kprimeindnoadj))/(kdense0(kprimeindnoadj+1) - kdense0(kprimeindnoadj));
        else if (kprimeindnoadj>=kdensenum) then
            kprimeindnoadj = kdensenum-1; weightnoadj = 1.0
        end if
        
        
        do zprimect=1,znum
            
            !transfer the weight to kprimeind that does adjust
            distkz(zprimect,kprimeind,t+1) = distkz(zprimect,kprimeind,t+1) &
                + pr_mat_z(zct,zprimect) * padjustp(zct,kct) * (1.0 - weight) * distkz(zct,kct,t)
            
            !transfer the weight to kprimeind+1 that does adjust
            distkz(zprimect,kprimeind+1,t+1) = distkz(zprimect,kprimeind+1,t+1) &
                + pr_mat_z(zct,zprimect) * padjustp(zct,kct) * weight * distkz(zct,kct,t)
    
            !transfer the weight to kprimeindnoadj that doesn't adjust
            distkz(zprimect,kprimeindnoadj,t+1) = distkz(zprimect,kprimeindnoadj,t+1) &
                + pr_mat_z(zct,zprimect) * (1.0 - padjustp(zct,kct)) * (1.0 - weightnoadj) * distkz(zct,kct,t)
            
            
            !transfer the weight to kprimeindnoadj+1 that doesn't adjust
            distkz(zprimect,kprimeindnoadj+1,t+1) = distkz(zprimect,kprimeindnoadj+1,t+1) &
                + pr_mat_z(zct,zprimect) * (1.0 - padjustp(zct,kct)) * weightnoadj * distkz(zct,kct,t)
            
        end do !zprimect
        end if !distkz tolerance
    end do !kct
    end do !zct
    
    !now that other simulation apparatus is complete, input correct data into moment storage
    !MICROsim(nummicro,numper)
    ! 1 = i/k
    ! 2 = stdev(i/k)
    ! 3 = P(inaction)
    ! 4 = P(i/k>=0.2)
    ! 5 = P(i/k<=-0.2)
    ! 6 = P(i/k > 0)
    ! 7 = P(i/k < 0)
    
    do zct=1,znum
    do kct=1,kdensenum
        if (distkz(zct,kct,t)>0.0) then
            zval = z0(zct);
            kval = kdense0(kct);
            kprimeval = kprimep(zct,kct)
            ival = kprimeval - (1.0-delta)*kval !investment conditional upon investment
            kprimevalnoadj = max(kmin,(1.0-delta)*kval)
            xival = xistarp(zct,kct)
            padjust = cdfxi(xival)
            
            !investment rate
            MICROsim(1,t) = MICROsim(1,t) + distkz(zct,kct,t) * padjust * (ival / kval)
            
            !investment rate squared - for stdev construction
            MICROsim(2,t) = MICROsim(2,t) + distkz(zct,kct,t) * padjust * ( (ival / kval)**2.0 )
            
            !P(inaction)
            MICROsim(3,t) = MICROsim(3,t) + distkz(zct,kct,t) * (1.0 - padjust)
            
            !P(pos. spike)
            if ((ival/kval)>=0.2) then
                MICROsim(4,t) = MICROsim(4,t) + distkz(zct,kct,t) * padjust 
            end if
            
            !P(neg. spike)
            if ((ival/kval)<=-0.2) then
                MICROsim(5,t) = MICROsim(5,t) + distkz(zct,kct,t) * padjust 
            end if
            
            !P(pos. invest)
            if ((ival/kval)>0.0) then
                MICROsim(6,t) = MICROsim(6,t) + distkz(zct,kct,t) * padjust 
            end if
            
            !P(neg. invest)
            if ((ival/kval)<0.0) then
                MICROsim(7,t) = MICROsim(7,t) + distkz(zct,kct,t) * padjust 
            end if
            
        end if
    end do !kct
    end do !zct

    !now, convert squared investment moment to stdev
    MICROsim(2,t) = sqrt(MICROsim(2,t) - (MICROsim(1,t)**2.0))
    
end do !t

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished simulation in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished simulation in ",finish-start," seconds."

!!!NOW THAT SIMULATION IS DONE, DO THE FCST RULE UPDATE

!now, from the simulated values of p and K', update the fcst arrays
call update_fcst()

if (doVFI==1) then; !store the arrays for future use
open(8,file="kprimeold.txt")
do zct=1,znum
do act=1,anum
do kct=1,knum
do kbarct=1,knum
write(8,*) kprimeold(zct,act,kct,kbarct)
end do
end do
end do
end do
close(8)
end if

open(8,file="kdense0.txt")
do kct=1,kdensenum; write(8,*) kdense0(kct); end do
close(8)

open(8,file="psim.txt")
do t=1,numper; write(8,*) psim(t); end do !t
close(8)

open(8,file="kbarsim.txt")
do t=1,numper; write(8,*) kbarsim(t); end do !t
close(8)

open(8,file="perrorsim.txt")
do t=1,numper; write(8,*) perrorsim(t); end do !t
close(8)

open(8,file="kbarfcstsim.txt")
do t=1,numper; write(8,*) kbarfcstsim(t); end do !t
close(8)

open(8,file="pfcstsim.txt")
do t=1,numper; write(8,*) pfcstsim(t); end do !t
close(8)

open(8,file="distkzsim.txt");
do zct=1,znum; do kct=1,kdensenum; do t=1,numper
write(8,*) distkz(zct,kct,t)
end do; end do; end do
close(8)

open(8,file="ysim.txt")
do t=1,numper; write(8,*) ysim(t); end do !t
close(8)

open(8,file="isim.txt")
do t=1,numper; write(8,*) isim(t); end do !t
close(8)

open(8,file="Nsim.txt")
do t=1,numper; write(8,*) Nsim(t); end do !t
close(8)


open(8,file="paggrule.txt")
do t=1,numper; write(8,*) paggrule(t); end do !t
close(8)

open(8,file="kbaraggrule.txt")
do t=1,numper; write(8,*) kbaraggrule(t); end do !t
close(8)

open(8,file="kprimep.txt")
do zct=1,znum; do kct=1,kdensenum; write(8,*) kprimep(zct,kct); end do; end do 
close(8)

open(8,file="xistarp.txt")
do zct=1,znum; do kct=1,kdensenum; write(8,*) xistarp(zct,kct); end do; end do 
close(8)

open(8,file="padjustp.txt")
do zct=1,znum; do kct=1,kdensenum; write(8,*) padjustp(zct,kct); end do; end do 
close(8)


open(8,file="MICROsim.txt")
do t=1,numper
do momct=1,nummicro
write(8,*) MICROsim(momct,t)
end do !momct
end do !t
close(8)

open(8,file="kfcststore.txt")
do ct=1,maxfcstit; do act=1,anum; write(8,*) kfcststore(act,:,ct); end do; end do
close(8)

open(8,file="pfcststore.txt")
do ct=1,maxfcstit; do act=1,anum; write(8,*) pfcststore(act,:,ct); end do; end do
close(8)


write(*,*) " "
if (doVFI==0) exit
if (pfcsterror<fcsterrortol.and.kfcsterror<fcsterrortol) exit    

!if forecast rule not yet converged, update
kfcstmatnew = kfcstmat + fcstgain * (kfcstmatnew - kfcstmat)
kfcststore(:,:,fcstiter+1) = kfcstmatnew
kfcstmat = kfcstmatnew

pfcstmatnew = pfcstmat + fcstgain * (pfcstmatnew - pfcstmat)
pfcststore(:,:,fcstiter+1) = pfcstmatnew
pfcstmat = pfcstmatnew

!now, based on the new coefficients, update the arrays used in model solution
call fcst_process()

write(*,*) " "
end do !fcstiter

write(*,*) "Done with GE at ",omp_get_wtime()-start," seconds."
write(13,*) "Done with GE at ",omp_get_wtime()-start," seconds."


do act=1,anum
    write(*,*) "kfcstmat(act,:) = ",kfcstmat(act,:)
end do !act

do act=1,anum
    write(*,*) "pfcstmat(act,:) = ",pfcstmat(act,:)
end do !act

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!IRF BLOCK
if (doIRF==1) then
    
!initialize distribution and agg capital
distIRF(:,:,:,:,:) = 0.0
do shockct=1,2
do simct=1,numsimIRF
    kbarsimIRF(1,simct,shockct) = 1.55
    distIRF(:,:,1,simct,shockct) = 1.0
    distIRF(:,:,1,simct,shockct) = distIRF(:,:,1,simct,shockct) / sum(distIRF(:,:,1,simct,shockct))
end do !simct
end do !shockct

do shockct=1,2
do simct=1,numsimIRF
do t=1,numperIRF-1

    kprimep(:,:) = 0.0
    xistarp(:,:) = 0.0
    padjustp(:,:) = 0.0
    
    act=asimposIRF(t,simct,shockct)
    aval = a0(act); avalstore = aval
    kbarval = kbarsimIRF(t,simct,shockct)
    kbarfcstval = kfcstmat(act,1) + kfcstmat(act,2)*log(kbarval); kbarfcstval = exp(kbarfcstval);
    pfcstval = pfcstmat(act,1) + pfcstmat(act,2)*log(kbarval); pfcstval = exp(pfcstval);


    kbarfcstind = kbarnum/2
    call hunt(kbar0,kbarnum,kbarfcstval,kbarfcstind)
    
    if (kbarfcstind<=0) then
        weight=0.0
        kbarfcstind=1
    else if (kbarfcstind>=1.and.kbarfcstind<=(kbarnum-1)) then
        weight = (kbarfcstval-kbar0(kbarfcstind))/(kbar0(kbarfcstind+1)-kbar0(kbarfcstind))
    else if (kbarfcstind>=kbarnum) then
        weight = 1.0
        kbarfcstind=kbarnum-1
    end if
    
    kbarfcstindstore = kbarfcstind
    weightstore = weight

    
    !note that we will in future iterate over p
    pvala = pfcstval - pwindow
    pvalb = pfcstval + pwindow
    pvalc = pvala + dble(0.67) * (pvalb-pvala) 
    
    do piter=1,maxpit
    
    !brent method for price optimization
    
    
    if (piter==1) pval = pvala
    if (piter==2) pval = pvalb
    if (piter==3) pval = pvalc
    
    if (piter==pcutoff) then
        !pvala = plb
        pvala = pfcstval - 4.0*pwindow
        pval = pvala
    else if (piter==pcutoff+1) then
        !pvalb = pub
        pvalb = pfcstval + 4.0*pwindow
        pval = pvalb
    else if (piter==pcutoff+2) then
        pvalc = pvala + dble(0.67) * (pvalb-pvala) 
            pval = pvalc
        end if
    
    
    if ((piter>3.and.piter<pcutoff).or.(piter>pcutoff+2)) then
            
                        !first, try inverse quadratic interpolation of the excess demand function
            pval = ( pvala * fb * fc ) / ( (fa - fb) * (fa - fc) ) &
                + ( pvalb * fa * fc ) / ( (fb - fa) * (fb - fc ) ) &
                + ( pvalc * fa * fb ) / ( (fc - fa) * (fc - fb ) )
       
            !if it lies within bounds, and isn't too close to the bounds, then done
            
            !o/w, take bisection step
            if ((minval( (/ abs(pvala - pval), abs(pvalb-pval) /) )<&
                    abs( (pvalb-pvala)/dble(9.0) ) ).or.(pval<pvala).or.(pval>pvalb))   then
                pval = (pvala + pvalb) / dble(2.0)
            
                    end if
        end if
    
    wval = phi/pval;

    pvalstore = pval
    wvalstore = wval
    

    !$omp parallel private(zct,kct,zval,kval,Vnaval,kprimeval,statect,stateprimect,&
    !$omp& aprimect,zprimect,Vnextval1,Vnextval2,Vnextval,Vaval,xival,Vval)
    !$omp do collapse(2)
    do zct=1,znum
	do kct=1,kdensenum
            
        aval = avalstore
        pval = pvalstore
        wval = wvalstore
        kbarfcstind = kbarfcstindstore
        weight = weightstore
                    
        if (distIRF(zct,kct,t,simct,shockct)>0.0) then
        
        zval = z0(zct)
        kval = kdense0(kct)
        
        !block to determine Vnaval
        Vnaval = pval * freduced(zval,kval)
        kprimeval = max(kmin,(1.0-delta)*kval)
        
        statect=(act-1)*znum+zct
        do aprimect=1,anum
        do zprimect=1,znum
            stateprimect=(aprimect-1)*znum+zprimect
            

            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
            knum,kprimeval,Vnextval1)
            call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
            knum,kprimeval,Vnextval2)
            
            Vnextval = Vnextval1 + weight * (Vnextval2-Vnextval1)
            
            Vnaval = Vnaval + beta * pr_mat(statect,stateprimect) * Vnextval
        
        end do !zprimect
        end do !aprimect
        
        !block to construct Vaval
        RHSzct =zct
        RHSact = act
        Vaval = brent(k0(1),k0(knum/2),k0(knum),fnRHS,brenttol,kprimeval) !note this sets correct kprimeval, not correct vaval
        Vaval = -1.0 * Vaval
        
        !at this point, Vaval = -p k' + beta E V'
        !now, construct current period return remaining portions
        Vaval = Vaval + pval * (freduced(zval,kval)+(1.0-delta)*kval)
        
        !now process this info
        xival = (Vaval - Vnaval)/phi
        Vval = -1.0 * phi * expecxi(xival) + cdfxi(xival)*Vaval &
            + (1.0-cdfxi(xival))*Vnaval
            
        kprimep(zct,kct)=kprimeval
        xistarp(zct,kct)=xival
        padjustp(zct,kct) = cdfxi(xival)
        

    end if !distkz tolerance
    
    end do !kct
    end do !zct
    !$omp end do nowait
    !$omp end parallel
    !write(*,*) "made it past policy calcs"
    Yactualp = 0.0
    Iactualp = 0.0
    Kprimeactualp = 0.0
    Cactualp = 0.0
    Nactualp = 0.0
    
    
    do zct=1,znum
    do kct=1,kdensenum
         if (distIRF(zct,kct,t,simct,shockct)>0.0) then
        zval = z0(zct); 
        kval = kdense0(kct)
        kprimeval=kprimep(zct,kct)
        xival=xistarp(zct,kct)
        
        
        yval = yreduced(zval,kval)
        
        Yactualp = Yactualp  + distIRF(zct,kct,t,simct,shockct) * yval
        Iactualp = Iactualp + distIRF(zct,kct,t,simct,shockct)* cdfxi(xival) &
            * (kprimeval - (1.0-delta)*kval)
    
        Kprimeactualp = Kprimeactualp + distIRF(zct,kct,t,simct,shockct)* cdfxi(xival) * kprimeval &
            +  distIRF(zct,kct,t,simct,shockct) * (1.0 - cdfxi(xival)) * (1.0-delta)*kval
            
        Nactualp = Nactualp +  distIRF(zct,kct,t,simct,shockct) * (nreduced(zval,kval) + expecxi(xival))              
            
        end if
    end do !kct
    end do  !zct
    
    Cactualp = Yactualp - Iactualp

    !are you initializing the brent?
    if (piter==1) fa = (1.0/pval) - Cactualp
    if (piter==2) fb = (1.0/pval) - Cactualp
    if (piter==3) fc = (1.0/pval) - Cactualp
    
    if (piter==pcutoff) fa = (1.0/pval) - Cactualp
    if (piter==pcutoff+1) fb = (1.0/pval) - Cactualp
    if (piter==pcutoff+2) fc = (1.0/pval) - Cactualp

    perror = 1.0/pval - Cactualp
    
    !if not restarting or initializing
   if ((piter>3.and.piter<pcutoff).or.(piter>pcutoff+2)) then 
        if (perror<0) then
            pvalc = pvalb; fc = fb;
            pvalb = pval; fb = perror;
            !pval a doesn't change
        else if (perror>=0) then
            pvalc = pvala; fc = fa;
            pvala = pval; fa = perror;
            !pval b doesn't change
        end if
    end if
    
    if (abs(perror)<perrortol) exit

    end do !piter

    if (mod(t,numperIRF/2)==1) then 
    
        write(13,"(A,I5,A,I5,A,I5,A,I5,A,F7.5,A,F7.5)") "t= ",t,", sim = ",simct,", ct = ",shockct,&
            ", pit = ",piter,", err = ",perror,", p = ",pval
        write(*,"(A,I5,A,I5,A,I5,A,I5,A,F7.5,A,F7.5)") "t= ",t,", sim = ",simct,", ct = ",shockct,&
            ", pit = ",piter,", err = ",perror,", p = ",pval
    end if 
    
    !now, the value of p has been determined, in a market clearing fashion
    kbarsimIRF(t+1,simct,shockct)=Kprimeactualp
    psimIRF(t,simct,shockct) = pval
    perrorsimIRF(t,simct,shockct) = perror
    kbarfcstsimIRF(t+1,simct,shockct) = kbarfcstval
    pfcstsimIRF(t,simct,shockct) = pfcstval
    ysimIRF(t,simct,shockct) = Yactualp
    isimIRF(t,simct,shockct) = Iactualp
    NsimIRF(t,simct,shockct) = Nactualp
    
    
    !now that market clearing has occured, insert distribution into next period's
    do zct=1,znum
    do kct=1,kdensenum
        if (distIRF(zct,kct,t,simct,shockct)>0.0) then
        zval = z0(zct);
        kval = kdense0(kct);
        kprimeval = kprimep(zct,kct)
        kprimevalnoadj = max(kmin,(1.0-delta)*kval)
        xival = xistarp(zct,kct)
        
        !bracket kprimeval (which is kprime given that you adjust)
        kprimeind=kct
        call hunt(kdense0,kdensenum,kprimeval,kprimeind)
        
        if (kprimeind<=0) then
            kprimeind = 1; weight = 0.0
        else if (kprimeind>=1.and.kprimeind<=(kdensenum-1)) then
            weight = (kprimeval - kdense0(kprimeind))/(kdense0(kprimeind+1) - kdense0(kprimeind));
        else if (kprimeind>=kdensenum) then
            kprimeind = kdensenum-1; weight = 1.0
        end if
        
        !bracket kprimevalnoadj (which is kprime given that you don't adjust)
        kprimeindnoadj=kct
        call hunt(kdense0,kdensenum,kprimevalnoadj,kprimeindnoadj)
        
        if (kprimeindnoadj<=0) then
            kprimeindnoadj = 1; weightnoadj = 0.0
        else if (kprimeindnoadj>=1.and.kprimeindnoadj<=(kdensenum-1)) then
            weightnoadj = (kprimevalnoadj - kdense0(kprimeindnoadj))/(kdense0(kprimeindnoadj+1) - kdense0(kprimeindnoadj));
        else if (kprimeindnoadj>=kdensenum) then
            kprimeindnoadj = kdensenum-1; weightnoadj = 1.0
        end if
        
        
        do zprimect=1,znum
            
            !transfer the weight to kprimeind that does adjust
            distIRF(zprimect,kprimeind,t+1,simct,shockct) = distIRF(zprimect,kprimeind,t+1,simct,shockct) &
                + pr_mat_z(zct,zprimect) * padjustp(zct,kct) * (1.0 - weight) * distIRF(zct,kct,t,simct,shockct)
            
            !transfer the weight to kprimeind+1 that does adjust
            distIRF(zprimect,kprimeind+1,t+1,simct,shockct) = distIRF(zprimect,kprimeind+1,t+1,simct,shockct) &
                + pr_mat_z(zct,zprimect) * padjustp(zct,kct) * weight * distIRF(zct,kct,t,simct,shockct)
    
            !transfer the weight to kprimeindnoadj that doesn't adjust
            distIRF(zprimect,kprimeindnoadj,t+1,simct,shockct) = distIRF(zprimect,kprimeindnoadj,t+1,simct,shockct) &
                + pr_mat_z(zct,zprimect) * (1.0 - padjustp(zct,kct)) * (1.0 - weightnoadj) * distIRF(zct,kct,t,simct,shockct)
            
            !transfer the weight to kprimeindnoadj+1 that doesn't adjust
            distIRF(zprimect,kprimeindnoadj+1,t+1,simct,shockct) = distIRF(zprimect,kprimeindnoadj+1,t+1,simct,shockct) &
                + pr_mat_z(zct,zprimect) * (1.0 - padjustp(zct,kct)) * weightnoadj * distIRF(zct,kct,t,simct,shockct)
            
        end do !zprimect
        end if !distkz tolerance
    end do !kct
    end do !zct
    
    
end do !t
end do !simct
end do !shockct    

write(*,*) "Done with IRF at ",omp_get_wtime()-start," seconds."
write(13,*) "Done with IRF at ",omp_get_wtime()-start," seconds."

!psimIRF 
open(8,file="psimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) psimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)

!kbarsimIRF 
open(8,file="kbarsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) kbarsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)    
    
!perrorsimIRF 
open(8,file="perrorsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) perrorsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)      
    
!pfcstsimIRF 
open(8,file="pfcstsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) pfcstsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)          
    
!kbarfcstsimIRF 
open(8,file="kbarfcstsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) kbarfcstsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)       
    
!ysimIRF 
open(8,file="ysimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) ysimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)     
    
!isimIRF 
open(8,file="isimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) isimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)    
    
!NsimIRF 
open(8,file="NsimIRF.txt")
do simct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) NsimIRF(t,simct,shockct)
end do !shockct
end do !t
end do !ct
close(8)    
    
end if !doIRF==1

finish = omp_get_wtime()
write(13,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
write(*,"(A,F15.1,A)") "Finished in ",finish-start," seconds."
close(13); !closing log file

contains

subroutine update_fcst()
implicit none

!takes simulated values of p, K', and runs OLS to get new coefficient matrices

integer :: act,perct,t
double precision :: xmean,ymean,x2mean,y2mean,xval,yval

do act=1,anum
    
    !do kprime first
    xmean = 0.0
    ymean = 0.0
    x2mean = 0.0
    xymean = 0.0
    perct=0
    do t = numdiscard,numper-1
        if (asimpos(t)==act) then
            perct=perct+1
            xval = log(kbarsim(t))
            yval = log(kbarsim(t+1))

            xmean = xmean + xval
            ymean = ymean + yval
            x2mean = x2mean + xval ** 2.0
            xymean = xymean + xval * yval
        end if
    
    end do !period
    
    xmean = xmean / dble(perct)
    ymean = ymean / dble(perct)
    x2mean = x2mean / dble(perct)
    xymean = xymean / dble(perct)
    
    kfcstmatnew(act,2) = ( xymean - xmean * ymean ) / ( x2mean - ( xmean ** 2.0 ) )
    kfcstmatnew(act,1) = ymean - kfcstmatnew(act,2) * xmean
    
    !then do price rule
    xmean = 0.0
    ymean = 0.0
    x2mean = 0.0
    xymean = 0.0
    perct=0
    do t = numdiscard,numper-1
        if (asimpos(t)==act) then
            perct=perct+1
            xval = log(kbarsim(t))
            yval = log(psim(t))

            xmean = xmean + xval
            ymean = ymean + yval
            x2mean = x2mean + xval ** 2.0
            xymean = xymean + xval * yval
        end if
    
    end do !period
    
    xmean = xmean / dble(perct)
    ymean = ymean / dble(perct)
    x2mean = x2mean / dble(perct)
    xymean = xymean / dble(perct)
    
    pfcstmatnew(act,2) = ( xymean - xmean * ymean ) / ( x2mean - ( xmean ** 2.0 ) )
    pfcstmatnew(act,1) = ymean - pfcstmatnew(act,2) * xmean
    
end do !act

write(13,*) " "; write(13,*) "kfcstmat (old, new) "
write(*,*) " "; write(*,*) "kfcstmat (old, new) "
do act=1,anum
write(13,*) kfcstmat(act,:),kfcstmatnew(act,:)
write(*,*) kfcstmat(act,:),kfcstmatnew(act,:)
end do !act
write(13,*) " "
write(*,*) " "

write(13,*) "pfcstmat (old, new) "
write(*,*) "pfcstmat (old, new) "
do act=1,anum
write(13,*) pfcstmat(act,:),pfcstmatnew(act,:)
write(*,*) pfcstmat(act,:),pfcstmatnew(act,:)
end do !act
write(13,*) " "
write(*,*) " "

kfcsterror = maxval(abs(kfcstmat-kfcstmatnew))
pfcsterror = maxval(abs(pfcstmat-pfcstmatnew))

write(13,*) " "; write(13,*) "Fcst Ct = ",fcstiter," K Fcst error = ",kfcsterror;
write(*,*) " "; write(*,*) "Fcst Ct = ",fcstiter," K Fcst error = ",kfcsterror;
write(13,*) "Fcst Ct = ",fcstiter," P Fcst error = ",pfcsterror; write(13,*) " "
write(*,*) "Fcst Ct = ",fcstiter," P Fcst error = ",pfcsterror; write(*,*) " "

end subroutine update_fcst

subroutine fcst_process()
implicit none

!this subroutine takes kfcstmat and pfcstmat and produces pfcstvals and kbarfcstvals, kbarfcstinds, kbarfcstwgts

integer :: act,kbarct,ind
double precision :: wgt,kbarfcstval,pfcstval,kbarval

!given the forecast rules above, store the implied forecast values, linear interpolation indexes, and wgts
!in predetermined arrays (which will be incremented each forecast rule iteration)
kbarfcstinds(:,:) =0
kbarfcstvals(:,:) = 0.0
kbarfcstwgts(:,:) = 0.0
pfcstvals(:,:) = 0.0
do act=1,anum
do kbarct=1,kbarnum
    
    !what is the agg capital value?
    kbarval = kbar0(kbarct)
    
    !what is the forecast value of capital and of price?
    kbarfcstval = exp(kfcstmat(act,1) + kfcstmat(act,2)*log(kbarval))

    pfcstval= exp(pfcstmat(act,1) + pfcstmat(act,2)*log(kbarval))
    
    !restrict agg cap to the grid
    kbarfcstval = minval( (/ maxval( (/ kbar0(1) , kbarfcstval /) ) , kbar0(kbarnum) /) )
    
    pfcstvals(act,kbarct) = pfcstval
    kbarfcstvals(act,kbarct) = kbarfcstval
    
    !what are the linear weights and indexes for kbarfcstval?
    ind = kbarct
    call hunt(kbar0,kbarnum,kbarfcstval,ind)
    
    if (ind>=1.and.ind<kbarnum) then
        wgt = ( kbarfcstval - kbar0(ind) ) / ( kbar0(ind+1) - kbar0(ind) )
    else if (ind<=0) then
        ind =  1
        wgt = 0.0
    else if (ind>=kbarnum) then
        ind = kbarnum-1
        wgt = 1.0
    end if
    
    !insert into arrays for forecast storage
    kbarfcstinds(act,kbarct) = ind
    kbarfcstwgts(act,kbarct) = wgt
    
end do !kbarct
end do !act

end subroutine fcst_process

subroutine discretize_simulate()
implicit none

double precision :: asimgrid(anum)
integer :: act,zct,statect,ct,t,aprimect,shockct

!BEGIN DISCRETIZATION PORTION
!discretize idio prod process
call tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)

!discretize agg prod process
call tauchen(anum,rhoa,sigmaa,nstdeva,pr_mat_a,a0)

open(8,file="pr_mat_z.txt")
do zct=1,znum
write(8,*) pr_mat_z(zct,:)
end do !zct
close(8)

open(8,file="z0.txt")
do zct=1,znum
write(8,*) z0(zct)
end do !zct
close(8)

open(8,file="pr_mat_a.txt")
do act=1,anum
write(8,*) pr_mat_a(act,:)
end do !azct
close(8)

open(8,file="a0.txt")
do act=1,anum
write(8,*) a0(act)
end do !act
close(8)

!create unified transition matrix
do zct=1,znum
do act=1,anum
    statect = (act-1)*znum + zct
    do zprimect=1,znum
    do aprimect=1,anum
        stateprimect = (aprimect-1)*znum+zprimect
        pr_mat(statect,stateprimect)=pr_mat_a(act,aprimect)*pr_mat_z(zct,zprimect)
    end do
    end do
end do
end do

do statect=1,znum*anum
    pr_mat(statect,:) = pr_mat(statect,:) / sum(pr_mat(statect,:))
end do !statect
!!!END DISCRETIZATION PORTION

!draw random seeds
call random_seed(size=seeddim)

!insert random seed into seedarray
allocate(seedarray(seeddim))
do ct=1,seeddim
    seedarray(ct) = seedint + ct
end do !ct
call random_seed(put=seedarray)

!draw random shocks
call random_number(asimshock)
asimpos(1) = ainit

do t=2,numper
    asimgrid(1)=pr_mat_a(asimpos(t-1),1)
    do act=2,anum
        asimgrid(act) = asimgrid(act-1) + pr_mat_a(asimpos(t-1),act)
    end do !act
    
      
    if (asimshock(t)<asimgrid(1)) then 
        aprimect=1; 
    end if
    do act=2,anum
        if (asimgrid(act-1)<=asimshock(t).and.asimshock(t)<asimgrid(act)) then
            aprimect=act
            
        end if
    end do !act
    asimpos(t) = aprimect
end do

open(8,file="asimpos.txt")
do t=1,numper
write(8,*) asimpos(t)
end do !t
close(8)

open(8,file="asimshock.txt")
do t=1,numper
write(8,*) asimshock(t)
end do !t
close(8)
!END UNCONDITIONAL SIMULATION

!START IRF SIMULATION
!if doing IRF, determine ergodic distribution of agg prod, then IRF shock prob, then perform IRF agg prod simulation
if (doIRF==1) then

call random_number(asimshockIRF)
call random_number(arriveshockIRF)

!initialize ergodic distribution
ergdistaold(:) = 0.0
ergdistaold(anum/2) = 1.0

do ct=1,maxergit
    do act=1,anum
    do aprimect=1,anum
        ergdista(aprimect) = ergdista(aprimect) + pr_mat_a(act,aprimect)*ergdistaold(act)
    end do !aprimect
    end do !act
    
    if (maxval(abs(ergdista-ergdistaold))<ergdistatol) exit
    
    ergdistaold = ergdista
    ergdista(:) = 0.0
    
end do !ct

!now that the ergodic distribution is done, compute size of IRF shock probability
shockprobdenom = 0.0
shockprobdenom = log(a0(shockanumIRF))
do act=1,anum
    shockprobdenom = shockprobdenom  - ergdista(act) * log(a0(act))
end do !act
shockprobIRF = shocksizeIRF / shockprobdenom

!now, perform the simulation
!asimshockIRF(numperIRF,numsimIRF),asimposIRF(numperIRF,numsimIRF,2)

asimposIRF(1,:,:) = ainit

do ct=1,numsimIRF !counting IRF simulations
do t=2,numperIRF !counting IRF
    
    if (t<shockperIRF) then
        
        !both transit as normal
        act = asimposIRF(t-1,ct,1)
        
        !bins to use for comparison
        asimgrid(:) = 0.0
        asimgrid(1) = pr_mat_a(act,1)
        do aprimect=2,anum
            asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
        end do !aprimect
        
        if (asimshockIRF(t,ct)<asimgrid(1)) then
            asimposIRF(t,ct,1) = 1
            asimposIRF(t,ct,2) = 1
        else if (asimshockIRF(t,ct)>=asimgrid(1)) then
            do aprimect=2,anum
                if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                    asimposIRF(t,ct,1) = aprimect
                    asimposIRF(t,ct,2) = aprimect
                end if
            end do !aprimect
        end if
        
    else if (t==shockperIRF) then
        
        !first, decide if you're transiting normally or not
        if (arriveshockIRF(ct)>shockprobIRF) then
                
                !in this case, transit normally as above
                act = asimposIRF(t-1,ct,1)
                
                !bins to use for comparison
                asimgrid(:) = 0.0
                asimgrid(1) = pr_mat_a(act,1)
                do aprimect=2,anum
                    asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
                end do !aprimect
                
                if (asimshockIRF(t,ct)<asimgrid(1)) then
                    asimposIRF(t,ct,1) = 1
                    asimposIRF(t,ct,2) = 1
                else if (asimshockIRF(t,ct)>=asimgrid(1)) then
                    do aprimect=2,anum
                        if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                            asimposIRF(t,ct,1) = aprimect
                            asimposIRF(t,ct,2) = aprimect
                        end if
                    end do !aprimect
                end if

        else if (arriveshockIRF(ct)<=shockprobIRF) then
            !in this case, transit normally only for version 2, not version 1
            
            !version 2 is shocked, goes to shockasnumIRF position
            asimposIRF(t,ct,2) = shockanumIRF
            
            !version 1 transits normally
            act = asimposIRF(t-1,ct,1)
            
            !bins to use for comparison
            asimgrid(:) = 0.0
            asimgrid(1) = pr_mat_a(act,1)
            do aprimect=2,anum
                asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
            end do !aprimect
            
            if (asimshockIRF(t,ct)<asimgrid(1)) then
                asimposIRF(t,ct,1) = 1
            else if (asimshockIRF(t,ct)>=asimgrid(1)) then
                do aprimect=2,anum
                    if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                        asimposIRF(t,ct,1) = aprimect
                    end if
                end do !aprimect
            end if
         
        end if 

    else if (t>shockperIRF) then
        
        !both transit as normal, but separately

        !version 1 transits normally
        act = asimposIRF(t-1,ct,1)
        
        !bins to use for comparison
        asimgrid(:) = 0.0
        asimgrid(1) = pr_mat_a(act,1)
        do aprimect=2,anum
            asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
        end do !aprimect
        
        if (asimshockIRF(t,ct)<asimgrid(1)) then
            asimposIRF(t,ct,1) = 1
        else if (asimshockIRF(t,ct)>=asimgrid(1)) then
            do aprimect=2,anum
                if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                    asimposIRF(t,ct,1) = aprimect
                end if
            end do !aprimect
        end if
            
        !version 2 transits normally
        act = asimposIRF(t-1,ct,2)
        
        !bins to use for comparison
        asimgrid(:) = 0.0
        asimgrid(1) = pr_mat_a(act,1)
        do aprimect=2,anum
            asimgrid(aprimect) = asimgrid(aprimect-1) + pr_mat_a(act,aprimect)
        end do !aprimect
        
        if (asimshockIRF(t,ct)<asimgrid(1)) then
            asimposIRF(t,ct,2) = 1
        else if (asimshockIRF(t,ct)>=asimgrid(1)) then
            do aprimect=2,anum
                if (asimshockIRF(t,ct)>=asimgrid(aprimect-1).and.asimshockIRF(t,ct)<asimgrid(aprimect)) then
                    asimposIRF(t,ct,2) = aprimect
                end if
            end do !aprimect
        end if
    
    end if
    
end do !t
end do !ct


open(8,file="asimposIRF.txt")
do ct=1,numsimIRF
do t=1,numperIRF
do shockct=1,2
    write(8,*) asimposIRF(t,ct,shockct)
end do !shockct
end do !t
end do !ct
close(8)

end if !doIRF==1
!!END IRF SIMULATION

end subroutine discretize_simulate

double precision function cdfxi(xi)
implicit none
double precision :: xi
if (xi<0.0) cdfxi = 0.0
if (xi>0.0.and.xi<=xibar) cdfxi = xi/xibar
if (xi>xibar) cdfxi = 1.0
end function cdfxi
    
double precision function expecxi(xi)
implicit none
double precision :: xi
if(xi<0.0) expecxi = 0.0
if(xi<=xibar.and.xi>=0.0)expecxi = ( xi ** 2.0 ) / (2.0*xibar)
if(xi>xibar)expecxi = ( xibar ** 2.0 ) / (2.0*xibar)
end function expecxi

double precision function nreduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)
nreduced = ( nu **  exponentnu ) * ( wval ** ( -1.0 * exponentnu ) )
nreduced = nreduced * (z ** exponentnu) * (aval ** exponentnu) * ( k ** (alpha * exponentnu) )

end function nreduced

double precision function freduced(z,k)
implicit none
double precision :: z,k
double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)
freduced = ( nu ** (nu * exponentnu) ) * (1.0 - nu) * ( wval ** ( -1.0 * nu * exponentnu ) )
freduced = freduced * (z ** exponentnu) * (aval ** exponentnu) * ( k ** (alpha * exponentnu) )

end function freduced

double precision function yreduced(z,k)
implicit none
double precision :: z,k

double precision :: exponentnu

exponentnu = 1.0 / (1.0 - nu)
yreduced = ( (nu/wval) ** (nu * exponentnu) ) * ( (aval * z ) ** exponentnu  )
yreduced = yreduced * ( k ** ( alpha * exponentnu ) )

end function yreduced

double precision function fnRHS(kprimeval)
implicit none

!returns -1 (-pk' + beta E V(z',k',a',K'^F) ) given k'

double precision :: kprimeval,Vnextval1,Vnextval2
integer :: zprimect,aprimect,stateprimect,statect

!current return    
fnRHS = -pval * kprimeval

!note that RHSzct is the current zct, and RHSact is the current act, with diff naming so not a loop index themselves
statect=(RHSact-1)*znum+RHSzct
do zprimect=1,znum
do aprimect=1,anum
    
    stateprimect=(aprimect-1)*znum+zprimect
    call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind),V2old(zprimect,aprimect,:,kbarfcstind),&
        knum,kprimeval,Vnextval1)
    call splint(k0,Vold(zprimect,aprimect,:,kbarfcstind+1),V2old(zprimect,aprimect,:,kbarfcstind+1),&
        knum,kprimeval,Vnextval2)
        
    fnRHS = fnRHS + beta*pr_mat(statect,stateprimect) * ( (1.0 - weight) * Vnextval1 + weight * Vnextval2)        
    
end do !aprimect 
end do !zprimect


!you want to MINIMIZE this
fnRHS = -1.0 * fnRHS

end function fnRHS

end program kt_ks