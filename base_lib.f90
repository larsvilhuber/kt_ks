module base_lib
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!base_lib.f90
!
!This version of the base_lib module is current as of 2/25/14. This
!file contains a number of utility functions and subroutines, some
!original and some taken from Numerical Recipes in Fortran.
!
! Stephen Terry
! This Version: 02/25/14
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

double precision function erfcc(x)
implicit none

!input/output declarations
double precision :: x

!other declarations
double precision :: t,z

z = abs(x)
t = 1.0 / (1.0 + 0.5*z)

erfcc = t * exp(-z * z-1.26551223+t*(1.00002368+t*(0.37409196+&
    t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+&
    t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))

if (x.lt.0.0) erfcc = 2.0 - erfcc
end function erfcc

double precision function normcdf(x,mu,sigma)
implicit none

!input/output declarations
!x: the input value for the normal cdf
!mu: the mean of the normal distribution
!sigma: the standard deviation of the normal distribution
double precision :: x,mu,sigma

!other declarations
double precision :: z

!standardized value ~ N(0,1)
z = (x-mu)/sigma

normcdf =  0.5 * erfcc( ( -1.0 * z ) / sqrt(2.0) )
    
end function normcdf

subroutine tauchen(znum,rhoz,sigmaz,nstdevz,pr_mat_z,z0)
implicit none

!input/output declarations
integer, intent(in) :: znum
double precision, intent(in) :: rhoz,sigmaz,nstdevz
double precision, intent(out) :: pr_mat_z(znum,znum),z0(znum)
    
!other declarations
integer :: zct,zprimect
double precision :: zmin,zmax,gridinc,stdev

!determine end points of the grid (log space)
stdev =  ((sigmaz**2.0)/(1-rhoz**2.0))**0.5
zmin = - nstdevz * stdev
zmax = nstdevz * stdev

!insert points into the grid (log space)
call linspace(z0,zmin,zmax,znum)
gridinc = z0(2)-z0(1)

!loop over z states
do zct=1,znum
    
    !insert transition matrix middle rows
    do zprimect=2,(znum-1)
        pr_mat_z(zct,zprimect) = &
            normcdf(z0(zprimect)+gridinc/2.0,rhoz*z0(zct),sigmaz) - &
            normcdf(z0(zprimect)-gridinc/2.0,rhoz*z0(zct),sigmaz)
    end do !zct
    
    !first interval and last interval take the rest of the weight
    pr_mat_z(zct,1) = normcdf(z0(1)+gridinc/2.0,rhoz*z0(zct),sigmaz)
    pr_mat_z(zct,znum) = 1.0 - normcdf(z0(znum)-gridinc/2.0,rhoz*z0(zct),sigmaz)
    
end do !zct

!round the transition matrix
do zct=1,znum
    pr_mat_z(zct,:) = pr_mat_z(zct,:)/sum(pr_mat_z(zct,:))
end do !zct

!convert grid back to z-space
z0 = exp(z0)
    
end subroutine tauchen

subroutine qsimpweightsnodes(a,b,n,weights,nodes)
implicit none

!a = lower limit of integration
!b = upper limit of integration
!n = number of intervals, must be even.  There will be n+1 weights/nodes
!weights = n+1 x 1 vector of integration weights
!nodes = n+1 x 1 vector integration nodes

!input/output declarations
integer :: n
double precision :: a,b,weights(n+1),nodes(n+1)

!other declarations
integer :: ct

weights(:) = 1.0

do ct=1,n+1
    nodes(ct) = a + ((b-a)/dble(n)) * (dble(ct)-1.0)
end do !ct

do ct=2,n
    if (mod(ct,2)==0) then 
        weights(ct) = 4.0
    else 
        weights(ct) = 2.0
    end if
end do !ct

weights = weights * ( (b-a) / ( 3.0 * dble(n) ) )

end subroutine qsimpweightsnodes

subroutine amoeba(p,y,ndim,ftol,iter,f,param,nparam)
implicit none
!input/output declarations
integer :: iter,ndim,NMAX,ITMAX,nparam
double precision :: ftol,p(ndim+1,ndim),y(ndim+1),f,TINY,param(nparam)
parameter(NMAX=20,ITMAX=5000,TINY=1.0e-10)
external :: f

!multidimensional minimization of the function f(x,param,nparam), where
!x is an ndim x 1 vector

!on input, p must be a ndim+1 x ndim matrix with each row an endpoint of the initial
!simplex, and y is an ndim+1 x 1 vector with the value of f at each of the vertices in p

!ftol is the fractional convergence tolerance to be achieved in the function value, and
!on output, p and y are ndim+1 points, all within ftol of the converged minimum function value,
!with iter yielding the number of function evaluations required

integer :: i,ihi,ilo,inhi,j,m,n
double precision :: rtol,sumval,swap,ysave,ytry,psum(ndim),ptry(ndim),fac1,fac2,fac
write(*,*) f((/0.0,0.0,0.0/),param,nparam)
write(*,*) "p = "
write(*,*) p(1,:)
write(*,*) p(2,:)
write(*,*) p(3,:)
write(*,*) p(4,:)

write(*,*) "y = ",y
write(*,*) "ndim = ",ndim
write(*,*) "ftol = ",ftol
write(*,*) "nparam = ",nparam
write(*,*) "param = ",param



iter = 0
1 do n=1,ndim
    sumval = 0.0
    do m=1,ndim+1
        sumval=sumval + p(m,n)
    end do !m
    psum(n) = sumval
end do !n
2 ilo=1
if (y(1).gt.y(2)) then
    ihi=1
    inhi=2
else
    ihi=2
    inhi=1
end if
do i=1,ndim+1
    if (y(i).le.y(ilo)) ilo=i
    if (y(i).gt.y(ihi)) then
        inhi=ihi
        ihi=i
    else if (y(i).gt.y(inhi)) then
        if(i.ne.ihi) inhi=i
    end if
end do !i
rtol = 2.0 * abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
if (rtol.lt.ftol) then
    swap = y(1)
    y(1)=y(ilo)
    y(ilo)=swap
    do n=1,ndim
        swap = p(1,n)
        p(1,n) = p(ilo,n)
        p(ilo,n) = swap
    end do !n
    return
end if
iter = iter+2

!ytry=amotry(p,y,psum,ndim,ihi,dble(-1.0),f,param,nparam)
fac=-1.0
fac1=(1.0-fac)/dble(ndim)
fac2=fac1-fac
do j=1,ndim
    ptry(j) = psum(j)*fac1 - p(ihi,j)*fac2
end do !j
write(*,*) "ptry = ",ptry
write(*,*) "param = ",param
write(*,*) "nparam = ",nparam
ytry = f(ptry,param,nparam)
write(*,*) "made it here"
if (ytry.lt.y(ihi)) then
    y(ihi) = ytry
    do j=1,ndim
        psum(j) = psum(j)-p(ihi,j)+ptry(j)
        p(ihi,j) = ptry(j)
    end do !j
end if

if (ytry.le.y(ilo)) then
        
!ytry=amotry(p,y,psum,ndim,ihi,dble(2.0),f,param,nparam)
fac=2.0
fac1=(1.0-fac)/dble(ndim)
fac2=fac1-fac
do j=1,ndim
    ptry(j) = psum(j)*fac1 - p(ihi,j)*fac2
end do !j

ytry = f(ptry,param,nparam)

if (ytry.lt.y(ihi)) then
    y(ihi) = ytry
    do j=1,ndim
        psum(j) = psum(j)-p(ihi,j)+ptry(j)
        p(ihi,j) = ptry(j)
    end do !j
end if
        
        
else if (ytry.ge.y(inhi)) then
    ysave =y(ihi)
    !ytry=amotry(p,y,psum,ndim,ihi,dble(0.5),f,param,nparam)
    fac=0.5
    fac1=(1.0-fac)/dble(ndim)
    fac2=fac1-fac
    do j=1,ndim
        ptry(j) = psum(j)*fac1 - p(ihi,j)*fac2
    end do !j

    ytry = f(ptry,param,nparam)

    if (ytry.lt.y(ihi)) then
        y(ihi) = ytry
        do j=1,ndim
            psum(j) = psum(j)-p(ihi,j)+ptry(j)
            p(ihi,j) = ptry(j)
        end do !j
    end if

    if (ytry.ge.ysave) then
        do i=1,ndim+1
            if(i.ne.ilo) then
                do j=1,ndim
                    psum(j)=0.5*(p(i,j)+p(ilo,j))
                    p(i,j) = psum(j)
                end do !j
                y(i)=f(psum,param,nparam)
            end if
        end do !i
        iter=iter+ndim
        goto 1
    end if
else
    iter=iter-1
end if
goto 2
end subroutine amoeba

subroutine qromb(a,b,ss,tol,k,jmax,f,param,nparam)
implicit none
integer :: jmax,nparam,k
double precision :: a,b,ss,tol,f,param(nparam)
external :: f

!note dependency on trapzd and polint

!a,b are the scalar limits of integeration of the function f(x,param)

!param is nparam x 1 real vector of parameters to pass to f

!the numerical integeration is performed using "Romberg's method of 
!order 2k", as described in Numerical Recipes in Fortran. k=2 is simpson's
!rule, and k=5 is typically fine. k>=2 is required.

!tol is the convergence criterion determining how many evaluations of the
!integral are required (tol = 1.0e-6 is fine)

!jmax is the maximum number of steps in the integration (jmax=20 is fine)

!ss is the returned value of the integral

!auxiliary declarations
integer :: jmaxp,km,j
double precision :: dss
double precision, allocatable :: h(:),s(:)

jmaxp=jmax+1
km=k-1

allocate(h(jmaxp),s(jmaxp))

h(1) = 1.0
do j=1,jmax
    call trapzd(a,b,s(j),j,f,param,nparam)
    if (j.ge.k) then
        call polint(h(j-km),s(j-km),k,dble(0.0),ss,dss)
        if(abs(dss).le.tol*abs(ss)) return
    end if
    s(j+1)=s(j)
    h(j+1)=0.25*h(j)
end do !j
end subroutine qromb

subroutine trapzd(a,b,s,n,f,param,nparam)
implicit none

!input/output declarations
integer :: n,nparam
double precision :: a,b,s,f,param(nparam)
external :: f

!this routine computes the nth state of refinement of an extended trapezoidal rule

!a,b are the limits of integration
!s is the returned value of the integral, but you need to call multiple times, with successively
!larger values of n and no change to s, to get the best evaluation of int_a^b f(x)dx in the returned
!value of s. s is both an input and an an output

integer :: it,j
double precision :: del,sumval,tnm,x
    
if (n.eq.1) then
    s=0.5*(b-a)*(f(a,param,nparam)+f(b,param,nparam))
else
    it=2**(n-2)
    tnm=it
    del=(b-a)/tnm
    x=a+0.5*del
    sumval=0.0
    do j=1,it
        sumval = sumval + f(x,param,nparam)
        x=x+del
    end do !j
    s = 0.5 * (s + (b-a)*sumval/tnm)
end if
return
end subroutine trapzd

subroutine polint(xa,ya,n,x,y,dy)
implicit none
!input/output declarations
integer :: n,NMAX
double precision :: dy,x,y,xa(n),ya(n)
parameter(NMAX=20)

!given arrays of length xa(n) and ya(n), the routine 
!returns a value y = P(x), where P(x) is the n-1 degree
!interpolating polynomial for the n points xa --> ya
!dy is an error estimate

integer :: i,m,ns
double precision :: den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
ns=1
dif=abs(x-xa(1))
do i=1,n
    dift=abs(x-xa(i))
    if (dift.lt.dif) then
        ns=i
        dif=dift
    end if
    c(i)=ya(i)
    d(i)=ya(i)
end do !i
y=ya(ns)
ns=ns-1
do m=1,n-1
    do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
    end do !i
    if (2*ns.lt.n-m) then
        dy=c(ns+1)
    else
        dy=d(ns)
        ns=ns-1
    end if
    y=y+dy
end do !m
return
end subroutine polint

subroutine heapsort(n,ra)
implicit none
integer :: n
double precision :: ra(n)

integer :: i,ir,j,l
double precision :: rra

if (n.lt.2) return

l = n/2+1
ir=n

10 continue
    if(l.gt.1) then
        l=l-1
        rra=ra(l)
    else 
        rra=ra(ir)
        ra(ir)=ra(1)
        ir=ir-1
        if(ir.eq.1) then
            ra(1)=rra
            return
        end if
    end if
    i=l
    j=l+1
20 if(j.le.ir)then
        if(j.lt.ir)then
            if(ra(j).lt.ra(j+1)) j = j+1
        end if
        if(rra.lt.ra(j)) then
            ra(i)=ra(j)
            i=j
            j=j+j
        else
            j=ir+1
        end if
    goto 20
end if
ra(i)=rra
goto 10
end subroutine heapsort

subroutine pso(x,fx,f,lb,ub,nvar,npart,xtol,xquicktol,xquicknum,ftol,maxit,phi)
implicit none

!input/output declarations
integer :: nvar,npart,maxit,xquicknum
double precision :: x(nvar),f,lb(nvar),ub(nvar),xtol,phi(2),fx,ftol,&
    xquicktol
external :: f

!local declarations
integer :: iter,seedint,varct,partct,n
integer, allocatable :: seed(:)
double precision :: xnorm,chi,phisum,xstore(nvar,npart),fxstore(npart),&
    pbeststore(nvar,npart),fpbeststore(npart),globbest(nvar),fglobbest,vstore(nvar,npart),&
    pbestshocks(nvar),globbestshocks(nvar),spacenorm,reltol,xquicknorm,quickthresh,fpbestsort(npart)

seedint=2501


call random_seed(size=n)
allocate(seed(n))
do iter=1,n
    seed(iter) = seedint+iter
end do !iter
call random_seed(put=seed)

!SET CONSTRAINT PARAMETER
!note that sum(phi) > 4 is required, with phi(1)=phi(2)=2.05 most common in lit
phisum = sum(phi)
chi = 2.0 / ( phisum - 2.0 + sqrt( ( phisum ** 2.0 ) - 4.0 * phisum )  )

!SET THE WIDTH OF THE SPACE TO INITIALIZE VELOCITIES
spacenorm = sqrt(sum((ub-lb)**2.0))



!!!INITIALIZE THE SWARM 
fglobbest = 1.0e30
do partct=1,npart
    call random_number(xstore(:,partct))
    xstore(:,partct) = lb * (1.0-xstore(:,partct)) + ub * xstore(:,partct)
    call random_number(vstore(:,partct))
    vstore(:,partct) = -1.0*spacenorm * (1.0 - vstore(:,partct)) + spacenorm * vstore(:,partct)
    fxstore(partct) = f(xstore(:,partct))
    
    if (fxstore(partct)<fglobbest) then
        globbest = xstore(:,partct)
        fglobbest = fxstore(partct)
    end if
    
end do !partct

pbeststore = xstore
fpbeststore = fxstore

!!!NOW, LOOP OVER THE ITERATIONS
do iter=1,maxit
    
    !!!UPDATE VELOCITIES AND POSITIONS
    do partct=1,npart
        call random_number(pbestshocks)
        call random_number(globbestshocks)
        vstore(:,partct) = vstore(:,partct) + phi(1) * pbestshocks * (pbeststore(:,partct)-xstore(:,partct)) &
            + phi(2)*globbestshocks * (globbest - xstore(:,partct))
        vstore(:,partct) = chi * vstore(:,partct)
        xstore(:,partct) = xstore(:,partct) + vstore(:,partct)
        do varct=1,nvar
            xstore(varct,partct) = min(max(lb(varct),xstore(varct,partct)),ub(varct))
        end do !varct
        fxstore(partct) = f(xstore(:,partct))
    end do !parct
    
    !!!UPDATE PERSONAL BEST CONDITIONS & GLOBAL BEST CONDITION
    do partct=1,npart
        if (fxstore(partct)<fpbeststore(partct)) then
            pbeststore(:,partct) = xstore(:,partct)
            fpbeststore(partct) = fxstore(partct)
        end if
        if (fxstore(partct)<fglobbest) then
            globbest = xstore(:,partct)
            fglobbest = fxstore(partct)
        end if
    end do !partct
    
    !!!!COMPUTE CONVERGENCE CRITERIA
    !sort to find the best points
    fpbestsort=fpbeststore
    call heapsort(npart,fpbestsort)
    quickthresh = fpbestsort(xquicknum)
    
    
    xnorm = 0.0
    reltol = 0.0
    xquicknorm = 0.0
    
    do partct=1,npart
        !max distance
        xnorm = max(xnorm,sqrt(sum((globbest-xstore(:,partct))**2.0)))
        
        !max difference between personal and global bests
        reltol = max(reltol,abs(fglobbest-fpbeststore(partct))/&
    (abs(fglobbest)+abs(fpbeststore(partct))+0.01))
    
        !quick max distance measure
        if (fpbeststore(partct)<quickthresh) then
            xquicknorm = max(xquicknorm,sqrt(sum((globbest-xstore(:,partct))**2.0)))
        end if
    end do !partct
    
    if (xnorm<xtol.or.reltol<ftol.or.xquicknorm<xquicktol) exit
    
end do !iter

x = globbest
fx = fglobbest

return
end subroutine pso



subroutine psorestart(x,fx,f,lb,ub,nvar,npart,xtol,xquicktol,xquicknum,ftol,maxit,phi,newrun,psoseed)
implicit none

!note: this is a much slower version of the PSO algorithm, suitable for 
!function evaluations which take a long time and are therefore costly to
!redo if the program is interrupted.  This version can restart the PSO algorithm 
!by rerunning the program with newrun=0.  The subroutine works via extensive text
!file input/output, which is why it is so slow.

!input/output declarations
integer :: nvar,npart,maxit,xquicknum,newrun,psoseed
double precision :: x(nvar),f,lb(nvar),ub(nvar),xtol,phi(2),fx,ftol,&
    xquicktol
external :: f

!local declarations
integer :: iter,varct,partct,n,ct,stiter,randct,oldpartct,olditer,oldrandct
integer, allocatable :: seed(:)
double precision :: xnorm,chi,phisum,xstore(nvar,npart),fxstore(npart),&
    pbeststore(nvar,npart),fpbeststore(npart),globbest(nvar),fglobbest,vstore(nvar,npart),&
    pbestshocks(nvar),globbestshocks(nvar),spacenorm(nvar),reltol,xquicknorm,quickthresh,fpbestsort(npart),&
    randshocks(2*(maxit+1)*npart*nvar)


!first, store the parameters if it's a new run - you can check an old version of this to make sure you match before restarting
if (newrun==1) then
open(9,file="PSOparam.txt")
write(9,*) "nvar = ",nvar
write(9,*) "npart = ",npart
write(9,*) "xtol = ",xtol
write(9,*) "xquicktol = ",xquicktol
write(9,*) "xquicknum = ",xquicknum
write(9,*) "ftol = ",ftol
write(9,*) "maxit = ",maxit
write(9,*) "phi = ",phi(1),phi(2)
write(9,*) "newrun = ",newrun
write(9,*) "psoseed = ",psoseed
write(9,*) "lb = ",lb
write(9,*) "ub = ",ub
close(9)
end if

!do the seeding based on psoseed
call random_seed(size=n)
allocate(seed(n))
do iter=1,n
    seed(iter) = psoseed+iter
end do !iter
call random_seed(put=seed)

!now, here's the part where draw the random numbers

!note that if you use the same psoseed, and the same compiler this will produce IDENTICAL U(0,1) shocks
call random_number(randshocks)

!set parameter governing velocity updates
!note that sum(phi) > 4 is required, with phi(1)=phi(2)=2.05 most common in lit
phisum = sum(phi)
chi = 2.0 / ( phisum - 2.0 + sqrt( ( phisum ** 2.0 ) - 4.0 * phisum )  )

!set the width of the space to initialize velocities
do ct=1,nvar
spacenorm(ct) = abs(ub(ct)-lb(ct))/dble(3.0)
end do !ct
!now, you need to determine whether this is a new run or an old run, and perform correct initializations
if (newrun==1) then

iter = 0; !initialize the iteration counter
randct = 0; !initialize the random number counter

!in this case, proceed as before for the swarm initialization
xstore(:,:) = 0.0
fxstore(:) = 0.0
vstore(:,:) = 0.0

do partct=1,npart
    
    !call random numbers for initial scatter shot positions of x
    xstore(:,partct) = randshocks((randct+1):(randct+nvar))
    randct = randct + nvar
    
    !distribution x across the space based on the U(0,1) from above
    xstore(:,partct) = lb * (1.0-xstore(:,partct)) + ub * xstore(:,partct)
    
    !call random numbers for initial velocity draws
    vstore(:,partct) = randshocks((randct+1):(randct+nvar))
    randct = randct + nvar
    
    !distribution of velocities based on the width of the space
    vstore(:,partct) = -1.0*spacenorm(:) * (1.0 - vstore(:,partct)) + spacenorm(:) * vstore(:,partct)
    
    !evaluate the function at the positions above
    fxstore(partct) = f(xstore(:,partct))
    
    !record progress so far
    open(9,file="PSOprogress.txt")
    write(9,*) iter
    write(9,*) partct
    write(9,*) randct
    close(9)
    
    open(9,file="PSOxstore.txt")
    do ct=1,nvar
        write(9,*) xstore(ct,:)
    end do !ct
    close(9)
    
    open(9,file="PSOvstore.txt")
    do ct=1,nvar
        write(9,*) vstore(ct,:)
    end do !ct
    close(9)
    
    open(9,file="PSOfxstore.txt")
    do ct=1,npart
        write(9,*) fxstore(ct)
    end do !ct
    close(9)
        
    
end do !partct


fglobbest = 1.0e30
do partct=1,npart
if (fxstore(partct)<fglobbest) then
    globbest = xstore(:,partct)
    fglobbest = fxstore(partct)
end if
end do !partct

pbeststore = xstore
fpbeststore = fxstore

!now, record personal and global bests from initial evaluations
open(9,file="PSOpbeststore.txt")
do ct=1,nvar
    write(9,*) pbeststore(ct,:)
end do !ct
close(9)

open(9,file="PSOfpbeststore.txt")
do ct=1,npart
    write(9,*) fpbeststore(ct)
end do !ct
close(9)

open(9,file="PSOglobbest.txt")
do ct=1,nvar
    write(9,*) globbest(ct)
end do !ct
close(9)

open(9,file="PSOfglobbest.txt")
write(9,*) fglobbest
close(9)

!note that at this point you have initialized 
!xstore(nvar,npart) - the current swarm position
!fxstore(npart) - the current swarm objective values
!vstore(nvar,npart) - the current swarm velocities
!pbeststore(nvar,npart) - the current swarm personal best positions
!fpbeststore(npart) - the current swarm personal best objective values
!globbest(nvar) - the global swarm best position
!fglobbest(1) - the globabl swarm best objective value

stiter = 1

else if (newrun==0) then

!in this case, you need to know which iteration you're on

!read progress so far
open(9,file="PSOprogress.txt")
read(9,*) olditer
read(9,*) oldpartct
read(9,*) oldrandct
close(9)

if (olditer==0) then
    !in this case, you need to finish initialization
    
    !first, read in everything so far
    open(9,file="PSOxstore.txt")
    do ct=1,nvar
        read(9,*) xstore(ct,:)
    end do !ct
    close(9)
    
    open(9,file="PSOvstore.txt")
    do ct=1,nvar
        read(9,*) vstore(ct,:)
    end do !ct
    close(9)
    
    open(9,file="PSOfxstore.txt")
    do ct=1,npart
        read(9,*) fxstore(ct)
    end do !ct
    close(9)
    
    randct = oldrandct
    iter = olditer
    
    !then, finish the initialization if there's any left to be done
    if (oldpartct<npart) then
    do partct = oldpartct+1,npart
        
        !call random numbers for initial scatter shot positions of x
        xstore(:,partct) = randshocks((randct+1):(randct+nvar))
        randct = randct + nvar
        
        !distribution x across the space based on the U(0,1) from above
        xstore(:,partct) = lb * (1.0-xstore(:,partct)) + ub * xstore(:,partct)
        
        !call random numbers for initial velocity draws
        vstore(:,partct) = randshocks((randct+1):(randct+nvar))
        randct = randct + nvar
        
        !distribution of velocities based on the width of the space
        vstore(:,partct) = -1.0*spacenorm(:) * (1.0 - vstore(:,partct)) + spacenorm(:) * vstore(:,partct)
        
        !evaluate the function at the positions above
        fxstore(partct) = f(xstore(:,partct))
        
        !record progress so far
        open(9,file="PSOprogress.txt")
        write(9,*) iter
        write(9,*) partct
        write(9,*) randct
        close(9)
        
        open(9,file="PSOxstore.txt")
        do ct=1,nvar
            write(9,*) xstore(ct,:)
        end do !ct
        close(9)
        
        open(9,file="PSOvstore.txt")
        do ct=1,nvar
            write(9,*) vstore(ct,:)
        end do !ct
        close(9)
        
        open(9,file="PSOfxstore.txt")
        do ct=1,npart
            write(9,*) fxstore(ct)
        end do !ct
        close(9)

        
    end do !partct
    end if !oldpartct<npart
    
    fglobbest = 1.0e30
    do partct=1,npart
    if (fxstore(partct)<fglobbest) then
        globbest = xstore(:,partct)
        fglobbest = fxstore(partct)
    end if
    end do 
    
    
pbeststore = xstore
fpbeststore = fxstore

!now, record personal and global bests from initial evaluations
open(9,file="PSOpbeststore.txt")
do ct=1,nvar
    write(9,*) pbeststore(ct,:)
end do !ct
close(9)

open(9,file="PSOfpbeststore.txt")
do ct=1,npart
    write(9,*) fpbeststore(ct)
end do !ct
close(9)

open(9,file="PSOglobbest.txt")
do ct=1,nvar
    write(9,*) globbest(ct)
end do !ct
close(9)

open(9,file="PSOfglobbest.txt")
write(9,*) fglobbest
close(9)

!note that at this point you have initialized 
!xstore(nvar,npart) - the current swarm position
!fxstore(npart) - the current swarm objective values
!vstore(nvar,npart) - the current swarm velocities
!pbeststore(nvar,npart) - the current swarm personal best positions
!fpbeststore(npart) - the current swarm personal best objective values
!globbest(nvar) - the global swarm best position
!fglobbest(1) - the globabl swarm best objective value

stiter = olditer + 1

else if (olditer>=1) then
    !in this case, you just need to finish off the previous iteration
    
    !first, read in everything that's been done
    open(9,file="PSOxstore.txt")
    do ct=1,nvar
        read(9,*) xstore(ct,:)
    end do !ct
    close(9)
    
    open(9,file="PSOvstore.txt")
    do ct=1,nvar
        read(9,*) vstore(ct,:)
    end do !ct
    close(9)
    
    open(9,file="PSOfxstore.txt")
    do ct=1,npart
        read(9,*) fxstore(ct)
    end do !ct
    close(9)
    
    open(9,file="PSOpbeststore.txt")
    do ct=1,nvar
        read(9,*) pbeststore(ct,:)
    end do !ct
    close(9)

    open(9,file="PSOfpbeststore.txt")
    do ct=1,npart
        read(9,*) fpbeststore(ct)
    end do !ct
    close(9)

    open(9,file="PSOglobbest.txt")
    do ct=1,nvar
        read(9,*) globbest(ct)
    end do !ct
    close(9)

    open(9,file="PSOfglobbest.txt")
    read(9,*) fglobbest
    close(9)
    
    
    randct = oldrandct
    iter = olditer
    
    !finish the remainder of the last iteration, if there's any to be done
    if (oldpartct<npart) then
    do partct = oldpartct+1,npart
        
        
        !!!UPDATE VELOCITIES AND POSITIONS
       
            !call random_number(pbestshocks)
            pbestshocks = randshocks((randct+1):(randct+nvar))
            randct = randct + nvar
            
            !call random_number(globbestshocks)
            globbestshocks = randshocks((randct+1):(randct+nvar))
            randct = randct + nvar
            
            
            vstore(:,partct) = vstore(:,partct) + phi(1) * pbestshocks * (pbeststore(:,partct)-xstore(:,partct)) &
                + phi(2)*globbestshocks * (globbest - xstore(:,partct))
            vstore(:,partct) = chi * vstore(:,partct)
            xstore(:,partct) = xstore(:,partct) + vstore(:,partct)
            do varct=1,nvar
                xstore(varct,partct) = min(max(lb(varct),xstore(varct,partct)),ub(varct))
            end do !varct
            fxstore(partct) = f(xstore(:,partct))
               
        !record progress so far
        open(9,file="PSOprogress.txt")
        write(9,*) iter
        write(9,*) partct
        write(9,*) randct
        close(9)
        
        open(9,file="PSOxstore.txt")
        do ct=1,nvar
            write(9,*) xstore(ct,:)
        end do !ct
        close(9)
        
        open(9,file="PSOvstore.txt")
        do ct=1,nvar
            write(9,*) vstore(ct,:)
        end do !ct
        close(9)
        
        open(9,file="PSOfxstore.txt")
        do ct=1,npart
            write(9,*) fxstore(ct)
        end do !ct
        close(9)
               

        
    end do !partct
    end if  !oldpartct<npart
    
    
    
    !!!NOW, UPDATE PERSONAL BEST CONDITIONS & GLOBAL BEST CONDITION
    do partct=1,npart
        if (fxstore(partct)<fpbeststore(partct)) then
            pbeststore(:,partct) = xstore(:,partct)
            fpbeststore(partct) = fxstore(partct)
        end if
        if (fxstore(partct)<fglobbest) then
            globbest = xstore(:,partct)
            fglobbest = fxstore(partct)
        end if
    end do !partct
    
    !store the results of personal and global best calculations
    open(9,file="PSOpbeststore.txt")
    do ct=1,nvar
        write(9,*) pbeststore(ct,:)
    end do !ct
    close(9)

    open(9,file="PSOfpbeststore.txt")
    do ct=1,npart
        write(9,*) fpbeststore(ct)
    end do !ct
    close(9)

    open(9,file="PSOglobbest.txt")
    do ct=1,nvar
        write(9,*) globbest(ct)
    end do !ct
    close(9)

    open(9,file="PSOfglobbest.txt")
    write(9,*) fglobbest
    close(9)
            
    
    
    stiter = olditer + 1
end if !olditer

end if !newrun

!at this point, you should have initialized everything, and you are ready to go from "stiter" to "maxit" as usual


!!!NOW, LOOP OVER THE ITERATIONS
do iter=stiter,maxit
    
        do partct=1,npart
        !!!UPDATE VELOCITIES AND POSITIONS
       
            !call random_number(pbestshocks)
            pbestshocks = randshocks((randct+1):(randct+nvar))
            randct = randct + nvar
            
            !call random_number(globbestshocks)
            globbestshocks = randshocks((randct+1):(randct+nvar))
            randct = randct + nvar
            
            
            vstore(:,partct) = vstore(:,partct) + phi(1) * pbestshocks * (pbeststore(:,partct)-xstore(:,partct)) &
                + phi(2)*globbestshocks * (globbest - xstore(:,partct))
            vstore(:,partct) = chi * vstore(:,partct)
            xstore(:,partct) = xstore(:,partct) + vstore(:,partct)
            do varct=1,nvar
                xstore(varct,partct) = min(max(lb(varct),xstore(varct,partct)),ub(varct))
            end do !varct
            fxstore(partct) = f(xstore(:,partct))
               
        !record progress so far
        open(9,file="PSOprogress.txt")
        write(9,*) iter
        write(9,*) partct
        write(9,*) randct
        close(9)
        
        open(9,file="PSOxstore.txt")
        do ct=1,nvar
            write(9,*) xstore(ct,:)
        end do !ct
        close(9)
        
        open(9,file="PSOvstore.txt")
        do ct=1,nvar
            write(9,*) vstore(ct,:)
        end do !ct
        close(9)
        
        open(9,file="PSOfxstore.txt")
        do ct=1,npart
            write(9,*) fxstore(ct)
        end do !ct
        close(9)
        
        
        
        
        end do !partct
    
    !!!UPDATE PERSONAL BEST CONDITIONS & GLOBAL BEST CONDITION
    do partct=1,npart
        if (fxstore(partct)<fpbeststore(partct)) then
            pbeststore(:,partct) = xstore(:,partct)
            fpbeststore(partct) = fxstore(partct)
        end if
        if (fxstore(partct)<fglobbest) then
            globbest = xstore(:,partct)
            fglobbest = fxstore(partct)
        end if
    end do !partct
    
    
    !store the results of personal and global best calculations
    open(9,file="PSOpbeststore.txt")
    do ct=1,nvar
        write(9,*) pbeststore(ct,:)
    end do !ct
    close(9)

    open(9,file="PSOfpbeststore.txt")
    do ct=1,npart
        write(9,*) fpbeststore(ct)
    end do !ct
    close(9)

    open(9,file="PSOglobbest.txt")
    do ct=1,nvar
        write(9,*) globbest(ct)
    end do !ct
    close(9)

    open(9,file="PSOfglobbest.txt")
    write(9,*) fglobbest
    close(9)
    
    
    !!!!COMPUTE CONVERGENCE CRITERIA
    !sort to find the best points
    fpbestsort=fpbeststore
    call heapsort(npart,fpbestsort)
    quickthresh = fpbestsort(xquicknum)
    
    
    xnorm = 0.0
    reltol = 0.0
    xquicknorm = 0.0
    
    do partct=1,npart
        !max distance
        xnorm = max(xnorm,sqrt(sum((globbest-xstore(:,partct))**2.0)))
        
        !max difference between personal and global bests
        reltol = max(reltol,abs(fglobbest-fpbeststore(partct))/&
    (abs(fglobbest)+abs(fpbeststore(partct))+0.01))
    
        !quick max distance measure
        if (fpbeststore(partct)<quickthresh) then
            xquicknorm = max(xquicknorm,sqrt(sum((globbest-xstore(:,partct))**2.0)))
        end if
    end do !partct
    
    if (xnorm<xtol.or.reltol<ftol.or.xquicknorm<xquicktol) exit
    
end do !iter
    
x = globbest
fx = fglobbest
    
end subroutine psorestart


subroutine psoparam(x,fx,f,lb,ub,nvar,npart,xtol,xquicktol,xquicknum,ftol,maxit,param,nparam,phi)
implicit none

!input/output declarations
integer :: nvar,npart,maxit,nparam,xquicknum
double precision :: x(nvar),f,lb(nvar),ub(nvar),xtol,param(nparam),phi(2),fx,ftol,&
    xquicktol
external :: f

!local declarations
integer :: iter,seedint,varct,partct,n
integer, allocatable :: seed(:)
double precision :: xnorm,chi,phisum,xstore(nvar,npart),fxstore(npart),&
    pbeststore(nvar,npart),fpbeststore(npart),globbest(nvar),fglobbest,vstore(nvar,npart),&
    pbestshocks(nvar),globbestshocks(nvar),spacenorm,reltol,xquicknorm,quickthresh,fpbestsort(npart)

seedint=2501
call random_seed(size=n)
allocate(seed(n))
do iter=1,n
    seed(iter) = seedint+iter
end do !iter
call random_seed(put=seed)

!SET CONSTRAINT PARAMETER
!note that sum(phi) > 4 is required, with phi(1)=phi(2)=2.05 most common in lit
phisum = sum(phi)
chi = 2.0 / ( phisum - 2.0 + sqrt( ( phisum ** 2.0 ) - 4.0 * phisum )  )

!SET THE WIDTH OF THE SPACE TO INITIALIZE VELOCITIES
spacenorm = sqrt(sum((ub-lb)**2.0))


!!!INITIALIZE THE SWARM 
fglobbest = 1.0e30
do partct=1,npart
    call random_number(xstore(:,partct))
    xstore(:,partct) = lb * (1.0-xstore(:,partct)) + ub * xstore(:,partct)
    call random_number(vstore(:,partct))
    vstore(:,partct) = -1.0*spacenorm * (1.0 - vstore(:,partct)) + spacenorm * vstore(:,partct)
    fxstore(partct) = f(xstore(:,partct),param,nparam)
    
    if (fxstore(partct)<fglobbest) then
        globbest = xstore(:,partct)
        fglobbest = fxstore(partct)
    end if
    
end do !partct

pbeststore = xstore
fpbeststore = fxstore

!!!NOW, LOOP OVER THE ITERATIONS
do iter=1,maxit
    
    !!!UPDATE VELOCITIES AND POSITIONS
    do partct=1,npart
        call random_number(pbestshocks)
        call random_number(globbestshocks)
        vstore(:,partct) = vstore(:,partct) + phi(1) * pbestshocks * (pbeststore(:,partct)-xstore(:,partct)) &
            + phi(2)*globbestshocks * (globbest - xstore(:,partct))
        vstore(:,partct) = chi * vstore(:,partct)
        xstore(:,partct) = xstore(:,partct) + vstore(:,partct)
        do varct=1,nvar
            xstore(varct,partct) = min(max(lb(varct),xstore(varct,partct)),ub(varct))
        end do !varct
        fxstore(partct) = f(xstore(:,partct),param,nparam)
    end do !parct
    
    !!!UPDATE PERSONAL BEST CONDITIONS & GLOBAL BEST CONDITION
    do partct=1,npart
        if (fxstore(partct)<fpbeststore(partct)) then
            pbeststore(:,partct) = xstore(:,partct)
            fpbeststore(partct) = fxstore(partct)
        end if
        if (fxstore(partct)<fglobbest) then
            globbest = xstore(:,partct)
            fglobbest = fxstore(partct)
        end if
    end do !partct
    
    !!!!COMPUTE CONVERGENCE CRITERIA
    !sort to find the best points
    fpbestsort=fpbeststore
    call heapsort(npart,fpbestsort)
    quickthresh = fpbestsort(xquicknum)
    
    
    xnorm = 0.0
    reltol = 0.0
    xquicknorm = 0.0
    
    do partct=1,npart
        !max distance
        xnorm = max(xnorm,sqrt(sum((globbest-xstore(:,partct))**2.0)))
        
        !max difference between personal and global bests
        reltol = max(reltol,abs(fglobbest-fpbeststore(partct))/&
    (abs(fglobbest)+abs(fpbeststore(partct))+0.01))
    
        !quick max distance measure
        if (fpbeststore(partct)<quickthresh) then
            xquicknorm = max(xquicknorm,sqrt(sum((globbest-xstore(:,partct))**2.0)))
        end if
    end do !partct
    
    if (xnorm<xtol.or.reltol<ftol.or.xquicknorm<xquicktol) exit
    
end do !iter

x = globbest
fx = fglobbest

return
end subroutine psoparam

subroutine neldermead2d(x,fx,f,lb,ub,ftol,xtol,maxit,param,nparam)
    implicit none
    
    !input/output declarations
    integer :: nparam,maxit
    double precision :: x(2),f,lb(2),ub(2),param(nparam),ftol,xtol,fx
    external :: f
    
    !note that f is a function taking as arguments x,param,nparam, where x is 2x1 and param(nparam) is parameter vector
    
    !local declarations
    integer :: iter
    double precision :: x1(2),x2(2),x3(2),xB(2),xG(2),xW(2),xR(2),xE(2),&
        xC(2),xS(2),xM(2),f1,f2,f3,fB,fG,fW,fR,fE,fC,fS,fM,&
        xnorm,reltol
        
    
    !!!SET UP INITIAL SIMPLEX X1,X2,X3
    
    x1 = x; !initial guess is one point in simplex
    x2 = ( lb + ub ) / 2.0; !another point is mean of bounded range
    x3 = ( x + lb ) / 2.0; !another point is mean of guess and lb

    
    f1 = f(x1,param,nparam)
    f2 = f(x2,param,nparam)
    f3 = f(x3,param,nparam)
    
    !!!DETERMINE BEST,GOOD,WORSE
    if (f1<=f2.and.f1<=f3) then
        xB=x1; fB=f1;
        if (f2<=f3) then 
            xG=x2; fG=f2;
            xW=x3; fW=f3;
        else 
            xG=x3; fG=f3;
            xW=x2; fW=f2;
        end if
    else if (f2<=f1.and.f2<=f3) then
        xB=x2; fB=f2;
        if (f1<=f3) then
            xG=x1; fG=f1;
            xW=x3; fW=f3;
        else
            xG=x3; fG=f3;
            xW=x1; fW=f1;
        end if
    else if (f3<=f1.and.f3<=f2) then
        xB=x3; fB=f3;
        if (f1<=f2) then
            xG=x1; fG=f1;
            xW=x2; fW=f2;
        else
            xG=x2; fG=f2;
            xW=x1; fW=f1;
        end if
    end if 
    
    !!!ITERATIVELY UPDATE THESE LOCATIONS
    do iter = 1,maxit
        
        !find midpoint, note that this is always within bounds
        xM = (xB+xG)/2.0; fM = f(xM,param,nparam);
        
        !find reflection point, adjust for bounds
        xR = 2.0 * xM - xW;
        xR(1) = min(max(lb(1),xR(1)),ub(1))
        xR(2) = min(max(lb(2),xR(2)),ub(2))
        fR = f(xR,param,nparam)
        
        !Case I in Kopecky notation
        if (fR<fG) then
            
            if (fB<fR) then 
                xW=xR; fW=fR;
            else 
                xE = 2*xR - xM;
                xE(1) = min(max(lb(1),xE(1)),ub(1))
                xE(2) = min(max(lb(2),xE(2)),ub(2))
                fE = f(xE,param,nparam)
                if (fE<fB) then
                    xW=xE; fW=fE;
                else 
                    xW=xR; fW=fR;
                end if
            end if
        
        !Case II in Kopecky notation
        else if (fR>=fG) then
            
            if (fR<fW) then
                xW=xR; fW=fR;
            else 
                !note that I'm ignoring C1/C2 distinction, and C1 is always in bounds
                xC = (xM+xW)/2.0
                fC = f(xC,param,nparam)
                if (fC<fW) then
                    xW=xC; fW=fC;
                else 
                    !note that S is always in bounds
                    xS = (xB+xW)/2.0;
                    fS = f(xS,param,nparam)
                    xW=xS; fW=fS;
                    xG=xM; fG=fM;
                end if
                
            end if
        
        
        end if
        
        reltol = abs(fB-fW)/(abs(fB)+abs(fW)+0.01)
        xnorm = sqrt(sum((xB-xW)**2.0)) + &
            sqrt(sum((xB-xG)**2.0)) + &
            sqrt(sum((xW-xG)**2.0))
        
        if (reltol<ftol.or.xnorm<xtol) exit
        
        !note, if no convergence has been obtained, do the reordering again
        x1=xB; f1=fB; x2=xG; f2=fG; x3=xW; f3=fW;
        if (f1<=f2.and.f1<=f3) then
            xB=x1; fB=f1;
            if (f2<=f3) then 
                xG=x2; fG=f2;
                xW=x3; fW=f3;
            else 
                xG=x3; fG=f3;
                xW=x2; fW=f2;
            end if
        else if (f2<=f1.and.f2<=f3) then
            xB=x2; fB=f2;
            if (f1<=f3) then
                xG=x1; fG=f1;
                xW=x3; fW=f3;
            else
                xG=x3; fG=f3;
                xW=x1; fW=f1;
            end if
        else if (f3<=f1.and.f3<=f2) then
            xB=x3; fB=f3;
            if (f1<=f2) then
                xG=x1; fG=f1;
                xW=x2; fW=f2;
            else
                xG=x2; fG=f2;
                xW=x1; fW=f1;
            end if
        end if 
    
    end do !iter
    
    x = xB; !return the best guess so far to the "x" variable
    fx = fB; !return best function value so far to the "fx" variable
end subroutine neldermead2d

subroutine linspace(z,x,y,n)
    implicit none
    
    !n = the dimension of the output vector
    !z = the n x 1 output vector, with equally spaced points between x and y
    !x = the minimum of the linear grid
    !y = the maximum of the linear grid
    
    
    !input/output declarations
    integer :: n
    double precision :: z(n),x,y
    
    !local declarations
    integer :: i
    double precision :: d
    
    d = (y-x)/dble(n-1)
    z(1) = x
    
    do i = 2,n-1
        z(i) = z(i-1) + d
    end do
    
    z(n) = y
    
    return

end subroutine linspace

subroutine hunt(xx,n,x,jlo)
    implicit none
    
    !xx = the n x 1 table of values that you're comparing x to
    !n = the dimension of xx
    !x = the value which you're interested in
    !jlo = (on input) the guess for the integer such that x is in between xx(jlo) and xx(jlo+1)
    !jlo = (on output) the value of the integer such that x is in between xx(jlo) and xx(jlo+1)
    
    !input/output declarations
    integer :: jlo,n
    double precision  :: x,xx(n)
    
    !local declarations
    integer :: inc,jhi,jm
    logical :: ascnd
    
    !determine if table is ascending
    ascnd = xx(n).ge.xx(1)
    
    !in case input guess isn't useful, for robustness
    if (jlo.le.0.or.jlo.gt.n) then 
        jlo = 0
        jhi = n+1
        goto 3
    endif
    
    inc=1 !initialize the hunting increment
    
    !hunt up
    if (x.ge.xx(jlo).eqv.ascnd) then
1       jhi = jlo+inc
        if (jhi.gt.n) then
            jhi = n+1
        else if (x.ge.xx(jhi).eqv.ascnd) then
            jlo = jhi
            inc = inc+inc
            goto 1
        end if    
    !hunt down        
    else
        jhi = jlo
2       jlo = jhi - inc
        if (jlo.lt.1) then
            jlo = 0
        else if (x.lt.xx(jlo).eqv.ascnd) then
            jhi = jlo
            inc = inc+inc
            goto 2
        end if
        
        
        
    endif
    
    !now, hunt is done, begin the bisection phase
3   if (jhi-jlo.eq.1) then
        if (x.eq.xx(n)) jlo=n-1
        if (x.eq.xx(1)) jlo=1
        return
    end if
    jm = (jhi + jlo)/2
    if (x.ge.xx(jm).eqv.ascnd) then
        jlo = jm
    else 
        jhi = jm
    end if
    goto 3
    
end subroutine hunt

subroutine linterp(xx,yy,n,x,xindguess,y)
    implicit none
    
    !input/output declarations
    integer :: n,xindguess
    double precision :: xx(n),yy(n),x,y
    
    !local declarations
    integer :: xindactual
    double precision :: d
    
    
    call hunt(xx,n,x,xindguess)
    
    xindactual = xindguess
    
    d = (x - xx(xindactual)) / (xx(xindactual + 1 ) - xx(xindactual))
    
    y = yy(xindactual) + d * ( yy(xindactual+1) - yy(xindactual)  )
    
    return

end subroutine linterp

subroutine linterp2(xx1,xx2,yy,n,x,xindguess,y)
    implicit none
    
    !input/output declarations
    integer :: n(2),xindguess(2)
    double precision :: xx1(n(1)),xx2(n(2)),yy(n(1),n(2)),x(2),y
    
    !local declarations
    integer :: xindactual(2)
    double precision :: d(2),yinterp(2)
    
    
    call hunt(xx1,n(1),x(1),xindguess(1))
    call hunt(xx2,n(2),x(2),xindguess(2))
    
    xindactual = xindguess
    
    d(1) = (x(1) - xx1(xindactual(1))) / (xx1(xindactual(1) + 1 ) - xx1(xindactual(1)))
    d(2) = (x(2) - xx2(xindactual(2))) / (xx2(xindactual(2) + 1 ) - xx2(xindactual(2)))
    
    yinterp(1) = yy(xindactual(1),xindactual(2)) + &
        d(1) * ( yy(xindactual(1)+1,xindactual(2)) - yy(xindactual(1),xindactual(2)) )
    yinterp(2) = yy(xindactual(1),xindactual(2)+1) + &
        d(1) * ( yy(xindactual(1)+1,xindactual(2)+1) - yy(xindactual(1),xindactual(2)+1) )
        
    y = yinterp(1) + d(2) * (yinterp(2) - yinterp(1) )
    
    return

end subroutine linterp2

subroutine linterp3(xx1,xx2,xx3,yy,n,x,xindguess,y)
    implicit none
    
    !input/output declarations
    integer :: n(3),xindguess(3)
    double precision :: xx1(n(1)),xx2(n(2)),xx3(n(3)),yy(n(1),n(2),n(3)),x(3),y
    
    !local declarations
    integer :: xind(3)
    double precision :: d(3)
    
    call hunt(xx1,n(1),x(1),xindguess(1))
    call hunt(xx2,n(2),x(2),xindguess(2))
    call hunt(xx3,n(3),x(3),xindguess(3))
    
    xind = xindguess
    
    d(1) = (x(1) - xx1(xind(1))) / (xx1(xind(1) + 1 ) - xx1(xind(1)))
    d(2) = (x(2) - xx2(xind(2))) / (xx2(xind(2) + 1 ) - xx2(xind(2)))
    d(3) = (x(3) - xx3(xind(3))) / (xx3(xind(3) + 1 ) - xx3(xind(3)))
    
    y = yy(xind(1),xind(2),xind(3)) * (1.0 - d(1)) * (1.0 - d(2)) * (1.0 - d(3))
    
    y = y + yy(xind(1)+1,xind(2),xind(3)) * d(1) * (1.0 - d(2)) * (1.0 - d(3))
    
    y = y + yy(xind(1),xind(2)+1,xind(3)) * (1.0 - d(1)) * d(2) * (1.0 - d(3))
    
    y = y + yy(xind(1),xind(2),xind(3)+1) * (1.0 - d(1)) * (1.0 - d(2)) * d(3)
    
    y = y + yy(xind(1)+1,xind(2)+1,xind(3)) * d(1) * d(2) * (1.0 - d(3))
    
    y = y + yy(xind(1),xind(2)+1,xind(3)+1) * (1.0 - d(1)) * d(2) * d(3)
    
    y = y + yy(xind(1)+1,xind(2),xind(3)+1) * d(1) * (1.0 - d(2)) * d(3)
    
    y = y + yy(xind(1)+1,xind(2)+1,xind(3)+1) * d(1) * d(2) * d(3)
    
end subroutine linterp3

subroutine spline(x,y,n,yp1,ypn,y2)
    implicit none
    
    !input/output declarations
    integer :: n
    double precision :: yp1,ypn,x(n),y(n),y2(n)
    
    integer :: i,k
    double precision :: p,qn,sig,un,u(n)
    
    if (yp1.gt.0.99e30) then
        y2(1) = 0.0
        u(1) = 0.0
    else
        y2(1) = -0.5
        u(1) = (3.0/(x(2) - x(1))) * ( (y(2) - y(1)) / (x(2) - x(1)) - yp1)
    end if
    
    do i=2,n-1
        sig = (x(i) - x(i-1))/(x(i+1)-x(i-1))
        p = sig * y2(i-1) + 2.0
        y2(i)=(sig-1.0)/p
        u(i) = (6.0 * ((y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))&
            /(x(i)-x(i-1)))/(x(i+1)-x(i-1)) - sig * u(i-1))/p
    end do !i
    
    if (ypn.gt.0.99e30) then
        qn=0.0
        un=0.0
    else
        qn = 0.5
        un = (3.0/(x(n)-x(n-1))) * (ypn - (y(n)-y(n-1))/(x(n)-x(n-1)))
    end if
    
    y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0)
    do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
    end do !k

    return
end subroutine spline

subroutine splint(xa,ya,y2a,n,x,y)
    implicit none
    
    !input/output declarations
    integer :: n
    double precision :: x,y,xa(n),y2a(n),ya(n)
    
    integer :: k,khi,klo
    double precision :: a,b,h
    
    klo = 1
    khi = n
1    if (khi-klo.gt.1) then
        k = (khi + klo)/2
        if (xa(k).gt.x) then
            khi=k
        else
            klo=k
        end if
        goto 1
    end if
    h = xa(khi)-xa(klo)
    if (h.eq.0) write(*,*) 'bad xa input in splint'
    a = (xa(khi)-x)/h
    b = (x - xa(klo))/h
    y = a * ya(klo) + b*ya(khi) + ((a ** 3.0 - a)*y2a(klo)+(b**3.0-b)*y2a(khi))*(h**2.0)/6.0

end subroutine splint

subroutine splie2(x1a,x2a,ya,m,n,y2a)
    implicit none
    
    !note x1a is included only for consistency with splin2
    
    !input/output declarations
    integer :: m,n,NN
    double precision :: x1a(m),x2a(n),y2a(m,n),ya(m,n)
    parameter(NN=2000); !max val of m or n allowed
    
    integer :: j,k
    double precision :: y2tmp(NN),ytmp(NN)
    
    do j=1,m
        do k=1,n
            ytmp(k) = ya(j,k)
        end do !k
        call spline(x2a,ytmp,n,dble(1.0e30),dble(1.0e30),y2tmp)
        do k=1,n
            y2a(j,k)=y2tmp(k)
        end do !k
    end do !j
    return
end subroutine splie2

subroutine splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
    implicit none
    
    !input/output declarations
    integer :: m,n,NN
    double precision :: x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
    parameter(NN=2000) !max val of m or n allowed
    
    
    integer :: j,k
    double precision :: y2tmp(NN),ytmp(NN),yytmp(NN)


    do j=1,m
        do k=1,n
            ytmp(k) = ya(j,k)
            y2tmp(k) = y2a(j,k)
        end do !k
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
    end do !j
    call spline(x1a,yytmp,m,dble(1.0e30),dble(1.0e30),y2tmp)
    call splint(x1a,yytmp,y2tmp,m,x1,y)
    
    return
end subroutine splin2

double precision function brent(ax,bx,cx,f,tol,xmin)
    implicit none
    
    !input/output declarations
    double precision :: ax,bx,cx,tol,xmin,f
    
    !parameter declarations
    integer :: ITMAX
    double precision :: CGOLD,ZEPS
    parameter(ITMAX=1000,CGOLD=0.3819660,ZEPS=1.0e-10)
    
    external f
    
    !other declarations
    integer :: iter
    double precision :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    
    !a and b are in ascending order
    a=min(ax,cx)
    b=max(ax,cx)
    
    !initializations
    v=bx
    w=v
    x=v
    e=0.0
    
    fx=f(x)
    fv=fx
    fw=fx
    
    do iter = 1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.0*tol1
        if(abs(x-xm).le.(tol2-0.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0*(q-r)
            if(q.gt.0.0) p = -p
            q=abs(q)
            etemp=e
            e=d
            if(abs(p).ge.abs(0.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
        end if
1       if (x.ge.xm) then
            e=a-x
        else
            e=b-x
        end if
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
            u=x+d
        else
            u=x+sign(tol1,d)
        end if
        fu=f(u)
        if(fu.le.fx) then
            if(u.ge.x) then
                a=x
            else
                b=x
            end if
            v=2
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
        else
            if(u.lt.x) then
                a=u
            else
                b=u
            end if
            if(fu.le.fw .or. w.eq.x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
                v=u
                fv=fu
            end if
        end if
    end do  !iter
    write(*,*) "Brent exceeded maximum function iterations"
3   xmin=x
    brent=fx   
    
    return
end function brent

double precision function brentparam(ax,bx,cx,f,tol,xmin,param,nparam)
implicit none
!input/output declarations
    integer :: nparam
    double precision :: ax,bx,cx,tol,xmin,f,param(nparam)
    
    external f
    
    !parameter declarations
    integer :: ITMAX
    double precision :: CGOLD,ZEPS
    parameter(ITMAX=1000,CGOLD=0.3819660,ZEPS=1.0e-10)
        
    !other declarations
    integer :: iter
    double precision :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    
    !a and b are in ascending order
    a=min(ax,cx)
    b=max(ax,cx)
    
    !initializations
    v=bx
    w=v
    x=v
    e=0.0
    
    fx = f(x,param,nparam)
    
    fv=fx
    fw=fx
    
    do iter = 1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.0*tol1
        if(abs(x-xm).le.(tol2-0.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0*(q-r)
            if(q.gt.0.0) p = -p
            q=abs(q)
            etemp=e
            e=d
            if(abs(p).ge.abs(0.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) goto 1
            d=p/q
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
        end if
1       if (x.ge.xm) then
            e=a-x
        else
            e=b-x
        end if
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
            u=x+d
        else
            u=x+sign(tol1,d)
        end if
        fu=f(u,param,nparam)
        if(fu.le.fx) then
            if(u.ge.x) then
                a=x
            else
                b=x
            end if
            v=2
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
        else
            if(u.lt.x) then
                a=u
            else
                b=u
            end if
            if(fu.le.fw .or. w.eq.x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
                v=u
                fv=fu
            end if
        end if
    end do  !iter
    write(*,*) "Brent exceeded maximum function iterations"
3   xmin=x
    brentparam=fx   
    
    return



end function brentparam

end module base_lib
