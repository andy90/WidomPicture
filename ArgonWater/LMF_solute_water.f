! this file is used to calculate the LMFforce 
! by brute force

      program forceLMF
      implicit none
      real(8),dimension(:),allocatable :: r,rho,Fr,ur,FrBA,urBA
      integer :: nbinr
      real(8) :: stepr,rhobulk
      real(8) :: theta,deltheta
      integer :: NTHETA
      real(8) :: PI
      real(8) :: Fx,u1,FxBA,u1BA
      PARAMETER(PI=3.141592653)

      real(8) :: junk,gr,factor
      integer :: i,j,k
      integer :: NBINgro,IBIN
      real(8) :: delrgro,rrr
      real(8) :: ri,rj

! the total bin number of r,rhi and Fr
      nbinr=770
! theta bin size
      deltheta=0.0010d0*PI
! theta bin number
      NTHETA=INT(PI/deltheta)


      allocate(r(nbinr))
      allocate(rho(nbinr))
      allocate(Fr(nbinr))
      allocate(ur(nbinr))
      allocate(FrBA(nbinr))
      allocate(urBA(nbinr))

      do i=1,nbinr
        Fr(i)=0.0d0
        ur(i)=0.0d0
        FrBA(i)=0.0d0
        urBA(i)=0.0d0
      enddo

      open(unit=8,file='rdf_OA.xvg')
      do i=1,nbinr
        read(8,*) r(i),rho(i)
      enddo
      close(8)

      factor=0.0d0
      do i=nbinr-50,nbinr
        factor=factor+rho(i)/51
      enddo
      write(*,*) factor

      do i=1,nbinr
        rho(i)=rho(i)/factor
      enddo

      rhobulk=1000/2.96635**3
      do i=1,nbinr
        rho(i)=rho(i)*rhobulk
      enddo

! the size of the bin
      stepr=r(2)-r(1)
      

      do i=1,nbinr
        ri=r(i)+0.5*stepr
        do j=1,nbinr
          rj=r(j)+0.5*stepr
          do k=1,NTHETA
            theta=(k-0.50d0)*deltheta
            Fr(i)=Fr(i)+Fx(rj,ri,theta)*2*PI*rj**2*stepr
     X       *(rho(j)-rhobulk)*sin(theta)*deltheta
          enddo
        enddo
      enddo

      do i=1,nbinr
        do j=i,nbinr
          ur(i)=ur(i)+Fr(j)*stepr
        enddo
      enddo


      delrgro=0.005
      NBINgro=INT(6/delrgro)+2



      open(unit=9,file='ulmf_ArO.txt')
      do i=1,NBINgro
        rrr=(i-1)*delrgro
        IBIN=int(rrr/stepr)+1
        if (IBIN.le.nbinr) then
          write(9,"(2f10.4)") rrr,ur(IBIN)
        else
          write(9,"(2f10.4)") rrr,0.0
        endif
      enddo
      close(9)


      end program

! Fx is the AA interaction
      function Fx(RR,r,theta)
      implicit none
      real(8) :: RR,r,theta,Fx
      
      real(8) :: dis,Fr
      real(8) :: eps,sigma,rwca,rcut
      
      eps=sqrt(0.997740d0*0.6501940d0)
      sigma=(0.340d0+0.3165570d0)/2.0d0
      rwca=sigma*(2.0d0**(1.0d0/6.0d0))
      rcut=1.0d0

      dis=sqrt(RR**2+r**2-2.0d0*RR*r*cos(theta))
      
      if (dis.le.rwca) then
        Fx=0.0d0
      elseif(dis.le.rcut) then
        Fr=24.0d0*eps/dis*(sigma/dis)**6
     X      *(2.0d0*(sigma/dis)**6-1.0d0)
        Fx=-Fr*(RR*cos(theta)-r)/dis
      else
        Fx=0.0d0
      endif
      
      return
      end function

