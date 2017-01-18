!/opt/local/bin/f2py-3.4 -c -m tcorr tcorr.f
!/opt/local/bin/f2py-3.4 tcorr.f -m tcorr -h tcorr.pyf
! Ahsan Zeb
! This subroutine takes the hamiltonian, initial state and other parameters, evolve the state and compute the two time correlation that is given as output.
! hope it would be faster than numpy

!input: double/int: H(:,:), psi0(:), dt, ntmax, kappa, gamma, n1, ntot
!output: double complex: correlation(:)
      subroutine tcorr(n1,ntot,nnz,nrp,prntstep,Val,Col,RowPtr,
     .        psi0,dt,ntmax,kappa2,gamma2,corr)
      implicit none
      integer :: n1,ntot,nnz,nrp,ntmax, prntstep
      double precision :: dt, kappa2, gamma2
      double precision, dimension(nnz):: Val
      integer, dimension(nnz):: Col
      integer, dimension(nrp):: RowPtr
      double complex, dimension(ntot):: psi0
      double complex, dimension(ntmax):: corr
Cf2py intent(in) n1,ntot,nnz,nrp,ntmax, prntstep, dt, kappa2, gamma2
Cf2py intent(in) Val, Col, RowPtr,psi0
Cf2py depend(nnz) Val, Col
Cf2py depend(nrp) RowPtr
Cf2py depend(ntot) psi0
Cf2py intent(out) corr
Cf2py depend(ntmax) corr

      ! aux
      double complex, dimension(ntot) :: psit,k1,k2,k3
      integer :: i
      logical :: dkappa
      double precision :: dth,dt6
      double complex :: iotam
      dth = dt/2.0d0;
      dt6 = dt/6.0d0;
      iotam = (0.0d0,-1.0d0);

      !write(6,*)'nnz, ntot = ',nnz,ntot
      dkappa = .false.
      if (abs(kappa2-gamma2) > 1e-6) dkappa = .true.

      psit(:) = psi0(:);
      corr(1) = (1.0d0,0.0d0) !complex(1.0d0,0.0d0);

      do i=2,ntmax
       if (mod(i,prntstep) == 0)
     .   write(6,'(a20,2x,i10,a10,i10)') ' time evolution step ',i,
     .   ' out of ',ntmax
       call rk4step(psit)
       corr(i) = DOT_PRODUCT(psi0, psit)
      end do
      contains
      !----------------------------------------
      subroutine yprime(vin,vout)
      implicit none
      double complex, intent(in) :: vin(ntot)
      double complex, intent(out) :: vout(ntot)
      integer :: i,j	
      vout = 0.0d0

      do i=1,nrp-1
       do j=RowPtr(i),RowPtr(i+1)-1
        vout(i) = vout(i) + Val(j)*vin(Col(j));
       end do
      end do
      if (dkappa .eqv. .true.) then
       k1(1:n1) = kappa2*vin(1:n1);
       k1(n1+1:ntot) = gamma2*vin(n1+1:ntot);
       vout = vout*iotam - k1
      else
       vout = vout*iotam - kappa2*vin
      endif
      end subroutine yprime
      !----------------------------------------
      subroutine rk4step(y)
      implicit none
      double complex, intent(inout) :: y(ntot)
      call yprime(y,k1)
      call yprime(y+k1*dth,k2)
      k1 = k1 + 2*k2
      call yprime(y+k2*dth,k3)
      k1 = k1 + 2*k3
      call yprime(y+k3*dt,k2)
      k1 = k1 + k2
      y = y + k1*dt6
      end subroutine rk4step
      !----------------------------------------
      end subroutine tcorr









      
