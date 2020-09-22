!! :: ==================================================
!! :: 層流境界層を求めるプログラム
!! :: title :: main_BoundaryLayer.f90
!! :: author :: Hirotoshi Uji
!! :: date :: 2020.02.28 Fri
!! :: latest :: 2020.09.22 Tue
!! :: ==================================================
program main 
	use mod_nonDimensionalNumber
	use mod_PhysicalProp_2d
	implicit none
	!! :: parameter --- --- --- --- --- --- --- --- --- --- --- ---
	integer,parameter :: Nt = 100
	integer,parameter :: Neta = 1000
	real(8),parameter :: ERR_FUNC = 1.d-10
	real(8),parameter :: Re_x = 10000d0
	real(8),parameter :: inf_cov = 0.99999d0
	real(8),parameter :: sup_cov = 1.00001d0
	!! :: variables --- --- --- --- --- --- --- --- --- --- --- ---
	integer :: i = 0 , ip1 = 0 , j = 0 , ctt = 0 , ch_conver = 0
	real(8) :: dt , deta_BL
	real(8),dimension(0:Neta) :: eta, velo_x, velo_y, Temp, dens
	real(8),dimension(1:4,0:2,0:Neta) :: f_k=0d0 , g_k=0d0
	real(8),dimension(0:3,0:2,0:Neta) :: funcF=1d0 , funcG=1d0
	real(8),dimension(0:2,0:Neta) :: k_f=0d0 , k_g=0d0
	!! :: treatment --- --- --- --- --- --- --- --- --- --- --- ---
	dt = 1.0d-5
	deta_BL = 1d-2
	gamma = 1.4d0
	Pr = 0.71d0
	Ma = 0.3d0
	!! --- --- --- --- --- --- --- --- --- --- --- ---
	!! :: shooting method
	!! :: f の二階の導関数の初期値（eta=0）を仮定して，
	!! :: f の一階の導関数の境界条件（at eta=inf : f'=1）を満たすものを探す．
	!! --- --- --- --- --- --- --- --- --- --- --- ---
	do ctt = 0,Nt
		funcF(0,:,:) = 0d0 
		!! :: Given initial values
		funcF(0,0,0) = 0d0
		funcF(0,1,0) = 0d0
		!! :: 2階微分の初期値の仮定
		funcF(0,2,0) = 0.33205d0 + dble(ctt) * 0.000001d0
		!! :: 4th-ordered Runge-Kutta loop
		do j = 0,Neta-1
			call cal_RKstep_f(1, j)
			call cal_RKstep_f(2, j)
			call cal_RKstep_f(3, j)
			call cal_RKstep_f(4, j)
			call cal_RKlastStep_f(j)
		end do
		write(6,*)ctt,funcF(0,1,Neta),funcF(0,2,0)
		if( funcF(0,1,Neta) > inf_cov .AND. funcF(0,1,Neta) < sup_cov) then
			write(6,*)'Find F !!'
			exit
		end if
	end do
	!! :: Given initial values
	funcF(0,0,0) = 0d0
	funcF(0,1,0) = 0d0
	!! --- --- --- --- --- --- --- --- --- --- --- ---
	!! :: shooting method
	!! :: g の初期値（eta=0）を仮定して，
	!! :: g の一階の導関数の境界条件（at eta=inf : f'=1）を満たすものを探す．
	!! --- --- --- --- --- --- --- --- --- --- --- ---
	do ctt = 0,Nt
		funcG(0,:,:) = 0d0
		!! :: Given initial values
		funcG(0,1,0) = 0d0
		!! :: g の初期値の仮定
		funcG(0,0,0) = 1.01514d0 + dble(ctt) * 0.000001d0
		!! :: 4th-ordered Runge-Kutta loop
		do j = 0,Neta-1
			call cal_RKstep_g(1, j)
			call cal_RKstep_g(2, j)
			call cal_RKstep_g(3, j)
			call cal_RKstep_g(4, j)
			call cal_RKlastStep_g(j)
		end do
		write(6,*)ctt,funcG(0,0,Neta),funcG(0,0,0)
		if( funcG(0,0,Neta) > inf_cov .AND. funcG(0,0,Neta) < sup_cov) then
			write(6,*)'Find G !!'
			exit
		end if
	end do
	!! --- --- --- --- --- --- --- --- --- --- --- ---
	!! :: 
	!! :: 物理量を変換する．
	!! :: 
	!! --- --- --- --- --- --- --- --- --- --- --- ---
	do j = 0,Neta
		eta(j) = deta_BL * dble(j) 
		velo_x(j) = funcF(0,1,j)
		velo_y(j) = 0.5d0 * ( eta(j) * funcF(0,1,j) - funcF(0,0,j) ) / ( sqrt(Re_x) * funcG(0,0,j) )
		Temp(j) = funcG(0,0,j)
		dens(j) = 1d0 / funcG(0,0,j)
	end do
	!! --- --- --- --- --- --- --- --- --- --- --- ---
	!! :: 
	!! :: 出力する．
	!! :: 
	!! --- --- --- --- --- --- --- --- --- --- --- ---
	open(unit=11,file="Rs_BL_Rex30000.dat",action="write",status="replace")
	do j = 0,Neta
		write(11,'(5(f10.6))')eta(j),velo_x(j),velo_y(j),Temp(j),dens(j)
	end do
	close(11)
	write(6,*)"end program"
	stop
	contains
	!!
	subroutine cal_RKlastStep_f(j)
		!! :: argument 
		integer,intent(in) :: j
		!! :: variables
		integer :: i = 0
		!! :: treatment
		do i = 0,2
			k_f(i,j) = (f_k(1,i,j) + 2d0*f_k(2,i,j) + 2d0*f_k(3,i,j) + f_k(4,i,j)) / 6d0
			funcF(0,i,j+1) = funcF(0,i,j) + k_f(i,j)*deta_BL
		end do
		!!
	end subroutine
	!!
	subroutine cal_RKstep_f(RKt, j)
		!! :: argument 
		integer,intent(in) :: RKt !! :: Runge-Kutta times
		integer,intent(in) :: j
		!! :: variables
		integer :: i = 0
		integer :: RKtm1 = 0
		real(8) :: dw = 0 !! :: deta_BL の重み delta-eta-weight
		!! :: treatment
		RKtm1 = RKt - 1
		!! :: :: RKt = 1 => dw = 0.5d0
		!! :: :: RKt = 2 => dw = 0.5d0
		!! :: :: RKt = 3 => dw = 1d0
		dw = (dble(RKt)**2d0 - 3d0*dble(RKt) + 4d0 ) * 0.25d0 
		!!
		f_k(RKt,2,j) = -0.5d0 * funcF(RKtm1,0,j) * funcF(RKtm1,2,j)
		f_k(RKt,1,j) = funcF(RKtm1,2,j)
		f_k(RKt,0,j) = funcF(RKtm1,1,j)
		!!
		if(Rkt < 4) then
			funcF(RKt,2,j) = funcF(0,2,j) + f_k(RKt,2,j)*deta_BL*dw
			funcF(RKt,1,j) = funcF(0,1,j) + f_k(RKt,1,j)*deta_BL*dw
			funcF(RKt,0,j) = funcF(0,0,j) + f_k(RKt,0,j)*deta_BL*dw
		end if
		!!
	end subroutine
	!!
	subroutine cal_RKlastStep_g(j)
		!! :: argument 
		integer,intent(in) :: j
		!! :: variables
		integer :: i = 0
		!! :: treatment
		do i = 0,1
			k_g(i,j) = (g_k(1,i,j) + 2d0*g_k(2,i,j) + 2d0*g_k(3,i,j) + g_k(4,i,j)) / 6d0
			funcG(0,i,j+1) = funcG(0,i,j) + k_g(i,j)*deta_BL
		end do
		!!
	end subroutine
	!!
	subroutine cal_RKstep_g(RKt, j)
		!! :: argument 
		integer,intent(in) :: RKt !! :: Runge-Kutta times
		integer,intent(in) :: j
		!! :: variables
		integer :: i = 0
		integer :: RKtm1 = 0
		real(8) :: dw = 0 !! :: deta_BL の重み delta-eta-weight
		!! :: treatment
		RKtm1 = RKt - 1
		!! :: :: RKt = 1 => dw = 0.5d0
		!! :: :: RKt = 2 => dw = 0.5d0
		!! :: :: RKt = 3 => dw = 1d0
		dw = (dble(RKt)**2d0 - 3d0*dble(RKt) + 4d0 ) * 0.25d0 
		!!
		g_k(RKt,1,j) = -0.5d0 * Pr * funcG(RKtm1,1,j) * funcF(RKtm1,0,j) - Pr * (gamma - 1d0) * Ma**2d0 * funcF(RKtm1,2,j)**2d0
		g_k(RKt,0,j) = funcG(RKtm1,1,j)	
		!!
		if(Rkt < 4) then
			funcG(RKt,1,j) = funcG(0,1,j) + g_k(RKt,1,j)*deta_BL*dw
			funcG(RKt,0,j) = funcG(0,0,j) + g_k(RKt,0,j)*deta_BL*dw
		end if
		!!
	end subroutine
	!!
end program 
