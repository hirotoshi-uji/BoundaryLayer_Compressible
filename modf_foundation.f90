!! --- --- --- --- --- --- --- --- --- ---
!! :: ./ CavityFlow_2d_1 / mod_foundation.f90
!! --- --- --- --- --- --- --- --- --- ---
!! --- --- --- --- --- --- --- --- --- ---
!! :: mod_Warning
!! --- --- --- --- --- --- --- --- --- ---
! 	call Warn_fileFaild2Open
! 	call Warn_minusLength
! 	call Warn_CourantCondition(r)
!! --- --- --- --- --- --- --- --- --- ---
module mod_Warning
	implicit none
	integer :: isNA = 0
	contains
	!! --- --- --- --- --- --- --- --- --- ---
	!! :: subroutine
	!! --- --- --- --- --- --- --- --- --- ---
	subroutine Warn_fileFaild2Open
		write(6,*) 'Failed to open'
		write(6,*)'program is stopped.'
		stop
	end subroutine
	!! --- --- --- --- --- --- --- --- --- ---
	!! :: subroutine
	!! --- --- --- --- --- --- --- --- --- ---
	subroutine Warn_minusLength
		write(6,*)"Length is unexpected in minus."
		write(6,*)'program is stopped.'
		stop
	end subroutine
	!! --- --- --- --- --- --- --- --- --- ---
	!! :: subroutine
	!! --- --- --- --- --- --- --- --- --- ---
	subroutine Warn_CourantCondition(r)
		real(8),intent(in) :: r
		write(6,*)"It is not satisfied Courant condition."
		write(6,*)'r = ',r
		write(6,*)'program is stopped.'
		stop
	end subroutine
end module
!! --- --- --- --- --- --- --- --- --- ---
!! :: 物性値，数値粘性，無次元数，グローバル変数
!! --- --- --- --- --- --- --- --- --- ---
module mod_PhysicalProp_2d
	implicit none
	real(8) :: Gamma = 0d0
	real(8) :: Cp = 0d0
	real(8),allocatable,dimension(:,:) :: mu
	real(8) :: Rgass = 0d0
	real(8) :: perRgass = 0d0
end module
!!
module mod_numericalProp
	implicit none
	real(8) :: rhoCou = 0d0							!! :: Courant 条件判別用
	real(8) :: rhoVisc = 0d0 							!! :: Viscous の条件判別用
	real(8) :: zeta = 0d0									!! :: 移流項流束分割用
	real(8) :: sigma = 0d0								!! :: DCS の係数
	real(8) :: sigma_NSCBC = 0.25d0			!! :: NSCBC 用の sigma
	real(8) :: K_NSCBC = 0d0						!! :: NSCBC 用の K
	real(8) :: preK_NSCBC = 0d0					!! :: NSCBC 用の K （preset）
end module
!!
module mod_nonDimensionalNumber
	implicit none
	real(8) :: Re =0d0			!! :: Reynolds 数
	real(8) :: perRe = 0d0 		!! :: Reynolds 数の逆数
	real(8) :: Pr = 0d0			!! :: Prandtl 数
	real(8) :: Ma = 0d0			!! :: Mach 数
	real(8) :: Sc = 0d0			!! :: Sutherland 数
	real(8) :: Ma_max = 0d0	!! :: 最大 Mach 数
	real(8),allocatable,dimension(:,:) :: Ma_local_x 
	!!
	real(8) :: Lref =1d0			!! :: 代表長さ
end module
!!
module mod_declarVar_2d
	implicit none
	!! :: 等温条件
	real(8) :: constTemp = 0d0
	real(8) :: constPres = 0d0
	real(8) :: checkChar = 0d0
	!! :: Nj の 1 次元ベクトル --- --- --- --- --- --- --- --- --- ---
	real(8),allocatable,dimension(:) :: x_Sec1 , x_Sec2 , y_Sec3 , y_Sec4 , x , y 
	!! :: Nj の 2 次元ベクトル --- --- --- --- --- --- --- --- --- ---
	real(8),allocatable,dimension(:,:) :: Temp,Temp_xpr,Temp_ypr
	real(8),allocatable,dimension(:,:) :: soundSpd
	real(8),allocatable,dimension(:,:) :: dilatation
	real(8),allocatable,dimension(:,:) :: ViscStr_xx , ViscStr_xy , ViscStr_yy
	real(8),allocatable,dimension(:,:) :: ViscStr_xx_pr , ViscStr_xy_xpr , ViscStr_xy_ypr , ViscStr_yy_pr
	real(8),allocatable,dimension(:,:) :: vorticity
	!! :: Nj の 3 次元ベクトル --- --- --- --- --- --- --- --- --- ---
	real(8),allocatable,dimension(:,:,:) :: heatFlux , heatFlux_pr
	real(8),allocatable,dimension(:,:,:) :: rhoup , rhoup_xpr , rhoup_ypr
	!! :: 時間変化項で用いるベクトル --- --- --- --- --- --- --- --- --- ---
	real(8),allocatable,dimension(:,:,:) :: Qn , Qn_pr , Qnp1		
	real(8),allocatable,dimension(:,:,:) :: Qn_pr_1st , Qn_pr_2nd , Qn_pr_3rd							
	!! :: 移流項で用いるベクトル --- --- --- --- --- --- --- --- --- ---
	real(8),allocatable,dimension(:,:,:) :: Fnx , Fnx_p , Fnx_p_pr ,  Fnx_m , Fnx_m_pr , Fnx_pr 
	real(8),allocatable,dimension(:,:,:) :: Fny , Fny_p , Fny_p_pr ,  Fny_m , Fny_m_pr , Fny_pr	
	!! :: 粘性項で用いるベクトル --- --- --- --- --- --- --- --- --- ---
	real(8),allocatable,dimension(:,:,:) :: Vnx , Vnp1x , Vnx_pr
	real(8),allocatable,dimension(:,:,:) :: Vny , Vnp1y , Vny_pr 
	!! :: NSCBC で用いるベクトル --- --- --- --- --- --- --- --- --- ---
	real(8),allocatable,dimension(:,:,:) :: Phix ,Phiy , dVecx , dVecy
	real(8),allocatable,dimension(:,:) :: cx , cy ,MaVecx , MaVecy
end module
module mod_BLvar
	implicit none
	!! :: parameter
	!! :: variable
	real(8) :: deta_BL
	real(8) :: f_0 , f_1 , f_2 , f_3
	real(8) :: g_0 , g_1 , g_2 

end module
module mod_Grid
	use mod_Warning
	implicit none
	integer,allocatable,dimension(:) :: Nx_Sec1 , Nx_Sec2 , Nx_Sec5 , Nx_Sec7 , Nx_Sec9 , Nx_Sec10
	integer,allocatable,dimension(:) :: Ny_Sec3 , Ny_Sec4 , Ny_Sec6 , Ny_Sec8 , Ny_Sec11
	contains
	subroutine cal_WholeGridNumber(Nj_max,j0_max,jEnd_max,jd)
		!! :: argument --- --- --- --- --- --- --- --- --- ---
		integer,intent(inout) :: Nj_max
		real(8),intent(in) :: j0_max
		real(8),intent(in) :: jEnd_max
		real(8),intent(in) :: jd
		!! :: variable --- --- --- --- --- --- --- --- --- ---
		real(8) :: Length
		!! :: treatment --- --- --- --- --- --- --- --- --- ---
		Length = jEnd_max - j0_max
		if (Length < 0) call Warn_minusLength																						!! :: mod_foundation.f90 / mod_Warning 
		Nj_max = ceiling(Length * jd)
	end subroutine
	!!
	subroutine cal_GridNumber(Nj_Sec,j0_Sec,jEnd_Sec,jd,j0_max)
		!! :: argument --- --- --- --- --- --- --- --- --- ---
		integer,intent(inout),dimension(0:2) :: Nj_Sec
		real(8),intent(in) :: j0_Sec
		real(8),intent(in) :: jEnd_Sec
		real(8),intent(in) :: jd
		real(8),intent(in) :: j0_max
		!! :: variable --- --- --- --- --- --- --- --- --- ---
		real(8) :: Length
		!! :: treatment --- --- --- --- --- --- --- --- --- ---
		!! :: Section の幅を計算
		Length = jEnd_Sec - j0_Sec
		if (Length < 0) call Warn_minusLength																						!! :: mod_foundation.f90 / mod_Warning 
		Nj_Sec(2) = ceiling(Length * jd)
		!! :: Section の始まりと終わりの座標
		Length = j0_Sec - j0_max
		if (Length < 0) call Warn_minusLength																						!! :: mod_foundation.f90 / mod_Warning 
		Nj_Sec(0) = ceiling(Length * jd)
		Nj_Sec(1) = Nj_Sec(0) + Nj_Sec(2)
	end subroutine
end module
!! --- --- --- --- --- --- --- --- --- ---
!! :: ファイル操作
!! --- --- --- --- --- --- --- --- --- ---
module mod_parOnFile
	use mod_Warning
	use mod_declarVar_2d
	use mod_Grid
	implicit none
	!! :: argument --- --- --- --- --- --- --- --- --- ---
	integer,parameter :: fnumPar = 11			!! :: ParameterFile のファイル番号
	integer,parameter :: fnumDAT = 12			!! :: DATfile のファイル番号
	integer,parameter :: fnumCSV = 13			!! :: CSVfile のファイル番号
	integer,parameter :: fnumCon = 14			!! :: Condition のファイル番号
	!! :: variable --- --- --- --- --- --- --- --- --- ---
	integer :: ios = 0 , ios1 = 0
	integer :: opStepT = 0 , opStepX = 0 ,opStepY = 0
	!!
	character(50) :: fnameCSV 
	character(50) :: fnameDAT_t0 			= './../Result/Rs_F_rhoup0000.00.dat' 															!! :: 出力ファイル名
	character(50) :: fnameDAT 					= './../Result/Rs_F_rhoup0000.00.dat' 															!! :: 出力ファイル名
	character(50) :: fnameInitDAT 			= './../Result/Rs_F_rhoup0_init_2.dat'															!! :: 出力ファイル名
	character(50) :: fnameDAT_RKF			= './../Result/Rs_F_RKF_rhoup0000.dat' 												!! :: 出力ファイル名
	character(50) :: fnameDAT_RKS			= './../Result/Rs_F_RKS_rhoup0000.dat' 												!! :: 出力ファイル名
	character(50) :: fnameCon					= './../Result/00_condition.txt'																	!! :: 出力ファイル名	
	character(50) :: fnameVortix				= './../Result/Rs_F_vorticity000.dat' 														!! :: 出力ファイル名
	character(50) :: fnameParSpace			= './Par_Space.dat' 																					!! :: 入力ファイル名
	character(50) :: fnameParTime			= './Par_Time.dat' 																					!! :: 入力ファイル名
	character(50) :: fnameParPhysProp	= './Par_PhysProp.dat' 																			!! :: 入力ファイル名
	contains
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		subroutine op_FieldOnDAT_AP(Nx,Ny,Mat,DfnameDAT)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: Nx
			integer,intent(in) :: Ny
			real(8),intent(in),dimension(0:Ny,0:Nx) :: Mat
			character(50) :: DfnameDAT
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: j = 0 , k = 0
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnumDAT ,iostat=ios1,file=DfnameDAT,action='write',form='formatted',status='replace')
			do j = 0,Ny
				do k = 0,Nx
					write(fnumDAT,'(3(f18.12))')x(k),y(j),Mat(j,k)
				end do
				write(fnumDAT,*)""
			end do
			close(fnumDAT)
		end subroutine
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		subroutine op_FieldOnDAT_rhoup(Nx,Ny,DfnameDAT)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: Nx
			integer,intent(in) :: Ny
			character(50) :: DfnameDAT
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: j = 0 , k = 0
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnumDAT ,iostat=ios1,file=DfnameDAT,action='write',form='formatted',status='replace') 	
			do j = Ny_Sec6(0),Ny_Sec6(1)
				do k = Nx_Sec2(0),Nx_Sec2(1)
					write(fnumDAT,'(6(f18.12))')x(k),y(j),rhoup(0,j,k),rhoup(1,j,k),rhoup(2,j,k),rhoup(3,j,k)
				end do
				write(fnumDAT,*)""
			end do
			do j = Ny_Sec4(0),Ny_Sec4(1)
				do k = 0,Nx
					write(fnumDAT,'(6(f18.12))')x(k),y(j),rhoup(0,j,k),rhoup(1,j,k),rhoup(2,j,k),rhoup(3,j,k)
				end do
				write(fnumDAT,*)""
			end do
			close(fnumDAT)
		end subroutine
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		subroutine op_FieldOnDAT_NaN_AP(Nx,Ny)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: Nx
			integer,intent(in) :: Ny
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: j = 0 , k = 0 , h = 0
			character(50) :: DfnameDAT = "./../Result/Rs_F_NaN_rhoup.dat"
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnumDAT ,iostat=ios1,file=DfnameDAT,action='write',form='formatted',status='replace')
			do j = Ny_Sec6(0),Ny_Sec6(1)
				do k = Nx_Sec2(0),Nx_Sec2(1)
					write(fnumDAT,'(6(f18.12))')x(k),y(j),rhoup(0,j,k),rhoup(1,j,k),rhoup(2,j,k),rhoup(3,j,k)
				end do
				write(fnumDAT,*)""
			end do
			do j = Ny_Sec4(0),Ny_Sec4(1)
				do k = 0,Nx
					write(fnumDAT,'(6(f18.12))')x(k),y(j),rhoup(0,j,k),rhoup(1,j,k),rhoup(2,j,k),rhoup(3,j,k)
				end do
				write(fnumDAT,*)""
			end do
			close(fnumDAT)
		end subroutine
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		subroutine op_FieldOnDATY_AP(Nx,Ny,Mat,k,DfnameDAT)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: Nx
			integer,intent(in) :: Ny
			integer,intent(in) :: k
			real(8),intent(in),dimension(0:Ny,0:Nx) :: Mat
			character(50) :: DfnameDAT
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: j = 0
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnumDAT ,iostat=ios1,file=DfnameDAT,action='write',form='formatted',status='replace') ; ios = ios + ios1
			if (ios /= 0) call Warn_fileFaild2Open	
			do j = 0,Ny
				write(fnumDAT,'(3(f18.12))')y(j),Mat(j,k)
			end do
			close(fnumDAT)
		end subroutine
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		subroutine op_FieldOnDATX_AP(Nx,Ny,Mat,j,DfnameDAT)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: Nx
			integer,intent(in) :: Ny
			integer,intent(in) :: j
			real(8),intent(in),dimension(0:Ny,0:Nx) :: Mat
			character(50) :: DfnameDAT
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer ::  k = 0
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnumDAT ,iostat=ios1,file=DfnameDAT,action='write',form='formatted',status='replace') ; ios = ios + ios1
			if (ios /= 0) call Warn_fileFaild2Open	
			do k = 0,Nx
				write(fnumDAT,'(3(f18.12))')x(k),Mat(j,k)
			end do
			close(fnumDAT)
		end subroutine
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		subroutine op_FieldOnDATY_xPyAP(Nx,Ny,Mat,k,fnameDatY)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: Nx
			integer,intent(in) :: Ny
			integer,intent(in) :: k
			real(8),intent(in),dimension(0:Ny,0:Nx-1) :: Mat
			character(40) :: fnameDatY
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: j = 0 
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnumDAT ,iostat=ios1,file=fnameDatY,action='write',form='formatted',status='replace') ; ios = ios + ios1
			if (ios /= 0) call Warn_fileFaild2Open	
			do j = 0,Ny
				write(fnumDAT,'(2(f18.12))')y(j),Mat(j,k)
			end do
			close(fnumDAT)
		end subroutine
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		subroutine op_FieldOnDATX_xPyAP(Nx,Ny,Mat,j,DfnameDatX)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: Nx
			integer,intent(in) :: Ny
			integer,intent(in) :: j
			real(8),intent(in),dimension(0:Ny,0:Nx-1) :: Mat
			character(40) :: DfnameDatX
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: k = 0 
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnumDAT ,iostat=ios1,file=DfnameDatX,action='write',form='formatted',status='replace') ; ios = ios + ios1
			if (ios /= 0) call Warn_fileFaild2Open	
			do k = 0,Nx-1
				write(fnumDAT,'(2(f18.12))')x(k),Mat(j,k)
			end do
			close(fnumDAT)
		end subroutine
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		subroutine op_FieldOnDAT_xPyAP(Nx,Ny,Mat,DfnameDAT)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: Nx
			integer,intent(in) :: Ny
			real(8),intent(in),dimension(0:Ny,0:Nx-1) :: Mat
			character(40) :: DfnameDat
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: j = 0 , k = 0
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnumDAT ,iostat=ios1,file=DfnameDat,action='write',form='formatted',status='replace') ; ios = ios + ios1
			if (ios /= 0) call Warn_fileFaild2Open	
			do j = 0,Ny
				write(fnumDAT,*)"# y = ",j
				do k = 0,Nx-1
					write(fnumDAT,'(3(f18.12))')x(k),y(j),Mat(j,k)
				end do
				write(fnumDAT,*)""
			end do
			close(fnumDAT)
		end subroutine
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		function openWriteFile(fnum,fname)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: fnum
			character(50),intent(in) :: fname
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: openWriteFile
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnum ,iostat=ios1,file=fname,action='write',form='formatted',status='replace')
			openWriteFile = ios + ios1
		end function
		!! --- --- --- --- --- --- --- --- --- ---
		!! :: function :: 
		!! --- --- --- --- --- --- --- --- --- ---
		function openReadFile(fnum,fname)
			!! :: argument --- --- --- --- --- --- --- --- --- ---
			integer,intent(in) :: fnum
			character(50),intent(in) :: fname
			!! :: variable --- --- --- --- --- --- --- --- --- ---
			integer :: openReadFile
			!! :: treatment --- --- --- --- --- --- --- --- --- ---
			open (unit = fnum ,iostat=ios1,file=fname,action='read',form='formatted',status='old')
			openReadFile = ios + ios1
		end function
end module

!! --- --- --- --- --- --- --- --- --- ---
!! :: mod_WriteAndDisplay.f90
!! --- --- --- --- --- --- --- --- --- ---
!	call writer_3VectorOnFile(fnum,Nj,opStepX,t,Vec1,Vec2,Vec3)
!	call writer_nVectorOnFile(fnum,Nj,opStepX,Nvec,t,Vec)
!	call writer_nVectorScaler(fnum,Nj,opStepX,Nvec,Vec)
!! --- --- --- --- --- --- --- --- --- ---
! 	call displayer_Vector(fnum,Nj,opStepX,Vec)
! 	call displayer_Matrix(fnum,Nj,opStepX,Mat)
!! --- --- --- --- --- --- --- --- --- ---
! 	call ewriter_Vector(fnum,Nj,opStepX,Vec)
! 	call ewriter_timeColumn(fnum,t)
!	call ewriter_newline(fnum)
! 	call ewriter_2newline(fnum)
!	call ewriter_blank(fnum)
!! --- --- --- --- --- --- --- --- --- ---
module mod_WriteAndDisplay
	implicit none
	contains
!! --- --- --- --- --- --- --- --- --- ---
!! :: writer 
!! --- --- --- --- --- --- --- --- --- ---
	!! --- --- --- --- ---
	!! :: subroutine ::
	!! --- --- --- --- ---
	subroutine writer_3VectorOnFile(fnum,Nj,opStepX,t,Vec1,Vec2,Vec3)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		integer,intent(in) :: Nj
		integer,intent(in) :: opStepX
		real(8),intent(in) :: t
		real(8),intent(in),dimension(0:Nj) :: Vec1 , Vec2 , Vec3
		!! :: treatment --- --- --- --- ---
		write(6,*)t, 'fnum : ',fnum
		call ewriter_timeColumn(fnum,t)
		call ewriter_Vector(fnum,Nj,opStepX,Vec1)
		call ewriter_blank(fnum)
		call ewriter_timeColumn(fnum,t)
		call ewriter_Vector(fnum,Nj,opStepX,Vec2)
		call ewriter_blank(fnum)
		call ewriter_timeColumn(fnum,t)
		call ewriter_Vector(fnum,Nj,opStepX,Vec3)
		call ewriter_newline(fnum)
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine ::
	!! --- --- --- --- ---
	subroutine writer_nVectorOnFile(fnum,Nj,opStepX,Nvec,t,Vec)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		integer,intent(in) :: Nj
		integer,intent(in) :: opStepX
		integer,intent(in) :: Nvec
		real(8),intent(in) :: t
		real(8),intent(in),dimension(0:(Nvec-1),0:Nj) :: Vec
		!! :: variable --- --- --- --- ---
		integer :: iVec = 0
		!! :: treatment --- --- --- --- ---
		write(6,*)t, 'fnum : ',fnum
		do iVec = 0 , Nvec-1
			call ewriter_timeColumn(fnum,t)
			call ewriter_Vector(fnum,Nj,opStepX,Vec(iVec,:))
			if ( iVec /= Nvec -1 ) call ewriter_blank(fnum)
		end do
		call ewriter_newline(fnum)
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine ::
	!! --- --- --- --- ---
	subroutine writer_nScaler(fnum,Nj,opStepX,Nvec,x)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		integer,intent(in) :: Nj
		integer,intent(in) :: opStepX
		integer,intent(in) :: Nvec
		real(8),intent(in),dimension(0:Nj) :: x
		!! :: variable --- --- --- --- ---
		integer :: iVec = 0
		!! :: treatment --- --- --- --- ---
		do iVec = 0 , Nvec-1
			call ewriter_partialBlank(fnum) !! :: time の列を開けるため
			call ewriter_Vector(fnum,Nj,opStepX,x)
			if ( iVec /= Nvec -1 ) call ewriter_blank(fnum)
		end do
		call ewriter_newline(fnum)
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine ::
	!! --- --- --- --- ---
	subroutine writer_nDataName(fnum,Nj,opStepX,Nvec,nameRec)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		integer,intent(in) :: Nj
		integer,intent(in) :: opStepX
		integer,intent(in) :: Nvec
		character(10),intent(in),dimension(Nvec) :: nameRec
		!! :: variable --- --- --- --- ---
		integer :: iVec = 0 , j = 0
		!! :: treatment --- --- --- --- ---
		do iVec = 1 , Nvec
			write(fnum,'(a)',advance="no")'#t' !! :: time の列を開けるため
			write(fnum,'(",",a)',advance="no")nameRec(iVec)
			do j = 1,Nj
				write(fnum,'(a)',advance="no")','
			end do
			if ( iVec /= Nvec ) call ewriter_blank(fnum)
		end do
		call ewriter_newline(fnum)
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine :: CSVFile の初期設定
	!! --- --- --- --- ---
	subroutine writer_csvInit(fnumCSV,Nj,opStepX,Nvec,nameRec,t,x,Vec)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnumCSV
		integer,intent(in) :: Nj
		integer,intent(in) :: opStepX
		integer,intent(in) :: Nvec
		character(10),intent(in),dimension(Nvec) :: nameRec
		real(8),intent(in) :: t
		real(8),intent(in),dimension(0:Nj) :: x
		real(8),intent(in),dimension(0:(Nvec-1),0:Nj) :: Vec
		!! :: treatment --- --- --- --- ---
		call writer_nDataName(fnumCSV,Nj,opStepX,Nvec,nameRec)
		call writer_nScaler(fnumCSV,Nj,opStepX,Nvec,x)
		call writer_nVectorOnFile(fnumCSV,Nj,opStepX,Nvec,t,Vec)
	end subroutine
!! --- --- --- --- --- --- --- --- --- ---
!! :: displayer 
!! --- --- --- --- --- --- --- --- --- ---
	!! --- --- --- --- --- --- --- --- --- ---
	!! :: subroutine :: 
	!! --- --- --- --- --- --- --- --- --- ---
	subroutine displayer_Vector(fnum,Nj,opStepX,Vec)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		integer,intent(in) :: Nj
		integer,intent(in) :: opStepX
		real(8),intent(in),dimension(0:Nj) :: Vec
		!! :: variable --- --- --- --- ---
		call ewriter_Vector(fnum,Nj,opStepX,Vec)
		call ewriter_2newline(fnum) 
	end subroutine!! --- --- --- --- --- --- --- --- --- ---
	!! :: subroutine 
	!! --- --- --- --- --- --- --- --- --- ---
	subroutine displayer_Matrix(fnum,Nj,opStepX,Mat)
	!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		integer,intent(in) :: Nj
		integer,intent(in) :: opStepX
		real(8),intent(in),dimension(0:Nj,0:Nj) :: Mat
	!! :: variable value --- --- --- --- ---
		integer :: i = 0 , j = 0
	!! :: treatment --- --- --- --- ---
		do i = 0,Nj
			write(fnum,'(1x,f14.8)',advance = 'no')Mat(i,0)
			do j = 1,Nj
				write(fnum,'(",",1x,f14.8)',advance = 'no')Mat(i,j)
			end do
			call ewriter_newline(fnum)
		end do
		call ewriter_2newline(fnum)
	end subroutine
!! --- --- --- --- --- --- --- --- --- ---
!! :: element Writer
!! --- --- --- --- --- --- --- --- --- ---
	!! --- --- --- --- --- --- --- --- --- ---
	!! :: subroutine :: 
	!! --- --- --- --- --- --- --- --- --- ---
	subroutine ewriter_Vector(fnum,Nj,opStepX,Vec)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		integer,intent(in) :: Nj
		integer,intent(in) :: opStepX
		real(8),intent(in),dimension(0:Nj) :: Vec
		!! :: variable --- --- --- --- ---
		integer :: i = 0
		!! :: treatment --- --- --- --- ---
		write (fnum,'(1x,f14.8)',advance='no')Vec(0) !! :: i = 0
		do i = 1, Nj
			if ( mod( i, Nj / opStepX ) == 0 ) then 
				write (fnum,'(",",1x,f14.8)',advance='no')Vec(i)
			end if 
		end do
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine :: 時間の列を追加する
	!! --- --- --- --- ---
	subroutine ewriter_timeColumn(fnum,t)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		real(8),intent(in) :: t
		!! :: treatment --- --- --- --- ---
		write (fnum,'(f14.8)',advance='no')t
		write (fnum,'(a)',advance='no')','
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine :: 改行する
	!! --- --- --- --- ---
	subroutine ewriter_newline(fnum)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		!! :: treatment --- --- --- --- ---
		write(fnum,'(a)')''
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine :: 2階改行する
	!! --- --- --- --- ---
	subroutine ewriter_2newline(fnum)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		!! :: treatment --- --- --- --- ---
		call ewriter_newline(fnum)
		call ewriter_newline(fnum)
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine :: 空白
	!! --- --- --- --- ---
	subroutine ewriter_blank(fnum)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		!! :: treatment --- --- --- --- ---
		write (fnum,'(a)',advance='no')',,'
	end subroutine
	!! --- --- --- --- ---
	!! :: subroutine :: （先頭）片側空白
	!! --- --- --- --- ---
	subroutine ewriter_partialBlank(fnum)
		!! :: argument --- --- --- --- ---
		integer,intent(in) :: fnum
		!! :: treatment --- --- --- --- ---
		write (fnum,'(a)',advance='no')','
	end subroutine
end module
