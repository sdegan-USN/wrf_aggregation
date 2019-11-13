Module module_volc_emiss_driver
CONTAINS

!subroutine volcemiss(emis_vol,ims,ime,kms,kme,jms,jme,num_emis_vol,ni,nj,julday_wrf,curr_secs,xlong,xlat,z,i,j,kts,dx,p_e_vash1,p_e_vash2,p_e_vash3,p_e_vash4,p_e_vash5,p_e_vash6,p_e_vash7,p_e_vash8,p_e_vash9,p_e_vash10,p_e_vso2,p_ke_h2o,z_at_w,gmt,config_flags) !Added h2o - SDE 1NOV18
 subroutine volcemiss(emis_vol,ims,ime,kms,kme,jms,jme,num_emis_vol,ni,nj,julday_wrf,curr_secs,xlong,xlat,z,i,j,kts,dx,p_e_vash1,p_e_vash2,p_e_vash3,p_e_vash4,p_e_vash5,p_e_vash6,p_e_vash7,p_e_vash8,p_e_vash9,p_e_vash10,p_e_vso2,p_ke_h2o,z_at_w,gmt,emissoptvol) !Added h2o - SDE 1NOV18
!USE module_configure ! Added for config_flags - SDE 1NOV18
! USE module_configure,only:  grid_config_rec_type
IMPLICIT NONE
! integer :: ims,ime,kms,kme,jms,jme,num_emis_vol,ni,nj,julday_wrf
integer :: ims,ime,kms,kme,jms,jme,num_emis_vol,ni,nj,julday_wrf,emissoptvol
real :: emis_vol(ims:ime,kms:kme,jms:jme,num_emis_vol),xlong(ims:ime,jms:jme),xlat(ims:ime,jms:jme)
real :: z(ims:ime,kms:kme,jms:jme)
real(kind=8) :: curr_secs,secs
REAL :: kcurr_time
integer :: k,i,j,kts
!TYPE(grid_config_rec_type), INTENT(IN  ) :: config_flags

INTEGER :: hgt,timestep,ntimestep
REAL,ALLOCATABLE :: ASH(:,:),SO2(:,:),ASH_EMIS(:),SO2_EMIS(:),H2O_EMIS(:)
INTEGER,ALLOCATABLE :: DATE(:),TIME(:),JULDAY(:),DURATION(:),HEIGHT(:),TIMEW(:)
integer :: year,month,day
integer :: timestep_curr

integer :: nlevels,ref_date
integer :: nr_abschnitte_top,nr_abschnitte_unten
real :: volc_top,volc_unten
REAL,ALLOCATABLE :: VOLC_HEIGHTS_TOP(:),VOLC_HEIGHTS_BOTTOM(:)
REAL,ALLOCATABLE :: hohe_wrf_top(:),hohe_wrf_middle(:),hohe_wrf_unten(:),volc_emiss_wrf(:)
REAL :: total_mass_extern,total_mass_wrf,mass_wrf_k
REAL :: total_mass_extern_so2,total_mass_wrf_so2,mass_wrf_k_so2
INTEGER :: kk
LOGICAL :: level_check
REAL :: HEIGHT_SURFACE
REAL :: area
REAL :: dx
INTEGER :: p_e_vash1,p_e_vash2,p_e_vash3,p_e_vash4,p_e_vash5,p_e_vash6,p_e_vash7,p_e_vash8,p_e_vash9,p_e_vash10,p_e_vso2,p_ke_h2o ! Added h20 - SDE 1NOV18
REAL,ALLOCATABLE :: BIN_DISTR(:)
INTEGER :: NBINS,bi

INTEGER :: time_before_volc,time_wrf,time_after_volc
REAL :: OFFSET_EXT_WRF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL :: percen_mass_umbrel,base_umbrel,curr_hours
INTEGER :: ivolcano
REAL :: gmt
INTEGER :: gmtm,gmtp
REAL :: begday,beghr,begmin,endday,endhr,endmin
REAL,ALLOCATABLE :: erup_beg(:,:,:), erup_end(:,:,:),erup_hgt(:,:,:),erup_ash(:,:,:)!,erup_so2(:,:,:),erup_h2o(:,:,:) !!!!!
INTEGER :: num_erup
INTEGER :: ki
REAL :: so2_mass,h2o_flux,emiss_ash_mass
REAL :: emiss_ash_height,eh,erup_so2,erup_h2o
REAL :: ashz_above_vent
REAL :: z_at_w_volc,z_at_w(ims:ime,kms:kme,jms:jme )
INTEGER :: k_final,k_initial
INTEGER :: kk4,x1,ko,kl
REAL :: vert_mass_dist(kts:kme)
INTEGER :: p_ksbin1,p_ksbin2,p_ksbin3,p_ksbin4,p_ksbin5,p_ksbin6,p_ksbin7,p_ksbin8,p_ksbin9,p_ksbin10 
REAL :: VOLC_HOUR,VOLC_MIN,VOLC_SEC,VOLC_TIME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

z_at_w_volc=z(ni,kts,nj)
num_erup=4

percen_mass_umbrel=.75
base_umbrel=.25    ! fraction


ivolcano=0
area=dx*dx

ntimestep=9
nlevels=46
emis_vol=0
NBINS=10
OFFSET_EXT_WRF=0

ALLOCATE(DATE(ntimestep),TIME(ntimestep),HEIGHT(ntimestep),DURATION(ntimestep),ASH_EMIS(ntimestep),SO2_EMIS(ntimestep),H2O_EMIS(ntimestep),TIMEW(ntimestep))
ALLOCATE(JULDAY(ntimestep))
ALLOCATE(BIN_DISTR(NBINS))

! READ IN volc.asc BEG
OPEN(14,FILE="volc_d01.asc")

DO timestep=1,6
  READ (14,*)
END DO

READ(14,*) (BIN_DISTR(bi),bi=1,NBINS)
DO timestep=1,7
  READ (14,*)
END DO

DO timestep=1,ntimestep
  READ(14,*) DATE(timestep),TIME(timestep),DURATION(timestep),HEIGHT(timestep)&
    &,ASH_EMIS(timestep),SO2_EMIS(timestep),H2O_EMIS(timestep)
  year=INT(DATE(timestep)/10000)
  ref_date=year*10000+0101
  JULDAY(timestep)=juldate(DATE(timestep),000000)-juldate(ref_date,000000)
  VOLC_HOUR=INT(TIME(timestep)/10000.)
  VOLC_MIN=INT((TIME(timestep)-VOLC_HOUR*10000.)/100.)
  VOLC_SEC=INT((TIME(timestep)-VOLC_HOUR*10000.-VOLC_MIN*100.))
  TIMEW(timestep)=VOLC_HOUR*3600.+VOLC_MIN*60.+VOLC_SEC
  write(*,*) "MARCUS",timestep,JULDAY(timestep),TIME(timestep),HEIGHT(timestep)
  write(*,*) "MARCUS",TIMEW(timestep)
END DO

CLOSE(14)
IF(SUM(BIN_DISTR(:)).ne.1.) write(*,*) "MARCUS MASS NOT 1"
! READ IN volc.asc END

secs=mod(gmt*3600.+curr_secs,86400.)
write(*,*) "PROBLEM2",julday_wrf,gmt,curr_secs,secs
write(*,*) "PROBLEM2"

! CHECK IF TIMES FIT BEG
timestep_curr=0
DO timestep=1,ntimestep-1 ! last time step in input file should have 0 emission as it should be the end time for the last intervall
  time_before_volc=JULDAY(timestep)*100000.+TIMEW(timestep)
  time_wrf=julday_wrf*100000.+secs ! year*10000000000 should be here as well
  time_after_volc=JULDAY(timestep+1)*100000.+TIMEW(timestep+1)

write(*,*) "MARCUS ddhh3",time_before_volc,time_wrf,time_after_volc

write(*,*) "TIMESTEP1",JULDAY(timestep),DATE(timestep),TIMEW(timestep)
write(*,*) "TIMESTEP1",julday_wrf,secs
write(*,*) "TIMESTEP1"

  IF((time_wrf.ge.time_before_volc).and.(time_wrf.lt.time_after_volc)) THEN
    timestep_curr=timestep
    write(*,*) "TIMESTEP2",time_before_volc,time_wrf,time_after_volc
    write(*,*) "TIMESTEP2",timestep,DATE(timestep),TIMEW(timestep)
    write(*,*) "TIMESTEP2",HEIGHT(timestep),ASH_EMIS(timestep)
    write(*,*) "TIMESTEP2"
  END IF
END DO
! CHECK IF TIMES FIT END

!Commented out Mastin's, using our caluclations for mass
emiss_ash_height = HEIGHT(timestep_curr)
!eh=2600.*(emiss_ash_height*.0005)**4.1494
eh=ASH_EMIS(timestep_curr)
emiss_ash_mass=eh*1.e9/area
ivolcano=1
!ading SO2 and water vapor - SDE 1NOV18
erup_so2=SO2_EMIS(timestep_curr)
erup_h2o=H2O_EMIS(timestep_curr)
so2_mass=erup_so2*3600.*1.e9/64./area !molecular weight SO2 is 64
h2o_flux=erup_h2o/area ! kg/s ---> kg/(m2s)


! UMBRELLA PART BEG
!ashz_above_vent=emiss_ash_height - z_at_w_volc
ashz_above_vent=emiss_ash_height - z_at_w(ni,kts,nj)
write(*,*) "MARCUS HEIGHT",kts,ashz_above_vent,emiss_ash_height,z_at_w(ni,kts,nj)

do k=kme-1,kts,-1
write(*,*) "LUCIA1",k,z_at_w(ni,k,nj),emiss_ash_height
  if(z_at_w(ni,k,nj) < emiss_ash_height) then
    k_final=k+1
write(*,*) "LUCIA1",z_at_w(ni,k,nj),emiss_ash_height,k_final
    exit
  endif
enddo
do k=kme-1,kts,-1
write(*,*) "LUCIA2",k,z_at_w(ni,k,nj),(1.-base_umbrel),ashz_above_vent,z_at_w(ni,kts,nj)
  if(z_at_w(ni,k,nj) < ((1.-base_umbrel)*ashz_above_vent)+z_at_w(ni,kts,nj))then
    k_initial=k
write(*,*) "LUCIA2",z_at_w(ni,k,nj),((1.-base_umbrel)*ashz_above_vent)+z_at_w(ni,kts,nj),k_initial
    exit
  endif
enddo
write(*,*) "MARCUS HEIGHT",k_initial,k_final


vert_mass_dist=0.
kk4 = k_final-k_initial+2
do ko=1,kk4-1
  kl=ko+k_initial-1
  vert_mass_dist(kl) = 6.*percen_mass_umbrel* float(ko)    &
  /float(kk4)**2 * (1. - float(ko)/float(kk4))
enddo
if(sum(vert_mass_dist(kts:kme)) .ne. percen_mass_umbrel) then
  x1= ( percen_mass_umbrel- sum(vert_mass_dist(kts:kme)) )   &
  /float(k_final-k_initial+1)
  do ko=k_initial,k_final
    vert_mass_dist(ko) = vert_mass_dist(ko)+ x1
  enddo
endif

do ko=1,k_initial-1   
  vert_mass_dist(ko)=float(ko)/float(k_initial-1)
enddo
x1=sum(vert_mass_dist(1:k_initial-1))
do ko=1,k_initial-1
  vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
enddo
!UMBRELLA PART END


! WRITE OUT TO ARREY
  do ko=1,k_final
write(*,*) "MARCUS VERT",ko,BIN_DISTR(1),vert_mass_dist(ko),emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash1)=BIN_DISTR(1)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash2)=BIN_DISTR(2)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash3)=BIN_DISTR(3)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash4)=BIN_DISTR(4)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash5)=BIN_DISTR(5)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash6)=BIN_DISTR(6)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash7)=BIN_DISTR(7)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash8)=BIN_DISTR(8)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash9)=BIN_DISTR(9)*vert_mass_dist(ko)*emiss_ash_mass
    emis_vol(ni,ko,nj,p_e_vash10)=BIN_DISTR(10)*vert_mass_dist(ko)*emiss_ash_mass
    !if(config_flags%emiss_opt_vol == 2)emis_vol(ni,ko,nj,p_e_vso2)=vert_mass_dist(ko)*so2_mass
    if(emissoptvol == 2)emis_vol(ni,ko,nj,p_e_vso2)=vert_mass_dist(ko)*so2_mass
    !if(config_flags%emiss_opt_vol == 2)emis_vol(ni,ko,nj,p_ke_h2o)=vert_mass_dist(ko)*h2o_flux
    if(emissoptvol == 2)emis_vol(ni,ko,nj,p_ke_h2o)=vert_mass_dist(ko)*h2o_flux
  enddo
  do ko=k_final+1,kme
    emis_vol(ni,ko,nj,p_e_vash1)=0.
    emis_vol(ni,ko,nj,p_e_vash2)=0.
    emis_vol(ni,ko,nj,p_e_vash3)=0.
    emis_vol(ni,ko,nj,p_e_vash4)=0.
    emis_vol(ni,ko,nj,p_e_vash5)=0.
    emis_vol(ni,ko,nj,p_e_vash6)=0.
    emis_vol(ni,ko,nj,p_e_vash7)=0.
    emis_vol(ni,ko,nj,p_e_vash8)=0.
    emis_vol(ni,ko,nj,p_e_vash9)=0.
    emis_vol(ni,ko,nj,p_e_vash10)=0.
    !if(config_flags%emiss_opt_vol == 2)emis_vol(ni,ko,nj,p_e_vso2)=0.
    if(emissoptvol == 2)emis_vol(ni,ko,nj,p_e_vso2)=0.
    !if(config_flags%emiss_opt_vol == 2)emis_vol(ni,ko,nj,p_ke_h2o)=0.
    if(emissoptvol == 2)emis_vol(ni,ko,nj,p_ke_h2o)=0.
  enddo

end subroutine volcemiss


FUNCTION juldate(YYYYMMDD,HHMISS)
IMPLICIT NONE
integer, parameter :: double_k = selected_real_kind(15,307)
integer, parameter :: dp=double_k
INTEGER :: YYYYMMDD,YYYY,MM,DD,HH,MI,SS,HHMISS
INTEGER :: JULDAY,JY,JM,JA,IGREG
REAL (kind=dp) ::  JULDATE
PARAMETER (IGREG=15+31*(10+12*1582))
YYYY=YYYYMMDD/10000
MM=(YYYYMMDD-10000*YYYY)/100
DD=YYYYMMDD-10000*YYYY-100*MM
HH=HHMISS/10000
MI=(HHMISS-10000*HH)/100
SS=HHMISS-10000*HH-100*MI
IF (YYYY.EQ.0) PAUSE 'There is no Year Zero.'
IF (YYYY.LT.0) YYYY=YYYY+1
IF (MM.GT.2) THEN
  JY=YYYY
  JM=MM+1
ELSE
  JY=YYYY-1
  JM=MM+13
ENDIF
JULDAY=INT(365.25*JY)+INT(30.6001*JM)+DD+1720995
IF (DD+31*(MM+12*YYYY).GE.IGREG) THEN
  JA=INT(0.01*JY)
  JULDAY=JULDAY+2-JA+INT(0.25*JA)
ENDIF
JULDATE=DBLE(REAL(JULDAY))+DBLE(REAL(HH)/24.)+&
&DBLE(REAL(MI)/1440.)+DBLE(REAL(SS)/86400.)
END FUNCTION juldate


!FUNCTION julday(dummydate)
!IMPLICIT NONE
!!integer, parameter :: double_k = selected_real_kind(15,307)
!!integer, parameter :: dp=double_k
!INTEGER,INTENT(IN) :: dummydate
!INTEGER :: ref_year
!REAL :: date_act,date_ref
!REAL :: juldate
!!REAL(KIND=dp) :: date_act,date_ref
!!REAL(KIND=dp) :: juldate
!CHARACTER*8 :: chdate
!INTEGER :: julday
!write(chdate,'(I8)') dummydate
!chdate(1:8)=chdate(1:4)//'0101'
!read(chdate(1:8),'(I8)') ref_year
!!date_ref=juldate(ref_year,000000)
!!date_act=juldate(dummydate,000000)
!julday=date_act-date_ref+1
!RETURN
!END FUNCTION julday





subroutine read_dyn_emiss(time_pass)
IMPLICIT NONE

INTEGER :: i,dd,hh,hgt
CHARACTER(20) :: dummy
REAL,ALLOCATABLE :: HEIGHT(:),ASH(:,:,:),SO2(:,:,:)
INTEGER,ALLOCATABLE :: DATE(:),TIME(:)
integer :: time_pass

ALLOCATE(DATE(3),TIME(8),HEIGHT(19))
ALLOCATE(ASH(3,8,19),SO2(3,8,19))

OPEN(14,FILE="emissions.txt")
DO dd=1,3
  DO hh=1,8
    DO hgt=1,19
       READ(14,*) DATE(dd),TIME(hh),HEIGHT(hgt)&
               &,ASH(dd,hh,hgt),SO2(dd,hh,hgt)
       write(*,*) "MARCUS",dd,hh,hgt,DATE(dd),TIME(hh),HEIGHT(hgt)
    END DO
  END DO
ENDDO

CLOSE(14)
time_pass=1

END subroutine


subroutine getij(ni,nj)
integer :: ni,nj
ni=154 !247
nj=305 !336
end subroutine getij

END Module module_volc_emiss_driver
