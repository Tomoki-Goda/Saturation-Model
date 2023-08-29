! integer         :: Nsize            ! size parameter, max number of final-states
! integer         :: Nfnst            ! total number of final states
! integer         :: Ngroups          ! number of groups (>1 for multi-parton scattering)
! integer         :: NfinalState(i)   ! number of final states in group i
! integer         :: process(i)       ! id of the process of group i
! real(kind(1d0)) :: Ecm              ! total center-off-mass energy
! real(kind(1d0)) :: eventWeight      ! event weight
! real(kind(1d0)) :: xsection         ! total cross section
! real(kind(1d0)) :: errest           ! statistical error estimate on total cross section
! real(kind(1d0)) :: pInst(0:3,j,i)   ! initial-state momentum j of group i
! real(kind(1d0)) :: sInst(j,i)       ! initial-state squared momentum j of group i
! real(kind(1d0)) :: pFnst(0:3,j)     ! final-state momentum j
! real(kind(1d0)) :: sFnst(j)         ! final-state squared momentum j
! real(kind(1d0)) :: matrixElement(i) ! value of matrix element for group i
! real(kind(1d0)) :: partonLumin(i)   ! value of parton luminosity for group i
! real(kind(1d0)) :: alphaStrong(i)   ! value of strong coupling for group i
! real(kind(1d0)) :: muScale(i)       ! value energy scale for group i
! character(256)  :: filename         ! name (including path) of event file
! real(kind(1d0) function mass(p(0:3))            ! mass = sqrt(abs(square))
! real(kind(1d0) function pTrans(p(0:3))          ! transverse momentum
! real(kind(1d0) function ETrans(p(0:3))          ! transverse mass
! real(kind(1d0) function rapidity(p(0:3))        ! rapidity
! real(kind(1d0) function pseudoRap(p(0:3))       ! pseudo rapidity
! real(kind(1d0) function phi(p(0:3))             ! azimuthal angle
! real(kind(1d0) function deltaPhi(p(0:3),q(0:3)) ! difference between azim.angles
! real(kind(1d0) function deltaR(p(0:3),q(0:3))   ! delta R
! integer        :: unity(i)   ! the number i 
! character(1)   :: numeral(i) ! the number i as a character
! real(kind(1d0) :: r1PI       ! the number pi
! real(kind(1d0) :: r2PI       ! the number 2*pi
! subroutine sort_big2small(l,a,n) ! sort a(1:n) from big to small or small to big
! subroutine sort_small2big(l,a,n) ! the array l contains the associated permutation
!
! type :: histo_1d_type ! 1-dimensional histogram
! type :: histo_2d_type ! 2-dimensional histogram
!
! Initialize 1-dim histogram with:
!   call hst%init( left=0d0 ,right=2d0 ,Nbins=4 )
! or alternatively:
!   call hst%init( [0.0,0.5 ,0.5,1.0 ,1.0,1.5 ,1.5,2.0] )
! Data collection:
!   call hst%collect( value ,weight )
! Write histogram file:
!   call hst%write('somename')
!
! Initialize 2-dim histogram with:
!   call hst%init( xLow=0d0 ,xUpp=2d0 ,NxBins=4 ,yLow=-1d0 ,yUpp=1d0 ,Nybins=2 )
! or alternatively:
!   call hst%init( [0.0,0.5 ,0.5,1.0 ,1.0,1.5 ,1.5,2.0] ,[-1.0,0.0 ,0.0,1.0] )

!==========================================================================
program read_event_file ;use read_events_mod ! == DO NOT TOUCH THIS LINE ==
!==========================================================================

!== DECLARATIONS ==
implicit none
!type(histo_1d_type) :: h_pT(1:Nsize),h_y(1:Nsize)
!real(kind(1d0)) :: pT(Nsize)
!integer :: Njet,ii,pTordered(Nsize)
!type(histo_1d_type) :: h_Phi1(1:1)
!type(histo_1d_type) :: h_Phi2(1:1)
!type(histo_1d_type) :: h_Phi3(1:1)
!real(kind(1d0)) :: Phi1(1)
!real(kind(1d0)) :: Phi2(1)
!real(kind(1d0)) :: Phi3(1)
type(histo_1d_type) :: h_pt2(1:1)
real(kind(1d0)) ::pt2(1)
!type(histo_1d_type) :: h_pt22(1:1)
!real(kind(1d0)) ::pt22(1)

real(kind(1d0)) :: p1Breit(0:3), p2Breit(0:3), EnergyRap, Qvirtual(0:3)
integer :: Njet,ii,Phi_ordered(Nsize)

 type(breit_type) :: obj ! Lorentz transformation to a Breit frame
!============================================
call open_file ! == DO NOT TOUCH THIS LINE ==
!============================================

!== INITIALIZATIONS ==
!Njet = 2
!do ii=1,Njet
!  call h_pT(ii)%init( left=0d0 ,right=200d0 ,Nbins=100 )
!  call h_y( ii)%init( left=-4.7d0 ,right=4.7d0 ,Nbins=100 )
!enddo
call h_pt2(1)%init( left=5d0, right=5d1 , Nbins=25 )
!call h_pt22(1)%init( left=5d0, right=5d1 , Nbins=50 )

!call obj%init(pInst(0:3,1,1)) ! prepares transformation to Breit frame of q
!===================================================================
do ;call read_event ;if(exitLoop)exit ! == DO NOT TOUCH THIS LINE ==
!===================================================================

Qvirtual = -pInst(0:3,1,1)-pFnst(0:3,3)
call obj%init(Qvirtual)

p1Breit(0:3) = obj%act(pFnst(0:3,1)) ! pBreit is the transformed version of p
p2Breit(0:3) = obj%act(pFnst(0:3,2)) ! pBreit is the transformed version of p
pt2(1)=(pTrans(p1Breit)+pTrans(p2Breit))/2
!pt22(1)=(pTrans(p1Breit+p2Breit))/2

call h_pt2(1)%collect( pt2(1) ,eventWeight )
!call h_pt22(1)%collect( pt22(1) ,eventWeight )
!call h_Phi2(1)%collect( Phi2(1) ,eventWeight )
!call h_Phi3(1)%collect( Phi3(1) ,eventWeight )

!====================================================
enddo ;call close_file ! == DO NOT TOUCH THIS LINE ==
!====================================================

!== WRITE HISTOGRAM FILES ==
!do ii=1,Njet
!  call h_pT(ii)%write('pT'//numeral(ii)//'.hst')
!  call h_y( ii)%write('y'//numeral(ii)//'.hst')
!enddo
call h_pt2(1)%write('pt21.hst')
!call h_pt22(1)%write('pt22.hst')
!call h_Phi2(1)%write('Phi2.hst')
!call h_Phi3(1)%write('Phi3.hst')

!=========================================
end program ! == DO NOT TOUCH THIS LINE ==
!=========================================

!# The 4 collumns in a histogram file represent
!#
!#   left-bin-border  right-bin-border  value  statistical-error-estimate
!#
!# Histogram files are best manipulated with standard unix commands.
!# For example, to put pT distributions without error estimates in one file:
!
!$ paste pT1.hst pT2.hst pT3.hst | awk '{print $1" "$2" "$3" "$7" "$11}' > result.hst
!
!# To get the ratio of pT1 and pT2:
!
!$ paste pT1.hst pT2.hst | awk '{print $1" "$2" "$3" "$7}' \
!$ | awk '{if($4+0!=0) print $1" "$2" "$3/$4; else print $1" "$2" ""0"}' \
!$ > result.hst
!
