module reconstruction

  use constant
  use types
  use efficiency

  implicit none

contains

  subroutine evaluate_cell(mesh,sol,k,xc,yc,dx,dy,gauss_weight,gauss_point,u)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k
    real(dp), intent(in) :: xc,yc,dx,dy
    real(dp), dimension(:), intent(in) :: gauss_weight,gauss_point
    real(dp), dimension(:), intent(inout) :: u
    integer :: p1,p2,order
    real(dp), dimension(:), allocatable :: utemp
    real(dp), dimension(:,:), allocatable :: pol2

    allocate(utemp(size(u)))

    call reconstruct1(mesh,sol,k,1,gauss_weight,pol2)
    order=size(gauss_point)
    
    u=0.0_dp
    do p1=1,order
       do p2=1,order
          call evaluate(mesh,sol,pol2,k,xc+gauss_point(p1)*dx/2.0_dp,yc+gauss_point(p2)*dy/2.0_dp,utemp)
          u=u+utemp*gauss_weight(p1)*gauss_weight(p2)/4.0_dp
       enddo
    enddo

    deallocate(utemp)

    return
  end subroutine evaluate_cell
  
  subroutine evaluate(mesh,sol,pol,k,x,y,u)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    real(dp), dimension(:,:), intent(in) :: pol
    integer, intent(in) :: k
    real(dp), intent(in) :: x,y
    real(dp), dimension(:), intent(inout) :: u
    integer :: i1,i2,i,d,isol,order
    real(dp) :: xc,yc,Kk,intk,temp
    integer, dimension(2) :: alpha
    
    xc=mesh%cell(k)%xc
    yc=mesh%cell(k)%yc
    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    u(:)=sol%val(k,:)

    select case (size(pol(:,1)))
    case (0)
       order=1
    case (2)
       order=2
    case (5)
       order=3
    case (9)
       order=4
    case (14)
       order=5
    case default
       print*,"Order not available in evaluate"
       call exit()
    end select
    d=order-1
    i=1

    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quad_monome2(mesh,alpha,k,intk)
          temp=(x-xc)**alpha(1)*(y-yc)**alpha(2)-intk
          do isol=1,sol%nvar
             u(isol)=u(isol)+pol(i,isol)*temp
          enddo
          i=i+1
       enddo
    enddo

    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quad_monome2(mesh,alpha,k,intk)
       temp=(x-xc)**alpha(1)*(y-yc)**alpha(2)-intk
       do isol=1,sol%nvar
          u(isol)=u(isol)+pol(i,isol)*temp
       enddo
       i=i+1
    enddo
    
    return
  end subroutine evaluate

  subroutine reconstruct(mesh,sol,k,order,scheme,gauss_weight,period,pol_OUT)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order,scheme
    real(dp), dimension(:), intent(in) :: gauss_weight
    real(dp), dimension(:), allocatable :: gauss_point
    logical, dimension(2), intent(in) :: period
    real(dp), dimension(:,:), allocatable, intent(inout) :: pol_OUT

    allocate (gauss_point(size(gauss_weight)))
    select case (order)
    case (1)
       gauss_point=gauss_point1
    case (2)
       gauss_point=gauss_point2
    case (3)
       gauss_point=gauss_point3
    case (4)
       gauss_point=gauss_point4
    case (5)
       gauss_point=gauss_point5
    end select
    
    select case (scheme)
    case (1)
       if (mesh%cell(k)%bound<=2) then
          if (mesh%cell(k)%bound==0.or.period(max(mesh%cell(k)%bound,1))) then
             call reconstruct1(mesh,sol,k,order,gauss_weight,pol_OUT)
          else
             if((ALL(mesh%cell(k)%stencil_type==5)).and.(size(mesh%cell(k)%stencil_type)==5.or. &
                                                       size(mesh%cell(k)%stencil_type)==14)) then
                call reconstruct1_bound(mesh,sol,k,order,period,pol_OUT)
             else
                call reconstruct1(mesh,sol,k,order,gauss_weight,pol_OUT)
             endif
          endif
       else
          if ((period(1).and.(.NOT.period(2))).or.(period(2).and.(.NOT.period(1)))) then
             call reconstruct1_bound(mesh,sol,k,order,period,pol_OUT)
          else
             call reconstruct1(mesh,sol,k,order,gauss_weight,pol_OUT)
          endif
       endif
    case (2)
       if ((mesh%cell(k)%bound<=2).and.(mesh%cell(k)%bound==0.or.period(min(2,max(mesh%cell(k)%bound,1))))) then
          call reconstruct2(mesh,sol,k,order,pol_OUT)
       else
          select case (mesh%cell(k)%bound)
          case (1)
             if (mesh%edge(mesh%cell(k)%edge4(1))%boundType=="NEUMANN".or. &
              mesh%edge(mesh%cell(k)%edge4(2))%boundType=="NEUMANN") then
                call reconstruct2(mesh,sol,k,order,pol_OUT)
             else
                call reconstruct1(mesh,sol,k,order,gauss_weight,pol_OUT)
             endif
          case (2)
             if (mesh%edge(mesh%cell(k)%edge4(3))%boundType=="NEUMANN".or. &
              mesh%edge(mesh%cell(k)%edge4(4))%boundType=="NEUMANN") then
                call reconstruct2(mesh,sol,k,order,pol_OUT)
             else
                call reconstruct1(mesh,sol,k,order,gauss_weight,pol_OUT)
             endif
          case (3)
             if (period(1).and.period(2)) then
                call reconstruct2(mesh,sol,k,order,pol_OUT)
             else
                call reconstruct1(mesh,sol,k,order,gauss_weight,pol_OUT)
             endif
          end select
       endif
    end select

    deallocate (gauss_point)

    return
  end subroutine reconstruct
  
  subroutine reconstruct1(mesh,sol,k,order,gauss_weight,pol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:), intent(in) :: gauss_weight
    real(dp), dimension(:,:), allocatable, intent(inout) :: pol
    real(dp), dimension(:,:), allocatable :: X,U
    integer :: d,Ni,Nj,i,isol
    real(dp) :: Kk,pond,dist
    real(dp), dimension(2) :: c
    character(len=99) :: id_char
    integer(16) :: id
    
    d=order-1
    
    Ni=size(mesh%cell(k)%stencil)
    Nj=d*(d+1)/2+d
    
    allocate(X(Ni,Nj),U(Ni,sol%nvar))

    Kk=mesh%cell(k)%dx*mesh%cell(k)%dy
    c(1)=mesh%cell(k)%xc
    c(2)=mesh%cell(k)%yc 
    pond=0.0_dp !1.5_dp

    write(id_char,'(99I1)')order,0,mesh%cell(k)%stencil_type
    read(id_char, '(I99)' )id

    select case (id)
    case (2055555555_dp)
       X=X_2_55555555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case (3055555555_dp)
       X=X_3_55555555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case (40555555555555555555555555_16)
       X=X_4_555555555555555555555555
       call adjust_X(mesh%cell(k)%dx,mesh%cell(k)%dy,d,X)
    case default
       call calculate_X(mesh,mesh%cell(k)%stencil,mesh%cell(k)%stencil_bound,gauss_weight,k,c,d,Ni,X)
    end select

    do i=1,Ni
       dist=((mesh%cell(mesh%cell(k)%stencil(i))%xc-c(1)+mesh%cell(k)%stencil_bound(i,1))**2+ &
            (mesh%cell(mesh%cell(k)%stencil(i))%yc-c(2)+mesh%cell(k)%stencil_bound(i,2))**2)**(pond/4.0_dp)
       do isol=1,sol%nvar
          U(i,isol)=(sol%val(mesh%cell(k)%stencil(i),isol)-sol%val(k,isol))/dist
       enddo
       X(i,:)=X(i,:)/dist
    enddo

    if (allocated(pol)) deallocate(pol)
    allocate(pol(Nj,sol%nvar))

    do isol=1,sol%nvar
       call solve(X,U(:,isol),pol(:,isol))
    enddo

    deallocate(X,U)

    return
  end subroutine reconstruct1

  subroutine reconstruct1_bound(mesh,sol,k,order,period,pol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    logical, dimension(2), intent(in) :: period
    real(dp), dimension(:,:), allocatable, intent(inout) :: pol
    integer :: d,Nj,isol
    real(dp), dimension(:), allocatable :: U0,U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,U12,U13,U14,U15,U16,U17,U18,U19
    
    d=order-1
    Nj=d*(d+1)/2+d

    if (mesh%cell(k)%neigh(3)<0.and.(.NOT.period(2))) then

       allocate (U0(sol%nvar),U1(sol%nvar),U2(sol%nvar),U3(sol%nvar),U4(sol%nvar),U5(sol%nvar))
       U0=sol%val(k,:)
       U1=sol%val(mesh%cell(k)%stencil(1),:)
       U2=sol%val(mesh%cell(k)%stencil(2),:)
       U3=sol%val(mesh%cell(k)%stencil(3),:)
       U4=sol%val(mesh%cell(k)%stencil(4),:)
       U5=sol%val(mesh%cell(k)%stencil(5),:)
       if (order>=3) then
          allocate (U6(sol%nvar),U7(sol%nvar),U8(sol%nvar),U9(sol%nvar),U10(sol%nvar))
          allocate (U11(sol%nvar),U12(sol%nvar),U13(sol%nvar),U14(sol%nvar))
          U6=sol%val(mesh%cell(k)%stencil(6),:)
          U7=sol%val(mesh%cell(k)%stencil(7),:)
          U8=sol%val(mesh%cell(k)%stencil(8),:)
          U9=sol%val(mesh%cell(k)%stencil(9),:)
          U10=sol%val(mesh%cell(k)%stencil(10),:)
          U11=sol%val(mesh%cell(k)%stencil(11),:)
          U12=sol%val(mesh%cell(k)%stencil(12),:)
          U13=sol%val(mesh%cell(k)%stencil(13),:)
          U14=sol%val(mesh%cell(k)%stencil(14),:)
          if (order>=4) then
             allocate (U15(sol%nvar),U16(sol%nvar),U17(sol%nvar),U18(sol%nvar),U19(sol%nvar))
             U15=sol%val(mesh%cell(mesh%cell(k)%stencil(13))%neigh(4),:)
             U16=sol%val(mesh%cell(mesh%cell(k)%stencil(11))%neigh(4),:)
             U17=sol%val(mesh%cell(mesh%cell(k)%stencil(10))%neigh(4),:)
             U18=sol%val(mesh%cell(mesh%cell(k)%stencil(12))%neigh(4),:)
             U19=sol%val(mesh%cell(mesh%cell(k)%stencil(14))%neigh(4),:)
          endif
       endif

       if (allocated(pol)) deallocate(pol)
       allocate(pol(Nj,sol%nvar))
    
       select case (order)
       case (1)
          do isol=1,sol%nvar
             pol=0.0_dp
          enddo
       case (2)
          do isol=1,sol%nvar
             pol(1,:)=(U4-U1+U5-U2)/(2.0_dp*mesh%cell(k)%dy)
             pol(2,:)=(U2-U1+U5-U4)/(4.0_dp*mesh%cell(k)%dx)
          enddo
       case (3)
          do isol=1,sol%nvar
             pol(1,:)=(-U13+4.0_dp*U7-3.0_dp*U6-U11+4.0_dp*U4-3.0_dp*U1-U12+4.0_dp*U5-3.0_dp*U2- &
                U14+4.0_dp*U9-3.0_dp*U8)/(8.0_dp*mesh%cell(k)%dy)
             pol(2,:)=(U9-U8+2.0_dp*U5-2.0_dp*U2+U10-U3-2.0_dp*U4+2.0_dp*U1-U7+U6)/(8.0_dp*mesh%cell(k)%dx*mesh%cell(k)%dy)
             pol(3,:)=(U6+U13-2.0_dp*U7+U1+U11-2.0_dp*U4+U12+U2-2.0_dp*U5+U8+U14-2.0_dp*U9)/(8.0_dp*mesh%cell(k)%dy**2)
             pol(4,:)=(0.0625_dp*U0-0.03125_dp*U3-0.21875_dp*U4+0.3125_dp*U2+0.09375_dp*U5+0.15625_dp*U9+ &
                0.0625_dp*U10+0.0625_dp*U11-0.125_dp*U14-0.375*U1)/mesh%cell(k)%dx
             pol(5,:)=(-0.0875_dp*U0-0.05625_dp*U3+0.00625_dp*U4-0.03125_dp*U5+0.0125_dp*U2+0.08125_dp*U9- &
                0.0875_dp*U10-0.0375_dp*U11-0.05_dp*U12+0.1_dp*U13+0.075_dp*U14+0.075_dp*U1)/(mesh%cell(k)%dx**2)
          enddo
       case (4)
          do isol=1,sol%nvar
             pol(1,:)=(-0.314743239651358*U2+0.172008002725080*U19+0.664105128006643e-1*U18-0.583606159380352e-1*U17- &
                0.664377955363000e-1*U16+0.178046562184026*U15-0.466412028185066*U14-0.267139565663571*U13-0.210996422604319*U12- &
                0.247227779419067*U11-0.183224202402931*U10+0.205311021980442*U9+0.815701163480125*U5+0.864009639241732*U4+ &
                0.989978173141377*U3-1.28295137285046*U0-0.193972050261253*U1)/mesh%cell(k)%dy
             pol(2,:)=(-0.407102582605936*U2-0.194507093464875*U19+0.681247727901954e-1*U18+0.400145507673599e-2*U17- &
                0.641233177108981e-1*U16+0.186504183308844*U15+0.229319752604949*U14-0.197308112002911*U13+0.260030920167430*U12- &
                0.276036740470970*U11-0.160058202984986e-1*U10-0.456165878341404e-1*U9+0.122162604437543*U5- &
                0.797471806615512e-1*U4+0.320116405815234e-2*U3-0.560203704003114e-1*U0+ &
                0.463122953632654*U1)/(mesh%cell(k)%dx*mesh%cell(k)%dy)
             pol(3,:)=(-0.600672971795136e-1*U2+0.656602400956103e-1*U19-0.523826846268949e-1*U18-0.647508184784925e-1*U17- & 
                0.123681338730380e-1*U16+0.638413968828155e-1*U15+0.210076391339866e-1*U14-0.390141870016351e-1*U13- &
                0.954892679027685e-3*U12+0.995816661614412e-2*U11+0.900327393053071e-2*U10-0.618406693484164e-1*U9+ &
                0.140960349309219e-1*U5-0.454710798570785e-3*U4+0.481993452160651e-1*U3+0.156511458622459*U0- &
                0.964441614990498e-1*U1)/(mesh%cell(k)%dx**2*mesh%cell(k)%dy)
             pol(4,:)=(0.222717351500523*U2-0.224226991453861*U19-0.120007275289077*U18+0.280465623264732e-1*U17+ &
                0.480538377325019e-1*U16-0.231866132941028*U15+0.538232084021817*U14+0.286140414469142*U13+0.370989450424273*U12+ &
                0.416824299424557*U11+0.387813750160175*U10-0.259730811039845*U9-0.690796652813276*U5-0.751909784824555*U4- &
                0.797562749447247*U3+0.707348123242087*U0+0.699345214899333e-1*U1)/(mesh%cell(k)%dy**2)
             pol(5,:)=(0.993543106036089e-1*U2+0.732266278486593e-1*U19+0.619316111015537e-2*U18+0.363768643584769e-3*U17- &
                0.582939246633762e-2*U16-0.739541651360616e-1*U15-0.700618406519754e-1*U14+0.729719897995759e-1*U13- &
                0.672699162807216e-1*U12+0.658148417066927e-1*U11-0.145507457357077e-2*U10-0.414696253316965e-2*U9- &
                0.343488540859204e-1*U5+0.382048017045157e-1*U4+0.291014914572674e-3*U3-0.509276115345818e-2*U0- &
                0.942615495977490e-1*U1)/(mesh%cell(k)%dx*mesh%cell(k)%dy**2)
             pol(6,:)=(-0.349975748252638e-1*U2+0.525888201387028e-1*U19+0.437795561882611e-1*U18+0.131441736525044e-1*U17+ &
                0.269795077631986e-2*U16+0.544561658358688e-1*U15-0.128234509443592*U14-0.666121013207166e-1*U13- &
                0.956863101178672e-1*U12-0.106890384319755*U11-0.102576694498062*U10+0.634897538200422e-1*U9+0.135528070708834*U5+ &
                0.150466836314015*U4+0.150515338782093*U3-0.134018430211774*U0+0.234933918338461e-2*U1)/(mesh%cell(k)%dy**3)
             pol(7,:)=(0.399622590067706*U2-0.133018066536490e-1*U19+0.847232326311024e-1*U18-0.389535588915452e-2*U17- &
                0.886185885227499e-1*U16+0.210925184344515e-1*U15-0.124754456221861*U14+0.935916091201611e-1*U13+ &
                0.484736874096435e-1*U12-0.328922638563393e-1*U11+0.155814235483962e-1*U10+0.444070570967883e-1*U9+ &
                0.133443979660482*U5-0.174734752049100*U4-0.311628470817188e-2*U3+0.545349820917557e-1*U0- &
                0.454157572469545*U1)/mesh%cell(k)%dx
             pol(8,:)=(0.143188432107778*U2-0.838486722604430e-1*U19+0.966715169420829e-1*U18+0.813022917466328e-1*U17- &
                0.153692251634171e-1*U16-0.787559112648585e-1*U15-0.882138958076972e-2*U14+0.159239723591021*U13- &
                0.223263004969039e-1*U12-0.528828665219398e-1*U11-0.752091669914036e-1*U10+0.173153874167376*U9- &
                0.394688978041894e-1*U5+0.127319023646277e-2*U4-0.134958166599649*U3-0.388232084252594*U0+ &
                0.245043652197661*U1)/(mesh%cell(k)%dx**2)
             pol(9,:)=(-0.159906633041169e-1*U2+0.341336243293083e-1*U19-0.672820419676518e-1*U18-0.272826482562946e-2*U17+ &
                0.645537771402763e-1*U16-0.286770946763028e-1*U15+0.254638050398206e-1*U14-0.472899236367768e-1*U13- &
                0.579756274707273e-1*U12+0.688886867709237e-1*U11+0.109130592967589e-1*U10+0.311022189844584e-1*U9- &
                0.548835940140250e-1*U5+0.259639868878629e-1*U4-0.218261185829528e-2*U3+0.381957074458601e-1*U0- &
                0.222050442222504e-1*U1)/(mesh%cell(k)%dx**3)
          enddo
       end select

       deallocate(U0,U1,U2,U3,U4,U5)
       if (order>=3) then 
          deallocate(U6,U7,U8,U9,U10,U11,U12,U13,U14)
          if (order>=4) deallocate(U15,U16,U17,U18,U19)
       endif

    elseif (mesh%cell(k)%neigh(4)<0.and.(.NOT.period(2))) then

       allocate (U0(sol%nvar),U1(sol%nvar),U2(sol%nvar),U3(sol%nvar),U4(sol%nvar),U5(sol%nvar))
       U0=sol%val(k,:)
       U1=sol%val(mesh%cell(k)%stencil(1),:)
       U2=sol%val(mesh%cell(k)%stencil(2),:)
       U3=sol%val(mesh%cell(k)%stencil(3),:)
       U4=sol%val(mesh%cell(k)%stencil(4),:)
       U5=sol%val(mesh%cell(k)%stencil(5),:)
       if (order>=3) then
          allocate (U6(sol%nvar),U7(sol%nvar),U8(sol%nvar),U9(sol%nvar),U10(sol%nvar))
          allocate (U11(sol%nvar),U12(sol%nvar),U13(sol%nvar),U14(sol%nvar))
          U6=sol%val(mesh%cell(k)%stencil(6),:)
          U7=sol%val(mesh%cell(k)%stencil(7),:)
          U8=sol%val(mesh%cell(k)%stencil(8),:)
          U9=sol%val(mesh%cell(k)%stencil(9),:)
          U10=sol%val(mesh%cell(k)%stencil(10),:)
          U11=sol%val(mesh%cell(k)%stencil(11),:)
          U12=sol%val(mesh%cell(k)%stencil(12),:)
          U13=sol%val(mesh%cell(k)%stencil(13),:)
          U14=sol%val(mesh%cell(k)%stencil(14),:)
          if (order>=4) then
             allocate (U15(sol%nvar),U16(sol%nvar),U17(sol%nvar),U18(sol%nvar),U19(sol%nvar))
             U15=sol%val(mesh%cell(mesh%cell(k)%stencil(13))%neigh(3),:)
             U16=sol%val(mesh%cell(mesh%cell(k)%stencil(11))%neigh(3),:)
             U17=sol%val(mesh%cell(mesh%cell(k)%stencil(10))%neigh(3),:)
             U18=sol%val(mesh%cell(mesh%cell(k)%stencil(12))%neigh(3),:)
             U19=sol%val(mesh%cell(mesh%cell(k)%stencil(14))%neigh(3),:)
          endif
       endif

       if (allocated(pol)) deallocate(pol)
       allocate(pol(Nj,sol%nvar))
    
       select case (order)
       case (1)
          do isol=1,sol%nvar
             pol=0.0_dp
          enddo
       case (2)
          do isol=1,sol%nvar
             pol(1,:)=-(U4-U1+U5-U2)/(2.0_dp*mesh%cell(k)%dy)
             pol(2,:)=(U2-U1+U5-U4)/(4.0_dp*mesh%cell(k)%dx)
          enddo
       case (3)
          do isol=1,sol%nvar
             pol(1,:)=-(-U13+4.0_dp*U7-3.0_dp*U6-U11+4.0_dp*U4-3.0_dp*U1-U12+4.0_dp*U5-3.0_dp*U2- &
                U14+4.0_dp*U9-3.0_dp*U8)/(8.0_dp*mesh%cell(k)%dy)
             pol(2,:)=-(U9-U8+2.0_dp*U5-2.0_dp*U2+U10-U3-2.0_dp*U4+2.0_dp*U1-U7+U6)/(8.0_dp*mesh%cell(k)%dx*mesh%cell(k)%dy)
             pol(3,:)=(U6+U13-2.0_dp*U7+U1+U11-2.0_dp*U4+U12+U2-2.0_dp*U5+U8+U14-2.0_dp*U9)/(8.0_dp*mesh%cell(k)%dy**2)
             pol(4,:)=(0.0625_dp*U0-0.03125_dp*U3-0.21875_dp*U4+0.3125_dp*U2+0.09375_dp*U5+0.15625_dp*U9+ &
                0.0625_dp*U10+0.0625_dp*U11-0.125_dp*U14-0.375*U1)/mesh%cell(k)%dx
             pol(5,:)=(-0.0875_dp*U0-0.05625_dp*U3+0.00625_dp*U4-0.03125_dp*U5+0.0125_dp*U2+0.08125_dp*U9- &
                0.0875_dp*U10-0.0375_dp*U11-0.05_dp*U12+0.1_dp*U13+0.075_dp*U14+0.075_dp*U1)/(mesh%cell(k)%dx**2)
          enddo
       case (4)
          do isol=1,sol%nvar
             pol(1,:)=-(-0.314743239651358*U2+0.172008002725080*U19+0.664105128006643e-1*U18-0.583606159380352e-1*U17- &
                0.664377955363000e-1*U16+0.178046562184026*U15-0.466412028185066*U14-0.267139565663571*U13-0.210996422604319*U12- &
                0.247227779419067*U11-0.183224202402931*U10+0.205311021980442*U9+0.815701163480125*U5+0.864009639241732*U4+ &
                0.989978173141377*U3-1.28295137285046*U0-0.193972050261253*U1)/mesh%cell(k)%dy
             pol(2,:)=-(-0.407102582605936*U2-0.194507093464875*U19+0.681247727901954e-1*U18+0.400145507673599e-2*U17- &
                0.641233177108981e-1*U16+0.186504183308844*U15+0.229319752604949*U14-0.197308112002911*U13+0.260030920167430*U12- &
                0.276036740470970*U11-0.160058202984986e-1*U10-0.456165878341404e-1*U9+0.122162604437543*U5- &
                0.797471806615512e-1*U4+0.320116405815234e-2*U3-0.560203704003114e-1*U0+ &
                0.463122953632654*U1)/(mesh%cell(k)%dx*mesh%cell(k)%dy)
             pol(3,:)=-(-0.600672971795136e-1*U2+0.656602400956103e-1*U19-0.523826846268949e-1*U18-0.647508184784925e-1*U17- & 
                0.123681338730380e-1*U16+0.638413968828155e-1*U15+0.210076391339866e-1*U14-0.390141870016351e-1*U13- &
                0.954892679027685e-3*U12+0.995816661614412e-2*U11+0.900327393053071e-2*U10-0.618406693484164e-1*U9+ &
                0.140960349309219e-1*U5-0.454710798570785e-3*U4+0.481993452160651e-1*U3+0.156511458622459*U0- &
                0.964441614990498e-1*U1)/(mesh%cell(k)%dx**2*mesh%cell(k)%dy)
             pol(4,:)=(0.222717351500523*U2-0.224226991453861*U19-0.120007275289077*U18+0.280465623264732e-1*U17+ &
                0.480538377325019e-1*U16-0.231866132941028*U15+0.538232084021817*U14+0.286140414469142*U13+0.370989450424273*U12+ &
                0.416824299424557*U11+0.387813750160175*U10-0.259730811039845*U9-0.690796652813276*U5-0.751909784824555*U4- &
                0.797562749447247*U3+0.707348123242087*U0+0.699345214899333e-1*U1)/(mesh%cell(k)%dy**2)
             pol(5,:)=(0.993543106036089e-1*U2+0.732266278486593e-1*U19+0.619316111015537e-2*U18+0.363768643584769e-3*U17- &
                0.582939246633762e-2*U16-0.739541651360616e-1*U15-0.700618406519754e-1*U14+0.729719897995759e-1*U13- &
                0.672699162807216e-1*U12+0.658148417066927e-1*U11-0.145507457357077e-2*U10-0.414696253316965e-2*U9- &
                0.343488540859204e-1*U5+0.382048017045157e-1*U4+0.291014914572674e-3*U3-0.509276115345818e-2*U0- &
                0.942615495977490e-1*U1)/(mesh%cell(k)%dx*mesh%cell(k)%dy**2)
             pol(6,:)=-(-0.349975748252638e-1*U2+0.525888201387028e-1*U19+0.437795561882611e-1*U18+0.131441736525044e-1*U17+ &
                0.269795077631986e-2*U16+0.544561658358688e-1*U15-0.128234509443592*U14-0.666121013207166e-1*U13- &
                0.956863101178672e-1*U12-0.106890384319755*U11-0.102576694498062*U10+0.634897538200422e-1*U9+0.135528070708834*U5+ &
                0.150466836314015*U4+0.150515338782093*U3-0.134018430211774*U0+0.234933918338461e-2*U1)/(mesh%cell(k)%dy**3)
             pol(7,:)=(0.399622590067706*U2-0.133018066536490e-1*U19+0.847232326311024e-1*U18-0.389535588915452e-2*U17- &
                0.886185885227499e-1*U16+0.210925184344515e-1*U15-0.124754456221861*U14+0.935916091201611e-1*U13+ &
                0.484736874096435e-1*U12-0.328922638563393e-1*U11+0.155814235483962e-1*U10+0.444070570967883e-1*U9+ &
                0.133443979660482*U5-0.174734752049100*U4-0.311628470817188e-2*U3+0.545349820917557e-1*U0- &
                0.454157572469545*U1)/mesh%cell(k)%dx
             pol(8,:)=(0.143188432107778*U2-0.838486722604430e-1*U19+0.966715169420829e-1*U18+0.813022917466328e-1*U17- &
                0.153692251634171e-1*U16-0.787559112648585e-1*U15-0.882138958076972e-2*U14+0.159239723591021*U13- &
                0.223263004969039e-1*U12-0.528828665219398e-1*U11-0.752091669914036e-1*U10+0.173153874167376*U9- &
                0.394688978041894e-1*U5+0.127319023646277e-2*U4-0.134958166599649*U3-0.388232084252594*U0+ &
                0.245043652197661*U1)/(mesh%cell(k)%dx**2)
             pol(9,:)=(-0.159906633041169e-1*U2+0.341336243293083e-1*U19-0.672820419676518e-1*U18-0.272826482562946e-2*U17+ &
                0.645537771402763e-1*U16-0.286770946763028e-1*U15+0.254638050398206e-1*U14-0.472899236367768e-1*U13- &
                0.579756274707273e-1*U12+0.688886867709237e-1*U11+0.109130592967589e-1*U10+0.311022189844584e-1*U9- &
                0.548835940140250e-1*U5+0.259639868878629e-1*U4-0.218261185829528e-2*U3+0.381957074458601e-1*U0- &
                0.222050442222504e-1*U1)/(mesh%cell(k)%dx**3)
          enddo
       end select

       deallocate(U0,U1,U2,U3,U4,U5)
       if (order>=3) then 
          deallocate(U6,U7,U8,U9,U10,U11,U12,U13,U14)
          if (order>=4) deallocate(U15,U16,U17,U18,U19)
       endif

    elseif (mesh%cell(k)%neigh(1)<0.and.(.NOT.period(1))) then

       allocate (U0(sol%nvar),U1(sol%nvar),U2(sol%nvar),U3(sol%nvar),U4(sol%nvar),U5(sol%nvar))
       U0=sol%val(k,:)
       U1=sol%val(mesh%cell(k)%stencil(1),:)
       U2=sol%val(mesh%cell(k)%stencil(2),:)
       U3=sol%val(mesh%cell(k)%stencil(3),:)
       U4=sol%val(mesh%cell(k)%stencil(4),:)
       U5=sol%val(mesh%cell(k)%stencil(5),:)
       if (order>=3) then
          allocate (U6(sol%nvar),U7(sol%nvar),U8(sol%nvar),U9(sol%nvar),U10(sol%nvar))
          allocate (U11(sol%nvar),U12(sol%nvar),U13(sol%nvar),U14(sol%nvar))
          U6=sol%val(mesh%cell(k)%stencil(6),:)
          U7=sol%val(mesh%cell(k)%stencil(7),:)
          U8=sol%val(mesh%cell(k)%stencil(8),:)
          U9=sol%val(mesh%cell(k)%stencil(9),:)
          U10=sol%val(mesh%cell(k)%stencil(10),:)
          U11=sol%val(mesh%cell(k)%stencil(11),:)
          U12=sol%val(mesh%cell(k)%stencil(12),:)
          U13=sol%val(mesh%cell(k)%stencil(13),:)
          U14=sol%val(mesh%cell(k)%stencil(14),:)
          if (order>=4) then
             allocate (U15(sol%nvar),U16(sol%nvar),U17(sol%nvar),U18(sol%nvar),U19(sol%nvar))
             U15=sol%val(mesh%cell(mesh%cell(k)%stencil(14))%neigh(2),:)
             U16=sol%val(mesh%cell(mesh%cell(k)%stencil(8))%neigh(2),:)
             U17=sol%val(mesh%cell(mesh%cell(k)%stencil(6))%neigh(2),:)
             U18=sol%val(mesh%cell(mesh%cell(k)%stencil(7))%neigh(2),:)
             U19=sol%val(mesh%cell(mesh%cell(k)%stencil(13))%neigh(2),:)
          endif
       endif

       if (allocated(pol)) deallocate(pol)
       allocate(pol(Nj,sol%nvar))
    
       select case (order)
       case (1)
          do isol=1,sol%nvar
             pol=0.0_dp
          enddo
       case (2)
          do isol=1,sol%nvar
             pol(1,:)=(U3-U2+U5-U4)/(4.0_dp*mesh%cell(k)%dy)
             pol(2,:)=(U4-U2+U5-U3)/(2.0_dp*mesh%cell(k)%dx)
          enddo
       case (3)
          do isol=1,sol%nvar
             pol(1,:)=(0.0625_dp*U0-0.03125_dp*U1-0.21875_dp*U5+0.3125_dp*U2+0.09375_dp*U4+0.15625_dp*U10+ &
                0.0625_dp*U6+0.0625_dp*U8-0.125_dp*U13-0.375*U3)/mesh%cell(k)%dy
             pol(2,:)=(U10-U9+2.0_dp*U4-2.0_dp*U2+U6-U1-2.0_dp*U5+2.0_dp*U3-U12+U11)/(8.0_dp*mesh%cell(k)%dx*mesh%cell(k)%dy)
             pol(3,:)=(-0.0875_dp*U0-0.05625_dp*U1+0.00625_dp*U5-0.03125_dp*U4+0.0125_dp*U2+0.08125_dp*U10- &
                0.0875_dp*U6-0.0375_dp*U8-0.05_dp*U7+0.1_dp*U14+0.075_dp*U13+0.075_dp*U3)/(mesh%cell(k)%dy**2)
             pol(4,:)=(-U14+4.0_dp*U12-3.0_dp*U11-U8+4.0_dp*U5-3.0_dp*U3-U7+4.0_dp*U4-3.0_dp*U2- &
                U13+4.0_dp*U10-3.0_dp*U9)/(8.0_dp*mesh%cell(k)%dx)
             pol(5,:)=(U11+U14-2.0_dp*U12+U3+U8-2.0_dp*U5+U2+U7-2.0_dp*U4+U9+U13-2.0_dp*U10)/(8.0_dp*mesh%cell(k)%dx**2)
          enddo
       case (4)
          do isol=1,sol%nvar
             pol(1,:)=(0.399622590067706*U2-0.133018066536490e-1*U19+0.847232326311024e-1*U18-0.389535588915452e-2*U17- &
                0.886185885227499e-1*U16+0.210925184344515e-1*U15-0.124754456221861*U13+0.935916091201611e-1*U14+ &
                0.484736874096435e-1*U7-0.328922638563393e-1*U8+0.155814235483962e-1*U6+0.444070570967883e-1*U10+ &
                0.133443979660482*U4-0.174734752049100*U5-0.311628470817188e-2*U1+0.545349820917557e-1*U0- &
                0.454157572469545*U3)/mesh%cell(k)%dx
             pol(2,:)=(-0.407102582605936*U2-0.194507093464875*U19+0.681247727901954e-1*U18+0.400145507673599e-2*U17- &
                0.641233177108981e-1*U16+0.186504183308844*U15+0.229319752604949*U13-0.197308112002911*U14+0.260030920167430*U7- &
                0.276036740470970*U8-0.160058202984986e-1*U6-0.456165878341404e-1*U10+0.122162604437543*U4- &
                0.797471806615512e-1*U5+0.320116405815234e-2*U1-0.560203704003114e-1*U0+ &
                0.463122953632654*U3)/(mesh%cell(k)%dx*mesh%cell(k)%dy)
             pol(3,:)=(0.993543106036089e-1*U2+0.732266278486593e-1*U19+0.619316111015537e-2*U18+0.363768643584769e-3*U17- &
                0.582939246633762e-2*U16-0.739541651360616e-1*U15-0.700618406519754e-1*U13+0.729719897995759e-1*U14- &
                0.672699162807216e-1*U7+0.658148417066927e-1*U8-0.145507457357077e-2*U6-0.414696253316965e-2*U10- &
                0.343488540859204e-1*U4+0.382048017045157e-1*U5+0.291014914572674e-3*U1-0.509276115345818e-2*U0- &
                0.942615495977490e-1*U3)/(mesh%cell(k)%dx*mesh%cell(k)%dy**2)
             pol(4,:)=(0.143188432107778*U2-0.838486722604430e-1*U19+0.966715169420829e-1*U18+0.813022917466328e-1*U17- &
                0.153692251634171e-1*U16-0.787559112648585e-1*U15-0.882138958076972e-2*U13+0.159239723591021*U14- &
                0.223263004969039e-1*U7-0.528828665219398e-1*U8-0.752091669914036e-1*U6+0.173153874167376*U10- &
                0.394688978041894e-1*U4+0.127319023646277e-2*U5-0.134958166599649*U1-0.388232084252594*U0+ &
                0.245043652197661*U3)/(mesh%cell(k)%dx**2)
             pol(5,:)=(-0.600672971795136e-1*U2+0.656602400956103e-1*U19-0.523826846268949e-1*U18-0.647508184784925e-1*U17- & 
                0.123681338730380e-1*U16+0.638413968828155e-1*U15+0.210076391339866e-1*U13-0.390141870016351e-1*U14- &
                0.954892679027685e-3*U7+0.995816661614412e-2*U8+0.900327393053071e-2*U6-0.618406693484164e-1*U10+ &
                0.140960349309219e-1*U4-0.454710798570785e-3*U5+0.481993452160651e-1*U1+0.156511458622459*U0- &
                0.964441614990498e-1*U3)/(mesh%cell(k)%dx**2*mesh%cell(k)%dy)
             pol(6,:)=(-0.159906633041169e-1*U2+0.341336243293083e-1*U19-0.672820419676518e-1*U18-0.272826482562946e-2*U17+ &
                0.645537771402763e-1*U16-0.286770946763028e-1*U15+0.254638050398206e-1*U13-0.472899236367768e-1*U14- &
                0.579756274707273e-1*U7+0.688886867709237e-1*U8+0.109130592967589e-1*U6+0.311022189844584e-1*U10- &
                0.548835940140250e-1*U4+0.259639868878629e-1*U5-0.218261185829528e-2*U1+0.381957074458601e-1*U0- &
                0.222050442222504e-1*U3)/(mesh%cell(k)%dx**3)
             pol(7,:)=(-0.314743239651358*U2+0.172008002725080*U19+0.664105128006643e-1*U18-0.583606159380352e-1*U17- &
                0.664377955363000e-1*U16+0.178046562184026*U15-0.466412028185066*U13-0.267139565663571*U14-0.210996422604319*U7- &
                0.247227779419067*U8-0.183224202402931*U6+0.205311021980442*U10+0.815701163480125*U4+0.864009639241732*U5+ &
                0.989978173141377*U1-1.28295137285046*U0-0.193972050261253*U3)/mesh%cell(k)%dy
             pol(8,:)=(0.222717351500523*U2-0.224226991453861*U19-0.120007275289077*U18+0.280465623264732e-1*U17+ &
                0.480538377325019e-1*U16-0.231866132941028*U15+0.538232084021817*U13+0.286140414469142*U14+0.370989450424273*U7+ &
                0.416824299424557*U8+0.387813750160175*U6-0.259730811039845*U10-0.690796652813276*U4-0.751909784824555*U5- &
                0.797562749447247*U1+0.707348123242087*U0+0.699345214899333e-1*U3)/(mesh%cell(k)%dy**2)
             pol(9,:)=(-0.349975748252638e-1*U2+0.525888201387028e-1*U19+0.437795561882611e-1*U18+0.131441736525044e-1*U17+ &
                0.269795077631986e-2*U16+0.544561658358688e-1*U15-0.128234509443592*U13-0.666121013207166e-1*U14- &
                0.956863101178672e-1*U7-0.106890384319755*U8-0.102576694498062*U6+0.634897538200422e-1*U10+0.135528070708834*U4+ &
                0.150466836314015*U5+0.150515338782093*U1-0.134018430211774*U0+0.234933918338461e-2*U3)/(mesh%cell(k)%dy**3)
          enddo
       end select

       deallocate(U0,U1,U2,U3,U4,U5)
       if (order>=3) then 
          deallocate(U6,U7,U8,U9,U10,U11,U12,U13,U14)
          if (order>=4) deallocate(U15,U16,U17,U18,U19)
       endif

    else

       allocate (U0(sol%nvar),U1(sol%nvar),U2(sol%nvar),U3(sol%nvar),U4(sol%nvar),U5(sol%nvar))
       U0=sol%val(k,:)
       U1=sol%val(mesh%cell(k)%stencil(1),:)
       U2=sol%val(mesh%cell(k)%stencil(2),:)
       U3=sol%val(mesh%cell(k)%stencil(3),:)
       U4=sol%val(mesh%cell(k)%stencil(4),:)
       U5=sol%val(mesh%cell(k)%stencil(5),:)
       if (order>=3) then
          allocate (U6(sol%nvar),U7(sol%nvar),U8(sol%nvar),U9(sol%nvar),U10(sol%nvar))
          allocate (U11(sol%nvar),U12(sol%nvar),U13(sol%nvar),U14(sol%nvar))
          U6=sol%val(mesh%cell(k)%stencil(6),:)
          U7=sol%val(mesh%cell(k)%stencil(7),:)
          U8=sol%val(mesh%cell(k)%stencil(8),:)
          U9=sol%val(mesh%cell(k)%stencil(9),:)
          U10=sol%val(mesh%cell(k)%stencil(10),:)
          U11=sol%val(mesh%cell(k)%stencil(11),:)
          U12=sol%val(mesh%cell(k)%stencil(12),:)
          U13=sol%val(mesh%cell(k)%stencil(13),:)
          U14=sol%val(mesh%cell(k)%stencil(14),:)
          if (order>=4) then
             allocate (U15(sol%nvar),U16(sol%nvar),U17(sol%nvar),U18(sol%nvar),U19(sol%nvar))
             U15=sol%val(mesh%cell(mesh%cell(k)%stencil(14))%neigh(1),:)
             U16=sol%val(mesh%cell(mesh%cell(k)%stencil(8))%neigh(1),:)
             U17=sol%val(mesh%cell(mesh%cell(k)%stencil(6))%neigh(1),:)
             U18=sol%val(mesh%cell(mesh%cell(k)%stencil(7))%neigh(1),:)
             U19=sol%val(mesh%cell(mesh%cell(k)%stencil(13))%neigh(1),:)
          endif
       endif

       if (allocated(pol)) deallocate(pol)
       allocate(pol(Nj,sol%nvar))
    
       select case (order)
       case (1)
          do isol=1,sol%nvar
             pol=0.0_dp
          enddo
       case (2)
          do isol=1,sol%nvar
             pol(1,:)=(U3-U2+U5-U4)/(4.0_dp*mesh%cell(k)%dy)
             pol(2,:)=-(U4-U2+U5-U3)/(2.0_dp*mesh%cell(k)%dx)
          enddo
       case (3)
          do isol=1,sol%nvar
             pol(1,:)=(0.0625_dp*U0-0.03125_dp*U1-0.21875_dp*U5+0.3125_dp*U2+0.09375_dp*U4+0.15625_dp*U10+ &
                0.0625_dp*U6+0.0625_dp*U8-0.125_dp*U13-0.375*U3)/mesh%cell(k)%dy
             pol(2,:)=-(U10-U9+2.0_dp*U4-2.0_dp*U2+U6-U1-2.0_dp*U5+2.0_dp*U3-U12+U11)/(8.0_dp*mesh%cell(k)%dx*mesh%cell(k)%dy)
             pol(3,:)=(-0.0875_dp*U0-0.05625_dp*U1+0.00625_dp*U5-0.03125_dp*U4+0.0125_dp*U2+0.08125_dp*U10- &
                0.0875_dp*U6-0.0375_dp*U8-0.05_dp*U7+0.1_dp*U14+0.075_dp*U13+0.075_dp*U3)/(mesh%cell(k)%dy**2)
             pol(4,:)=-(-U14+4.0_dp*U12-3.0_dp*U11-U8+4.0_dp*U5-3.0_dp*U3-U7+4.0_dp*U4-3.0_dp*U2- &
                U13+4.0_dp*U10-3.0_dp*U9)/(8.0_dp*mesh%cell(k)%dx)
             pol(5,:)=(U11+U14-2.0_dp*U12+U3+U8-2.0_dp*U5+U2+U7-2.0_dp*U4+U9+U13-2.0_dp*U10)/(8.0_dp*mesh%cell(k)%dx**2)
          enddo
       case (4)
          do isol=1,sol%nvar
             pol(1,:)=-(0.399622590067706*U2-0.133018066536490e-1*U19+0.847232326311024e-1*U18-0.389535588915452e-2*U17- &
                0.886185885227499e-1*U16+0.210925184344515e-1*U15-0.124754456221861*U13+0.935916091201611e-1*U14+ &
                0.484736874096435e-1*U7-0.328922638563393e-1*U8+0.155814235483962e-1*U6+0.444070570967883e-1*U10+ &
                0.133443979660482*U4-0.174734752049100*U5-0.311628470817188e-2*U1+0.545349820917557e-1*U0- &
                0.454157572469545*U3)/mesh%cell(k)%dx
             pol(2,:)=-(-0.407102582605936*U2-0.194507093464875*U19+0.681247727901954e-1*U18+0.400145507673599e-2*U17- &
                0.641233177108981e-1*U16+0.186504183308844*U15+0.229319752604949*U13-0.197308112002911*U14+0.260030920167430*U7- &
                0.276036740470970*U8-0.160058202984986e-1*U6-0.456165878341404e-1*U10+0.122162604437543*U4- &
                0.797471806615512e-1*U5+0.320116405815234e-2*U1-0.560203704003114e-1*U0+ &
                0.463122953632654*U3)/(mesh%cell(k)%dx*mesh%cell(k)%dy)
             pol(3,:)=-(0.993543106036089e-1*U2+0.732266278486593e-1*U19+0.619316111015537e-2*U18+0.363768643584769e-3*U17- &
                0.582939246633762e-2*U16-0.739541651360616e-1*U15-0.700618406519754e-1*U13+0.729719897995759e-1*U14- &
                0.672699162807216e-1*U7+0.658148417066927e-1*U8-0.145507457357077e-2*U6-0.414696253316965e-2*U10- &
                0.343488540859204e-1*U4+0.382048017045157e-1*U5+0.291014914572674e-3*U1-0.509276115345818e-2*U0- &
                0.942615495977490e-1*U3)/(mesh%cell(k)%dx*mesh%cell(k)%dy**2)
             pol(4,:)=(0.143188432107778*U2-0.838486722604430e-1*U19+0.966715169420829e-1*U18+0.813022917466328e-1*U17- &
                0.153692251634171e-1*U16-0.787559112648585e-1*U15-0.882138958076972e-2*U13+0.159239723591021*U14- &
                0.223263004969039e-1*U7-0.528828665219398e-1*U8-0.752091669914036e-1*U6+0.173153874167376*U10- &
                0.394688978041894e-1*U4+0.127319023646277e-2*U5-0.134958166599649*U1-0.388232084252594*U0+ &
                0.245043652197661*U3)/(mesh%cell(k)%dx**2)
             pol(5,:)=(-0.600672971795136e-1*U2+0.656602400956103e-1*U19-0.523826846268949e-1*U18-0.647508184784925e-1*U17- & 
                0.123681338730380e-1*U16+0.638413968828155e-1*U15+0.210076391339866e-1*U13-0.390141870016351e-1*U14- &
                0.954892679027685e-3*U7+0.995816661614412e-2*U8+0.900327393053071e-2*U6-0.618406693484164e-1*U10+ &
                0.140960349309219e-1*U4-0.454710798570785e-3*U5+0.481993452160651e-1*U1+0.156511458622459*U0- &
                0.964441614990498e-1*U3)/(mesh%cell(k)%dx**2*mesh%cell(k)%dy)
             pol(6,:)=-(-0.159906633041169e-1*U2+0.341336243293083e-1*U19-0.672820419676518e-1*U18-0.272826482562946e-2*U17+ &
                0.645537771402763e-1*U16-0.286770946763028e-1*U15+0.254638050398206e-1*U13-0.472899236367768e-1*U14- &
                0.579756274707273e-1*U7+0.688886867709237e-1*U8+0.109130592967589e-1*U6+0.311022189844584e-1*U10- &
                0.548835940140250e-1*U4+0.259639868878629e-1*U5-0.218261185829528e-2*U1+0.381957074458601e-1*U0- &
                0.222050442222504e-1*U3)/(mesh%cell(k)%dx**3)
             pol(7,:)=(-0.314743239651358*U2+0.172008002725080*U19+0.664105128006643e-1*U18-0.583606159380352e-1*U17- &
                0.664377955363000e-1*U16+0.178046562184026*U15-0.466412028185066*U13-0.267139565663571*U14-0.210996422604319*U7- &
                0.247227779419067*U8-0.183224202402931*U6+0.205311021980442*U10+0.815701163480125*U4+0.864009639241732*U5+ &
                0.989978173141377*U1-1.28295137285046*U0-0.193972050261253*U3)/mesh%cell(k)%dy
             pol(8,:)=(0.222717351500523*U2-0.224226991453861*U19-0.120007275289077*U18+0.280465623264732e-1*U17+ &
                0.480538377325019e-1*U16-0.231866132941028*U15+0.538232084021817*U13+0.286140414469142*U14+0.370989450424273*U7+ &
                0.416824299424557*U8+0.387813750160175*U6-0.259730811039845*U10-0.690796652813276*U4-0.751909784824555*U5- &
                0.797562749447247*U1+0.707348123242087*U0+0.699345214899333e-1*U3)/(mesh%cell(k)%dy**2)
             pol(9,:)=(-0.349975748252638e-1*U2+0.525888201387028e-1*U19+0.437795561882611e-1*U18+0.131441736525044e-1*U17+ &
                0.269795077631986e-2*U16+0.544561658358688e-1*U15-0.128234509443592*U13-0.666121013207166e-1*U14- &
                0.956863101178672e-1*U7-0.106890384319755*U8-0.102576694498062*U6+0.634897538200422e-1*U10+0.135528070708834*U4+ &
                0.150466836314015*U5+0.150515338782093*U1-0.134018430211774*U0+0.234933918338461e-2*U3)/(mesh%cell(k)%dy**3)
          enddo
       end select

       deallocate(U0,U1,U2,U3,U4,U5)
       if (order>=3) then 
          deallocate(U6,U7,U8,U9,U10,U11,U12,U13,U14)
          if (order>=4) deallocate(U15,U16,U17,U18,U19)
       endif

    endif

    return
  end subroutine reconstruct1_bound

  subroutine reconstruct1_corner(mesh,sol,k,order,pol)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:,:), allocatable, intent(inout) :: pol
    integer :: d,Nj,isol
    real(dp), dimension(:), allocatable :: U0,U1,U2,U3,U4,U5,U6,U7,U8,U9,U10,U11,U12,U13,U14,U15
    
    d=order-1
    Nj=d*(d+1)/2+d

    if (mesh%cell(k)%neigh(7)>0) then

       allocate (U0(sol%nvar),U1(sol%nvar),U2(sol%nvar),U3(sol%nvar))
       U0=sol%val(k,:)
       U1=sol%val(mesh%cell(k)%stencil(1),:)
       U2=sol%val(mesh%cell(k)%stencil(2),:)
       U3=sol%val(mesh%cell(k)%stencil(3),:)
       if (order>2) then
          allocate (U4(sol%nvar),U5(sol%nvar))
          U4=sol%val(mesh%cell(k)%stencil(4),:)
          U5=sol%val(mesh%cell(k)%stencil(5),:)
       endif
       if (order>3) then
          allocate (U6(sol%nvar),U7(sol%nvar),U8(sol%nvar),U9(sol%nvar),U10(sol%nvar))
          allocate (U11(sol%nvar),U12(sol%nvar),U13(sol%nvar),U14(sol%nvar),U15(sol%nvar))
          U6=sol%val(mesh%cell(k)%stencil(6),:)
          U7=sol%val(mesh%cell(k)%stencil(7),:)
          U8=sol%val(mesh%cell(k)%stencil(8),:)
          U9=sol%val(mesh%cell(k)%stencil(9),:)
          U10=sol%val(mesh%cell(k)%stencil(10),:)
          U11=sol%val(mesh%cell(k)%stencil(11),:)
          U12=sol%val(mesh%cell(k)%stencil(12),:)
          U13=sol%val(mesh%cell(k)%stencil(13),:)
          U14=sol%val(mesh%cell(k)%stencil(14),:)
          U15=sol%val(mesh%cell(k)%stencil(15),:)
       endif

       if (allocated(pol)) deallocate(pol)
       allocate(pol(Nj,sol%nvar))
    
       select case (order)
       case (1)
          do isol=1,sol%nvar
             pol=0.0_dp
          enddo
       case (2)
          do isol=1,sol%nvar
             pol(1,:)=(2.0_dp*(U2-U0)+U3-U1)/(3.0_dp*mesh%cell(k)%dy)
             pol(2,:)=(2.0_dp*(U0-U1)+U2-U3)/(3.0_dp*mesh%cell(k)%dx)
          enddo
       case (3)
       case (4)
       end select

       deallocate(U0,U1,U2,U3)
       if (order>2) deallocate (U4,U5)
       if (order>3) deallocate(U6,U7,U8,U9,U10,U11,U12,U13,U14,U15)

    elseif (mesh%cell(k)%neigh(8)>0) then

       allocate (U0(sol%nvar),U1(sol%nvar),U2(sol%nvar),U3(sol%nvar))
       U0=sol%val(k,:)
       U1=sol%val(mesh%cell(k)%stencil(1),:)
       U2=sol%val(mesh%cell(k)%stencil(2),:)
       U3=sol%val(mesh%cell(k)%stencil(3),:)
       if (order>2) then
          allocate (U4(sol%nvar),U5(sol%nvar))
          U4=sol%val(mesh%cell(k)%stencil(4),:)
          U5=sol%val(mesh%cell(k)%stencil(5),:)
       endif
       if (order>3) then
          allocate (U6(sol%nvar),U7(sol%nvar),U8(sol%nvar),U9(sol%nvar),U10(sol%nvar))
          allocate (U11(sol%nvar),U12(sol%nvar),U13(sol%nvar),U14(sol%nvar),U15(sol%nvar))
          U6=sol%val(mesh%cell(k)%stencil(6),:)
          U7=sol%val(mesh%cell(k)%stencil(7),:)
          U8=sol%val(mesh%cell(k)%stencil(8),:)
          U9=sol%val(mesh%cell(k)%stencil(9),:)
          U10=sol%val(mesh%cell(k)%stencil(10),:)
          U11=sol%val(mesh%cell(k)%stencil(11),:)
          U12=sol%val(mesh%cell(k)%stencil(12),:)
          U13=sol%val(mesh%cell(k)%stencil(13),:)
          U14=sol%val(mesh%cell(k)%stencil(14),:)
          U15=sol%val(mesh%cell(k)%stencil(15),:)
       endif

       if (allocated(pol)) deallocate(pol)
       allocate(pol(Nj,sol%nvar))
    
       select case (order)
       case (1)
          do isol=1,sol%nvar
             pol=0.0_dp
          enddo
       case (2)
          do isol=1,sol%nvar
             pol(1,:)=(2.0_dp*(U2-U0)+U3-U1)/(3.0_dp*mesh%cell(k)%dy)
             pol(2,:)=(2.0_dp*(U1-U0)+U3-U2)/(3.0_dp*mesh%cell(k)%dx)
          enddo
       case (3)
       case (4)
       end select

       deallocate(U0,U1,U2,U3)
       if (order>2) deallocate (U4,U5)
       if (order>3) deallocate(U6,U7,U8,U9,U10,U11,U12,U13,U14,U15)

    elseif (mesh%cell(k)%neigh(5)>0) then

       allocate (U0(sol%nvar),U1(sol%nvar),U2(sol%nvar),U3(sol%nvar))
       U0=sol%val(k,:)
       U1=sol%val(mesh%cell(k)%stencil(1),:)
       U2=sol%val(mesh%cell(k)%stencil(2),:)
       U3=sol%val(mesh%cell(k)%stencil(3),:)
       if (order>2) then
          allocate (U4(sol%nvar),U5(sol%nvar))
          U4=sol%val(mesh%cell(k)%stencil(4),:)
          U5=sol%val(mesh%cell(k)%stencil(5),:)
       endif
       if (order>3) then
          allocate (U6(sol%nvar),U7(sol%nvar),U8(sol%nvar),U9(sol%nvar),U10(sol%nvar))
          allocate (U11(sol%nvar),U12(sol%nvar),U13(sol%nvar),U14(sol%nvar),U15(sol%nvar))
          U6=sol%val(mesh%cell(k)%stencil(6),:)
          U7=sol%val(mesh%cell(k)%stencil(7),:)
          U8=sol%val(mesh%cell(k)%stencil(8),:)
          U9=sol%val(mesh%cell(k)%stencil(9),:)
          U10=sol%val(mesh%cell(k)%stencil(10),:)
          U11=sol%val(mesh%cell(k)%stencil(11),:)
          U12=sol%val(mesh%cell(k)%stencil(12),:)
          U13=sol%val(mesh%cell(k)%stencil(13),:)
          U14=sol%val(mesh%cell(k)%stencil(14),:)
          U15=sol%val(mesh%cell(k)%stencil(15),:)
       endif

       if (allocated(pol)) deallocate(pol)
       allocate(pol(Nj,sol%nvar))
    
       select case (order)
       case (1)
          do isol=1,sol%nvar
             pol=0.0_dp
          enddo
       case (2)
          do isol=1,sol%nvar
             pol(1,:)=(2.0_dp*(U0-U2)+U1-U3)/(3.0_dp*mesh%cell(k)%dy)
             pol(2,:)=(2.0_dp*(U0-U1)+U2-U3)/(3.0_dp*mesh%cell(k)%dx)
          enddo
       case (3)
       case (4)
       end select

       deallocate(U0,U1,U2,U3)
       if (order>2) deallocate (U4,U5)
       if (order>3) deallocate(U6,U7,U8,U9,U10,U11,U12,U13,U14,U15)

    else

       allocate (U0(sol%nvar),U1(sol%nvar),U2(sol%nvar),U3(sol%nvar))
       U0=sol%val(k,:)
       U1=sol%val(mesh%cell(k)%stencil(1),:)
       U2=sol%val(mesh%cell(k)%stencil(2),:)
       U3=sol%val(mesh%cell(k)%stencil(3),:)
       if (order>2) then
          allocate (U4(sol%nvar),U5(sol%nvar))
          U4=sol%val(mesh%cell(k)%stencil(4),:)
          U5=sol%val(mesh%cell(k)%stencil(5),:)
       endif
       if (order>3) then
          allocate (U6(sol%nvar),U7(sol%nvar),U8(sol%nvar),U9(sol%nvar),U10(sol%nvar))
          allocate (U11(sol%nvar),U12(sol%nvar),U13(sol%nvar),U14(sol%nvar),U15(sol%nvar))
          U6=sol%val(mesh%cell(k)%stencil(6),:)
          U7=sol%val(mesh%cell(k)%stencil(7),:)
          U8=sol%val(mesh%cell(k)%stencil(8),:)
          U9=sol%val(mesh%cell(k)%stencil(9),:)
          U10=sol%val(mesh%cell(k)%stencil(10),:)
          U11=sol%val(mesh%cell(k)%stencil(11),:)
          U12=sol%val(mesh%cell(k)%stencil(12),:)
          U13=sol%val(mesh%cell(k)%stencil(13),:)
          U14=sol%val(mesh%cell(k)%stencil(14),:)
          U15=sol%val(mesh%cell(k)%stencil(15),:)
       endif

       if (allocated(pol)) deallocate(pol)
       allocate(pol(Nj,sol%nvar))
    
       select case (order)
       case (1)
          do isol=1,sol%nvar
             pol=0.0_dp
          enddo
       case (2)
          do isol=1,sol%nvar
             pol(1,:)=(2.0_dp*(U0-U2)+U1-U3)/(3.0_dp*mesh%cell(k)%dy)
             pol(2,:)=(2.0_dp*(U1-U0)+U3-U2)/(3.0_dp*mesh%cell(k)%dx)
          enddo
       case (3)
       case (4)
       end select

       deallocate(U0,U1,U2,U3)
       if (order>2) deallocate (U4,U5)
       if (order>3) deallocate(U6,U7,U8,U9,U10,U11,U12,U13,U14,U15)

    endif

    return
  end subroutine reconstruct1_corner

  subroutine reconstruct2(mesh,sol,k,order,pol_OUT)
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(in) :: sol
    integer, intent(in) :: k,order
    real(dp), dimension(:,:), allocatable, intent(inout) :: pol_OUT
    integer, dimension(:), allocatable :: stencil
    real(dp), dimension(:,:), allocatable :: X,U,pol
    integer :: d,Ni,Nj,isol
    real(dp) :: pond,dist,dx,dy,dx2,dy2
    real(dp), dimension(2) :: c
    real(dp), dimension(16) :: Xs,Ys

    Ni=16
    d=order-1
    Nj=d*(d+1)/2+d
    
    allocate(X(Ni,Nj),U(Ni,sol%nvar),stencil(12),pol(0,0))

    stencil=abs(mesh%cell(k)%stencil2)
    c(1)=mesh%cell(k)%xc
    c(2)=mesh%cell(k)%yc
    dx=mesh%cell(k)%dx
    dy=mesh%cell(k)%dy
    pond=2.0_dp

    Xs(1)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(1,1)
    Xs(2)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(2,1)
    Xs(3)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(2,1)
    Xs(4)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(3,1)
    Xs(5)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(4,1)
    Xs(6)=c(1)-0.5_dp*dx-mesh%cell(k)%stencil2_bound(5,1)
    Xs(7)=c(1)+0.0_dp*dx-mesh%cell(k)%stencil2_bound(5,1)
    Xs(8)=c(1)+0.5_dp*dx-mesh%cell(k)%stencil2_bound(6,1)
    Xs(9)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(7,1)
    Xs(10)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(8,1)
    Xs(11)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(8,1)
    Xs(12)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(9,1)
    Xs(13)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(10,1)
    Xs(14)=c(1)+0.5_dp*dx-mesh%cell(k)%stencil2_bound(11,1)
    Xs(15)=c(1)+0.0_dp*dx-mesh%cell(k)%stencil2_bound(11,1)
    Xs(16)=c(1)-0.5_dp*dx-mesh%cell(k)%stencil2_bound(12,1)

    Ys(1)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(1,2)
    Ys(2)=c(2)+0.5_dp*dy-mesh%cell(k)%stencil2_bound(2,2)
    Ys(3)=c(2)+0.0_dp*dy-mesh%cell(k)%stencil2_bound(2,2)
    Ys(4)=c(2)-0.5_dp*dy-mesh%cell(k)%stencil2_bound(3,2)
    Ys(5)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(4,2)
    Ys(6)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(5,2)
    Ys(7)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(5,2)
    Ys(8)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(6,2)
    Ys(9)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(7,2)
    Ys(10)=c(2)-0.5_dp*dy-mesh%cell(k)%stencil2_bound(8,2)
    Ys(11)=c(2)+0.0_dp*dy-mesh%cell(k)%stencil2_bound(8,2)
    Ys(12)=c(2)+0.5_dp*dy-mesh%cell(k)%stencil2_bound(9,2)
    Ys(13)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(10,2)
    Ys(14)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(11,2)
    Ys(15)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(11,2)
    Ys(16)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(12,2)

    if (mesh%cell(k)%stencil2(2)<0.and.mesh%edge(mesh%cell(k)%edge4(1))%boundType=="NEUMANN") then
       Xs(1)=c(1)
       Xs(2)=c(1)
       Xs(3)=c(1)
       Xs(4)=c(1)
       Xs(5)=c(1)
       stencil(1)=abs(mesh%cell(k)%stencil2(12))
       stencil(2)=k
       stencil(3)=k
       stencil(4)=abs(mesh%cell(k)%stencil2(5))
    endif

    if (mesh%cell(k)%stencil2(5)<0.and.mesh%edge(mesh%cell(k)%edge4(3))%boundType=="NEUMANN") then
       Ys(5)=c(2)
       Ys(6)=c(2)
       Ys(7)=c(2)
       Ys(8)=c(2)
       Ys(9)=c(2)
       stencil(4)=abs(mesh%cell(k)%stencil2(3))
       stencil(5)=k
       stencil(6)=k
       stencil(7)=abs(mesh%cell(k)%stencil2(8))
    endif

    if (mesh%cell(k)%stencil2(8)<0.and.mesh%edge(mesh%cell(k)%edge4(2))%boundType=="NEUMANN") then
       Xs(9)=c(1)
       Xs(10)=c(1)
       Xs(11)=c(1)
       Xs(12)=c(1)
       Xs(13)=c(1)
       stencil(7)=abs(mesh%cell(k)%stencil2(6))
       stencil(8)=k
       stencil(9)=k
       stencil(10)=abs(mesh%cell(k)%stencil2(11))
    endif

    if (mesh%cell(k)%stencil2(12)<0.and.mesh%edge(mesh%cell(k)%edge4(4))%boundType=="NEUMANN") then
       Ys(13)=c(2)
       Ys(14)=c(2)
       Ys(15)=c(2)
       Ys(16)=c(2)
       Ys(1)=c(2)
       stencil(10)=abs(mesh%cell(k)%stencil2(9))
       stencil(11)=k
       stencil(12)=k
       stencil(1)=abs(mesh%cell(k)%stencil2(2))
    endif

       !Xs(1)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(1,1)
       !Xs(2)=c(1)-0.5_dp*dx-mesh%cell(k)%stencil2_bound(5,1)
       !Xs(3)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(2,1)
       !Xs(4)=c(1)+0.5_dp*dx-mesh%cell(k)%stencil2_bound(6,1)
       !Xs(5)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(2,1)
       !Xs(6)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(9,1)
       !Xs(7)=c(1)+0.0_dp*dx-mesh%cell(k)%stencil2_bound(5,1)
       !Xs(8)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(10,1)
       !Xs(9)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(3,1)
       !Xs(10)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(7,1)
       !Xs(11)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(8,1)
       !Xs(12)=c(1)+1.0_dp*dx-mesh%cell(k)%stencil2_bound(8,1)
       !Xs(13)=c(1)-1.0_dp*dx-mesh%cell(k)%stencil2_bound(4,1)
       !Xs(14)=c(1)+0.5_dp*dx-mesh%cell(k)%stencil2_bound(11,1)
       !Xs(15)=c(1)-0.0_dp*dx-mesh%cell(k)%stencil2_bound(11,1)
       !Xs(16)=c(1)-0.5_dp*dx-mesh%cell(k)%stencil2_bound(12,1)

       !Ys(1)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(1,2)
       !Ys(2)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(5,2)
       !Ys(3)=c(2)-0.0_dp*dy-mesh%cell(k)%stencil2_bound(2,2)
       !Ys(4)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(6,2)
       !Ys(5)=c(2)+0.5_dp*dy-mesh%cell(k)%stencil2_bound(2,2)
       !Ys(6)=c(2)+0.5_dp*dy-mesh%cell(k)%stencil2_bound(9,2)
       !Ys(7)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(5,2)
       !Ys(8)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(10,2)
       !Ys(9)=c(2)-0.5_dp*dy-mesh%cell(k)%stencil2_bound(3,2)
       !Ys(10)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(7,2)
       !Ys(11)=c(2)+0.0_dp*dy-mesh%cell(k)%stencil2_bound(8,2)
       !Ys(12)=c(2)-0.5_dp*dy-mesh%cell(k)%stencil2_bound(8,2)
       !Ys(13)=c(2)-1.0_dp*dy-mesh%cell(k)%stencil2_bound(4,2)
       !Ys(14)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(11,2)
       !Ys(15)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(11,2)
       !Ys(16)=c(2)+1.0_dp*dy-mesh%cell(k)%stencil2_bound(12,2)

    select case (order)
    !case (4)
       !X=X_4
       !call adjust_X(dx,dy,d,X)
    case (3)
       X=X_3
       call adjust_X(dx,dy,d,X)
    case (2)
       X=X_2
       call adjust_X(dx,dx,d,X)
    case default
       !call calculate_X(mesh,stencil,mesh%cell(k)%stencil2_bound,gauss_weight,k,c,d,Ni,X)
       if(allocated(pol_OUT))deallocate(pol_OUT)
       allocate(pol_OUT(Nj,sol%nvar))
       pol_OUT(:,:)=0.0_dp
       return
    end select
    
    dx2=0.5*dx
    dy2=0.5*dy
    !dist=((dx)**2+(dy)**2)**(pond/2.0_dp)
    dist=0.1_dp
    call evaluate(mesh,sol,pol,stencil(1),Xs(1),Ys(1),U(1,:))
    U(1,:)=(U(1,:)-sol%val(k,:))/dist
    X(1,:)=X(1,:)/dist
    call evaluate(mesh,sol,pol,stencil(4),Xs(5),Ys(5),U(5,:))
    U(5,:)=(U(5,:)-sol%val(k,:))/dist
    X(5,:)=X(5,:)/dist
    call evaluate(mesh,sol,pol,stencil(7),Xs(9),Ys(9),U(9,:))
    U(9,:)=(U(9,:)-sol%val(k,:))/dist
    X(9,:)=X(9,:)/dist
    call evaluate(mesh,sol,pol,stencil(10),Xs(13),Ys(13),U(13,:))
    U(13,:)=(U(13,:)-sol%val(k,:))/dist
    X(13,:)=X(13,:)/dist

    !dist=((1.0_dp*dx)**2+(0.5*dy)**2)**(pond/2.0_dp)
    dist=10.0_dp
    call evaluate(mesh,sol,pol,stencil(2),Xs(2),Ys(2),U(2,:))
    U(2,:)=(U(2,:)-sol%val(k,:))/dist
    X(2,:)=X(2,:)/dist
    call evaluate(mesh,sol,pol,stencil(3),Xs(4),Ys(4),U(4,:))
    U(4,:)=(U(4,:)-sol%val(k,:))/dist
    X(4,:)=X(4,:)/dist
    call evaluate(mesh,sol,pol,stencil(8),Xs(10),Ys(10),U(10,:))
    U(10,:)=(U(10,:)-sol%val(k,:))/dist
    X(10,:)=X(10,:)/dist
    call evaluate(mesh,sol,pol,stencil(9),Xs(12),Ys(12),U(12,:))
    U(12,:)=(U(12,:)-sol%val(k,:))/dist
    X(12,:)=X(12,:)/dist

    !dist=((0.5_dp*dx)**2+(1.0*dy)**2)**(pond/2.0_dp)
    dist=10.0_dp
    call evaluate(mesh,sol,pol,stencil(5),Xs(6),Ys(6),U(6,:))
    U(6,:)=(U(6,:)-sol%val(k,:))/dist
    X(6,:)=X(6,:)/dist
    call evaluate(mesh,sol,pol,stencil(6),Xs(8),Ys(8),U(8,:))
    U(8,:)=(U(8,:)-sol%val(k,:))/dist
    X(8,:)=X(8,:)/dist
    call evaluate(mesh,sol,pol,stencil(11),Xs(14),Ys(14),U(14,:))
    U(14,:)=(U(14,:)-sol%val(k,:))/dist
    X(14,:)=X(14,:)/dist
    call evaluate(mesh,sol,pol,stencil(12),Xs(16),Ys(16),U(16,:))
    U(16,:)=(U(16,:)-sol%val(k,:))/dist
    X(16,:)=X(16,:)/dist

    dist=1.0_dp
    call evaluate(mesh,sol,pol,stencil(2),Xs(3),Ys(3),U(3,:))
    U(3,:)=(U(3,:)-sol%val(k,:))/dist
    X(3,:)=X(3,:)/dist
    call evaluate(mesh,sol,pol,stencil(5),Xs(7),Ys(7),U(7,:))
    U(7,:)=(U(7,:)-sol%val(k,:))/dist
    X(7,:)=X(7,:)/dist
    call evaluate(mesh,sol,pol,stencil(8),Xs(11),Ys(11),U(11,:))
    U(11,:)=(U(11,:)-sol%val(k,:))/dist
    X(11,:)=X(11,:)/dist
    call evaluate(mesh,sol,pol,stencil(11),Xs(15),Ys(15),U(15,:))
    U(15,:)=(U(15,:)-sol%val(k,:))/dist
    X(15,:)=X(15,:)/dist

    if (allocated(pol_OUT)) deallocate(pol_OUT)
    allocate(pol_OUT(Nj,sol%nvar))
    
    do isol=1,sol%nvar
       call solve(X,U(:,isol),pol_OUT(:,isol))
    enddo

    deallocate(stencil,X,U)
    
    return
  end subroutine reconstruct2

  subroutine adjust_X(dx,dy,d,X)
    real(dp), intent(in) :: dx,dy
    integer, intent(in) :: d
    real(dp), dimension(:,:), intent(inout) :: X
    integer :: i1,i2,i

    i=1
    do i2=1,d
       do i1=0,d-i2
          X(:,i)=X(:,i)*(dx**i1)*(dy**i2)
          i=i+1
       enddo
    enddo
    do i1=1,d
       X(:,i)=X(:,i)*(dx**i1)
       i=i+1
    enddo

  end subroutine adjust_X

  subroutine calculate_X(mesh,stencil,stencil_bound,gauss_weight,k,c,d,Ni,X)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(:), intent(in) :: stencil
    real(dp), dimension(:,:), intent(in) :: stencil_bound
    real(dp), dimension(:), intent(in) :: gauss_weight
    integer, intent(in) :: k,d,Ni
    real(dp), dimension(2), intent(in) :: c
    real(dp), dimension(:,:), intent(inout) :: X
    integer :: i1,i2,i,j
    integer, dimension(2) :: alpha
    real(dp) :: intj,intk,Kj

    i=1
    do i2=1,d
       do i1=0,d-i2
          alpha(1)=i1
          alpha(2)=i2
          call quad_monome2(mesh,alpha,k,intk)
          do j=1,Ni
             Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
             call quad_monome(mesh,c-stencil_bound(j,:),alpha,stencil(j),gauss_weight,intj)
             X(j,i)=intj-intk
          enddo
          i=i+1
       enddo
    enddo
    do i1=1,d
       alpha(1)=i1
       alpha(2)=0
       call quad_monome2(mesh,alpha,k,intk)
       do j=1,Ni
          Kj=mesh%cell(stencil(j))%dx*mesh%cell(stencil(j))%dy
          call quad_monome(mesh,c-stencil_bound(j,:),alpha,stencil(j),gauss_weight,intj)
          X(j,i)=intj-intk
       enddo
       i=i+1
    enddo

  end subroutine calculate_X
  
  subroutine Nequa(d,N)
    integer, intent(in) :: d
    integer, intent(out) :: N

    N=(d+1)*(d+2)/2-1

    return
  end subroutine Nequa

  subroutine couronne(mesh,period,stencil,stencil_type,n)
    type(meshStruct), intent(in) :: mesh
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), allocatable, intent(inout) :: stencil,stencil_type
    integer, intent(in) :: n
    type(cellStruct) :: cell
    integer, dimension(:), allocatable :: stencil2,stencil2_type
    integer :: i,j,k,neigh,s

    s=3*(2*n+1)**2
    allocate(stencil2(s),stencil2_type(s))
    stencil2=-1
    stencil2(1:size(stencil))=stencil
    stencil2_type=0
    stencil2_type(1:size(stencil_type))=stencil_type
    k=size(stencil)

    do i=1,size(stencil)
       cell=mesh%cell(stencil(i))
       do j=1,size(cell%edge)
          neigh=cell%neigh(j)
          if (neigh>0.or.period(mesh%edge(cell%edge(j))%dir)) then
             if (all(stencil2/=abs(neigh))) then
                k=k+1
                stencil2(k)=abs(neigh)
                stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
             endif
          endif
       enddo

       neigh=cell%corner_cell(1)
       if (neigh>0) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(1))%dir).and.mesh%cell(abs(neigh))%yc<mesh%cell(stencil(i))%yc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(3))%dir).and.mesh%cell(abs(neigh))%xc<mesh%cell(stencil(i))%xc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(1))%dir).and.period(mesh%edge(cell%edge(3))%dir)) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       endif

       neigh=cell%corner_cell(2)
       if (neigh>0) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(2))%dir).and.mesh%cell(abs(neigh))%yc<mesh%cell(stencil(i))%yc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(3))%dir).and.mesh%cell(abs(neigh))%xc>mesh%cell(stencil(i))%xc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(2))%dir).and.period(mesh%edge(cell%edge(3))%dir)) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       endif

       neigh=cell%corner_cell(3)
       if (neigh>0) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(1))%dir).and.mesh%cell(abs(neigh))%yc>mesh%cell(stencil(i))%yc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(4))%dir).and.mesh%cell(abs(neigh))%xc<mesh%cell(stencil(i))%xc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(1))%dir).and.period(mesh%edge(cell%edge(4))%dir)) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       endif

       neigh=cell%corner_cell(4)
       if (neigh>0) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(2))%dir).and.mesh%cell(abs(neigh))%yc>mesh%cell(stencil(i))%yc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(4))%dir).and.mesh%cell(abs(neigh))%xc>mesh%cell(stencil(i))%xc) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       else if (period(mesh%edge(cell%edge(2))%dir).and.period(mesh%edge(cell%edge(4))%dir)) then
          if (all(stencil2/=abs(neigh))) then
             k=k+1
             stencil2(k)=abs(neigh)
             stencil2_type(k)=mesh%cell(stencil(1))%level-mesh%cell(abs(neigh))%level+5
          endif
       endif

    enddo

    deallocate(stencil,stencil_type)
    allocate(stencil(k),stencil_type(k))

    stencil=stencil2(1:k)
    stencil_type=stencil2_type(1:k)

    deallocate(stencil2,stencil2_type)
    
    return
  end subroutine couronne

  subroutine buildStencil(mesh,k,N,order,period)
    type(meshStruct), intent(inout) :: mesh
    integer, intent(in) :: k,N,order
    logical, dimension(2), intent(in) :: period
    integer, dimension(:), allocatable :: stencil2,stencil2_type
    integer :: i,j,s

    allocate(stencil2(1),stencil2_type(1))
    stencil2(1)=k
    stencil2_type(1)=0
    i=1
    s=ceiling(size_stencil(order)*N)
    
    do while (size(stencil2)<s)
       call couronne(mesh,period,stencil2,stencil2_type,i)
       i=i+1
    enddo

    s=size(stencil2)-1
    allocate(mesh%cell(k)%stencil(s),mesh%cell(k)%stencil_type(s),mesh%cell(k)%stencil_bound(s,2))
    j=1
    do i=1,s+1
       if (stencil2(i)/=k) then
          mesh%cell(k)%stencil(j)=stencil2(i)
          mesh%cell(k)%stencil_type(j)=stencil2_type(i)
          j=j+1
       endif
    enddo
    mesh%cell(k)%stencil_bound=0.0_dp
    do i=1,s
       if (mesh%cell(k)%xc-mesh%cell(abs(mesh%cell(k)%stencil(i)))%xc>mesh%Lx/2.0_dp) then
          mesh%cell(k)%stencil_bound(i,1)=mesh%cell(k)%stencil_bound(i,1)+mesh%Lx
       endif
       if (mesh%cell(k)%xc-mesh%cell(abs(mesh%cell(k)%stencil(i)))%xc<-mesh%Lx/2.0_dp) then
          mesh%cell(k)%stencil_bound(i,1)=mesh%cell(k)%stencil_bound(i,1)-mesh%Lx
       endif
       if (mesh%cell(k)%yc-mesh%cell(abs(mesh%cell(k)%stencil(i)))%yc>mesh%Ly/2.0_dp) then
          mesh%cell(k)%stencil_bound(i,2)=mesh%cell(k)%stencil_bound(i,2)+mesh%Ly
       endif
       if (mesh%cell(k)%yc-mesh%cell(abs(mesh%cell(k)%stencil(i)))%yc<-mesh%Ly/2.0_dp) then
          mesh%cell(k)%stencil_bound(i,2)=mesh%cell(k)%stencil_bound(i,2)-mesh%Ly
       endif
    enddo
    
    deallocate(stencil2,stencil2_type)

    return
  end subroutine buildStencil

  subroutine buildStencil2(mesh,k)
    type(meshStruct), intent(inout) :: mesh
    integer, intent(in) :: k
    integer :: i

    i=1

    if (mesh%cell(k)%level<mesh%cell(abs(mesh%cell(k)%neigh(i)))%level) then
       mesh%cell(k)%stencil2(2)=mesh%cell(k)%neigh(i+1)
       mesh%cell(k)%stencil2(3)=mesh%cell(k)%neigh(i)
       i=i+2
    else
       mesh%cell(k)%stencil2(2)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(3)=mesh%cell(k)%neigh(i)
       i=i+1
    endif

    if (mesh%cell(k)%level<mesh%cell(abs(mesh%cell(k)%neigh(i)))%level) then
       mesh%cell(k)%stencil2(8)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(9)=mesh%cell(k)%neigh(i+1)
       i=i+2
    else
       mesh%cell(k)%stencil2(8)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(9)=mesh%cell(k)%neigh(i)
       i=i+1
    endif

    if (mesh%cell(k)%level<mesh%cell(abs(mesh%cell(k)%neigh(i)))%level) then
       mesh%cell(k)%stencil2(5)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(6)=mesh%cell(k)%neigh(i+1)
       i=i+2
    else
       mesh%cell(k)%stencil2(5)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(6)=mesh%cell(k)%neigh(i)
       i=i+1
    endif

    if (mesh%cell(k)%level<mesh%cell(abs(mesh%cell(k)%neigh(i)))%level) then
       mesh%cell(k)%stencil2(11)=mesh%cell(k)%neigh(i+1)
       mesh%cell(k)%stencil2(12)=mesh%cell(k)%neigh(i)
       i=i+2
    else
       mesh%cell(k)%stencil2(11)=mesh%cell(k)%neigh(i)
       mesh%cell(k)%stencil2(12)=mesh%cell(k)%neigh(i)
       i=i+1
    endif

    mesh%cell(k)%stencil2(1)=mesh%cell(k)%corner_cell(3)
    mesh%cell(k)%stencil2(4)=mesh%cell(k)%corner_cell(1)
    mesh%cell(k)%stencil2(7)=mesh%cell(k)%corner_cell(2)
    mesh%cell(k)%stencil2(10)=mesh%cell(k)%corner_cell(4)

    mesh%cell(k)%stencil2_bound=0
    do i=1,12
       if (mesh%cell(k)%xc-mesh%cell(abs(mesh%cell(k)%stencil2(i)))%xc>mesh%Lx/2.0_dp) then
          mesh%cell(k)%stencil2_bound(i,1)=mesh%cell(k)%stencil2_bound(i,1)+mesh%Lx
       endif
       if (mesh%cell(k)%xc-mesh%cell(abs(mesh%cell(k)%stencil2(i)))%xc<-mesh%Lx/2.0_dp)  then
          mesh%cell(k)%stencil2_bound(i,1)=mesh%cell(k)%stencil2_bound(i,1)-mesh%Lx
       endif
       if (mesh%cell(k)%yc-mesh%cell(abs(mesh%cell(k)%stencil2(i)))%yc>mesh%Ly/2.0_dp) then
          mesh%cell(k)%stencil2_bound(i,2)=mesh%cell(k)%stencil2_bound(i,2)+mesh%Ly
       endif
       if (mesh%cell(k)%yc-mesh%cell(abs(mesh%cell(k)%stencil2(i)))%yc<-mesh%Ly/2.0_dp) then
          mesh%cell(k)%stencil2_bound(i,2)=mesh%cell(k)%stencil2_bound(i,2)-mesh%Ly
       endif
    enddo

    return
  end subroutine buildStencil2

  subroutine polynomialProduct(x,y,c,alpha,s)
    real(dp), intent(in) :: x,y
    real(dp), dimension(2), intent(in) :: c
    integer, dimension(2), intent(in) :: alpha
    real(dp), intent(out) :: s

    s=(x-c(1))**alpha(1)*(y-c(2))**alpha(2)

    return
  end subroutine polynomialProduct

  subroutine solve(X,U,R)
    real(dp), dimension(:,:), intent(in) :: X
    real(dp), dimension(:), intent(in) :: U
    real(dp), dimension(:), intent(inout) :: R
    real(dp), dimension(:,:), allocatable :: A
    real(dp), dimension(:), allocatable :: b

    allocate(A(size(R),size(R)),b(size(R)))

    A=matmul(transpose(X),X)
    b=matmul(transpose(X),U)
    
    call cholesky(A,b,R)
    !call QR(X,b,R)

    deallocate(A,b)

    return
  end subroutine solve

  subroutine cholesky(A,b,x)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(inout) :: x
    real(dp), dimension(:,:), allocatable :: L
    integer :: i,j,k,n
    real(dp), dimension(:), allocatable :: x1

    n=size(b)
    allocate(L(n,n),x1(n))
    
    L=0.0_dp
    do j=1,n
       do i=1,j-1
          L(i,j)=A(i,j)
          do k=1,i-1
             L(i,j)=L(i,j)-L(k,i)*L(k,j)
          enddo
          L(i,j)=L(i,j)/L(i,i)
       enddo
       L(j,j)=A(j,j)
       do k=1,j-1
          L(j,j)=L(j,j)-L(k,j)*L(k,j)
       enddo
       L(j,j)=sqrt(L(j,j))
    enddo

    do j=1,n
       x1(j)=b(j)
       do k=1,j-1
          x1(j)=x1(j)-x1(k)*L(k,j)
       enddo
       x1(j)=x1(j)/L(j,j)
    enddo

    do j=1,n
       x(n-j+1)=x1(n-j+1)
       do k=1,j-1
          x(n-j+1)=x(n-j+1)-x(n-k+1)*L(n-j+1,n-k+1)
       enddo
       x(n-j+1)=x(n-j+1)/L(n-j+1,n-j+1)
    enddo

    deallocate(L,x1)
    
    return
  end subroutine cholesky

  subroutine QR(A,b,x)
    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(:), intent(in) :: b
    real(dp), dimension(:), intent(inout) :: x
    real(dp), dimension(:,:), allocatable :: R
    real(dp), dimension(:), allocatable :: w,w2
    integer :: k,k2,j,m,n
    real(dp) :: norm,s,u1,tau

    m=size(A(:,1))
    n=size(A(1,:))
    allocate(R(m,n),w(m),w2(n))
    R=A

    do j=1,n
       norm=0
       do k=j,m
          norm=norm+R(k,j)**2
       enddo
       norm=sqrt(norm)
       s=sign(1.0_dp,-R(j,j))
       u1=R(j,j)-s*norm
       w(j:m)=R(j:m,j)/u1
       w(j)=1
       tau=-s*u1/norm
       w2=0.0_dp
       do k2=1,n
          do k=j,m
             w2(k2)=w2(k2)+w(k)*R(k,k2)
          enddo
       enddo
       do k=j,m
          do k2=1,n
             R(k,k2)=R(k,k2)-tau*w(k)*w2(k2)
          enddo
       enddo
    enddo
    
    do j=1,n
       w2(j)=b(j)
       do k=1,j-1
          w2(j)=w2(j)-w2(k)*R(k,j)
       enddo
       w2(j)=w2(j)/R(j,j)
    enddo
    
    do j=1,n
       x(n-j+1)=w2(n-j+1)
       do k=1,j-1
          x(n-j+1)=x(n-j+1)-x(n-k+1)*R(n-j+1,n-k+1)
       enddo
       x(n-j+1)=x(n-j+1)/R(n-j+1,n-j+1)
    enddo

    deallocate(R,w,w2)

    return
  end subroutine QR

  subroutine quad_monome(mesh,c,alpha,k,gauss_weight,int)
    type(meshStruct), intent(in) :: mesh
    real(dp), dimension(2), intent(in) :: c
    integer, dimension(2), intent(in) :: alpha
    integer, intent(in) :: k
    real(dp), dimension(:), intent(in) :: gauss_weight
    real(dp), intent(out) :: int
    integer :: p1,p2
    real(dp) :: s,s2

    if (size(gauss_weight)==2) then
       int=0.0_dp
       do p1=1,size(gauss_weight)
          s2=0.0_dp
          do p2=1,size(gauss_weight)
             call polynomialProduct(mesh%cell(k)%X_gauss2(p1),mesh%cell(k)%Y_gauss2(p2),c,alpha,s)
             s2=s2+s*gauss_weight(p2)
          enddo
          int=int+gauss_weight(p1)*s2
       enddo
       int=int*0.25_dp
    else
       int=0.0_dp
       do p1=1,size(gauss_weight)
          s2=0.0_dp
          do p2=1,size(gauss_weight)
             call polynomialProduct(mesh%cell(k)%X_gauss(p1),mesh%cell(k)%Y_gauss(p2),c,alpha,s)
             s2=s2+s*gauss_weight(p2)
          enddo
          int=int+gauss_weight(p1)*s2
       enddo
       int=int*0.25_dp
    endif

    return
  end subroutine quad_monome

   subroutine quad_monome2(mesh,alpha,k,int)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(2), intent(in) :: alpha
    integer, intent(in) :: k
    real(dp), intent(out) :: int

    if (mod(alpha(1),2)==0.and.mod(alpha(2),2)==0) then
       int=mesh%cell(k)%dx**(alpha(1))*mesh%cell(k)%dy**(alpha(2))/((alpha(1)+1)*(alpha(2)+1)*2**(alpha(1)+alpha(2)))
    else
       int=0.0_dp
    endif

    return
  end subroutine quad_monome2
  
end module reconstruction
