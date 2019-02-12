module inout

  use constant
  use types
  use efficiency
  use phys
  use ISO_C_BINDING
  
  implicit none

contains

  subroutine get_config(config_file)
    character(len=20), intent(out) :: config_file
    logical :: file_exists

    call getarg(1,config_file)
    INQUIRE(file="config/"//trim(config_file), EXIST=file_exists)
    if (.not.file_exists) then
       print*,"Configuration file ",trim(config_file)," doesn't exist"
       call exit()
    endif

    return
  end subroutine get_config
  
  subroutine init(config_file,test_case,xL,xR,yL,yR,level,nvar,cfl,tf,fs,namefile,verbosity,sol, &
       str_equa,str_flux,str_time_scheme,order,L_str_criteria,L_var_criteria,L_eps, &
       gauss_point,gauss_weight,bool_AMR,fn_adapt,f_adapt,recursivity)
    character(len=20), intent(out) :: config_file,test_case,namefile,str_equa,str_flux,str_time_scheme
    real(dp), intent(out) :: xL,xR,yL,yR,cfl,tf
    integer, intent(out) :: level,nvar,fs,verbosity,order,f_adapt,recursivity
    type(solStruct), intent(out) :: sol
    character(len=20), dimension(:), allocatable, intent(out) :: L_str_criteria
    integer, dimension(:), allocatable, intent(out) :: L_var_criteria
    real(dp), dimension(:), allocatable, intent(out) :: L_eps
    real(dp), dimension(:), allocatable, intent(out) :: gauss_point,gauss_weight
    character(len=20), intent(out) :: fn_adapt
    logical, intent(out) :: bool_AMR
    character(len=20) :: blank
    integer :: i,ncriteria

    config_file="config/"//trim(config_file)
    open(11,file=config_file,form="formatted")
    
    read(11,*)blank
    read(11,*)test_case
    read(11,*)blank
    read(11,*)xL
    read(11,*)xR
    read(11,*)yL
    read(11,*)yR
    read(11,*)blank
    read(11,*)level
    read(11,*)cfl
    read(11,*)tf
    read(11,*)str_equa
    read(11,*)str_flux
    read(11,*)str_time_scheme
    read(11,*)order
    read(11,*)blank
    read(11,*)ncriteria
    allocate(L_str_criteria(ncriteria),L_var_criteria(ncriteria),L_eps(ncriteria))
    do i=1,ncriteria
       read(11,*)L_str_criteria(i),L_var_criteria(i),L_eps(i)
    enddo
    read(11,*)blank
    read(11,*)fs
    read(11,*)namefile
    read(11,*)verbosity
    read(11,*)nvar
    allocate(sol%name(nvar))
    do i=1,nvar
       read(11,*)sol%name(i)
    enddo
    read(11,*)sol%nsolUser
    allocate(sol%var_user(sol%nsolUser),sol%name_user(sol%nsolUser))
    allocate(sol%conserv_var(nvar,2))
    do i=1,sol%nsolUser
       read(11,*)sol%var_user(i)
    enddo
    read(11,*)blank
    read(11,*)bool_AMR
    read(11,*)fn_adapt
    read(11,*)f_adapt
    read(11,*)recursivity
    
    close(11)
    
    sol%nvar=nvar

    allocate(gauss_point(order),gauss_weight(order))
    select case (order)
    case (1)
       gauss_point=gauss_point1
       gauss_weight=gauss_weight1
    case(2)
       gauss_point=gauss_point2
       gauss_weight=gauss_weight2
    case(3)
       gauss_point=gauss_point3
       gauss_weight=gauss_weight3
    case(4)
       gauss_point=gauss_point4
       gauss_weight=gauss_weight4
    case(5)
       gauss_point=gauss_point5
       gauss_weight=gauss_weight5
    case(6)
       gauss_point=gauss_point6
       gauss_weight=gauss_weight6
    case default
       print*,"Space order too high, no good enough quadrature implemented"
       call exit()
    end select

    return
  end subroutine init
  
  subroutine IC(IC_func,mesh,sol,order)
    procedure (sub_IC), pointer, intent(in) :: IC_func
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(inout) :: sol
    integer, intent(in) :: order
    real(dp), dimension(:), allocatable :: U0,S
    integer :: k,p1,p2
    real(dp) :: Xg,Yg

    allocate(U0(sol%nvar),S(sol%nvar))
    
    do k=1,mesh%nc
       U0=0.0_dp
       do p1=1,6
          do p2=1,6
             Xg=mesh%cell(k)%xc+gauss_point6(p1)*mesh%cell(k)%dx/2.0_dp
             Yg=mesh%cell(k)%yc+gauss_point6(p2)*mesh%cell(k)%dx/2.0_dp
             call IC_func(Xg,Yg,S)
             U0=U0+S*gauss_weight6(p1)*gauss_weight6(p2)/4.0_dp
          enddo
       enddo
       sol%val(k,:)=U0
    enddo

    if (order>6) then
       print*,"The order of initial condition's quadrature is to low"
       call exit()
    endif

    deallocate(U0,S)
    
    return
  end subroutine IC
  
  subroutine buildMesh(xL,xR,yL,yR,nx,ny,gauss_point,order,mesh,sol)
    real(dp), intent(in) :: xL,xR,yL,yR
    integer, intent(in) :: nx,ny,order
    real(dp), dimension(:), intent(in) :: gauss_point
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    integer :: i,j,k,p,dir
    real(dp) :: dx,dy,a,b,c,center,diff,xc,yc
    type(edgeStruct) :: edge

    dx=(xR-xL)/nx
    dy=(yR-yL)/ny
    mesh%np=(nx+1)*(ny+1)
    mesh%ne=ny*(nx+1)+nx*(ny+1)
    mesh%nc=nx*ny

    allocate(mesh%node(mesh%np),mesh%cell(mesh%nc),mesh%edge(mesh%ne))
    allocate(sol%val(mesh%nc,sol%nvar),sol%user(mesh%nc,sol%nsolUser))

    !Initialisation des points
    
    do j=0,ny
       do i=0,nx
          k=j*(nx+1)+i+1
          mesh%node(k)%x=xL+i*dx
          mesh%node(k)%y=yL+j*dy
       enddo
    enddo

    !Initialisation des cells
    
    do j=1,ny
       do i=1,nx
          k=(j-1)*nx+i
          mesh%cell(k)%dx=dx
          mesh%cell(k)%dy=dy
          mesh%cell(k)%xc=xL+i*dx-dx/2.0_dp
          mesh%cell(k)%yc=yL+j*dy-dy/2.0_dp

          allocate(mesh%cell(k)%node(4),mesh%cell(k)%neigh(8),mesh%cell(k)%edge(4))
          allocate(mesh%cell(k)%X_gauss(size(gauss_point)),mesh%cell(k)%Y_gauss(size(gauss_point)))
          allocate(mesh%cell(k)%polMax(order*(order-1)/2+order-1,sol%nvar))

          mesh%cell(k)%corner(1)=(j-1)*(nx+1)+i
          mesh%cell(k)%corner(2)=(j-1)*(nx+1)+i+1
          mesh%cell(k)%corner(3)=j*(nx+1)+i
          mesh%cell(k)%corner(4)=j*(nx+1)+i+1
          mesh%cell(k)%node(1)=mesh%cell(k)%corner(1)
          mesh%cell(k)%node(2)=mesh%cell(k)%corner(2)
          mesh%cell(k)%node(3)=mesh%cell(k)%corner(4)
          mesh%cell(k)%node(4)=mesh%cell(k)%corner(3)

          mesh%cell(k)%edge(1)=k+j-1
          mesh%cell(k)%edge(2)=(nx+1)*ny+j+(i-1)*(ny+1)
          mesh%cell(k)%edge(3)=k+j
          mesh%cell(k)%edge(4)=(nx+1)*ny+j+(i-1)*(ny+1)+1
          
          mesh%cell(k)%neigh(1)=k-nx-1
          mesh%cell(k)%neigh(2)=k-1
          mesh%cell(k)%neigh(3)=k+nx-1
          mesh%cell(k)%neigh(4)=k+nx
          mesh%cell(k)%neigh(5)=k+nx+1
          mesh%cell(k)%neigh(6)=k+1
          mesh%cell(k)%neigh(7)=k-nx+1
          mesh%cell(k)%neigh(8)=k-nx

          xc=mesh%cell(k)%xc
          yc=mesh%cell(k)%yc
          do p=1,size(gauss_point)
             mesh%cell(k)%X_gauss(p)=xc+dx*gauss_point(p)/2.0_dp
             mesh%cell(k)%Y_gauss(p)=yc+dy*gauss_point(p)/2.0_dp
          enddo
          do p=1,size(gauss_point2)
             mesh%cell(k)%X_gauss2(p)=mesh%cell(k)%xc+mesh%cell(k)%dx*gauss_point2(p)/2.0_dp
             mesh%cell(k)%Y_gauss2(p)=mesh%cell(k)%yc+mesh%cell(k)%dy*gauss_point2(p)/2.0_dp
          enddo
          
          if (i==1) then
             mesh%cell(k)%neigh(1)=-((j-1)*nx)
             mesh%cell(k)%neigh(2)=-(j*nx)
             mesh%cell(k)%neigh(3)=-((j+1)*nx)
          endif

          if (j==1) then
             mesh%cell(k)%neigh(7)=-(i+(ny-1)*nx+1)
             mesh%cell(k)%neigh(8)=-(i+(ny-1)*nx)
             mesh%cell(k)%neigh(1)=-(i+(ny-1)*nx-1)
          endif
          
          if (i==nx) then
             mesh%cell(k)%neigh(5)=-(j*nx+1)
             mesh%cell(k)%neigh(6)=-((j-1)*nx+1)
             mesh%cell(k)%neigh(7)=-((j-2)*nx+1)
          endif

          if (j==ny) then
             mesh%cell(k)%neigh(3)=-(i-1)
             mesh%cell(k)%neigh(4)=-i
             mesh%cell(k)%neigh(5)=-(i+1)
          endif
       enddo
    enddo
    mesh%cell(1)%neigh(1)=-nx*ny
    mesh%cell(nx*ny)%neigh(5)=-1
    mesh%cell(nx)%neigh(7)=-((ny-1)*nx+1)
    mesh%cell((ny-1)*nx+1)%neigh(3)=-nx

    !Initialisation des edges
    
    do j=1,ny
       do i=1,nx+1
          k=(j-1)*(nx+1)+i
          mesh%edge(k)%cell1=(j-1)*nx+i-1
          mesh%edge(k)%cell2=(j-1)*nx+i
          mesh%edge(k)%dir=1
          mesh%edge(k)%length=dy
          mesh%edge(k)%period=k
          mesh%edge(k)%sub=0
       enddo
       mesh%edge((j-1)*(nx+1)+1)%cell1=-j*nx
       mesh%edge((j-1)*(nx+1)+1)%period=(j-1)*(nx+1)+nx+1
       mesh%edge((j-1)*(nx+1)+nx+1)%cell2=-((j-1)*nx+1)
       mesh%edge((j-1)*(nx+1)+nx+1)%period=(j-1)*(nx+1)+1
    enddo
    do i=1,nx
       do j=1,ny+1
          k=(i-1)*(ny+1)+j+(nx+1)*ny
          mesh%edge(k)%cell1=i+(j-1)*nx-nx
          mesh%edge(k)%cell2=i+(j-1)*nx
          mesh%edge(k)%dir=2
          mesh%edge(k)%length=dx
          mesh%edge(k)%period=k
          mesh%edge(k)%sub=0
       enddo
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%cell1=-(nx*(ny-1)+i)
       mesh%edge((i-1)*(ny+1)+1+(nx+1)*ny)%period=(i-1)*(ny+1)+ny+1+(nx+1)*ny
       mesh%edge((i-1)*(ny+1)+ny+1+(nx+1)*ny)%cell2=-i
       mesh%edge((i-1)*(ny+1)+ny+1+(nx+1)*ny)%period=(i-1)*(ny+1)+1+(nx+1)*ny
    enddo

    do i=1,mesh%ne
       edge=mesh%edge(i)
       dir=edge%dir
       allocate(mesh%edge(i)%X_gauss(size(gauss_point)),mesh%edge(i)%Y_gauss(size(gauss_point)))
       allocate(mesh%edge(i)%flux_acc(size(gauss_point)))
       allocate(mesh%edge(i)%flux(order,sol%nvar))
       mesh%edge(i)%flux_acc=.false.
       select case (dir)
       case(1)
          if (edge%cell1>0) then
             a=mesh%node(mesh%cell(edge%cell1)%corner(2))%y
             b=mesh%node(mesh%cell(edge%cell1)%corner(4))%y
             c=mesh%node(mesh%cell(edge%cell1)%corner(2))%x
          else
             a=mesh%node(mesh%cell(edge%cell2)%corner(1))%y
             b=mesh%node(mesh%cell(edge%cell2)%corner(3))%y
             c=mesh%node(mesh%cell(edge%cell2)%corner(1))%x
          endif
          center=(a+b)/2.0_dp
          diff=(b-a)/2.0_dp
          do p=1,size(gauss_point)
             mesh%edge(i)%X_gauss(p)=c
             mesh%edge(i)%Y_gauss(p)=center+diff*gauss_point(p)
          enddo
       case(2)
          if (edge%cell1>0) then
             a=mesh%node(mesh%cell(edge%cell1)%corner(3))%x
             b=mesh%node(mesh%cell(edge%cell1)%corner(4))%x
             c=mesh%node(mesh%cell(edge%cell1)%corner(3))%y
          else
             a=mesh%node(mesh%cell(edge%cell2)%corner(1))%x
             b=mesh%node(mesh%cell(edge%cell2)%corner(2))%x
             c=mesh%node(mesh%cell(edge%cell2)%corner(1))%y
          endif
          center=(a+b)/2.0_dp
          diff=(b-a)/2.0_dp
          do p=1,size(gauss_point)
             mesh%edge(i)%X_gauss(p)=center+diff*gauss_point(p)
             mesh%edge(i)%Y_gauss(p)=c
          enddo
       end select
    enddo
    
    return
  end subroutine buildMesh

  subroutine buildP4EST(level,connectivity,p4est)
    integer, intent(in) :: level
    type(c_ptr), intent(out) :: connectivity,p4est

    call p4_new(level,connectivity,p4est)

    return
  end subroutine buildP4EST
  
  subroutine buildMesh_P4EST(p4est,xL,xR,yL,yR,gauss_point,order,mesh,sol,quadrants)
    type(c_ptr), intent(in) :: p4est
    real(dp), intent(in) :: xL,xR,yL,yR
    integer, intent(in) :: order
    real(dp), dimension(:), intent(in) :: gauss_point
    type(meshStruct), intent(inout) :: mesh
    type(solStruct), intent(inout) :: sol
    type(c_ptr), intent(out) :: quadrants
    integer(c_int) :: tt
    type(c_ptr) :: p4_mesh,nodes,edges
    type(c_ptr) :: C_corners,C_neighbors,C_sub,C_cell1,C_cell2,C_iedge,C_period,C_nodes
    integer, dimension(:), pointer :: F_corners,F_neighbors,F_sub,F_cell1,F_cell2
    integer, dimension(:), pointer :: F_iedge,F_period,F_nodes
    integer :: k,i,lev,p,Nneigh,N_edge,Nnodes,ie,nedge,i1,iloc
    real(dp) :: a,b,c,center,diff

    call p4_build_mesh(p4est,tt,p4_mesh,quadrants,nodes,edges,mesh%np,mesh%nc,mesh%ne)
    
    if (allocated(mesh%node)) deallocate(mesh%node,mesh%cell,mesh%edge)
    allocate(mesh%node(mesh%np),mesh%cell(mesh%nc),mesh%edge(mesh%ne))
    if (allocated(sol%val)) deallocate(sol%val,sol%user)
    allocate(sol%val(mesh%nc,sol%nvar),sol%user(mesh%nc,sol%nsolUser))

    do k=1,mesh%np
       call p4_get_node(p4est,tt,nodes,k-1,mesh%node(k)%x,mesh%node(k)%y)
       mesh%node(k)%x=(xR-xL)*mesh%node(k)%x+xL
       mesh%node(k)%y=(yR-yL)*mesh%node(k)%y+yL
    enddo

    ie=0
    do k=1,mesh%nc
       call p4_get_cell(p4est,p4_mesh,tt,quadrants,nodes,edges,k-1,mesh%cell(k)%xc,mesh%cell(k)%yc, &
            mesh%cell(k)%dx,mesh%cell(k)%dy,C_corners,Nneigh,C_neighbors,Nnodes,C_nodes,N_edge,lev)

       allocate(mesh%cell(k)%node(Nnodes))
       allocate(mesh%cell(k)%neigh(Nneigh))
       allocate(mesh%cell(k)%X_gauss(max(size(gauss_point),2)))
       allocate(mesh%cell(k)%Y_gauss(max(size(gauss_point),2)))
       allocate(mesh%cell(k)%polMax(order*(order-1)/2+order-1,sol%nvar))
       
       mesh%cell(k)%xc=(xR-xL)*mesh%cell(k)%xc+xL
       mesh%cell(k)%yc=(yR-yL)*mesh%cell(k)%yc+yL
       mesh%cell(k)%dx=(xR-xL)*mesh%cell(k)%dx
       mesh%cell(k)%dy=(yR-yL)*mesh%cell(k)%dy
       mesh%cell(k)%level=lev
       
       call c_f_pointer(C_corners,F_corners,(/4/))
       call c_f_pointer(C_nodes,F_nodes,(/Nnodes/))
       call c_f_pointer(C_neighbors,F_neighbors,(/12/))
       
       mesh%cell(k)%corner=F_corners
       mesh%cell(k)%node=F_nodes
       mesh%cell(k)%neigh(1:Nneigh)=F_neighbors(1:Nneigh)

       do p=1,size(gauss_point)
          mesh%cell(k)%X_gauss(p)=mesh%cell(k)%xc+mesh%cell(k)%dx*gauss_point(p)/2.0_dp
          mesh%cell(k)%Y_gauss(p)=mesh%cell(k)%yc+mesh%cell(k)%dy*gauss_point(p)/2.0_dp
       enddo
       do p=1,size(gauss_point2)
          mesh%cell(k)%X_gauss2(p)=mesh%cell(k)%xc+mesh%cell(k)%dx*gauss_point2(p)/2.0_dp
          mesh%cell(k)%Y_gauss2(p)=mesh%cell(k)%yc+mesh%cell(k)%dy*gauss_point2(p)/2.0_dp
       enddo

       allocate(mesh%cell(k)%edge(N_edge))
       iloc=1
       do i=1,4

          call p4_get_edge(p4est,p4_mesh,quadrants,edges,k-1,i-1,C_iedge,nedge,C_cell1,C_cell2,C_sub,C_period)
          if (abs(nedge)==2) then
             call c_f_pointer(C_iedge,F_iedge,(/2/))
             call c_f_pointer(C_period,F_period,(/2/))
          else
             call c_f_pointer(C_iedge,F_iedge,(/1/))
             call c_f_pointer(C_period,F_period,(/1/))
          endif

          do i1=1,abs(nedge)
             mesh%cell(k)%edge(iloc)=abs(F_iedge(i1))
             iloc=iloc+1
             if (F_iedge(i1)>0) then
                allocate(mesh%edge(F_iedge(i1))%flux_acc(size(gauss_point)))
                allocate(mesh%edge(F_iedge(i1))%flux(order,sol%nvar))
                allocate(mesh%edge(F_iedge(i1))%X_gauss(size(gauss_point)),mesh%edge(F_iedge(i1))%Y_gauss(size(gauss_point)))
                
                call c_f_pointer(C_sub,F_sub,(/2/))

                select case (nedge)
                case (1)
                   call c_f_pointer(C_cell1,F_cell1,(/1/))
                   call c_f_pointer(C_cell2,F_cell2,(/1/))
                   ie=ie+1
                case default
                   call c_f_pointer(C_cell1,F_cell1,(/2/))
                   call c_f_pointer(C_cell2,F_cell2,(/2/))
                   ie=ie+2
                end select

                mesh%edge(F_iedge(i1))%flux_acc=.false.
                mesh%edge(F_iedge(i1))%sub=F_sub
                mesh%edge(F_iedge(i1))%cell1=F_cell1(i1)
                mesh%edge(F_iedge(i1))%cell2=F_cell2(i1)
                mesh%edge(F_iedge(i1))%period=F_period(i1)
                   
                select case (i)
                case (1,2)
                   mesh%edge(F_iedge(i1))%dir=1
                   mesh%edge(F_iedge(i1))%length=(yR-yL)/(nedge*2**lev)
                case (3,4)
                   mesh%edge(F_iedge(i1))%dir=2
                   mesh%edge(F_iedge(i1))%length=(xR-xL)/(nedge*2**lev)
                end select

                select case (mesh%edge(F_iedge(i1))%dir)
                case(1)
                   if (F_cell1(i1)>0) then
                      a=mesh%node(mesh%cell(F_cell1(i1))%corner(2))%y
                      b=mesh%node(mesh%cell(F_cell1(i1))%corner(3))%y
                      c=mesh%node(mesh%cell(F_cell1(i1))%corner(2))%x
                   else
                      a=mesh%node(mesh%cell(F_cell2(i1))%corner(1))%y
                      b=mesh%node(mesh%cell(F_cell2(i1))%corner(4))%y
                      c=mesh%node(mesh%cell(F_cell2(i1))%corner(1))%x
                   endif
                   center=(a+b)/2.0_dp
                   diff=(b-a)/2.0_dp
                   if (nedge/=1) then
                      center=center+(i1-1.5_dp)*diff
                      diff=diff/2.0_dp
                   endif
                   do p=1,size(gauss_point)
                      mesh%edge(F_iedge(i1))%X_gauss(p)=c
                      mesh%edge(F_iedge(i1))%Y_gauss(p)=center+diff*gauss_point(p)
                   enddo
                case(2)
                   if (F_cell1(i1)>0) then
                      a=mesh%node(mesh%cell(F_cell1(i1))%corner(4))%x
                      b=mesh%node(mesh%cell(F_cell1(i1))%corner(3))%x
                      c=mesh%node(mesh%cell(F_cell1(i1))%corner(4))%y
                   else
                      a=mesh%node(mesh%cell(F_cell2(i1))%corner(1))%x
                      b=mesh%node(mesh%cell(F_cell2(i1))%corner(2))%x
                      c=mesh%node(mesh%cell(F_cell2(i1))%corner(1))%y
                   endif
                   center=(a+b)/2.0_dp
                   diff=(b-a)/2.0_dp
                   if (nedge/=1) then
                      center=center-(i1-1.5_dp)*diff
                      diff=diff/2.0_dp
                   endif
                   do p=1,size(gauss_point)
                      mesh%edge(F_iedge(i1))%X_gauss(p)=center+diff*gauss_point(p)
                      mesh%edge(F_iedge(i1))%Y_gauss(p)=c
                   enddo
                end select
             endif
          enddo
          call p4_free(C_iedge)
          call p4_free(C_cell1)
          call p4_free(C_cell2)
          call p4_free(C_sub)
          call p4_free(C_period)
       enddo
       call p4_free(C_corners)
       call p4_free(C_neighbors)
       call p4_free(C_nodes)
    enddo
    call p4_free(edges)
    call p4_destroy_mesh(p4_mesh,nodes)
    print*,"Mesh created with ",mesh%nc," quadrants"
    
    return
  end subroutine buildMesh_P4EST
  
  subroutine writeSol(mesh,sol,namefile,nfile)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    character(len=20), intent(in) :: namefile
    integer, intent(in) :: nfile
    integer :: k,n
    integer :: i,nnodes
    character(len=34) :: completenamefile
    
    write(completenamefile,'(A,A,I3.3,A)')'./results/',trim(namefile),nfile,'.vtk'
    open(11,file=completenamefile,form="formatted")
    
    write(11,'(a)')"# vtk DataFile Version 2.0"
    write(11,'(a)')"Results of the calculation"
    write(11,'(a)')"ASCII"
    write(11,'(a)')"DATASET UNSTRUCTURED_GRID"
    
    write(11,'(a,i8,a)')"POINTS ",mesh%np," float"
    do k=1,mesh%np
       write(11,'(e15.8,a,e15.8,a,e15.8)')mesh%node(k)%x," ",mesh%node(k)%y," ",0.0_dp
    enddo

    n=0
    do k=1,mesh%nc
       n=n+size(mesh%cell(k)%node)+1
    enddo
    write(11,'(a,i8,i9)')"CELLS ",mesh%nc,n
    do k=1,mesh%nc
       nnodes=size(mesh%cell(k)%node)
       write(11,'(i1)',advance='no')nnodes
       do i=1,nnodes-1
          write(11,'(i8)',advance='no')mesh%cell(k)%node(i)-1
       enddo
       write(11,'(i8)')mesh%cell(k)%node(nnodes)-1
    enddo

    write(11,'(a,i8)')"CELL_TYPES ",mesh%nc
    do k=1,mesh%nc
       write(11,'(i1)')9
    enddo
    
    write(11,'(a,i8)')"CELL_DATA ",mesh%nc
    do n=1,sol%nvar
       write(11,'(a,a,a)')"SCALARS ",sol%name(n)," float 1"
       write(11,'(a)')"LOOKUP_TABLE default"
       do k=1,mesh%nc
          if (abs(sol%val(k,n))>10.0**(-dp)) then
             write(11,'(e15.8)')sol%val(k,n)
          else
             write(11,'(e15.8)')0.0_dp
          endif
       enddo
    enddo

    do n=1,sol%nsolUser
       write(11,'(a,a,a)')"SCALARS ",sol%name_user(n)," float 1"
       write(11,'(a)')"LOOKUP_TABLE default"
       do k=1,mesh%nc
             write(11,'(e15.8)')sol%user(k,n)
       enddo
    enddo
    
    close(11)
    
    return
  end subroutine writeSol

  subroutine print(mesh,sol,t,n,total)
    type(meshStruct), intent(in) :: mesh
    type(solStruct), intent(in) :: sol
    real(dp), intent(in) :: t
    integer, intent(in) :: n
    real(dp), dimension(:), intent(inout) :: total

    print*,"-----------------------------------------"
    print*,"t=",t,"itération ",n
    call check_conservativity(mesh,sol,total)
    print*,"-----------------------------------------"
    
    return
  end subroutine print

  subroutine write_accept(mesh,NAC_cycle,NAC_reason,n)
    type(meshStruct), intent(in) :: mesh
    integer, dimension(:), intent(in) :: NAC_cycle,NAC_reason
    integer, intent(in) :: n
    integer :: k,i,nnodes
    character(len=34) :: completenamefile

    write(completenamefile,'(A,A,I3.3,A)')'./results/','accept',n,'.vtk'
    open(11,file=completenamefile,form="formatted")
    
    write(11,'(a)')"# vtk DataFile Version 2.0"
    write(11,'(a)')"Results of the calculation"
    write(11,'(a)')"ASCII"
    write(11,'(a)')"DATASET UNSTRUCTURED_GRID"
    
    write(11,'(a,i8,a)')"POINTS ",mesh%np," float"
    do k=1,mesh%np
       write(11,'(e15.8,a,e15.8,a,e15.8)')mesh%node(k)%x," ",mesh%node(k)%y," ",0.0_dp
    enddo

    nnodes=0
    do k=1,mesh%nc
       nnodes=nnodes+size(mesh%cell(k)%node)+1
    enddo
    write(11,'(a,i8,i9)')"CELLS ",mesh%nc,nnodes
    do k=1,mesh%nc
       nnodes=size(mesh%cell(k)%node)
       write(11,'(i1)',advance='no')nnodes
       do i=1,nnodes-1
          write(11,'(i8)',advance='no')mesh%cell(k)%node(i)-1
       enddo
       write(11,'(i8)')mesh%cell(k)%node(nnodes)-1
    enddo

    write(11,'(a,i8)')"CELL_TYPES ",mesh%nc
    do k=1,mesh%nc
       write(11,'(i1)')9
    enddo
    
    write(11,'(a,i8)')"CELL_DATA ",mesh%nc
    write(11,'(a,a,a)')"SCALARS ","cycle"," float 1"
    write(11,'(a)')"LOOKUP_TABLE default"
    do k=1,mesh%nc
       write(11,'(e15.8)')real(NAC_cycle(k))
    enddo
    write(11,'(a,a,a)')"SCALARS ","reason"," float 1"
    write(11,'(a)')"LOOKUP_TABLE default"
    do k=1,mesh%nc
       write(11,'(e15.8)')real(NAC_reason(k))
    enddo
    
    close(11)
    
    return
  end subroutine write_accept

  subroutine write_output_header(test_case,xL,xR,yL,yR,level,cfl,tf,namefile, &
       str_equa,str_flux,str_time_scheme,order,L_str_criteria, &
       bool_AMR,str_fn_adapt,f_adapt,recursivity,levelmin,levelmax)
    character(len=20), intent(in) :: test_case,namefile,str_equa,str_flux,str_time_scheme
    real(dp), intent(in) :: xL,xR,yL,yR,cfl,tf
    integer, intent(in) :: level,order,f_adapt,recursivity,levelmin,levelmax
    character(len=20), dimension(:), allocatable, intent(in) :: L_str_criteria
    character(len=20), intent(in) :: str_fn_adapt
    logical, intent(in) :: bool_AMR
    character(len=30) :: completenamefile
    
    write(completenamefile,'(a,a,a)')'./results/',trim(namefile),'.out'
    open(12,file=trim(completenamefile),form="formatted")

    write(12,'(a)')"-------------------- Configuration --------------------"
    write(12,'(a,a)')"ICBC : ",trim(test_case)
    write(12,'(a,es10.2,a,es10.2,a,es10.2,a,es10.2,a)')"Domain [",xl,";",xR,"] x [",yl,";",yR,"]"
    write(12,'(a,i2)')"Initial level of mesh refinement : ",level
    write(12,'(a,i2,a,i2)')"Level min : ",levelmin,"   Level max : ",levelmax
    write(12,'(a,a)')"Equation : ",trim(str_equa)
    write(12,'(a,f5.2)')"CFL : ",cfl
    write(12,'(a,es10.2)')"Final time : ",tf
    write(12,'(a,a)')"Flux : ",trim(str_flux)
    write(12,'(a,a)')"Time scheme : ",trim(str_time_scheme)
    write(12,'(a,i1)')"Order : ",order
    write(12,'(a,i1)')"Number of detection criterias : ",size(L_str_criteria)
    write(12,'(a)')"-------------------- AMR ------------------------------"
    if (.NOT.bool_AMR) then
       write(12,'(a)')"AMR : NO"
    else
       write(12,'(a)')"AMR : YES"
       write(12,'(a,a)')"Adaptation function : ",trim(str_fn_adapt)
       write(12,'(a,i5)')"Adaptation frequency : ",f_adapt
       write(12,'(a,i1)')"Recursivity : ",recursivity
    endif
    write(12,'(a)')"-------------------- Calculation ----------------------"

    return
  end subroutine write_output_header

  subroutine write_output_calculation(t,n,nc,total)
    real(dp), intent(in) :: t
    integer, intent(in) :: n,nc
    real(dp), dimension(:), intent(in) :: total
    integer :: i

    write(12,'(a,es10.2,a,i10)')"t = ",t,"          itération ",n
    write(12,'(a,i10)')"Number of cells : ",nc
    write(12,'(a)',advance='no')"Total quantities : "
    do i=1,size(total)-1
       write(12,'(es22.15,a)',advance='no')total(i),"   "
    enddo
    write(12,'(es22.15)')total(size(total))

    return
  end subroutine write_output_calculation

  subroutine write_output_summary(cpu,eL1,eL2,total_cell,average_cell)
    real(dp), intent(in) :: cpu,eL1,eL2
    integer :: total_cell,average_cell
    
    write(12,'(a)')"-------------------- Results --------------------------"
    write(12,'(a,es10.3,a)')"CPU time : ",cpu," seconds"
    write(12,'(a,i8)')"Average number of cells : ",average_cell
    write(12,'(a,es10.3,a)')"CPU time per cell : ",cpu/total_cell," seconds"
    if (eL1/=-1.0_dp) then
       write(12,'(a,es10.3)')"L1 error : ",eL1
       write(12,'(a,es10.3)')"L2 error : ",eL2
    else
       write(12,'(a)')"No analytical solution for this configuration"
    endif

    close(12)

    return
  end subroutine write_output_summary
    

end module inout
