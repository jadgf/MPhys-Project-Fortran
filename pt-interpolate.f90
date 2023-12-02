module parameters
    Implicit None
!--------to be midified by the usere
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.1,kzmax=0.2,alpha=0.77474747,tolerance=0.03
    integer,parameter::meshres=30,zmeshres=30,nkpoints=(2*meshres+1),nkzpoints=(2*zmeshres+1),nbmin=12,nbmax=13
    integer nb
    
end module parameters

Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dxy,dz
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,lwork,info,ikx,iky,ikz,ik,count
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,a,b
    real*8 avec(3,3),bvec(3,3),rvec(3),kpoint(3)
    real*8,allocatable:: rvec_data(:,:),rvec_data_t(:,:),ene(:),rwork(:)
	real*8, allocatable:: k_ene(:),k_ene_data(:,:),sam(:,:),oam(:,:),kmesh(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),top_Hr(:,:,:),triv_Hr(:,:,:),work(:)
!------------------------------------------------------
    write(top_file,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
    write(triv_file,'(a,a)')trim(adjustl(prefix)),"_hr_trivial.dat"
    write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

    pi2=4.0d0*atan(1.0d0)*2.0d0

!---------------  reciprocal vectors
    open(98,file=trim(adjustl(nnkp)))
111 read(98,'(a)')line
    if(trim(adjustl(line)).ne."begin real_lattice") goto 111
    read(98,*)avec
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,*)bvec

!------read H(R)
    open(99,file=trim(adjustl(top_file)))
    open(97,file=trim(adjustl(triv_file)))
    open(100,file='3denergy.csv')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_data(3,nr),rvec_data_t(3,nr),Hk(nb,nb),top_Hr(nb,nb,nr),triv_Hr(nb,nb,nr),ndeg(nr),ene(nb))
    read(99,*)ndeg
    read(97,'(80(a))')line
    do k=1,nr
       do i=1,nb
          do j=1,nb
             read(99,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),i1,i2,a,b
             top_Hr(i1,i2,k)=dcmplx(a,b)
             read(97,*)rvec_data_t(1,k),rvec_data_t(2,k),rvec_data_t(3,k),i1,i2,a,b
             triv_Hr(i1,i2,k)=dcmplx(a,b)
          enddo
       enddo
    enddo

   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
	
	dxy=kmax/meshres
	dz=kzmax/zmeshres

!----- Create header of dx files
!   write(100, '(a,3(1x,i8))') 'object 1 class gridpositions counts',nkpoints,nkpoints,nkzpoints
!   write(100, '(a,3(1x,f12.6))') 'origin',-kmax,-kmax,-kzmax
!   write(100, '(a,3(1x,f12.6))') 'delta',dxy,0d0,0d0
!   write(100, '(a,3(1x,f12.6))') 'delta',0d0,dxy,0d0
!   write(100, '(a,3(1x,f12.6))') 'delta',0d0,0d0,dz
!   write(100, '(a,3(1x,i8))') 'object 2 class gridconnections counts',nkpoints,nkpoints,nkzpoints
!   write(100, '(a,i8,a,i10,a)') 'object 3 class array type float rank 1 shape',nbmax-nbmin+1,&
!                                    ' item', nkpoints,' data follows'

!----- Perform fourier transform
    allocate(sam(3,nbmin:nbmax), oam(3,nbmin:nbmax), k_ene(nb), k_ene_data(nkpoints*nkpoints*nkzpoints,2))
	allocate(kmesh(nkpoints*nkpoints*nkzpoints,3))

	count=0
	do ikx=-meshres,meshres
		do iky=-meshres,meshres
			do ikz=-zmeshres,zmeshres
				count=count+1
				kpoint(1)= ikx*dxy
				kpoint(2)= iky*dxy
				kpoint(3)= ikz*dz

				do i=1,nr
					rvec = rvec_data(1,i)*avec(:,1) + rvec_data(2,i)*avec(:,2) + rvec_data(3,i)*avec(:,3)
					
					phase = dot_product(kpoint,rvec)
					HK=HK+((1-alpha)*(triv_Hr(:,:,i))+(alpha)*(top_Hr(:,:,i)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(i))
					!HK=HK+(top_Hr(:,:,i)*dcmplx(cos(phase),-sin(phase))/float(ndeg(i)))
				enddo
				call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
				if(k_ene(nbmax) - k_ene(nbmin) < tolerance) then
					write(100, '((1x,f12.6),a,(1x,f12.6),a,(1x,f12.6))') kpoint(1), ',', kpoint(2), ',', kpoint(3) 
				endif
			enddo
		enddo
	enddo
!   write(100,'(A,/,A,/,A,/,A)') &
!   'object "regular positions regular connections" class field', &
!   'component "positions" value 1', &
!   'component "connections" value 2', &
!   'component "data" value 3', &
!   'end'
end Program Projected_band_structure
