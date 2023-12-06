module parameters
    Implicit None
!--------to be midified by the usere
    character(len=80):: prefix="BiTeI"
    real*8,parameter::ef= 4.18903772,kmax=0.07, alpha = 0.7747474
    integer,parameter::meshres=30, nkpoints=(2*meshres+1),nbmin=12,nbmax=13
    integer nb
    
end module parameters
Program Projected_band_structure
    use parameters
    Implicit None
!------------------------------------------------------
    real*8 dx, dy
    character(len=80) top_file,triv_file,nnkp,line
    integer*4 i,j,k,nr,i1,i2,lwork,info,ikx,iky,j1,j2,count
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,a,b,x1,y1, min_eg
    real*8 avec(3,3),bvec(3,3),rvec(3),kpoint(3)
    real*8,allocatable:: rvec_data(:,:),rvec_data_t(:,:),ene(:),rwork(:),k_ene(:), sam(:,:), oam(:,:), e_g_data(:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),Top_hr(:,:,:),Triv_hr(:,:,:),work(:)
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
    open(100,file='interpolate.dx')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_data(3,nr),rvec_data_t(3,nr),Hk(nb,nb),Top_hr(nb,nb,nr),Triv_hr(nb,nb,nr),ndeg(nr),ene(nb))
    read(99,*)ndeg
    do i = 1, 80
        read(97, *)! Read and discard 80 lines
    end do
    do k=1,nr
       do i=1,nb
          do j=1,nb
             read(99,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),i1,i2,a,b
             top_hr(i1,i2,k)=dcmplx(a,b)
             read(97,*)rvec_data_t(1,k),rvec_data_t(2,k),rvec_data_t(3,k),j1,j2,x1,y1
             triv_hr(j1,j2,k)=dcmplx(x1,y1)
          enddo
       enddo
    enddo
   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))
!----- Create K-mesh
    dx = kmax / meshres
    dy = kmax / meshres
    
!----- Perform fourier transform
    allocate(sam(3,nbmin:nbmax), oam(3,nbmin:nbmax), k_ene(nb), e_g_data(nkpoints**2))
    count = 1
    do ikx=-meshres,meshres
        do iky=-meshres,meshres
            kpoint(1)= ikx*dx
            kpoint(2)= iky*dy
            kpoint(3)= 0.5d0*bvec(3,3)
            HK=(0d0,0d0)
            do i=1,nr
                rvec = rvec_data(1,i)*avec(:,1) + rvec_data(2,i)*avec(:,2) + rvec_data(3,i)*avec(:,3)
                phase = dot_product(kpoint,rvec)
                HK=HK+((1-alpha)*(triv_hr(:,:,i))+alpha*(top_hr(:,:,i)))*dcmplx(cos(phase),-sin(phase))/float(ndeg(i))
            enddo
            call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)
            e_g_data(count) = k_ene(13)-k_ene(12)
            count = count +1
            write(100, '(2(1x,f12.6))') k_ene(nbmin), k_ene(nbmax)
        enddo
    enddo

    min_eg = MINVAL(e_g_data(:))
    print *, "Min E_G: ", min_eg

end program Projected_band_structure