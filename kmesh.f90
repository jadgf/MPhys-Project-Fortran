Program Wannier_band_structure
    Implicit None
!--------to be midified by the usere
    character(len=80):: prefix="BiTeI"
    integer,parameter::nkpath=3,np=100,meshres=1000
    real*8,parameter::ef= 4.18903772,ikmax=0.1, tolerance = 0.05
!------------------------------------------------------
    real*8 kz,kmesh(3,meshres**2), CBM
    character(len=30)::klabel(nkpath)
    character(len=80) hamil_file,nnkp,line
    integer*4,parameter::nk=(nkpath-1)*np+1
    integer*4 i,j,k,nr,i1,i2,nb,lwork,info,count
    real*8,parameter::third=1d0/3d0!,kz=0d0
    real*8 phase,pi2,jk,a,b
    real*8 klist(3,1:nk),xk(nk),kpath(3,np),avec(3,3),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath),rvec(3)
    real*8,allocatable:: rvec_data(:,:),ene(:,:),rwork(:),k_ene(:),kpoints(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),Hamr(:,:,:),work(:)
    complex*16 temp1,temp2
!------------------------------------------------------
    write(hamil_file,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
    write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

    pi2=4.0d0*atan(1.0d0)*2.0d0

!---------------  reciprocal vectors
    open(98,file=trim(adjustl(nnkp)),err=333)
111   read(98,'(a)')line
    if(trim(adjustl(line)).ne."begin real_lattice") goto 111
    read(98,*)avec
    do i=1,3
        read(98, '(a)')line
    enddo
    read(98,*)bvec
!---------------kpath
    data kpath(:,1) /     0.5d0,      0.0d0,    0.5d0/  !L
    data kpath(:,2) /     0.0d0,      0.0d0,    0.5d0/  !A
    data kpath(:,3) /     third,      third,    0.5d0/  !H

    data klabel     /'L','A','H'/

    ktemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+(kpath(2,1)-kpath(2,2))*bvec(:,2)+(kpath(3,1)-kpath(3,2))*bvec(:,3)

    xk(1)= -sqrt(dot_product(ktemp1,ktemp1))
    xkl(1)=xk(1)

    k=0
    ktemp1=0d0
    do i=1,nkpath-1
     do j=1,np
      k=k+1
      jk=dfloat(j-1)/dfloat(np)
      klist(:,k)=kpath(:,i)+jk*(kpath(:,i+1)-kpath(:,i))
      ktemp2=klist(1,k)*bvec(:,1)+klist(2,k)*bvec(:,2)+klist(3,k)*bvec(:,3)
      if(k.gt.1) xk(k)=xk(k-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
      if(j.eq.1) xkl(i)=xk(k)
      ktemp1=ktemp2
     enddo
    enddo
    klist(:,nk)=kpath(:,nkpath)
    ktemp2=klist(1,nk)*bvec(:,1)+klist(2,nk)*bvec(:,2)+klist(3,nk)*bvec(:,3)
    xk(nk)=xk(nk-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
    xkl(nkpath)=xk(nk)
!      write(*,*)klist
    klist=klist*pi2

!------read H(R)
    open(99,file=trim(adjustl(hamil_file)),err=444)
    open(100,file='kmesh.dat')
    read(99,*)
    read(99,*)nb,nr
    allocate(rvec_data(3,nr),Hk(nb,nb),Hamr(nb,nb,nr),ndeg(nr),ene(nb,nk))
    read(99,*)ndeg
    do k=1,nr
       do i=1,nb
          do j=1,nb
             read(99,*)rvec_data(1,k),rvec_data(2,k),rvec_data(3,k),i1,i2,a,b
             hamr(i1,i2,k)=dcmplx(a,b)
          enddo
       enddo
    enddo

   lwork=max(1,2*nb-1)
   allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

!---- Fourier transform H(R) to H(k)
    ene=0d0
    do k=1,nk
       HK=(0d0,0d0)
       do j=1,nr

          phase=0.0d0
          do i=1,3
             phase=phase+klist(i,k)*rvec_data(i,j)
          enddo

          HK=HK+Hamr(:,:,j)*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))

       enddo

       call zheev('V','U',nb,Hk,nb,ene(:,k),work,lwork,rwork,info)
       
    enddo
    CBM = MINVAL(ene(13, :))

    deallocate(Hk,ene,work)
    allocate(Hk(nb,nb),k_ene(nb))
!----- Create K-mesh
    do i=1,meshres
        do j=1, meshres
            kmesh(1,i) = pi2*ikmax * (i-1) / meshres
            kmesh(2,i) = pi2*ikmax * (j-1) / meshres
            kmesh(3,i) = 0.5
        enddo 
    enddo

!----- Perform cartesian fourier transform
    count = 0
    do k=1,meshres**2
        HK=(0d0,0d0)
        do i=1,nr
            rvec = rvec_data(1,i) * avec(:,1) + rvec_data(2,i) * avec(:,2) + rvec_data(3,i) * avec(:,3)
            phase = rvec(1)*kmesh(1,k) + rvec(2)*kmesh(2,k) + rvec(3)*kmesh(3,k)
            HK=HK+Hamr(:,:,j)*dcmplx(cos(phase),-sin(phase))/float(ndeg(j))
        enddo
        call zheev('V','U',nb,Hk,nb,k_ene,work,lwork,rwork,info)
        if((ABS(k_ene(12) - CBM).lt.tolerance).or.(ABS(k_ene(13) - CBM).lt.tolerance)) then
            count = count + 1
            allocate(kpoints(3,count))
            kpoints(:,count) = kmesh(:,k)
        endif
    enddo

    do i=1,SIZE(kpoints)
        print *, kpoints(1,i)
    enddo
    do i=1, count
       do k=1, 2
         write(100,'(2(f12.6))') kpoints(k,i)
       enddo
         write(100,*)
         write(100,*)
    enddo
    call write_plt(nkpath,xkl,klabel,ef)
    stop
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
    stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file)),' not found'
    stop

    end

   subroutine write_plt(nkp,xkl,kl,ef)
   implicit none
   integer nkp,i
   real*8 xkl(nkp),ef
   character(len=30)kl(nkp)
   
   open(99,file='kmesh.plt')
   write(99,'(a,f12.6,a,f12.6,a)') 'set xrange [-0.8 : 0.8]'
   write(99,'(a)') &
        'set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in'
   write(99,'(a,f4.2,a)')'set output "kmesh.pdf"'
   write(99,'(9(a,/),a)') &
        'set encoding iso_8859_1',&
        'set size ratio 0 1.0,1.0',&
        'set ylabel "$k_{y}$ coordinate"',&
        'set yrange [ -0.8 : 0.8 ]',&
        'unset key',&
        'set ytics 1.0 scale 1 nomirror out',&
        'set mytics 2',&
        'set parametric',&
        'set trange [-10:10]',&
        'plot "kmesh.dat" u 1:2 with points lt 1 lw 3,\'
  end subroutine write_plt

! ------ gfortran -o bandplot wanr2k.f90 -lblas -llapack