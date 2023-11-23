Program Wannier_band_structure
    Implicit None
!--------to be midified by the usere
    character(len=80):: prefix="data/BiTeI"
    real*8,parameter::ef= 4.18903772,ikmax=0.14, tolerance = 0.001, energy_diff = 0.32
    integer,parameter::nkpath=3,np=100,meshres=100, pre_size = int(10*100*(meshres**2*tolerance))
!------------------------------------------------------
    real*8 kmesh(3,meshres**2), dx, dy, CBM, TVB, m_x, m_y, m_z, L_x, L_y, L_z
    character(len=30)::klabel(nkpath)
    character(len=80) hamil_file,nnkp,line
    integer*4,parameter::nk=(nkpath-1)*np+1
    integer*4 i,j,k,nr,i1,i2,nb,lwork,info,count
    real*8,parameter::third=1d0/3d0, two = 2.0d0, sqrt2 = sqrt(two)
    real*8 phase,pi2,jk,a,b
    real*8 klist(3,1:nk),xk(nk),kpath(3,np),avec(3,3),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath),rvec(3),kvec(3)
    real*8,allocatable:: rvec_data(:,:),ene(:,:),rwork(:),k_ene(:),kpoints(:,:), m_k_data(:,:),  L_k_data(:,:)
    integer*4,allocatable:: ndeg(:)
    complex*16,allocatable:: Hk(:,:),Hamr(:,:,:),work(:),kp_eivec(:,:,:),H_col(:)
    complex*16 chi_k(2,1),m_x_comp(1,1), m_y_comp(1,1), m_z_comp(1,1), L_x_comp(1,1), L_y_comp(1,1),&
                L_z_comp(1,1), phi_k(3,1)
    complex*8 pauli_x(2, 2), pauli_y(2, 2), pauli_z(2, 2), Lx(3,3), Ly(3,3), Lz(3,3), Y_lm(3,1)
    complex*8, parameter:: one = complex(1.d0,0.d0),im = complex(0.d0,1.d0), zero = complex(0.d0,0.d0)
!------------------------------------------------------
    write(hamil_file,'(a,a)')trim(adjustl(prefix)),"_hr_topological.dat"
    write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

    pi2=4.0d0*atan(1.0d0)*2.0d0

!---------------  reciprocal vectors
    open(98,file=trim(adjustl(nnkp)),err=333)
111 read(98,'(a)')line
    if(trim(adjustl(line)).ne."begin real_lattice") goto 111
    read(98,*)avec
    read(98,'(a)')line
    read(98,'(a)')line
    read(98,'(a)')line
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
    open(100,file='dat_files/kmesh.dat')
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

!----- Fourier transform H(R) to H(k)
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
    print *, "CBM: ", CBM
    TVB = MAXVAL(ene(12, :))
    print *, "TVB: ", TVB

    deallocate(Hk,ene,work,rwork)
    allocate(Hk(nb,nb),k_ene(nb), H_col(nb))
    lwork=max(1,2*nb-1)
    allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

!----- Create K-mesh
    dx = 2.0 * ikmax / meshres
    dy = 2.0 * ikmax / meshres
    kmesh(3,:) = 0.5
    count=1
    do i=1,meshres ! +1 to include +0.1 point
        do j=1,meshres ! +1 to include +0.1 point
            kmesh(1, count) = -ikmax + (dx * (i - 1))
            kmesh(2, count) = -ikmax + (dy * (j - 1))
            count = count + 1
        enddo 
    enddo
    
!----- Perform cartesian fourier transform
    allocate(kpoints(2,pre_size), kp_eivec(pre_size,nb,nb), m_k_data(3,pre_size),  L_k_data(3,pre_size))

    count=1
    k_ene = 0d0
    do k=1,meshres**2
        HK=(0d0,0d0)
        kvec = kmesh(1,k)*bvec(:,1) + kmesh(2,k)*bvec(:,2) + kmesh(3,k)*bvec(:,3)
        do i=1,nr
            rvec = rvec_data(1,i)*avec(:,1) + rvec_data(2,i)*avec(:,2) + rvec_data(3,i)*avec(:,3)
            phase = rvec(1)*kvec(1) + rvec(2)*kvec(2) + rvec(3)*kvec(3)
            HK=HK+Hamr(:,:,i)*dcmplx(cos(phase),-sin(phase))/float(ndeg(i))
        enddo
        call zheev('V','U',nb,HK,nb,k_ene,work,lwork,rwork,info)

        if(k_ene(13).lt.(CBM + energy_diff)) then
            kpoints(:,count) = kvec(:)
            kp_eivec(count,:,:) = Hk(:,:)
            count = count + 1
            !print *, kvec(:)
        endif
    enddo
    print *, "Number of intersections: ", count-1

!-----Spin projection
   !-Define Pauli matrices
    pauli_x(1, 1) = complex(0.d0, 0.d0); pauli_x(1, 2) = complex(1.d0, 0.d0)
        pauli_x(2, 1) = complex(1.d0, 0.d0); pauli_x(2, 2) = complex(0.d0, 0.d0)
    pauli_y(1, 1) = complex(0.d0, 0.d0); pauli_y(1, 2) = complex(0.d0, -1.d0)
        pauli_y(2, 1) = complex(0.d0, 1.d0); pauli_y(2, 2) = complex(0.d0, 0.d0)
    pauli_z(1, 1) = complex(1.d0, 0.d0); pauli_z(1, 2) = complex(0.d0, 0.d0)
        pauli_z(2, 1) = complex(0.d0, 0.d0); pauli_z(2, 2) = complex(-1.d0, 0.d0)
   !-Defined
    do i=1, count-1
            m_x = 0d0
            m_y = 0d0
            m_z = 0d0
            H_col = kp_eivec(i,:,13) !----CHECK THIS '13'
            do k=1, 9
                chi_k = reshape([H_col(k), H_col(k+9)], [2, 1])
                m_x_comp = matmul(conjg(transpose(chi_k)),matmul(pauli_x, chi_k))
                m_y_comp = matmul(conjg(transpose(chi_k)),matmul(pauli_y, chi_k))
                m_z_comp = matmul(conjg(transpose(chi_k)),matmul(pauli_z, chi_k))
                
                m_x = m_x + real(m_x_comp(1,1))
                m_y = m_y + real(m_y_comp(1,1))
                m_z = m_z + real(m_z_comp(1,1))
            enddo
        m_k_data(:,i) = [m_x, m_y, m_z]
    enddo

!-----L projection
   !-Define angular momentum Operators 
    Lx(1, 1) = zero; Lx(1, 2) = one; Lx(1,3) = zero
    Lx(2, 1) = one; Lx(2, 2) = zero; Lx(2,3) = one
    Lx(3, 1) = zero; Lx(3, 2) = one; Lx(3,3) = zero

    Ly(1, 1) = zero; Ly(1, 2) = -im; Ly(1,3) = zero
    Ly(2, 1) = im; Ly(2, 2) = zero; Ly(2,3) = -im
    Ly(3, 1) = zero; Ly(3, 2) = im; Ly(3,3) = zero

    Lz(1, 1) = one; Lz(1, 2) = zero; Lz(1,3) = zero
    Lz(2, 1) = zero; Lz(2, 2) = zero; Lz(2,3) = zero
    Lz(3, 1) = zero; Lz(3, 2) = zero; Lz(3,3) = -one

    Lx = Lx*1/sqrt2
    Ly = Ly*1/sqrt2

   !-Defined
    do i=1, count-1
        L_x = 0d0
        L_y = 0d0
        L_z = 0d0
        H_col = kp_eivec(i,:,13)
        do k=1, 6
            phi_k = reshape([H_col((k-1)*3+1), H_col((k-1)*3+2),H_col((k-1)*3+3)], [3, 1])

            Y_lm(1,1) = (-sqrt2/2) * (phi_k(1,1) + im*phi_k(2,1));
            Y_lm(2,1) = (sqrt2/2) * (phi_k(1,1) - im*phi_k(2,1));
            Y_lm(3,1) = phi_k(3,1); 

            ! Y_lm(1,1) = (sqrt2/2) * (phi_k(2,1) - phi_k(1,1));
            ! Y_lm(2,1) = im * (sqrt2/2) * (phi_k(2,1) + phi_k(1,1));
            ! Y_lm(3,1) = phi_k(3,1);

            L_x_comp = matmul(conjg(transpose(Y_lm)),matmul(Lx, Y_lm))
            L_y_comp = matmul(conjg(transpose(Y_lm)),matmul(Ly, Y_lm))
            L_z_comp = matmul(conjg(transpose(Y_lm)),matmul(Lz, Y_lm))
            
            L_x = L_x + real(L_x_comp(1,1))
            L_y = L_y + real(L_y_comp(1,1))
            L_z = L_z + real(L_z_comp(1,1))
        enddo
        L_k_data(:,i) = [L_x, L_y, L_z]
    enddo
    
!-----Write to kmesh.dat
    do i=1, count-1
        write(100,'(8(f12.6))') kpoints(1,i), kpoints(2,i), m_k_data(1,i)/30, m_k_data(2,i)/30, m_k_data(3,i),&
                                L_k_data(1,i)/20, L_k_data(2,i)/20, L_k_data(3,i)
    enddo
    call write_plt()
    stop
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
    stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file)),' not found'
    stop

end

!-----kmesh.plt subroutine
   subroutine write_plt()
    implicit none
    open(99, file='plt/kmesh.plt')
    
    write(99, '(6(a,/),a)') &
        '#set title "Spin Projection"', &
        'set title "Orbital Angular Momentum Projection"', &
        'set terminal pdfcairo enhanced font "DejaVu" transparent fontscale 1 size 8.00in, 8.00in', &
        'set output "pdfs/kmesh.pdf"', &
        'set encoding iso_8859_1', &
        'set size ratio 0 1.0,1.0', &
        ' '
    write(99, '(9(a,/),a)') &
        'FILE = "dat_files/kmesh.dat"', &
        'set xlabel "k_x"', &
        'set ylabel "k_y"', &
        'set xrange [-0.15 : 0.15]', &
        'set yrange [ -0.15 : 0.15 ]', &
        'unset key', &
        'set ytics 0.05 scale 1 nomirror out', &
        'set mytics 2', &
        'set multiplot', &
        ' '
!-----Arrow Colour
     write(99, '(6(a,/),a)') &
        '#Palette for arrows', &
        'set palette defined ( 0 "blue", 0.5 "white", 1 "red" )', &
        'set cblabel "m_z"', &
        'set style arrow 1 head filled size screen 0.02,10,45 lt 1 lc palette', &
        ' '
!------Sort by angle and radius
    write(99, '(5(a,/),a)') &
        '#Connect two rings separately', &
        'theta(x,y) = atan2(y,x)', &
        'RingX(colX,colY,b) = (x1=column(colX), y1=column(colY), sqrt(x1**2 + y1**2) < 0.04)^b ? \', &
        '             (f ? (x2=x1,y2=y1):(f=1,x0=x1,y0=y1),$1) : NaN', &
        'angle(colX,colY)   = theta(column(colY),column(colX))', &
        ' '
    write(99, '(5(a,/),a)') &
        'set style data linespoints', &
        'set datafile missing NaN', &
        'set table $Sorted', &
        '  plot FILE u 1:2:(theta($2,$1))  smooth zsort', &
        'unset table', &
        ' '
!-----Ring plot
    write(99, '(7(a,/),a)') &
        '#Plot two rings', &
        'plot f=0 $Sorted u (RingX(1,2,0)):2:(angle(1,2)) w l ls 1 lt 5 lw 7, \', &
        '     f=0  "" u (x0):(y0):(x2-x0):(y2-y0) w vec ls 1 lt 5 lw 7 nohead, \', &
        '          "" u (RingX(1,2,1)):2:(angle(1,2)) w l ls 1 lt 5 lw 7, \', &
        '     f=0  "" u (x0):(y0):(x2-x0):(y2-y0) w vec ls 1 lt 5 lw 7 nohead, \', &
        '     FILE u ($1-$3/2):($2-$4/2):3:4:5 with vectors arrowstyle 1  #Spin Projection', &
        '#FILE u 1:2:6:7:8 with vectors arrowstyle 1  #Angular Momentum Projection', &
        'unset multiplot'
    close(99)
   end subroutine write_plt

! ! ------ gfortran -o kmesh kmesh.f90 -lblas -llapack
! gfortran -o kplot  kmesh.f90 -LC:/msys64/ucrt64/lib/lapack-3.11.0/build/lib -llapack -lblas