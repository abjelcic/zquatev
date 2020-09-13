c kako kompajlirati:
c  export OPENBLAS_NUM_THREADS=1 && gfortran -o run -cpp -O3 -mcmodel=large -march=native ProvjeraJelRadi.f HFB_blocked.f -lblas -llapack
c jasno, treba imat OpenBLAS, LAPACK i GFortran instalirane....

      program HFB
      implicit none;

      integer n,i,j,k,l;
      parameter( n = 8000 );

      double precision A(n,n);
      double precision B(n,n);
      double precision U(n,n);
      double precision V(n,n);
      double precision D(n);

      double precision Aoriginal(n,n);
      double precision Boriginal(n,n);

      double precision tol;


      double precision H(2*n,2*n);

      double precision W(2*n);
      integer LWORK, INFO;
      parameter( LWORK = 900000000 );
      double precision WORK(LWORK);
      integer LIWORK;
      parameter( LIWORK = 1000000 );
      integer IWORK(LIWORK);



      integer it1,it2,iclock_rate,iclock_max;
      double precision vrijeme;
      double precision e,ee;

      integer iseed;



      write(6,*)'Unesi iseed';
      read(*,*)iseed;
      call srand(iseed);


      write(6,*) 'Velicina problema: n = ',n;


      A = 0.D0;
      B = 0.D0;
      do i = 1 , n
          do j = 1 , i
              A(i,j) = rand();
              B(i,j) = rand();
          enddo
          B(i,i) = 0.D0;
      enddo
      do i = 1 , n
          do j = i+1 , n
              A(i,j) = + A(j,i);
              B(i,j) = - B(j,i);
          enddo
      enddo

      do j = 1 , n
          do i = 1 , n
              Aoriginal(i,j) = A(i,j);
              Boriginal(i,j) = B(i,j);
          enddo
      enddo

      do i = 1 , n
          do j = 1 , n
              H(  i  ,   j )  = + A( i , j );
              H( n+i ,   j )  = + B( i , j );
              H(   i , n+j )  = - B( i , j );
              H( n+i , n+j )  = + A( i , j );
          enddo
      enddo






      write(6,*);
      write(6,*)'------------------------------------------------';
      call system_clock( it1, iclock_rate, iclock_max );
      call SymplecticHamiltonianEig( 'V', n,  A,n,B,n, D, U,n,V,n,
     &                                 WORK,LWORK,IWORK,LIWORK,INFO);
      if(INFO.ne.0) stop 'INFO nije jednak nuli...';
      call system_clock( it2, iclock_rate, iclock_max );
      vrijeme = DBLE(it2-it1)/DBLE(iclock_rate);
      write(6,'(a,f10.4,a)')'vrijeme(moje)= ',vrijeme,'s; ';
      write(6,*)'------------------------------------------------';
      write(6,*);



      write(6,*);
      write(6,*)'------------------------------------------------';
      call system_clock( it1, iclock_rate, iclock_max );
      call dsyevd('V','l',2*n,H,2*n,W,WORK,LWORK,IWORK,LIWORK,INFO);
      if(INFO.ne.0) stop 'INFO nije jednak nuli...';
      call system_clock( it2, iclock_rate, iclock_max );
      vrijeme = DBLE(it2-it1)/DBLE(iclock_rate);
      write(6,'(a,f10.4,2a)')'vrijeme(dsyevd)= ',vrijeme,'s ';
      write(6,*)'------------------------------------------------';
      write(6,*);









      call CheckOrthogonality(n,U,n,V,n);
      call CheckEigen(n, Aoriginal,n, Boriginal,n, U,n, V,n, D);

      write(6,*);
      write(6,*)'Usporedba mojih i dsyev eigenvalues:';
      tol = 1.d-12;
      do i = 1 , n
          e  = D(i);
          ee = W(2*i);

          if( DABS(e-ee)/DABS(e) .gt. tol )
     &        write(6,*)i,ee,e,DABS(e-ee)/DABS(e);

c          write(6,*)i,ee,e,DABS(e-ee)/DABS(e);
      enddo






      return;
      END;



      subroutine CheckOrthogonality( n, U,LDU, V,LDV )

      implicit none;

      integer n,LDU,LDV;
      double precision U(LDU,n);
      double precision V(LDV,n);

      double precision Q(2*n,2*n);
      double precision QQ(2*n,2*n);
      double precision Fnorm;
      integer i,j;

      do j = 1 , n
          do i = 1 , n
              Q(i  ,j  ) = +U(i,j);
              Q(i  ,j+n) = -V(i,j);
              Q(i+n,j  ) = +V(i,j);
              Q(i+n,j+n) = +U(i,j);
          enddo
      enddo




      call dsyrk('U','N',2*n,2*n,1.d0,Q,2*n,0.d0,QQ,2*n);

      do i = 1 , 2*n
          QQ(i,i) = QQ(i,i) - 1.d0;
      enddo

      Fnorm = 0.d0;
      do j = 1 , 2*n
          do i = 1 , 2*n
              if( i .gt. j ) then
                  Fnorm = Fnorm + QQ(j,i)**2.d0;
              else
                  Fnorm = Fnorm + QQ(i,j)**2.d0;
              endif
          enddo
      enddo
      Fnorm = DSQRT(Fnorm);

      write(6,*)'||QQ^T-I||_F = ',Fnorm;





      call dsyrk('U','T',2*n,2*n,1.d0,Q,2*n,0.d0,QQ,2*n);

      do i = 1 , 2*n
          QQ(i,i) = QQ(i,i) - 1.d0;
      enddo

      Fnorm = 0.d0;
      do j = 1 , 2*n
          do i = 1 , 2*n
              if( i .gt. j ) then
                  Fnorm = Fnorm + QQ(j,i)**2.d0;
              else
                  Fnorm = Fnorm + QQ(i,j)**2.d0;
              endif
          enddo
      enddo
      Fnorm = DSQRT(Fnorm);

      write(6,*)'||Q^TQ-I||_F = ',Fnorm;


      return;
      end

      subroutine CheckEigen(n, A,LDA, B,LDB, U,LDU,V,LDV, D )
      implicit none;

      integer n,LDA,LDB,LDU,LDV;
      double precision A(LDA,n);
      double precision B(LDB,n);
      double precision U(LDU,n);
      double precision V(LDV,n);
      double precision D(n);

      double precision Q(2*n,2*n);
      double precision H(2*n,2*n);
      double precision HQ(2*n,2*n);
      double precision QtHQ(2*n,2*n);
      double precision Fnorm;
      integer i,j;

      ! [U,-V;V,U]^T * [A,-B;B,A] * [U,-V;V,U] = [D,0;0,D] provjeravam

      do j = 1 , n
          do i = 1 , n
              Q(i  ,j  ) = +U(i,j);
              Q(i  ,j+n) = -V(i,j);
              Q(i+n,j  ) = +V(i,j);
              Q(i+n,j+n) = +U(i,j);
          enddo
      enddo
      do j = 1 , n
          do i = 1 , n
              H(i  ,j  ) = +A(i,j);
              H(i  ,j+n) = -B(i,j);
              H(i+n,j  ) = +B(i,j);
              H(i+n,j+n) = +A(i,j);
          enddo
      enddo

      call dgemm('N','N',2*n,2*n,2*n,1.d0,H,2*n,Q,2*n,0.d0,HQ,2*n);
      call dgemm('T','N',2*n,2*n,2*n,1.d0,Q,2*n,HQ,2*n,0.d0,QtHQ,2*n);

      do i = 1 , n
          QtHQ(  i,  i) = QtHQ(  i,  i) - D(i);
          QtHQ(n+i,n+i) = QtHQ(n+i,n+i) - D(i);
      enddo

      Fnorm = 0.d0;
      do j = 1 , 2*n
          do i = 1 , 2*n
              Fnorm = Fnorm + QtHQ(i,j)**2.d0;
          enddo
      enddo
      Fnorm = DSQRT(Fnorm);

      write(6,*)'||Q^T*H*Q-I||_F = ',Fnorm;

      return;
      end














