
        !abjelcic@phy.hr za bugove
        !provjeri sa flagovima jel se negdje gazi po memoriji
        !provjeri sa puno nxn random matirca tako da n bude svakakav broj, i za svakave bsize testiraj isto
        !napisat lapack-style input/output argumente sta su
      subroutine SymplecticHamiltonianEig( JOBZ , n  ,
     &                                     A , LDA , B , LDB ,
     &                                     W ,
     &                                     U , LDU , V , LDV ,
     &                                     WORK,LWORK, IWORK,LIWORK,
     &                                     INFO                     )
      ! A je simetricna matrica, B antisim., samo donji trokut
      ! U i V ortogonalne
      ! JOBZ je 'V' ili 'N' ovismo hoces/neces sv. vektore
      ! u U i V spremi fino matrice t.d. je za Q=[U,-V;V,U]
      ! vrijedi Qt*[A,-B;B,A]*Q = [D,0;0,D], gdje je D dijagonalna
      ! i da je Q'*Q=Q*Q'=I.


      implicit none;

      integer Aticks,Bticks;
      common /measureAB/ Aticks,Bticks;

      character JOBZ;
      integer n;
      integer LDA;
      integer LDB;
      integer LDU;
      integer LDV;
      integer LWORK;
      integer LIWORK;
      integer INFO;
      double precision A( LDA , n );
      double precision B( LDB , n );
      double precision W(       n );
      double precision U( LDU , n );
      double precision V( LDV , n );
      double precision  WORK( LWORK);
      integer          IWORK(LIWORK);


      integer i, j, ibA, ibB, lmt;
      double precision c(n-1);
      double precision s(n-1);
      double precision tau1(n); !vjv moze n-1
      double precision tau2(n);
      double precision D(n);
      double precision E(n-1);

      integer bsizeA, bsizeB;
      parameter(bsizeA=32);
      parameter(bsizeB=32);
      double precision xA(n,bsizeA);
      double precision yA(n,bsizeA);
      double precision xB(n,bsizeB);
      double precision yB(n,bsizeB);

      double precision Ae1(n);
      double precision Av1(n);
      double precision Av2(n);
      double precision Be1(n);
      double precision Bv1(n);
      double precision Bv2(n);!tu svuda valjda moze i n-1 al ajde...
      double precision WORKA(bsizeA);
      double precision WORKB(bsizeB);

      !prebaciti ovaj dio u posebnu rutinu da ako je JOBZ='V' da se ova memorija ne trosi
      integer bsizeQ;
      parameter(bsizeQ=64);
      double precision V1V2(n,bsizeQ+bsizeQ);
      double precision V1V2E_U(bsizeQ+bsizeQ+bsizeQ,n);
      double precision V1V2E_V(bsizeQ+bsizeQ+bsizeQ,n);
      double precision T11(bsizeQ,bsizeQ);
      double precision T12(bsizeQ,bsizeQ);
      double precision T13(bsizeQ,bsizeQ);
      double precision T21(bsizeQ,bsizeQ);
      double precision T22(bsizeQ,bsizeQ);
      double precision T23(bsizeQ,bsizeQ);
      double precision T31(bsizeQ,bsizeQ);
      double precision T32(bsizeQ,bsizeQ);
      double precision T33(bsizeQ,bsizeQ);
      double precision T(3*bsizeQ,3*bsizeQ);
      double precision R1(bsizeQ,bsizeQ);
      double precision R2(bsizeQ,bsizeQ);
      double precision R3(bsizeQ,bsizeQ);
      double precision R(3*bsizeQ,bsizeQ);
      double precision S1(bsizeQ,bsizeQ);
      double precision S2(bsizeQ,bsizeQ);
      double precision S3(bsizeQ,bsizeQ);
      double precision S123(bsizeQ,3*bsizeQ);
      integer i1,i2,k,nV1V2,kV1V2;
      double precision f1(bsizeQ);
      double precision f2(bsizeQ);
      double precision g1(bsizeQ);
      double precision g2(bsizeQ);
      double precision g3(bsizeQ);
      double precision g4(bsizeQ);
      double precision h1(bsizeQ);
      double precision h2(bsizeQ);
      double precision h3(bsizeQ);
      double precision X1(n,3*bsizeQ);
      double precision X2(n,1*bsizeQ);
      double precision Y1(bsizeQ,n);
      double precision Y2(bsizeQ,n);

      integer it1,it2,iclock_rate,iclock_max;
      double precision vrijeme;
      integer ticks1, ticks2, ticks3, ticks4, ticks5, ticks6, ticks7;
      ticks1 = 0;
      ticks2 = 0;
      ticks3 = 0;
      ticks4 = 0;
      ticks5 = 0;
      ticks6 = 0;
      ticks7 = 0;

      write(6,*)'bsizeA =',bsizeA,', bsizeB =',bsizeB;
      write(6,*)'bsizeQ =',bsizeQ;

      !tu provjere n>=0 i slicno...takodjer testirati na tim slucajevima n=1 i slicno
      !quick return, vidi kak je to u lapacku napravljeno
      if( JOBZ.ne.'V' .and. JOBZ.ne.'N' ) then
          INFO = -1;
          stop 'JOBZ nije V niti N';
      endif


      ibA = 1;
      ibB = 1;
      lmt = min( ((n-2)/bsizeA)*bsizeA , ((n-2)/bsizeB)*bsizeB );
      do i = 1 , lmt

          if( i-ibA.eq.bsizeA ) then
              call system_clock( it1, iclock_rate, iclock_max );
              call RefreshA( n-i, i-ibA,
     &                       A(i+1,i+1),LDA,
     &                       B(i+1,ibA),LDB,  A(i+1,ibA),LDA,
     &                       xA(i+1,1),n,   yA(i+1,1),n
     &                      );
              ibA = i;
              call system_clock( it2, iclock_rate, iclock_max );
              ticks3 = ticks3 + (it2-it1);
          endif
          if( i-ibB.eq.bsizeB ) then
              call system_clock( it1, iclock_rate, iclock_max );
              call RefreshB( n-i, i-ibB,
     &                       B(i+1,i+1),LDB,
     &                       B(i+1,ibB),LDB,  A(i+1,ibB),LDA,
     &                       xB(i+1,1),n,   yB(i+1,1),n
     &                      );
              ibB = i;
              call system_clock( it2, iclock_rate, iclock_max );
              ticks4 = ticks4 + (it2-it1);
          endif



          call GenerateHGH( n-i, A(i,i), B(i,i),
     &                      tau1(i),tau2(i),
     &                      c(i),s(i),
     &                      D(i),E(i) );





      call system_clock( it1, iclock_rate, iclock_max );
          call Av1v2e1( n-i, i-ibA,
     &                  A(i+1,i+1),LDA,
     &                  B(i+1,i), A(i+1,i),
     &                  B(i+1,ibA),LDB,  A(i+1,ibA),LDA,
     &                  xA(i+1,1),n,  yA(i+1,1),n,
     &                  WORKA,
     &                  Av1(i+1), Av2(i+1), Ae1(i+1)
     &                 );
      call system_clock( it2, iclock_rate, iclock_max );
      ticks1 = ticks1 + (it2-it1);




      call system_clock( it1, iclock_rate, iclock_max );
          call Bv1v2e1( n-i, i-ibB,
     &                  B(i+1,i+1),LDB,
     &                  B(i+1,i), A(i+1,i),
     &                  B(i+1,ibB),LDB,  A(i+1,ibB),LDB,
     &                  xB(i+1,1),n,  yB(i+1,1),n,
     &                  WORKB,
     &                  Bv1(i+1), Bv2(i+1), Be1(i+1)
     &                 );
      call system_clock( it2, iclock_rate, iclock_max );
      ticks2 = ticks2 + (it2-it1);




          call UpdateXY( n-i,
     &                   Av1(i+1), Av2(i+1), Ae1(i+1),
     &                   Bv1(i+1), Bv2(i+1), Be1(i+1),
     &                   B(i+1,i),  A(i+1,i),
     &                   tau1(i), tau2(i),
     &                   c(i), s(i),
     &                   xA(i+1,i-ibA+1), yA(i+1,i-ibA+1),
     &                   xB(i+1,i-ibB+1), yB(i+1,i-ibB+1),
     &                   A(i+1,i+1), B(i+1,i+1)
     &                  );

      enddo

      vrijeme = DBLE(ticks1)/DBLE(iclock_rate);
      write(6,'(a,f10.4,2a)')'vrijeme1 Av1v2e1= ',vrijeme,'s ';
      vrijeme = DBLE(ticks2)/DBLE(iclock_rate);
      write(6,'(a,f10.4,2a)')'vrijeme2 Bv1v2e1= ',vrijeme,'s ';
      vrijeme = DBLE(ticks3)/DBLE(iclock_rate);
      write(6,'(a,f10.4,2a)')'vrijeme3 RefreshA= ',vrijeme,'s ';
      vrijeme = DBLE(ticks4)/DBLE(iclock_rate);
      write(6,'(a,f10.4,2a)')'vrijeme4 RefreshB= ',vrijeme,'s ';

      call RefreshA( n-i, i-ibA,
     &               A(i+1,i+1),LDA,
     &               B(i+1,ibA),LDB,  A(i+1,ibA),LDA,
     &               xA(i+1,1),n,   yA(i+1,1),n
     &              );
      call RefreshB( n-i, i-ibB,
     &               B(i+1,i+1),LDB,
     &               B(i+1,ibB),LDB,  A(i+1,ibB),LDA,
     &               xB(i+1,1),n,   yB(i+1,1),n
     &              );




      !dalje nastavljamo neblock varijantu...
      do i = lmt+1 , n-2
          !iz vektora a=A(i:end,i) i b=B(i:end,i) generiram v1,v2,tau1,tau2,c,s,d,e...
          call GenerateHGH( n-i, A(i,i),B(i,i),
     &          tau1(i),tau2(i), c(i),s(i), D(i),E(i) );

          call AsymmHAH(n-i, B(i+1,i+1),LDB, B(i+1,i),tau1(i) );
          call  SymmHAH(n-i ,A(i+1,i+1),LDA, B(i+1,i),tau1(i) );
          call drot( n-(i+1) , B(i+2,i+1),1 , A(i+2,i+1),1 , c(i),s(i));
          call  SymmHAH(n-i ,A(i+1,i+1),LDA, A(i+1,i),tau2(i) );
          call AsymmHAH(n-i ,B(i+1,i+1),LDB, A(i+1,i),tau2(i) );

      enddo
      call drotg( A(n,n-1) , B(n,n-1) , c(n-1) , s(n-1) );
      D(n-1) = A( n-1 , n-1 );
      E(n-1) = A(   n , n-1 );
      D(n  ) = A(   n , n   );
      tau1(n-1) = 0.d0;
      tau2(n-1) = 0.d0;
      A(n,n-1)  = 1.d0;
      B(n,n-1)  = 1.d0;



      if( JOBZ .eq. 'N' ) then

          !vidi jel ima sta bolje od dsterf...
          call dsterf( n, D, E, INFO );
          if(INFO.ne.0) stop 'dsterf nije radio...';
          do i = 1 , n
              W(i) = D(i);
          enddo

      else



          U = 0.d0;
          V = 0.d0;

          !tridiag eigensolve
          call system_clock( it1, iclock_rate, iclock_max );
          call dstedc('I',n,D,E,U,LDU,WORK,-1,IWORK,-1,INFO);
          if( INT(WORK(1)+0.5d0) .gt. LWORK ) then
              INFO = -1;
              return;
          endif
          if( IWORK(1) .gt. LIWORK ) then
              INFO = -1;
              return;
          endif
          if( INFO .ne. 0 ) then
              return;
          endif
          call dstedc('I',n,D,E,U,LDU,WORK,LWORK,IWORK,LIWORK,INFO);
          if( INFO .ne. 0 ) then
              stop 'failo dstedc';
              return;
          endif
          do i = 1 , n
              W(i) = D(i);
          enddo
          call system_clock( it2, iclock_rate, iclock_max );
          vrijeme = DBLE(it2-it1)/DBLE(iclock_rate);
          write(6,'(a,f10.4,2a)')'vrijeme tridiag = ',vrijeme,'s ';




          !Kressner compact WY-like representation, paper referencu daj
          lmt = ((n-1)/bsizeQ)*bsizeQ + 1;
          do i1 = lmt , 1 , -bsizeQ
              i2 = min( n-1 , i1+bsizeQ-1 );


              T11 = 0.d0; T12 = 0.d0; T13 = 0.d0;
              T21 = 0.d0; T22 = 0.d0; T23 = 0.d0;
              T31 = 0.d0; T32 = 0.d0; T33 = 0.d0;
              R1  = 0.d0; R2  = 0.d0; R3  = 0.d0;
              S1  = 0.d0; S2  = 0.d0; S3  = 0.d0;

              nV1V2 = n-i1;
              kV1V2 = i2-i1+1;
              do j = 1 , kV1V2
                  do i = 1 , j-1
                      V1V2(i,j      ) = 0.d0;
                      V1V2(i,j+kV1V2) = 0.d0;
                  enddo
                  do i = j , nV1V2
                      V1V2(i,j      ) = B( i1-1 + i+1 , i1-1 + j );
                      V1V2(i,j+kV1V2) = A( i1-1 + i+1 , i1-1 + j );
                  enddo
              enddo


              call system_clock( it1, iclock_rate, iclock_max );
              do i = i1 , i2
                  k = i-i1;
                  !v1 = V1V2( k+1:end , k+1         )
                  !v2 = V1V2( k+1:end , k+1 + kV1V2 )


                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! ** I - tau1 * v1*v1'
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  f1 = 0.d0; f2 = 0.d0;
                  g1 = 0.d0; g2 = 0.d0; g3 = 0.d0;
                  h1 = 0.d0; h2 = 0.d0;
                  !f1 = V1'*v1
                  call dgemv('T', nV1V2-k,k, 1.d0,V1V2(k+1,1),n,
     &                        V1V2(k+1,k+1),1,  0.d0,f1,1  );
                  !f2 = V2'*v1
                  call dgemv('T', nV1V2-k,k, 1.d0,V1V2(k+1,kV1V2+1),n,
     &                        V1V2(k+1,k+1),1,  0.d0,f2,1  );

                  !g1 = T11*f1 + T13*f2
                  call dgemv('N',k,k, 1.d0,T11,bsizeQ, f1,1, 0.d0,g1,1); !ovo bi moglo sa dtrmv
                  call dgemv('N',k,k, 1.d0,T13,bsizeQ, f2,1, 1.d0,g1,1);
                  !g2 = T21*f1 + T23*f2
                  call dgemv('N',k,k, 1.d0,T21,bsizeQ, f1,1, 0.d0,g2,1); !ovo bi moglo sa dtrmv
                  call dgemv('N',k,k, 1.d0,T23,bsizeQ, f2,1, 1.d0,g2,1);
                  !g3 = T31*f1 + T33*f2
                  call dgemv('N',k,k, 1.d0,T31,bsizeQ, f1,1, 0.d0,g3,1); !ovo bi moglo sa dtrmv
                  call dgemv('N',k,k, 1.d0,T33,bsizeQ, f2,1, 1.d0,g3,1);

                  T11(k+1,k+1) = tau1(i);
                  do j = 1 , k
                      T11(j,k+1) = -tau1(i) * g1(j);
                      T21(j,k+1) = -tau1(i) * g2(j);
                      T31(j,k+1) = -tau1(i) * g3(j);
                  enddo

                  !h1 = S1*f1
                  call dgemv('N',k,k, 1.d0,S1,bsizeQ, f1,1, 0.d0,h1,1); !ovo bi moglo sa dtrmv
                  !h2 = S3*f2
                  call dgemv('N',k,k, 1.d0,S3,bsizeQ, f2,1, 0.d0,h2,1); !ovo bi moglo sa dtrmv
                  do j = 1 , k
                      S1(j,k+1) = -tau1(i) * ( h1(j) + h2(j) );
                  enddo
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! ** [C,-S;S,C], C = I + (c-1)*e1*e1', S = s*e1*e1'
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  f1 = 0.d0; f2 = 0.d0;
                  g1 = 0.d0; g2 = 0.d0; g3 = 0.d0; g4 = 0.d0;
                  h1 = 0.d0; h2 = 0.d0; h3 = 0.d0;
                  !f1 = V1'*e1
                  do j = 1 , k+1
                      f1(j) = V1V2( k+1 ,       j );
                  enddo
                  !f2 = V2'*e1
                  do j = 1 , k
                      f2(j) = V1V2( k+1 , kV1V2+j );
                  enddo

                 !g1 = T11*f1 + T13*f2
                 call dgemv('N',k+1,k+1,1.d0,T11,bsizeQ,f1,1,0.d0,g1,1); !ovo bi moglo sa dtrmv
                 call dgemv('N',k+1,k,  1.d0,T13,bsizeQ,f2,1,1.d0,g1,1);
                 !g2 = T21*f1 + T23*f2
                 call dgemv('N',k,k+1,1.d0,T21,bsizeQ, f1,1, 0.d0,g2,1); !ovo bi moglo sa dtrmv
                 call dgemv('N',k,k,  1.d0,T23,bsizeQ, f2,1, 1.d0,g2,1);
                 !g3 = T31*f1 + T33*f2
                 call dgemv('N',k,k+1,1.d0,T31,bsizeQ, f1,1, 0.d0,g3,1); !ovo bi moglo sa dtrmv
                 call dgemv('N',k,k,  1.d0,T33,bsizeQ, f2,1, 1.d0,g3,1);
                 !g4 = S1*f1 + S3*f2
                 call dgemv('N',k,k+1, 1.d0,S1,bsizeQ, f1,1, 0.d0,g4,1); !ovo bi moglo sa dtrmv
                 call dgemv('N',k,k,   1.d0,S3,bsizeQ, f2,1, 1.d0,g4,1);
                 !h1 = R1*g4
                 call dgemv('N',k+1,k, 1.d0,R1,bsizeQ, g4,1, 0.d0,h1,1); !ovo bi moglo sa dtrmv
                 !h2 = R2*g4
                 call dgemv('N',k,k, 1.d0,R2,bsizeQ, g4,1, 0.d0,h2,1); !ovo bi moglo sa dtrmv
                 !h3 = R3*g4
                 call dgemv('N',k,k, 1.d0,R3,bsizeQ, g4,1, 0.d0,h3,1); !ovo bi moglo sa dtrmv


                  T22(k+1,k+1) = -(c(i)-1.d0);
                  do j = 1 , k+1
                      T12(j,k+1) = (c(i)-1.d0)*g1(j) + s(i)*h1(j);
                  enddo
                  do j = 1 , k
                      T22(j,k+1) = (c(i)-1.d0)*g2(j) + s(i)*h2(j);
                      T32(j,k+1) = (c(i)-1.d0)*g3(j) + s(i)*h3(j);
                  enddo

                  R2(k+1,k+1) = 1.d0;
                  do j = 1 , k+1
                      R1(j,k+1) = - g1(j);
                  enddo
                  do j = 1 , k
                      R2(j,k+1) = - g2(j);
                      R3(j,k+1) = - g3(j);
                  enddo

                  S2(k+1,k+1) = s(i);
                  do j = 1 , k
                      S2(j,k+1) = (c(i)-1.d0)*g4(j);
                  enddo
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! ** I - tau2 * v2*v2'
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  f1 = 0.d0; f2 = 0.d0;
                  g1 = 0.d0; g2 = 0.d0; g3 = 0.d0;
                  h1 = 0.d0; h2 = 0.d0; h3 = 0.d0;
                  !f1 = V1'*v2
                  call dgemv('T', nV1V2-k,k+1, 1.d0,V1V2(k+1,1),n,
     &                        V1V2(k+1,kV1V2+k+1),1,  0.d0,f1,1  );
                  !f2 = V2'*v2
                  call dgemv('T', nV1V2-k,k, 1.d0,V1V2(k+1,kV1V2+1),n,
     &                        V1V2(k+1,kV1V2+k+1),1,  0.d0,f2,1  );

                 !g1 = T11*f1 + T12*e + T13*f2
                 do j = 1 , k+1
                     g1(j) = T12(j,k+1);
                 enddo
                 call dgemv('N',k+1,k+1,1.d0,T11,bsizeQ,f1,1,1.d0,g1,1); !ovo bi moglo sa dtrmv
                 call dgemv('N',k+1,k  ,1.d0,T13,bsizeQ,f2,1,1.d0,g1,1);
                 !g2 = T21*f1 + T22*e + T23*f2
                 do j = 1 , k+1
                     g2(j) = T22(j,k+1);
                 enddo
                 call dgemv('N',k+1,k+1,1.d0,T21,bsizeQ,f1,1,1.d0,g2,1); !ovo bi moglo sa dtrmv
                 call dgemv('N',k+1,k,  1.d0,T23,bsizeQ,f2,1,1.d0,g2,1);
                 !g3 = T31*f1 + T32*e + T33*f2
                 do j = 1 , k
                     g3(j) = T32(j,k+1);
                 enddo
                 call dgemv('N',k,k+1,1.d0,T31,bsizeQ,f1,1,1.d0,g3,1); !ovo bi moglo sa dtrmv
                 call dgemv('N',k,k,  1.d0,T33,bsizeQ,f2,1,1.d0,g3,1);

                  T33(k+1,k+1) = tau2(i);
                  do j = 1 , k+1
                      T13(j,k+1) = -tau2(i) * g1(j);
                      T23(j,k+1) = -tau2(i) * g2(j);
                  enddo
                  do j = 1 , k
                      T33(j,k+1) = -tau2(i) * g3(j);
                  enddo


                  !h1 = S1*f1
                  call dgemv('N',k+1,k+1,1.d0,S1,bsizeQ,f1,1,0.d0,h1,1); !ovo bi moglo sa dtrmv
                  !h2 = S2*e
                  do j = 1 , k+1
                      h2(j) = S2(j,k+1);
                  enddo
                  !h3 = S3*f2
                  call dgemv('N',k+1,k,1.d0,S3,bsizeQ,f2,1,0.d0,h3,1); !ovo bi moglo sa dtrmv

                  do j = 1 , k+1
                      S3(j,k+1) = -tau2(i) * ( h1(j) + h2(j) + h3(j) );
                  enddo
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


              enddo
              call system_clock( it2, iclock_rate, iclock_max );
              ticks5 = ticks5 + (it2-it1);



              call system_clock( it1, iclock_rate, iclock_max );
              ! V1V2E_U = [V1',V2',E'] * U
              call dgemm( 'T','N', 2*kV1V2,n,nV1V2,
     &                     1.d0, V1V2,n,
     &                           U(i1+1,1),LDU,
     &                     0.d0, V1V2E_U,bsizeQ+bsizeQ+bsizeQ );
              do j = 1 , n
                  do i = 1 , kV1V2
                      V1V2E_U( 2*kV1V2+i , j ) = U( i1+i , j );
                  enddo
              enddo

              ! V1V2E_V = [V1',V2',E'] * V
              call dgemm( 'T','N', 2*kV1V2,n,nV1V2,
     &                     1.d0, V1V2,n,
     &                           V(i1+1,1),LDV,
     &                     0.d0, V1V2E_V,bsizeQ+bsizeQ+bsizeQ );
              do j = 1 , n
                  do i = 1 , kV1V2
                      V1V2E_V( 2*kV1V2+i , j ) = V( i1+i , j );
                  enddo
              enddo

              ! X1 = [V1,V2,E] * [ T11,T13,T12; T31,T33,T32; T21,T23,T22 ]
              do j = 1 , kV1V2
                  do i = 1 , kV1V2
                      T(i      ,j        ) = T11(i,j);
                      T(i      ,j+  kV1V2) = T13(i,j);
                      T(i      ,j+2*kV1V2) = T12(i,j);
                      T(i+kV1V2,j        ) = T31(i,j);
                      T(i+kV1V2,j+  kV1V2) = T33(i,j);
                      T(i+kV1V2,j+2*kV1V2) = T32(i,j);
                  enddo
              enddo
              call dgemm('N','N', nV1V2,3*kV1V2, 2*kV1V2,
     &                      1.d0,V1V2,n,
     &                           T,3*bsizeQ,
     &                     0.d0,X1,n                     );
              do j = 1 , kV1V2
                  do i = 1 , kV1V2
                      X1(i,j        ) = X1(i,j        ) + T21(i,j);
                      X1(i,j  +kV1V2) = X1(i,j  +kV1V2) + T23(i,j);
                      X1(i,j+2*kV1V2) = X1(i,j+2*kV1V2) + T22(i,j);
                  enddo
              enddo

              ! X2 = [V1,V2,E] * [ R1 ; R3 ; R2 ]
              do j = 1 , kV1V2
                  do i = 1 , kV1V2
                      R(i      ,j) = R1(i,j);
                      R(i+kV1V2,j) = R3(i,j);
                  enddo
              enddo
              call dgemm('N','N', nV1V2,kV1V2, 2*kV1V2,
     &                      1.d0,V1V2,n,
     &                           R,3*bsizeQ,
     &                     0.d0,X2,n                     );
              do j = 1 , kV1V2
                  do i = 1 , kV1V2
                      X2(i,j) = X2(i,j) + R2(i,j);
                  enddo
              enddo

              ! Y1 = [S1,S3,S2] * V1V2E_V
              ! Y2 = [S1,S3,S2] * V1V2E_U
              do j = 1 , kV1V2
                  do i = 1 , kV1V2
                      S123(i,j        ) = S1(i,j);
                      S123(i,j  +kV1V2) = S3(i,j);
                      S123(i,j+2*kV1V2) = S2(i,j);
                  enddo
              enddo
              call dgemm('N','N',kV1V2,n, 3*kV1V2,
     &                     1.d0,S123,bsizeQ,
     &                      V1V2E_V,bsizeQ+bsizeQ+bsizeQ,
     &                     0.d0,Y1,bsizeQ );
              call dgemm('N','N',kV1V2,n, 3*kV1V2,
     &                     1.d0,S123,bsizeQ,
     &                      V1V2E_U,bsizeQ+bsizeQ+bsizeQ,
     &                     0.d0,Y2,bsizeQ );

              ! U = U - X1*V1V2E_U - X2*Y1
              call dgemm('N','N', nV1V2,n, 3*kV1V2,
     &                     -1.d0,X1,n,
     &                     V1V2E_U,bsizeQ+bsizeQ+bsizeQ,
     &                     1.d0,U(i1+1,1),LDU            );
              call dgemm('N','N', nV1V2,n, kV1V2,
     &                     -1.d0,X2,n,
     &                     Y1,bsizeQ,
     &                     1.d0,U(i1+1,1),LDU            );

              ! V = V - X1*V1V2E_V + X2*Y2
              call dgemm('N','N', nV1V2,n, 3*kV1V2,
     &                     -1.d0,X1,n,
     &                     V1V2E_V,bsizeQ+bsizeQ+bsizeQ,
     &                     1.d0,V(i1+1,1),LDV            );
              call dgemm('N','N', nV1V2,n, kV1V2,
     &                     +1.d0,X2,n,
     &                     Y2,bsizeQ,
     &                     1.d0,V(i1+1,1),LDV            );


              call system_clock( it2, iclock_rate, iclock_max );
              ticks6 = ticks6 + (it2-it1);

          enddo


      vrijeme = DBLE(ticks5)/DBLE(iclock_rate);
      write(6,'(a,f10.4,2a)')'TSR               = ',vrijeme,'s ';
      vrijeme = DBLE(ticks6)/DBLE(iclock_rate);
      write(6,'(a,f10.4,2a)')'UV updates=       =',vrijeme,'s ';




      endif


      INFO = 0;

      return;
      end;







      !Probaj prije same redukcije na tridijagonalnu reskalirat matricu
      !vidi kako je to napravljeno u dsyevd rutini na primjer, prije same tridiag
      !redukcije se reskalira matrica, a nakon se reskaliraju sv. vrijednosti

      subroutine SymmHAH( n , A , LDA , v , tau )
      ! A simetricna, samo donji trokut
      ! H = I - tau * vv',
      implicit none;

      integer n,i;
      integer LDA;

      double precision A(LDA,n);
      double precision v(n);
      double precision tau;

      double precision y(n);
      double precision alpha;
      double precision ddot;

      call dsymv( 'l' , n , tau , A,LDA , v,1 , 0.D0 , y,1 );
      alpha = -0.5D0 * tau * ddot(n,v,1,y,1);
      call daxpy( n , alpha , v , 1 , y , 1 );
      call dsyr2( 'l' , n , -1.D0 , y,1 , v,1 , A,LDA );

      return;
      end;

      subroutine AsymmHAH( n , A , LDA , v , tau )
      ! A antisimetricna, samo donji trokut
      ! H = I - tau * vv',
      implicit none;
      integer SIZE;
      parameter(SIZE=10000);
      integer n;
      integer LDA;

      double precision A(LDA,n);
      double precision v(n);
      double precision tau;

      integer i,j;
      double precision temp1, temp2;
      double precision y1(SIZE);
      double precision y2(SIZE);
      double precision y(SIZE);
      double precision alpha;
      double precision ddot;



      call dcopy(n, v,1, y1,1 );
      call dcopy(n, v,1, y2,1 );
      call dtrmv('l','n', 'u', n,A,LDA, y1,1);
      call dtrmv('l','t', 'u', n,A,LDA, y2,1);
      call daxpy(n,-1.D0, y2,1, y1,1);
      call dscal(n,tau,y1,1);
      do j = 1 , n
          do i = j+1 , n
              A(i,j) = A(i,j) + v(i)*y1(j) - y1(i)*v(j);
          enddo
      enddo



      return;
      end;

      !testirati kolko je uopce ova moja rutina brza/spora u usporedbi sa dsymv...
      subroutine dasymv( n ,  A,LDA ,  x,  Ax)
      ! A antisimetricna, samo donji trokut
      implicit none;

      integer LDA,n;
      double precision A(LDA,n);
      double precision x(n);
      double precision Ax(n);

      double precision tmp(n);

!napravit cu Ax=(L+L')*x = L*x - L'*x, gdje je L donji trokut
!mogu staviti da L ima unit dijagonalu...svejedno je
      call dcopy(n,  x,1,  Ax,1);
      call dcopy(n,  x,1,  tmp,1);
      call dtrmv('l','n', 'u', n , A,LDA, Ax,1 ); !primijeti da dijagonala nije
      call dtrmv('l','t', 'u', n , A,LDA, tmp,1 ); !bitna jer se skrati...
      call daxpy(n, -1.D0, tmp,1, Ax,1);

      return;
      end;

      subroutine Av1v2e1(  n , nb  ,
     &                     A,LDA,
     &                     vv1, vv2,
     &                     V1,LDV1,   V2,LDV2,
     &                     XA,LDXA,   YA,LDYA,
     &                     WORK,
     &                     Av1, Av2, Ae1
     &                   )

      implicit none;

      integer n,nb,LDA,LDV1,LDV2,LDXA,LDYA;
      double precision A(LDA,n);
      double precision V1(LDV1,nb);
      double precision V2(LDV1,nb);
      double precision XA(LDXA,nb);
      double precision YA(LDYA,nb);

      double precision vv1(n);
      double precision vv2(n);
      double precision Av1(n);
      double precision Av2(n);
      double precision Ae1(n);

      double precision WORK(nb);

      double precision one;
      double precision zero;
      parameter( one  = 1.0d+0 );
      parameter( zero = 0.0d+0 );


      call dsymv( 'l' , n,  one  , A,LDA,  vv1,1,  zero,  Av1,1 ); !this is memory and computationally
      call dsymv( 'l' , n,  one  , A,LDA,  vv2,1,  zero,  Av2,1 ); !expensive part
      call dcopy( n , A(1,1),1,  Ae1,1 );

      if(nb.eq.0) return;

      !.......................tu radim Av1
      call dgemv('t', n,nb, +one,  YA,LDYA,   vv1,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  V2,LDV2,  WORK,1,  one,  Av1,1 );
      call dgemv('t', n,nb, +one,  V2,LDV2,   vv1,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  YA,LDYA,  WORK,1,  one,  Av1,1 );
      call dgemv('t', n,nb, +one,  XA,LDXA,   vv1,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  V1,LDV1,  WORK,1,  one,  Av1,1 );
      call dgemv('t', n,nb, +one,  V1,LDV1,   vv1,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  XA,LDXA,  WORK,1,  one,  Av1,1 );


      !.......................tu radim Av2
      call dgemv('t', n,nb, +one,  YA,LDYA,   vv2,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  V2,LDV2,  WORK,1,  one,  Av2,1 );
      call dgemv('t', n,nb, +one,  V2,LDV2,   vv2,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  YA,LDYA,  WORK,1,  one,  Av2,1 );
      call dgemv('t', n,nb, +one,  XA,LDXA,   vv2,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  V1,LDV1,  WORK,1,  one,  Av2,1 );
      call dgemv('t', n,nb, +one,  V1,LDV1,   vv2,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  XA,LDXA,  WORK,1,  one,  Av2,1 );


      !.......................tu radim Ae1
      call dgemv('n', n,nb, -one,  V2,LDV2,  YA(1,1),LDYA, one, Ae1,1 );
      call dgemv('n', n,nb, -one,  YA,LDYA,  V2(1,1),LDV2, one, Ae1,1 );
      call dgemv('n', n,nb, -one,  V1,LDV1,  XA(1,1),LDXA, one, Ae1,1 );
      call dgemv('n', n,nb, -one,  XA,LDXA,  V1(1,1),LDV1, one, Ae1,1 );



      return;
      end

      subroutine Bv1v2e1(  n , nb  ,
     &                     B,LDB,
     &                     vv1, vv2,
     &                     V1,LDV1,   V2,LDV2,
     &                     XB,LDXB,   YB,LDYB,
     &                     WORK,
     &                     Bv1, Bv2, Be1
     &                   )
      implicit none;

      integer n,nb,LDB,LDV1,LDV2,LDXB,LDYB;
      double precision B(LDB,n);
      double precision V1(LDV1,nb);
      double precision V2(LDV1,nb);
      double precision XB(LDXB,nb);
      double precision YB(LDYB,nb);

      double precision vv1(n);
      double precision vv2(n);
      double precision Bv1(n);
      double precision Bv2(n);
      double precision Be1(n);

      double precision WORK(nb);

      double precision one;
      double precision zero;
      parameter( one  = 1.0d+0 );
      parameter( zero = 0.0d+0 );


      call dasymv2(n, B,LDB, vv1,vv2, Bv1,Bv2);
!      call dasymv( n , B,LDB,  vv1,  Bv1 ); !this is memory and computationally
!      call dasymv( n , B,LDB,  vv2,  Bv2 ); !expensive part

      call dcopy( n , B(1,1),1,  Be1,1   );

      if(nb.eq.0) return;

      !.......................tu radim Bv1
      call dgemv('t', n,nb, +one,  YB,LDYB,   vv1,1, zero, WORK,1 );
      call dgemv('n', n,nb, +one,  V2,LDV2,  WORK,1,  one,  Bv1,1 );
      call dgemv('t', n,nb, +one,  V2,LDV2,   vv1,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  YB,LDYB,  WORK,1,  one,  Bv1,1 );
      call dgemv('t', n,nb, +one,  XB,LDXB,   vv1,1, zero, WORK,1 );
      call dgemv('n', n,nb, +one,  V1,LDV1,  WORK,1,  one,  Bv1,1 );
      call dgemv('t', n,nb, +one,  V1,LDV1,   vv1,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  XB,LDXB,  WORK,1,  one,  Bv1,1 );


      !.......................tu radim Bv2
      call dgemv('t', n,nb, +one,  YB,LDYB,   vv2,1, zero, WORK,1 );
      call dgemv('n', n,nb, +one,  V2,LDV2,  WORK,1,  one,  Bv2,1 );
      call dgemv('t', n,nb, +one,  V2,LDV2,   vv2,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  YB,LDYB,  WORK,1,  one,  Bv2,1 );
      call dgemv('t', n,nb, +one,  XB,LDXB,   vv2,1, zero, WORK,1 );
      call dgemv('n', n,nb, +one,  V1,LDV1,  WORK,1,  one,  Bv2,1 );
      call dgemv('t', n,nb, +one,  V1,LDV1,   vv2,1, zero, WORK,1 );
      call dgemv('n', n,nb, -one,  XB,LDXB,  WORK,1,  one,  Bv2,1 );


      !.......................tu radim Be1
      call dgemv('n', n,nb, +one,  V2,LDV2,  YB(1,1),LDYB, one, Be1,1 );
      call dgemv('n', n,nb, -one,  YB,LDYB,  V2(1,1),LDV2, one, Be1,1 );
      call dgemv('n', n,nb, +one,  V1,LDV1,  XB(1,1),LDXB, one, Be1,1 );
      call dgemv('n', n,nb, -one,  XB,LDXB,  V1(1,1),LDV1, one, Be1,1 );


      if( abs( Be1(1) ) .gt. 1.D-12 ) then
          write(6,*) Be1(1);
          stop 'Be1 bi trebo imat prvi element nula...';
      endif



      return;
      end

      subroutine GenerateHGH( n, a,b, tau1,tau2, c,s, d,e )
      implicit none;

      integer n;
      double precision a(n+1);
      double precision b(n+1);
      double precision tau1;
      double precision tau2;
      double precision c;
      double precision s;
      double precision d;
      double precision e;

      double precision alpha;
      double precision beta;

      double precision one;
      parameter( one = 1.0d+0 );

      double precision ddot;



      call dlarfg( n , b(2),  b(3),1,  tau1 );
      beta = b(2);
      b(2) = one;

      alpha = - tau1 * ddot( n , b(2),1, a(2),1 );
      call daxpy( n , alpha , b(2),1,  a(2),1 );

      call drotg( a(2), beta, c , s );
      call dlarfg( n , a(2),  a(3),1,  tau2 );
      e    = a(2);
      a(2) = one;


      d = a(1);


      return;
      end

      subroutine UpdateXY( n,
     &                     Av1, Av2, Ae1,
     &                     Bv1, Bv2, Be1,
     &                     v1,  v2,
     &                     tau1, tau2,
     &                     c, s,
     &                     xA, yA,
     &                     xB, yB,
     &                     a,  b )

      implicit none;

      integer n;
      double precision Av1(n);
      double precision Av2(n);
      double precision Ae1(n);
      double precision Bv1(n);
      double precision Bv2(n);
      double precision Be1(n);
      double precision v1(n);
      double precision v2(n);
      double precision xA(n);
      double precision yA(n);
      double precision xB(n);
      double precision yB(n);
      double precision a(n);
      double precision b(n);

      double precision tau1;
      double precision tau2;
      double precision c;
      double precision s;

      integer i;
      double precision alpha;
      double precision half;
      parameter( half = 0.5d+0 );
      double precision one;
      parameter( one = 1.0d+0 );
      double precision ddot;

      double precision g(n);
      double precision h(n);
      double precision wA(n);
      double precision wB(n);
      double precision v1v2;
      double precision sp1;
      double precision sp2;


      !racunam xA
      call dscal( n, tau1, Av1,1 );
      alpha = -half * tau1 * ddot( n , Av1,1,  v1,1 );
      do i = 1 , n
          xA(i) = Av1(i) + alpha*v1(i);
      enddo

      !racunam xB
      do i = 1 , n
          xB(i) = tau1*Bv1(i);
      enddo


      !racunam g,h,wA,wB
      do i = 1 , n
          g(i) = Ae1(i) - xA(i) - xA(1)*v1(i);
          h(i) = Be1(i) - xB(i) + xB(1)*v1(i);
      enddo
      do i = 1 , n  !tu mozda moze drot al ajde...
          wA(i) = (c-one)*g(i) - s*h(i);
          wB(i) = (c-one)*h(i) + s*g(i);
      enddo
      wA(1) = wA(1) - (c-one)*g(1);


      !racunam yA
      v1v2 = ddot(n, v1,1, v2,1);
      sp1  = ddot(n, xA,1, v2,1 );
      sp2  = ddot(n, wA,1, v2,1 );
      do i = 1 , n
          yA(i) = tau2 * ( Av2(i) - v1v2*xA(i) - sp1*v1(i) + wA(i) );
      enddo
      yA(1) = yA(1) + tau2*sp2;
      alpha = - half * tau2 * ddot( n , yA,1, v2,1 );
      call daxpy( n, alpha, v2,1,  yA,1 );


      !racunam yB
      sp1 = ddot(n, xB,1,  v2,1 );
      sp2 = ddot(n, wB,1,  v2,1 );
      do i = 1 , n
          yB(i) = tau2 * ( Bv2(i) - v1v2*xB(i) + sp1*v1(i) + wB(i) );
      enddo
      yB(1) = yB(1) - tau2*sp2;



      !racunam a i b
      do i = 1 , n
          a(i) = g(i) - yA(i) - yA(1)*v2(i) + wA(i);
          b(i) = h(i) - yB(i) + yB(1)*v2(i) + wB(i);
      enddo
      a(1) = a(1) + wA(1); !wA(1)=0 pa ne treba...
      b(1) = b(1) - wB(1); !b(1) bi ionako trebo bit nula...



      if( abs( b(1) ) .gt. 1.D-12 ) then
           write(6,*) b(1);
          stop 'b bi trebo imat prvi element nula...';
      endif

      return;
      end

      subroutine RefreshA( n , bsize,
     &                     A,LDA,
     &                     V1,LDV1,  V2,LDV2,
     &                     XA,LDXA,  YA,LDYA
     &                    )
      implicit none;

      integer n, bsize, LDA, LDV1, LDV2, LDXA, LDYA;
      double precision A(LDA,n);
      double precision V1(LDV1,bsize);
      double precision V2(LDV2,bsize);
      double precision XA(LDXA,bsize);
      double precision YA(LDYA,bsize);

      double precision one;
      parameter( one = 1.0d+0 );

      call dsyr2k('l','n', n,bsize, -one, XA,LDXA, V1,LDV1, one,A,LDA);
      call dsyr2k('l','n', n,bsize, -one, YA,LDYA, V2,LDV2, one,A,LDA);

      return;
      end

      subroutine RefreshB( n , bsize,
     &                     B,LDB,
     &                     V1,LDV1,  V2,LDV2,
     &                     XB,LDXB,  YB,LDYB
     &                    )
      implicit none;

      integer n, bsize, LDB, LDV1, LDV2, LDXB, LDYB;
      double precision B(LDB,n);
      double precision V1(LDV1,bsize);
      double precision V2(LDV2,bsize);
      double precision XB(LDXB,bsize);
      double precision YB(LDYB,bsize);

      double precision one;
      parameter( one = 1.0d+0 );

      call dasyr2k( n,bsize,  XB,LDXB, V1,LDV1,  B,LDB );
      call dasyr2k( n,bsize,  YB,LDYB, V2,LDV2,  B,LDB );

      return;
      end





!mozda bi se dalo umjesto dasyr2k C = C-A*B'+B*A' napraviti C = C - A1*B1' + B1*A1' - A2*B2' + B2*A2
!jer mi zapravo to treba pa mozda malo ubrza stvar...iako vjerojatno jako malo ako uopce
#define BSIZE 128
#define TRANSPBLOCK 16
!mora biti BSIZE multipla od TRANSPBLOCK!!!!!!!!
!ako je n jako velik, tipa 4000, onda mozda bolje 256 BSIZE

      subroutine dasyr2k( n,k,  A,LDA,  B,LDB,  C,LDC )
      implicit none;
      ! C = C-A*B'+B*A', C donji trokut samo kao input
      integer n,k,LDA,LDB,LDC;
      double precision A(LDA,*);
      double precision B(LDB,*);
      double precision C(LDC,*);

      double precision one;
      parameter( one = 1.0d+0 );
      double precision zero;
      parameter( zero = 0.0d+0 );

      integer bsize;
      parameter(bsize=BSIZE);
      double precision tmp(bsize,bsize);

      integer i,j,jb,k1,k2;


      do jb = 1 , n/bsize
          j = (jb-1)*bsize+1;

          i = j;
          call dgemm('n','t',bsize,bsize,k,one,A(i,1),LDA,B(j,1),LDB,
     &                zero,tmp,bsize);
          call add(C(i,j),LDC,tmp);

          i = j + bsize;
          call dgemm('n','t',n-i+1,bsize,k,-one,A(i,1),LDA,B(j,1),LDB,
     &                one,C(i,j),LDC);
          call dgemm('n','t',n-i+1,bsize,k,+one,B(i,1),LDB,A(j,1),LDA,
     &                one,C(i,j),LDC);

      enddo

      j = (n/bsize)*bsize+1;
      if( j.le.n ) then

          i = j;
          call dgemm('n','t',n-i+1,n-j+1,k,one,A(i,1),LDA,B(j,1),LDB,
     &                zero,tmp,bsize);
          do k2 = 0 , (n-j+1)-1
              do k1 = 0 , (n-i+1)-1
                  C(i+k1,j+k2) = C(i+k1,j+k2) - tmp(k1+1,k2+1)
     &                                        + tmp(k2+1,k1+1);
              enddo
          enddo

      endif


      end
      subroutine add( C,LDC, tmp )
      implicit none;
      integer LDC,bl,i,j;
      integer bsize;
      parameter(bsize=TRANSPBLOCK);
      integer n;
      parameter(n=BSIZE);
      double precision C(LDC,n);
      double precision tmp(n,n);

      do j = 1 , n , bsize

         do i = j , j+bsize-1
             do bl = 0 , i-j
                 C(i,j+bl) = C(i,j+bl) - tmp(i,j+bl) + tmp(j+bl,i);
             enddo
         enddo

          do i = j+bsize , n
              do bl = 0 , bsize-1
                  C(i,j+bl) = C(i,j+bl) - tmp(i,j+bl) + tmp(j+bl,i);
              enddo
          enddo

      enddo

      end







#define BSIZE2 4
      !racuna y1=A*x1, y2=A*x2, gdje je A antisimetricna, samo donji trokut spemljen
      subroutine dasymv2(n, A,LDA, x1,x2, y1,y2)
      implicit none;

      integer n,LDA;
      double precision A(LDA,n);
      double precision x1(n);
      double precision x2(n);
      double precision y1(n);
      double precision y2(n);

      double precision one;
      parameter( one = 1.0d+0 );

      integer bsize;
      parameter( bsize = BSIZE2 );

      integer i,j,ib,jb,k1,k2;
      integer nblocks;


      y1 = 0.D0;
      y2 = 0.D0;


      nblocks = n/bsize;
      do jb = 1 , nblocks
          j = (jb-1)*bsize+1;

          i = j;
          do k2 = 0 , bsize-1
              do k1 = 0 , k2-1
                  A(i+k1,j+k2) = - A(j+k2,i+k1);
              enddo
          enddo
          call dgemv('n', bsize,bsize,
     &               one,  A(i,j),LDA,  x1(j),1,  one,  y1(i),1  );
          call dgemv('n', bsize,bsize,
     &               one,  A(i,j),LDA,  x2(j),1,  one,  y2(i),1  );

          i = jb*bsize + 1;
          if( n-i+1 .gt. 0 ) then
              call dgemv('n', n-i+1,bsize,
     &                   one,  A(i,j),LDA,  x1(j),1,  one,  y1(i),1  );
              call dgemv('n', n-i+1,bsize,
     &                   one,  A(i,j),LDA,  x2(j),1,  one,  y2(i),1  );
              call dgemv('t', n-i+1,bsize,
     &                   -one,  A(i,j),LDA,  x1(i),1,  one,  y1(j),1  );
              call dgemv('t', n-i+1,bsize,
     &                   -one,  A(i,j),LDA,  x2(i),1,  one,  y2(j),1  );
          endif
      enddo


      j = nblocks*bsize+1;
      if( j.le.n ) then
          i = j;
          do k2 = 0 , (n-j+1)-1
              do k1 = 0 , k2-1
                  A(i+k1,j+k2) = - A(j+k2,i+k1);
              enddo
          enddo
          call dgemv('n', n-i+1,n-j+1,
     &                one,  A(i,j),LDA,  x1(j),1,   one,  y1(i),1 );
          call dgemv('n', n-i+1,n-j+1,
     &                one,  A(i,j),LDA,  x2(j),1,   one,  y2(i),1 );
      endif




      end













