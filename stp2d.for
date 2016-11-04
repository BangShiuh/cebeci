C     MAIN PROGRAM
      COMPLEX ALFA,OMEGA,DELA
      COMMON/BLC0/ NP,NX,NXT,NX0,IX,IXT,RL,C,
     +             UINF,ETA(201),DETA(201),A(201)
      COMMON/BLCC/ X(201),S(201),UE(201),RX(201)
      COMMON/BLC5/ ALFA,OMEGA,DELA,DELR,REY(201)
      COMMON/BLC8/ QSUM(201),RLAST(201),LAST(201)
      COMMON/BLCP/ U(201),V(201),UUDP(201),UUB(201),UUDPB(201)
      COMMON/BLCD/ IEND
      COMMON/SAVE2/WS0(20), NXX(20),SX(201,20),SRX(201,20),
     +  SUE(201,20),SRALFA(201,20),SIALFA(201,20),SCLOGA(201,20)
	CHARACTER*80 input_name, output_name, summary_name      
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  read input data
c
 	write(6,*) "Enter input file name (include extension name)"
      read(5, *) input_name
	DO I = 1, 80
	   if ( input_name(I:I) .EQ. " " ) goto 32
	END DO
   32 Iend = I-1
C      print*,Iend    
	DO I = Iend,1,-1
	   if ( input_name(I:I) .EQ. "." ) go to 33
	ENDDO
	if (I .LE. 1 ) I = Iend
   33 IDOT = I
C     print*, "IDOT=", IDOT
      if (IDOT .GE. 4 )  IDOT = IDOT-3
      DO I = 1, IDOT
         OUTPUT_NAME(I:I) = input_name(I:I)
	   summary_name(I:I) = input_name(I:I)
	ENDDO
	if ( IDOT .EQ. Iend ) THEN
         OUTPUT_NAME(IDOT+1:IDOT+8)   = "_out.txt"
	   summary_name(IDOT+1:IDOT+14) = "SummaryOut.txt"
	ELSE
         OUTPUT_NAME(IDOT:IDOT+7)   = "_out.txt"
	   summary_name(IDOT:IDOT+13) = "SummaryOut.txt"
	ENDIF
	print*, "Your input file is ", input_name
      print*,"This program stp2d will generate several output files:"
	print*, "1. output file name ",   OUTPUT_NAME
	print*, "2. Summary output file name",  summary_name 

	OPEN(unit=55,file=input_name,STATUS="OLD")
	OPEN(unit=66,file=output_name) 
	open(unit=77, file=summary_name)              

      REWIND 1
	READ(55,*)
      READ(55,*)NXT,NX0,IXT,ICONT 
	READ(55,*)
      READ(55,8101) UINF,C,RL
	READ(55,*)
      READ(55,8101) ALFA,OMEGA
      WRITE(66,9000) NXT,NX0,IXT
      WRITE(66,9100) UINF,RL,C
	READ(55,*)
      READ (55,8101) (X(I),I=1,NXT)
	READ(55,*)
      READ (55,8101) (S(I),I=1,NXT )
	READ(55,*)
      READ (55,8101) (UE(I),I=1,NXT)
C
      DO I = 1,20
         NXX(I) = 0
         do j = 1,nxt
            SCLOGA(j,I) = 0.0 
         enddo
      ENDDO
C
C     CALCULATION OF STABILITY REYNOLDS NUMBER
C
      DO 400 I=1,NXT
      QSUM(I) = 0.0
C     SEE EQ. (8.5.5)
      RX(I) = UE(I)*S(I)*RL
C     SEE. EQ. (8.5.2)
      REY(I) = SQRT(RX(I))
  400 CONTINUE
      WRITE(66,9200) (I,X(I),UE(I),RX(I),REY(I),I=1,NXT)
      IX      = 0
      NX      = 1
      IEND    = 0
C
C     CALCULATION OF VELOCITY PROFILE
C
 100  CALL VELPRO(INDEX)
      IF(INDEX.EQ.1) THEN
C
C  END OF INPUT VELOCITY PROFILES, CAL. STOPS
C
         GOTO 1000
      ENDIF
C
      IF ( NX .LT. NX0 ) GO TO 500
c
      IX      = IX + 1
      LAST(IX) = 0
      IF ( IEND .GT. 0 ) LAST(IX) = 1
      IF (IXT.NE.0.AND.IX.GT.IXT) GO TO 350
c
c  added by kcc
c  get initial eigenvalue using continuation method
c
      if(nx .eq. nx0 .and. icont.eq.1) call getraw
C
C     CALCULATION OF NEUTRAL STABILITY CURVE
C
      IF ( IEND .EQ. 0 ) CALL NEWTON(0) 
      IF ( IXT.EQ.0.OR.IX .EQ. 1 ) GO TO 500
C
C     CALCULATION OF AMPLIFICATION FACTOR
C
  350 CALL NEWTONI(INDEX)
      IF (INDEX.EQ.1) GOTO 1000
  500 NX      = NX + 1
      IF (  NX.GT.NXT ) GO TO 1000
      GO TO 100
 8000 FORMAT(20I3)
 8101 FORMAT(7F10.0)
 9000 FORMAT(1H0,5X,10HNXT     = ,I10,2X,10HNX0     = ,I10,2X,
     +     10HIXT     = ,I10 )
 9100 FORMAT(1H0,5X,10HUINF    = ,F10.3,2X,10HRL      = ,F15.1,
     +       2X,10HC       = ,F10.5,///)
 9200 FORMAT(1H0,6X,2HNX,8X,1HX,14X,2HUE,13X,2HRX,12X,3HREY /
     +       (1H ,5X,I3,4E15.6) )
 9300 FORMAT(///,20X,'SUMMARY OF STABILITY RESULTS'/
     +           20X,'----------------------------'/)
 9400 FORMAT(/,5X,'INSTABILITY LINE NO = ',I3,2X,
     +       'PHYSICAL FREQUENCY = ',F10.3,' KH',
     +       /,5X,'----------------------------',
     +       '--------------------------------- ',//
     +       12X,'X',10X,'UE',10X,' RX ',9X,'PHYSICAL ALFA ',
     +       9X,'CLOGA',/)
 9500 FORMAT(6X,F10.4,2X,F10.4,2X,F12.2,2X,F10.5,2X,F10.5,2X,F10.5)
c
 1000 IF(IXT.GT.0) THEN
C
C     PRINT OUT SUMMARY OF STABILITY RESULTS
C
c        WRITE(66,9300)
	   WRITE(77,9300)
         TWOPI = 8.0 * ATAN(1.)
         sxmin = 9999.0
         sxmax = -9999.0
         DO IX = 1,IXT
c           WRITE(66,9400) IX, WS0(IX)/TWOPI/1000.
            WRITE(77,9400) IX, WS0(IX)/TWOPI/1000.
            DO I = 1,NXX(IX)
c              WRITE(66,9500) SX(I,IX),SUE(I,IX),SRX(I,IX),SRALFA(I,IX),
c    +             SIALFA(I,IX),SCLOGA(I,IX)
               WRITE(77,9500) SX(I,IX),SUE(I,IX),SRX(I,IX),SRALFA(I,IX),
     +             SIALFA(I,IX),SCLOGA(I,IX)
            ENDDO
         ENDDO
c        write(77,1777) (WS0(IX)/TWOPI/1000., IX=1, IXT)  
c1777    FORMAT(//' Frequency ', 8(2X, E13.6))
c        write(77,1778)
c1778 format('   X             N1             N2             N3      ', 
c    +       '      ......')
c       DO IX = 1, IXT
c          WRITE (77,111) SX(1,IX),(SCLOGA(I,IX),I=1,NXX(IX))
c  111       FORMAT(F10.5, 8(2X, F13.6))
c 	   ENDDO

      ENDIF
C
      close(66)
	close(55)
	close(77)
      CLOSE(6)
	CLOSE(5)
	PRINT*," "
	PRINT*,"Calculations are successfully completed."
C      PRINT*,"The output is saved in ", OUTPUT_NAME
	PRINT*," "
 	PRINT*,"Hit any key to close this DOS-window."
 	READ(5,*)
	 	 
      STOP
      END

      SUBROUTINE VELPRO(INDEX)
      COMMON/BLC0/ NP,NX,NXT,NX0,IX,IXT,RL,C,
     +             UINF,ETA(201),DETA(201),A(201)
      COMMON/BLCC/ X(201),S(201),UE(201),RX(201)
      COMMON/BLCP/ U(201),V(201),UUDP(201),UUB(201),UUDPB(201)
      COMMON/BLCD/ IEND
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      INDEX = 0
      READ (1,50,END=600) KX,NP
      IF ( KX .EQ. NX ) GO TO 150
      WRITE ( 66,3000 ) NX,KX,NP
      STOP 2
  150 READ (1,60) ( ETA(J),J=1,NP )
      READ (1,60) ( U(J),J=1,NP )
      READ (1,60) ( V(J),J=1,NP )
      DO 400 J=2,NP
      DETA(J-1)=ETA(J)-ETA(J-1)
      A(J)=0.5*DETA(J-1)
 400  CONTINUE
C
C     COMPUTE UUDP(J)
C
      NP1     = NP-1
      DETAS   = ETA(2)-ETA(1)
      DV      = V(2)-V(1)
      ANG2    = ATAN2(DV,DETAS)
      DETAI   = DETAS
      DVI     = DV
      DETA2   = DETAS
      DO 10 I = 2,NP1
      ANG1    = ANG2
      DETA1   = DETA2
      I1      = I+1
      DETAS   = ETA(I1)-ETA(I)
      DV      = V(I1)-V(I)
      ANG2    = ATAN2(DV,DETAS)
      DETA2   = DETAS
      ANG     = (DETA2*ANG1+DETA1*ANG2)/(DETA1+DETA2)
      UUDP(I) = TAN(ANG)
   10 CONTINUE
      DVF     = DV
      DETAF   = DETA2
C
C     EXTRAPOLATE FOR END VALUES
C
      UUDP(1) = 2.*DVI/DETAI - UUDP(2)
      UUDP(NP)= 2.*DVF/DETAF - UUDP(NP1)
      DO 180 J=2,NP
      UUB(J)   = 0.5*(U(J)+U(J-1))
      UUDPB(J) = 0.5*(UUDP(J)+UUDP(J-1))
 180  CONTINUE
 50   FORMAT ( 2I5 )
 60   FORMAT ( 6E14.6 )
 3000 FORMAT(1H0,2X,'INCONSISTENT DATA.  NX = ',I3,2X,'KX = ',I3,2X,
     +       'NP = ',I4 / 3X,'PROGRAM STOPS.' )
      RETURN
 600  INDEX = 1
      END

      SUBROUTINE CSAVE
      COMPLEX  ALFA,OMEGA,SQM1,ALSQ,REYI,GAMMA,GAMMAS,AIREY,
     1         SQUIG,SSIV,R,C1,C2,C3,C4,BC1,BC2,BC3,BC4,
     2         SS,FS,WS,GS,WSAR,SSAR,WSOM,SSOM,WSR,SSR,WB,SB,
     3         DELS,DELF,DELW,DELG,AA,AJ,AA1,AA2,AA3,AA4,DET,
     4         GG,W1,W2,W3,W4,BB1,BB2,BB3,BB4,BB5,BB6,DETT,CC1,CC2,DD1,

     5         DD2,DD3,DD4,DD5,DD6,DENO,UM,UMA,UMO,UMR,DELA
      COMMON/BLC0/ NP,NX,NXT,NX0,IX,IXT,RL,C,
     +             UINF,ETA(201),DETA(201),A(201)
      COMMON/BLCC/ X(201),S(201),UE(201),RX(201)
      COMMON/BLC5/ ALFA,OMEGA,DELA,DELR,REY(201)
      COMMON/BLC3/ C1(201),C2(201),C3(201),C4(201)
      COMMON/BLC4/ UM(4,1),UMA(4,1),UMO(4,1),UMR(4,1)
      COMMON/BLCP/ U(201),V(201),UUDP(201),UUB(201),UUDPB(201)
      DIMENSION AA(4,4,201),W1(201),W2(201),W3(201),W4(201),R(4,201),
     1          GG(4,4,201),DELW(201),DELS(201),DELF(201),DELG(201),
     2          SS(201),FS(201),WS(201),GS(201),WSAR(201),
     3          SSAR(201),WSOM(201),SSOM(201),WSR(201),SSR(201)
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      SQM1     = (0.0,1.0)
      GAMMAS   = ALFA**2
      GAMMA    = CSQRT(GAMMAS)
      REYI     = SQM1*REY(NX)
      SQUIG    = CSQRT(GAMMAS+REYI*(ALFA*U(NP)-OMEGA))
      II       = 0
      II       = II + 1
      DO 150 J = 2,NP
C     SEE EQ. (8.3.7)
      C3(J)    = 0.5*DETA(J-1)
      C1(J)    = C3(J)*GAMMAS
      C2(J)    = -C3(J)*REYI*ALFA*UUDPB(J)
      C4(J)    = C3(J)*(GAMMAS+REYI*(ALFA*UUB(J)-OMEGA))
C
C     STANDARD PROBLEM
C
      DO 150 K = 1,4
  150 R(K,J)   = (0.0,0.0)
      R(2,1)   = (1.0,0.0)
      R(1,1)   = (0.0,0.0)
      R(3,1)   = (0.0,0.0)
      R(4,1)   = (0.0,0.0)
C
C     EDGE BOUNDARY CONDITIONS
C
      ALSQ     = GAMMA+SQUIG
C     SEE EQ. (8.3.9)
      BC1      = GAMMA
      BC2      = (0.0,0.0)
      BC3      = 1.0/ALSQ
      BC4      = SQUIG
      GO TO 2000
C
C     VARIATIONAL EQUATIONS
C
  200 II       = II + 1
C
C     WRT  ALFA
C
      AIREY    = REYI*ALFA
      DO 250 J = 2,NP
      WB       = 0.5*(WS(J)+WS(J-1))
      SB       = 0.5*(SS(J)+SS(J-1))
C     SEE EQ. (8.3.19A)
      R(3,J-1) = -2.0*WB*DETA(J-1)*ALFA
C     SEE EQ. (8.3.19B)
      R(4,J-1) = -SB*DETA(J-1)*(2.0*ALFA+REYI*UUB(J))+
     1           WB*DETA(J-1)*REYI*UUDPB(J)
  250 CONTINUE
C     SEE EQ. (8.3.19C)
      R(3,NP)  = SS(NP)*(ALFA/GAMMA+ALFA/SQUIG+0.5*REYI*U(NP)/
     1           SQUIG)/(GAMMA+SQUIG)**2-WS(NP)*ALFA/GAMMA
C     SEE EQ. (8.3.19D)
      R(4,NP)  = -SS(NP)*(2.0*ALFA+REYI*U(NP))/SQUIG*0.5
      R(2,1)   = (0.0,0.0)
      GO TO 2000
  400 II       = II + 1
C
C     WRT OMEGA
C
      DO 425 J = 2,NP
      SB       = 0.5*(SS(J)+SS(J-1))
C     SEE EQ. (8.3.20A)
      R(3,J-1) = (0.0,0.0)
C     SEE EQ. (8.3.20B)
      R(4,J-1) = SB*DETA(J-1)*REYI
  425 CONTINUE
C     SEE EQ. (8.3.20D)
      R(4,NP)  = SS(NP)*REYI/SQUIG*0.5
C     SEE EQ. (8.3.20C)
      R(3,NP)  = -R(4,NP)/(GAMMA+SQUIG)**2
      R(2,1)   = (0.0,0.0)
      GO TO 2000
  500 II       = II + 1
C
C     WRT REY
C
      DO 525 J = 2,NP
      SB       = 0.5*(SS(J)+SS(J-1))
      WB       = 0.5*(WS(J)+WS(J-1))
C     SEE EQ. (8.3.24A)
      R(3,J-1) = (0.0,0.0)
C     SEE EQ. (8.3.24B)
      R(4,J-1) = -SQM1*DETA(J-1)*(SB*(ALFA*UUB(J)-OMEGA)
     1           -WB*ALFA*UUDPB(J))
  525 CONTINUE
      SSIV     = ALFA*U(NP)-OMEGA
C     SEE EQ. (8.3.24D)
      R(4,NP)  = -SS(NP)*SQM1*SSIV*0.5/SQUIG
C     SEE EQ. (8.3.24C)
      R(3,NP)  = -R(4,NP)/(GAMMA+SQUIG)**2
      R(2,1)   = (0.0,0.0)
C
C     BLOCK ELIMINATION ALGORITHM
C
 2000 DO 2050 I = 1,4
      DO 2050 K = 1,4
      AA(I,K,1) = (0.0,0.0)
 2050 CONTINUE
C
C     WALL BOUNDARY CONDITIONS
C
      AA(1,1,1) = (1.0,0.0)
      AA(2,2,1) = (1.0,0.0)
      W1(1)    = R(1,1)
      W2(1)    = R(2,1)
      W3(1)    = R(3,1)
      W4(1)    = R(4,1)
C
C     FORWARD SWEEP
C
      DO 2100 J = 2,NP
      AJ       = C3(J)
      AA1      = AA(1,1,J-1)-C1(J  )*AA(1,3,J-1)-C2(J  )*AA(1,4,J-1)
      AA2      = AA(2,1,J-1)-C1(J  )*AA(2,3,J-1)-C2(J  )*AA(2,4,J-1)
      AA3      = AA(1,2,J-1)-C3(J  )*AA(1,3,J-1)-C4(J  )*AA(1,4,J-1)
      AA4      = AA(2,2,J-1)-C3(J  )*AA(2,3,J-1)-C4(J  )*AA(2,4,J-1)
      DET      = AA1*AA4-AA2*AA3
      GG(1,1,J) = (AA4*(C1(J  )*AJ-1.0)-AA2*C3(J  )*AJ)/DET
      GG(1,2,J) = (AA1*C3(J  )*AJ-AA3*(C1(J  )*AJ-1.0))/DET
      GG(1,3,J) = -AJ-AA(1,3,J-1)*GG(1,1,J)-AA(2,3,J-1)*GG(1,2,J)
      GG(1,4,J) = -AA(1,4,J-1)*GG(1,1,J)-AA(2,4,J-1)*GG(1,2,J)
      GG(2,1,J) = (C2(J  )*AJ*AA4-AA2*(C4(J  )*AJ-1.0))/DET
      GG(2,2,J) = (AA1*(C4(J  )*AJ-1.0)-AA3*C2(J  )*AJ)/DET
      GG(2,3,J) = -AA(1,3,J-1)*GG(2,1,J)-AA(2,3,J-1)*GG(2,2,J)
      GG(2,4,J) = -AJ-AA(1,4,J-1)*GG(2,1,J)-AA(2,4,J-1)*GG(2,2,J)
      AA(1,1,J) = (1.0,0.0)-C1(J  )*GG(1,3,J)-C2(J  )*GG(1,4,J)
      AA(1,2,J) = -GG(1,3,J)*C3(J  )-GG(1,4,J)*C4(J  )
      AA(1,3,J) = -AJ+GG(1,3,J)
      AA(1,4,J) = GG(1,4,J)
      AA(2,1,J) = -C1(J  )*GG(2,3,J)-C2(J  )*GG(2,4,J)
      AA(2,2,J) = (1.0,0.0)-GG(2,3,J)*C3(J  )-GG(2,4,J)*C4(J  )
      AA(2,3,J) = GG(2,3,J)
      AA(2,4,J) = -AJ+GG(2,4,J)
      W1(J)    = R(1,J)-GG(1,1,J)*W1(J-1)-GG(1,2,J)*W2(J-1)
     1           -GG(1,3,J)*W3(J-1)-GG(1,4,J)*W4(J-1)
      W2(J)    = R(2,J)-GG(2,1,J)*W1(J-1)-GG(2,2,J)*W2(J-1)
     1           -GG(2,3,J)*W3(J-1)-GG(2,4,J)*W4(J-1)
      W3(J)    = R(3,J)
      W4(J)    = R(4,J)
 2100 CONTINUE
C
C     BACKWARD SWEEP
C
      J = NP
      BB1      = AA(1,1,J)-AA(1,3,J)*BC1  -AA(1,4,J)*BC2
      BB2      = AA(1,2,J)-AA(1,3,J)*BC3  -AA(1,4,J)*BC4
      BB3      = AA(2,1,J)-AA(2,3,J)*BC1  -AA(2,4,J)*BC2
      BB4      = AA(2,2,J)-AA(2,3,J)*BC3  -AA(2,4,J)*BC4
      BB5      = W1(J)-AA(1,3,J)*W3(J)-AA(1,4,J)*W4(J)
      BB6      = W2(J)-AA(2,3,J)*W3(J)-AA(2,4,J)*W4(J)
      DETT     = BB1*BB4-BB2*BB3
      DELW(J)  = (BB4*BB5-BB2*BB6)/DETT
      DELS(J)  = (BB1*BB6-BB3*BB5)/DETT
      DELF(J)  = W3(J)-BC1  *DELW(J)-BC3  *DELS(J)
      DELG(J)  = W4(J)-BC2  *DELW(J)-BC4  *DELS(J)
 2150 J        = J-1
      CC1      = W3(J)-C1(J+1)*DELW(J+1)-C3(J+1)*DELS(J+1)+DELF(J+1)
      CC2      = W4(J)-C2(J+1)*DELW(J+1)-C4(J+1)*DELS(J+1)+DELG(J+1)
      DD1      = AA(1,1,J)-AA(1,3,J)*C1(J+1)-AA(1,4,J)*C2(J+1)
      DD2      = AA(1,2,J)-AA(1,3,J)*C3(J+1)-AA(1,4,J)*C4(J+1)
      DD3      = AA(2,1,J)-AA(2,3,J)*C1(J+1)-AA(2,4,J)*C2(J+1)
      DD4      = AA(2,2,J)-AA(2,3,J)*C3(J+1)-AA(2,4,J)*C4(J+1)
      DD5      = W1(J)-AA(1,3,J)*CC1-AA(1,4,J)*CC2
      DD6      = W2(J)-AA(2,3,J)*CC1-AA(2,4,J)*CC2
      DENO     = DD1*DD4-DD2*DD3
      DELW(J)  = (DD5*DD4-DD2*DD6)/DENO
      DELS(J)  = (DD1*DD6-DD3*DD5)/DENO
      DELF(J)  = CC1-C1(J+1)*DELW(J)-C3(J+1)*DELS(J)
      DELG(J)  = CC2-C2(J+1)*DELW(J)-C4(J+1)*DELS(J)
      IF (J .GT. 1) GO TO 2150
      GO TO ( 2500, 2600, 2800,2900 ), II
C
C     SOLUTION OF THE STANDARD PROBLEM
C
 2500 DO 2510 J = 1,NP
      WS(J)    = DELW(J)
      SS(J)    = DELS(J)
      FS(J)    = DELF(J)
      GS(J)    = DELG(J)
 2510 CONTINUE
      UM(3,1)  = DELF(1)
      GO TO 200
C
C     SOLUTION OF THE VARIATIONAL EQUATION -- WRT ALFA
C
 2600 DO 2610 J = 1,NP
      WSAR(J)  = DELW(J)
      SSAR(J)  = DELS(J)
 2610 CONTINUE
      UMA(3,1) = DELF(1)
      GO TO 400
C
C     SOLUTION OF THE VARIATIONAL EQUATION --WRT OMEGA
C
 2800 DO 2810 J = 1,NP
      WSOM(J)  = DELW(J)
      SSOM(J)  = DELS(J)
 2810 CONTINUE
      UMO(3,1) = DELF(1)
      GO TO 500
C
C     SOLUTION OF THE VARIATIONAL EQUATION --WRT REY
C
 2900 DO 2910 J=1,NP
      WSR(J)=DELW(J)
      SSR(J)=DELS(J)
 2910 CONTINUE
      UMR(3,1)=DELF(1)
      RETURN
      END


      SUBROUTINE NEWTON(icont) 
      COMPLEX  ALSAVE,OMSAVE,FSSAVE,FASAVE,FOSAVE,FRSAVE,ALFA1,OMEGA1,
     1         UM,UMA,UMO,UMR,ALFA,OMEGA,FSTD,FALFA,FOMGA,DELA
      COMMON/SAVE1/ ALFA1,OMEGA1
      COMMON/SAVE2/WS0(20), NXX(20),SX(201,20),SRX(201,20),
     +  SUE(201,20),SRALFA(201,20),SIALFA(201,20),SCLOGA(201,20)
      COMMON/BLC0/ NP,NX,NXT,NX0,IX,IXT,RL,C,
     +             UINF,ETA(201),DETA(201),A(201)
      COMMON/BLCC/ X(201),S(201),UE(201),RX(201)
      COMMON/BLCP/ U(201),V(201),UUDP(201),UUB(201),UUDPB(201)
      COMMON/BLC5/ ALFA,OMEGA,DELA,DELR,REY(201)
      COMMON/BLC8/ QSUM(201),RLAST(201),LAST(201)
      COMMON/SAVE/ ALSAVE(2,20),OMSAVE(2,20),FSSAVE(2,20),
     1             FASAVE(2,20),FOSAVE(2,20),FRSAVE(2,20),
     2             RSAVE(2,20),UESAVE(2,20),QSAVE(2,20)
      COMMON/BLC4/ UM(4,1),UMA(4,1),UMO(4,1),UMR(4,1)
      COMMON/BLCD/ IEND
      DATA     EPS / 1.0E-04 /
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      if(icont .eq. 0)
     1   WRITE ( 66, 10 ) NX, X(NX), S(NX), REY(NX), NP
      ISOLV = 0
      ITIME    = 0
      DALPHA   = 0.
      DOMEGA   = 0.
      WRITE ( 66, 20 )
      IF ( IX .EQ. 1 ) GO TO 90
      ALFA     = ALFA1
      OMEGA    = OMEGA1
   90 AR       = REAL ( ALFA )
      OR       = REAL ( OMEGA )
      AI       = AIMAG ( ALFA )
      OI       = AIMAG ( OMEGA )
      WRITE ( 66, 30 ) ITIME,AR,OR,DALPHA,DOMEGA
      CALL CSAVE
      ITIME    = ITIME + 1
      FSTD     = UM (3,1)
      FALFA    = UMA(3,1)
      FOMGA    = UMO(3,1)
C     SEE EQ. (8.3.17C)
      DELTA0   = REAL(FALFA)*AIMAG(FOMGA) - REAL(FOMGA)*AIMAG(FALFA)
      IF ( ABS(DELTA0) .LE. 1.0E-06 ) GO TO 500
      D1       = 1.0 / DELTA0
C     SEE EQ. (8.3.17A)
      DALPHA   = D1 *(AIMAG(FSTD)*REAL(FOMGA)-REAL(FSTD)*AIMAG(FOMGA) )
C     SEE EQ. (8.3.17B)
      DOMEGA   = D1 *(REAL(FSTD)*AIMAG(FALFA)-AIMAG(FSTD)*REAL(FALFA) )
      GO TO 114
  500 WRITE (66,100) DELTA0
      IEND   = 1
      RETURN
C     SEE EQ. (8.3.15A)
  114 AR1      = AR + DALPHA
C     SEE EQ. (8.3.15B)
      OM1      = OR + DOMEGA
      ALFA     = CMPLX ( AR1,AI )
      OMEGA    = CMPLX ( OM1,OI )
      IF ( (ABS(DALPHA) .LE. EPS) .AND. (ABS(DOMEGA) .LE. EPS) )
     +      GO TO 115
      IF ( ITIME .LE. 10 ) GO TO 90
      WRITE ( 66, 50 )
      STOP 5
C
C     SAVE DATA FOR DEFAULT VALUE
C
  115 WRITE ( 66, 30 ) ITIME,AR1,OM1,DALPHA,DOMEGA
      if(icont .eq. 1) return
C     SEE EQ. (8.4.1)
      WS0(IX)  = REAL(OMEGA)*UE(NX)*UINF*REY(NX)/S(NX)/C
      IF ( IX .GT. 20 ) GO TO 130
      ALSAVE(1,IX) = ALFA
      OMSAVE(1,IX) = OMEGA
      FSSAVE(1,IX) = UM(3,1)
      FASAVE(1,IX) = UMA(3,1)
      FOSAVE(1,IX) = UMO(3,1)
      FRSAVE(1,IX) = UMR(3,1)
      RSAVE (1,IX) = REY(NX)
      UESAVE(1,IX) = UE(NX)
      QSAVE (1,IX) = 0.0
C
C  SAVE DATA FOR POST-PROCESSING
C
      SL       = C * S(NX) / REY(NX)
      NXX(IX)   = NXX(IX ) + 1
      SX(NXX(IX),IX) = X(NX)
      SRX(NXX(IX),IX) = RX(NX)
      SUE(NXX(IX),IX) = UE(NX)
      SRALFA(NXX(IX),IX) = REAL(ALFA)/SL
      SIALFA(NXX(IX),IX) = AIMAG(ALFA)/SL
      SCLOGA(NXX(IX),IX) = 0.
      WRITE ( 66, 60 ) IX, WS0(IX)
  130 CONTINUE
      ALFA1    = ALFA
      OMEGA1   = OMEGA
   10 FORMAT ( 5X,5HNX = ,I3,2X, 8HX(NX) = , F12.6, 2X,8HS(NX) = ,
     +         F12.6,2X,6HREY = , F12.3, 2X, 5HNP = , I3 )
   20 FORMAT ( 5X, 2HIT, 6X, 4HALFA, 11X,
     +         5HOMEGA,10X, 5HDALFA,  10X, 5HDOMGA / )
   30 FORMAT ( 5X, I2, 6(2X, E13.6) )
   50 FORMAT ( 5X,42HITERATION EXCEEDED MAXIMUM. PROGRAM STOPS. )
   60 FORMAT ( 5X, 9HLINE NO. , I3, 2X, 21HPHYSICAL FREQUENCY = ,
     +         2E15.6 )
  100 FORMAT ( 5X, 15H+++++ DELTA0 = , E13.6, 14H TOO SMALL TO ,
     +         25HCALCULATE DALFA AND DOMGA  )
      RETURN
      END


      SUBROUTINE NEWTONI(INDEX)
      COMPLEX  ALFA,OMEGA,DELA,UM,UMA,UMO,UMR,FSTD,FALFA,FOMGA,
     1         FREY,DODR,DADR,ALSAVE,OMSAVE,FSSAVE,FASAVE,FOSAVE,
     2         FRSAVE
      COMMON/BLC0/ NP,NX,NXT,NX0,IX,IXT,RL,C,
     +             UINF,ETA(201),DETA(201),A(201)
      COMMON/BLCP/ U(201),V(201),UUDP(201),UUB(201),UUDPB(201)
      COMMON/SAVE/ ALSAVE(2,20),OMSAVE(2,20),FSSAVE(2,20),
     1             FASAVE(2,20),FOSAVE(2,20),FRSAVE(2,20),
     2             RSAVE(2,20),UESAVE(2,20),QSAVE(2,20)
      COMMON/BLCC/ X(201),S(201),UE(201),RX(201)
      COMMON/BLC4/ UM(4,1),UMA(4,1),UMO(4,1),UMR(4,1)
      COMMON/BLC5/ ALFA,OMEGA,DELA,DELR,REY(201)
      COMMON/BLC8/ QSUM(201),RLAST(201),LAST(201)
      COMMON/SAVE2/WS0(20), NXX(20),SX(201,20),SRX(201,20),
     +  SUE(201,20),SRALFA(201,20),SIALFA(201,20),SCLOGA(201,20)
      DATA     EPS / 1.0E-04 /
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      WRITE ( 66, 10 ) NX, X(NX), S(NX), REY(NX), NP
      INDEX    = 0
      IXM      = IX-1
      DO 1000 I = 1,IXM
      DELAR    = 0.
      DELAI    = 0.
      IT       = 0
      K        = IX-I
      IF ( K .GT. IXT ) GO TO 1000
      J        = I+1
      IF ( LAST(K) .GT. 0 ) GO TO 1000
      WRITE ( 66, 83 ) K, J
      AR       = REAL ( ALSAVE(1,K) )
      AI       = AIMAG( ALSAVE(1,K) )
      OR       = REAL ( OMSAVE(1,K) )
      OI       = AIMAG( OMSAVE(1,K) )
      FSTD     = FSSAVE ( 1,K )
      FALFA    = FASAVE ( 1,K )
      FOMGA    = FOSAVE ( 1,K )
      FREY     = FRSAVE ( 1,K )
      UEO      = UESAVE ( 1,K )
      R0       = RSAVE  ( 1,K )
      DR       = REY(NX) - R0
      OMEGA    = OMSAVE(1,K)*(REY(NX)/R0)*(UEO**2/UE(NX)**2)
      CC       = OR*UEO**2 / R0
      DODR     = (OMEGA-OMSAVE(1,K))/DR
      DADR     = -(1.0/FALFA) * ( FOMGA*DODR+FREY )
      DELAR    = REAL( DADR ) * DR
      DELAI    = AIMAG( DADR ) * DR
      DELA     = CMPLX ( DELAR, DELAI )
      AR1      = AR + DELAR
      AI1      = AI + DELAI
      WRITE (66,84)
  520 ALFA     = CMPLX( AR1,AI1 )
      WRITE (66,85) IT,ALFA,OMEGA,DELA
      CALL CSAVE
      FSTD     = UM(3,1)
      FALFA    = UMA(3,1)
      FOMGA    = UMO(3,1)
      FREY     = UMR(3,1)
      DELA     = -FSTD / FALFA
      IF ((ABS(REAL(DELA)) .LE. EPS) .AND. (ABS(AIMAG(DELA)) .LE. EPS))
     +     GO TO 700
      AR1      = REAL( ALFA ) + REAL( DELA )
      AI1      = AIMAG( ALFA ) + AIMAG( DELA )
      IT       = IT + 1
      IF ( IT .LE. 10 ) GO TO 520
      WRITE ( 66, 70 )
      INDEX    = 1
      RETURN
CC    STOP 7
  700 IT       = IT + 1
      ALFA = ALFA +DELA
      WRITE (66,85) IT,ALFA,OMEGA,DELA
      ALSAVE(2,K) = ALFA
      OMSAVE(2,K) = OMEGA
      FSSAVE(2,K) = UM(3,1)
      FASAVE(2,K) = UMA(3,1)
      FOSAVE(2,K) = UMO(3,1)
      FRSAVE(2,K) = UMR(3,1)
      RSAVE (2,K) = REY(NX)
      UESAVE(2,K) = UE(NX)
C     SEE EQ. (8.5.6)
      QSAVE (2,K) = ( REY(NX) * AIMAG( ALFA ) ) / S(NX)
      QSUM(K) = QSUM(K)+((QSAVE(2,K)+QSAVE(1,K))/2.0)*(S(NX)-S(NX-1))
      CLOGA   = -QSUM(K)
      WRITE ( 66,50 ) QSAVE(2,K),CLOGA
C
C  SAVE RESULTS FOR POST-PROCESS
C
      SL       = C * S(NX) / REY(NX)
      NXX(K)   = NXX(K) + 1
      SX(NXX(K),K) = X(NX)
      SRX(NXX(K),K) = RX(NX)
      SUE(NXX(K),K) = UE(NX)
      SRALFA(NXX(K),K) = REAL(ALFA)/SL
      SIALFA(NXX(K),K) = AIMAG(ALFA)/SL
      SCLOGA(NXX(K),K) = CLOGA
      PROD     = QSAVE(2,K)*QSAVE(1,K)
      IF ( PROD .GE. 0.0 ) GO TO 960
      LAST(K)  = 1
      WRITE (66,88) K, J
      GO TO 1000
  960 IF ( ABS(QSAVE(2,K)).GT.1.0E-05 ) GO TO 1000
      LAST(K) = 1
      WRITE (66,87) K, J
 1000 CONTINUE
      DO 1400 I=1,IXM
      K      = IX-I
      IF ( K .GT. IXT ) GO TO 1400
      ALSAVE(1,K) = ALSAVE(2,K)
      OMSAVE(1,K) = OMSAVE(2,K)
      FSSAVE(1,K) = FSSAVE(2,K)
      FASAVE(1,K) = FASAVE(2,K)
      FOSAVE(1,K) = FOSAVE(2,K)
      FRSAVE(1,K) = FRSAVE(2,K)
      RSAVE (1,K) = RSAVE (2,K)
      UESAVE(1,K) = UESAVE(2,K)
      QSAVE (1,K) = QSAVE (2,K)
 1400 CONTINUE
   10 FORMAT ( 5X,5HNX = ,I3,2X, 8HX(NX) = , F12.6, 2X,8HS(NX) = ,
     +         F12.6,2X,6HREY = , F12.3, 2X, 5HNP = , I3 )
   50 FORMAT ( 5X, 8HQSAVE = ,F15.6,8HCLOGA = , F15.6 )
   70 FORMAT ( 5X,39HITERATION EXCEEDS 10 ( FOR AR AND AI ) , 2X,
     +         13HPROGRAM STOPS  )
   83 FORMAT ( 5X, 9HLINE NO. , I3, 2X, 10HPOINT NO. , I3 /
     +         5X, 27(1H-) / )
   84 FORMAT ( 5X, 2HIT, 6X, 5HALFAR, 10X, 5HALFAI, 10X, 5HOMGAR,
     +         10X, 5HOMGAI, 10X, 5HDELAR, 10X, 5HDELAI / )
   85 FORMAT ( 5X, I2, 6( 2X, E13.6 ) )
   87 FORMAT ( 5X, 43H***** ALFA(IMAG) IS LE 1.0E-05 AT LINE NO. ,
     +         I2, 16H  AND POINT NO. , I2 / )
   88 FORMAT ( 5X, 36H***** ALFAI IS POSITIVE AT LINE NO. , I2,
     +         16H  AND POINT NO. , I2 / )
      RETURN
      END

      subroutine getraw
c
c  this subroutine is added by kcc for getting initial guess of
c  alfa and omega
c
      COMPLEX ALFA,OMEGA,DELA
      COMMON/BLC0/ NP,NX,NXT,NX0,IX,IXT,RL,C,
     +             UINF,ETA(201),DETA(201),A(201)
      COMMON/BLCC/ X(201),S(201),UE(201),RX(201)
      COMMON/BLC5/ ALFA,OMEGA,DELA,DELR,REY(201)
      COMMON/BLC8/ QSUM(201),RLAST(201),LAST(201)
      COMMON/BLCP/ U(201),V(201),UUDP(201),UUB(201),UUDPB(201)
c
      dimension etarf(41),urf(41),wrf(41),rn(11)
      dimension u0(201),uudp0(201),urfint(201),wrfint(201)   
      data afa0,omg0,rey0,nmax/0.0854,0.0255,774.597,11/
      data etarf/0.000000,  0.050000,  0.103000,  0.159180,  0.218731,
     1           0.281855,  0.348766,  0.419692,  0.494873,  0.574566,
     1           0.659040,  0.748582,  0.843497,  0.944107,  1.050753,
     1           1.163798,  1.283626,  1.410643,  1.545282,  1.687999,
     1           1.839279,  1.999635,  2.169613,  2.349790,  2.540777,
     1           2.743224,  2.957817,  3.185286,  3.426403,  3.681987,
     1           3.952906,  4.240080,  4.544485,  4.867154,  5.209183,
     1           5.571734,  5.956038,  6.363400,  6.795203,  7.252915,
     1           7.738090/
      data   urf/0.000000,  0.016580,  0.034155,  0.052783,  0.072526,
     1           0.093449,  0.115617,  0.139098,  0.163961,  0.190273,
     1           0.218100,  0.247504,  0.278536,  0.311240,  0.345639,
     1           0.381734,  0.419496,  0.458853,  0.499682,  0.541793,
     1           0.584919,  0.628704,  0.672689,  0.716311,  0.758907,
     1           0.799735,  0.838006,  0.872943,  0.903852,  0.930207,
     1           0.951727,  0.968432,  0.980658,  0.989014,  0.994290,
     1           0.997325,  0.998893,  0.999607,  0.999886,  0.999977,
     1           1.000000/
      data   wrf/0.000037, -0.000105, -0.000333, -0.000742, -0.001366,
     1          -0.002240, -0.003404, -0.004903, -0.006787, -0.009110,
     1          -0.011930, -0.015310, -0.019311, -0.023992, -0.029406,
     1          -0.035591, -0.042560, -0.050292, -0.058713, -0.067677,
     1          -0.076950, -0.086187, -0.094923, -0.102571, -0.108447,
     1          -0.111825, -0.112029, -0.108559, -0.101238, -0.090334,
     1          -0.076623, -0.061330, -0.045937, -0.031894, -0.020304,
     1          -0.011700, -0.006010, -0.002700, -0.001036, -0.000329,
     1           0.000058/
      data nmax/11/
      data rn/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/
c
c  ------------------------------------------------------------------
c
c  save profiles at nx = nx0
c
      do j = 1,np
         u0(j) = u(j)
         uudp0(j) = uudp(j)
      enddo
      reyo = rey(nx)
c
c  interpolate urf, & wrf into eta grids
c
      call lntp(41,etarf,urf,np,eta,urfint)
      call lntp(41,etarf,wrf,np,eta,wrfint)
      urfint(np) = 1.0
c
      alfar = afa0
      omegar= omg0
      alfai = 0.0
      omegai= 0.0
      alfa   = cmplx ( alfar,alfai )
      omega  = cmplx ( omegar,omegai )
      iter  = 0
c
c  start continuation method now
c
      write(66,11)
11    format(//25x,' Results for continuation method '/
     1       25x,' ------------------------------- '/)
100   iter  = iter + 1
      cn    = rn(iter)
      do j=1,np
         u(j)=urfint(j)+cn*(u0(j)-urfint(j))
         uudp(j)=wrfint(j)+cn*(uudp0(j)-wrfint(j))
      enddo
      DO  J=2,NP
         UUB(J)   = 0.5*(U(J)+U(J-1))
         UUDPB(J) = 0.5*(UUDP(J)+UUDP(J-1))
      enddo
      rey(nx) = rey0  + cn * ( reyo - rey0  )
c
      write (66,1 ) nx, iter,alfar,omegar,rey(nx)
c     write (66,* )
c     write (66,2 ) (j,eta(j),u(j),uudp(j),j=1,np)
2     format(i5,3e13.5)
c
      call newton(1)
      alfar  = real ( alfa )
      omgar  = real ( omega )
      if(iter .lt. nmax) goto 100
      write(66,12)
12    format(25x,' End of continuation method calculation',/
     1       25x,' --------------------------------------'//)
1     format('   nx,iter,alfar,omgar,rey(nx)=',2i5,3e13.5)
      return
      end
c
      subroutine lntp(nn,xn,yn,no,xo,yo)
      dimension xn(nn),yn(nn),xo(no),yo(no)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(nn.le.1) return
      do i = 1,no
         ii = i
         if(xo(i) .ge. xn(1)) goto 10
         yo(i) = yn(1)
      enddo
10    continue
      do i = no,1,-1
         in = i
         if(xo(i) .lt. xn(nn) ) goto 20
         yo(i) = yn(nn)
      enddo
20    continue
      js    = 2
      do 100 i=ii,in
         do 50  j=js,nn
            jj      = j
            if(xn(j).ge.xo(i)) goto 75
   50    continue
   75    yo(i)   = yn(jj-1) + (yn(jj)-yn(jj-1))/(xn(jj)-xn(jj-1))*
     1             (xo(i)-xn(jj-1))
         js      = jj
  100 continue
      return
      end

