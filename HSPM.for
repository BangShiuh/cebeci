C     MAIN (With Wake)
      COMMON /BOD/ NODTOT,X(201),Y(201),          
     +             XMID(200),YMID(200),COSTHE(200),SINTHE(200)
      COMMON /NUM/ PI,PI2INV, ICYCLE                      
      COMMON/COUPL/DLSP(201),VNP(201),QW(51),DELW(51,2)                  
      COMMON /OFFBOD/ XOFF(201),YOFF(201),SOFF(200),COFF(200)   
      COMMON /OBKU/ STE(2),CTE(2),VNTE(2),VTTE(2),VTOT(2)  
      COMMON /COMP / FMACH,C2,C6,C7,C8                 
      COMMON /WWAKE/ NW,XW(51),YW(51),SINW(51),COSW(51),QWM(51) 
      DIMENSION TITLE(20),DBODY(200),SBGEO(201),CBGEO(201)
C	CHARACTER*80 input_name, output_name  
	                    
      PI      = 3.1415926585                                            
      PI2INV  = .5/PI 
     
C	   write(6,*) "Enter input file name (include extension name)"
C      read(5, *) input_name
C      write(6,*) "Enter output file name"
C      read(5, *) output_name 
      
	ICYCLE = 1

C 233  FORMAT(I1)
C 234  FORMAT(I2)
      print*, "Your input file is ",  input_name
	  print*, "Your output file is ",  output_name

      open(unit=5, file="naca.txt", status="OLD")
      open(unit=6, file="hspm.txt")

      
	  READ (5,*)                                            
      READ (5,*) NODTOT, NW                                           
 501  FORMAT(2I5)
      READ (5,*)                                                       
      READ (5,*)(X(I),I=1,NODTOT+1)
	  READ (5,*)                                   
      READ (5,*)(Y(I),I=1,NODTOT+1)
	
	DO I =1, NODTOT+1
	   DLSP(I) = 0.0
	   VNP(I)  = 0.0
	END DO
	DO I =1, NW+1
	   XW(I)  = 0.0
	   YW(I)  = 0.0
         DELW(I,1) = 0.0
	   DELW(I,2) = 0.0
         QW(I) = 0.0
	END DO
	
 504  FORMAT(6F10.0)
      READ (5,*)                                                    
      READ (5,504) ALPHA                                                 
      WRITE (6,1030) ALPHA                                              
 1030 FORMAT (//,'SOLUTION AT ALPHA = ',F10.5,/)           
      COSALF  = COS(ALPHA*PI/180.)                                      
      SINALF  = SIN(ALPHA*PI/180.)
	  READ (5,*)                                      
      READ (5,*) FMACH
      GAM    = 1.4                                                       
      C2     = 0.5 * (GAM-1.0)       
      C6     = C2 * FMACH**2   
      C7     = GAM/(GAM-1.0) 
      C8     = 0.5 * GAM * FMACH**2   
      DO 100 I=1,NODTOT  
C     XMI AND YMI, SEE EQ. (5.3.12)
      XMID(I) = .5*(X(I) + X(I+1)) 
      YMID(I) = .5*(Y(I) + Y(I+1)) 
      DX=X(I+1)-X(I)                                                     
      DY=Y(I+1)-Y(I)                                                     
      DIST=SQRT(DX*DX+DY*DY)         
C     SEE EQ. (5.3.2)
      SINTHE(I)=DY/DIST           
      COSTHE(I)=DX/DIST     
      DBODY(I)=DIST       
100   CONTINUE
      if (ICYCLE .ne. 1 ) then
        DO 150 I=1,NW                                                    
        DX=XW(I+1)-XW(I)                                                   
        DY=YW(I+1)-YW(I)                                                   
        DIST=SQRT(DX*DX+DY*DY)                                             
        SINW(I)=DY/DIST                                                    
        COSW(I)=DX/DIST                                                    
        QWM(I)=0.5*(QW(I)+QW(I+1))                                        
150     CONTINUE
      endif                                                          
      DO 200 I=2,NODTOT                                                 00006590
      SBGEO(I)=(SINTHE(I)*DBODY(I-1)+SINTHE(I-1)*DBODY(I))/(DBODY(I-1)+         
     +          DBODY(I))                                               00006610
200   CBGEO(I)=(COSTHE(I)*DBODY(I-1)+COSTHE(I-1)*DBODY(I))/(DBODY(I-1)+
     +          DBODY(I))                                               00006630
      SBGEO(1)=2.0*SINTHE(1)-SBGEO(2)                                   00006640
      SBGEO(NODTOT+1)=2.0*SINTHE(NODTOT)-SBGEO(NODTOT)                  00006650
      CBGEO(1)=2.0*COSTHE(1)-CBGEO(2)                                   00006660
      CBGEO(NODTOT+1)=2.0*COSTHE(NODTOT)-CBGEO(NODTOT)                  00006670
      DO 220 J=1,NODTOT+1                                               00003120
      XOFF(J)=X(J)-SBGEO(J)*DLSP(J)                                     00003130
      YOFF(J)=Y(J)+CBGEO(J)*DLSP(J)                                     00003140
220   CONTINUE
      DO 210 I=1,NODTOT                                                 00003150
      DX=XOFF(I+1)-XOFF(I)                                              00003160
      DY=YOFF(I+1)-YOFF(I)                                              00003170
      DIST=SQRT(DX*DX+DY*DY)                                            00003180
C     SEE EQ. (5.3.2)
      SOFF(I)=DY/DIST                                                   00003190
      COFF(I)=DX/DIST                                                   00003200
210   CONTINUE
      STE(1)=SOFF(1)                                                    00003230
      STE(2)=SOFF(NODTOT)                                               00003240
      CTE(1)=COFF(1)                                                    00003250
      CTE(2)=COFF(NODTOT)                                               00003260
      IT=1                                                              00003270
 240  CONTINUE                                                          00003280
      CALL COEF(SINALF,COSALF)                                          00003300
      CALL GAUSS(1)                                                     00003310
      CALL OBKUTA(SINALF,COSALF,INDEX)
	call WAKEDS(SINALF,COSALF)                                  
      IF(INDEX.EQ.1) GO TO 230                                          00003340
      IF(IT.GT.10) STOP 10                                              00003350
      IT=IT+1                                                           00003370
      DO 250 II=1,2                                                     00003460
      STEO=STE(II)                                                      00003470
      CTEO=CTE(II)                                                      00003480
      STE(II)=(STEO*VTTE(II)+CTEO*VNTE(II))/VTOT(II)                    00003500
      CTE(II)=(CTEO*VTTE(II)-STEO*VNTE(II))/VTOT(II)                    00003510
250   CONTINUE
      STE(1)=-STE(1)
      CTE(1)=-CTE(1)
      GO TO 240                                                         00003520
230   CONTINUE                                                          00003530
      CALL VPDIS(SINALF,COSALF)                                         00003550
      CALL VPDWK(SINALF,COSALF)                                         00010910
      CALL CLCM(SINALF,COSALF) 
	write(6,*)
	write(6,*) input_name
C	if (ICYCLE .lt. 10 ) then
C	write(6,*) input_name(1:10)
C	else
C      write(6,*) input_name(1:11)
C	endif
    CLOSE(6)
	CLOSE(5)
	PRINT*," "
	PRINT*,"Calculations are successfully completed."
C      PRINT*,"The output is saved in ", OUTPUT_NAME
	PRINT*," "
 	PRINT*,"Hit any key to close this DOS-window."
 	READ(5,*) 	 
	                                               
      STOP                                                              00003570
      END                                                               00003580
      SUBROUTINE COEF(SINALF,COSALF)                                    00007190
      COMMON /BOD/ NODTOT,X(201),Y(201),                  
     +             XMID(200),YMID(200),COSTHE(200),SINTHE(200)                  
      COMMON /COF/ A(201,201),BV(201),KUTTA                                    
      COMMON /NUM/ PI,PI2INV, ICYCLE                                       
      COMMON/COUPL/DLSP(201),VNP(201),QW(51),DELW(51,2)                  
      COMMON /OFFBOD/ XOFF(201),YOFF(201),SOFF(200),COFF(200)            
      COMMON /OBKU/ STE(2),CTE(2),VNTE(2),VTTE(2),VTOT(2)                
      COMMON /WWAKE/ NW,XW(51),YW(51),SINW(51),COSW(51),QWM(51)         00007230
      KUTTA   = NODTOT + 1                                              00007360
      DO 120  I = 1,NODTOT                                              00007380
      A(I,KUTTA) = 0.0                                                  00007410
      DO 110  J = 1,NODTOT                                              00007460
      FLOG    = 0.0                                                     00007470
      FTAN    = PI                                                      00007480
      IF (J .EQ. I)     GO TO 100                                       00007490
      DXJ     = XMID(I) - X(J)                                          00007500
      DXJP    = XMID(I) - X(J+1)                                        00007510
      DYJ     = YMID(I) - Y(J)                                          00007520
      DYJP    = YMID(I) - Y(J+1)                                        00007530
C     FLOG IS LN(R(I,J+1)/R(I,J)), SEE EQ. (5.3.12)
      FLOG    = .5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))        00007540
C     FTAN IS BETA(I,J), SEE EQ. (5.3.12)
      FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)              00007550
C     CTIMTJ IS COS(THETA(I)-THETA(J))
 100  CTIMTJ  = COSTHE(I)*COSTHE(J) + SINTHE(I)*SINTHE(J)               00007560
C     STIMTJ IS SIN(THETA(I)-THETA(J))
      STIMTJ  = SINTHE(I)*COSTHE(J) - COSTHE(I)*SINTHE(J)               00007570
C     ELEMENTS OF THE COEFFICIENT MATRIX, A(I,J), SEE EQ. (5.4.1A)
      A(I,J)  = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)                      00007580
      B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)                      00007590
C     ELEMENTS OF THE COEFFICIENT MATRIX, A(I,N+1), SEE EQ. (5.4.1B)
      A(I,KUTTA) = A(I,KUTTA) + B                                       00007600
 110  CONTINUE                                                          00007610
      WAQ=0.0
	if (ICYCLE .ne. 1 ) then                  
        DO 111 J=1,NW                                                   
        DXJ     = XMID(I) - XW(J)                                        
        DXJP    = XMID(I) - XW(J+1)                                     
        DYJ     = YMID(I) - YW(J)                                      
        DYJP    = YMID(I) - YW(J+1)                                   
        FLOG    =.5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))     
        FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)           
        CTIMTJ  = COSTHE(I)*COSW(J) + SINTHE(I)*SINW(J)                  
        STIMTJ  = SINTHE(I)*COSW(J) - COSTHE(I)*SINW(J)               
        AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)                   
 111    WAQ     = WAQ+AA*QWM(J)
      endif                                           
C     ELEMENTS OF VECTOR B FOR I=1,...,N, SEE EQ. (5.5.1)
      BV(I) = SINTHE(I)*COSALF - COSTHE(I)*SINALF                       00007780
     +        +0.5*(VNP(I)+VNP(I+1))-WAQ                                00007790
 120  CONTINUE                                                          00007810
200   DO 90   J = 1,KUTTA                                               00007850
 90   A(KUTTA,J)   = 0.0                                                00007860
      WAQ=0.0
      I=1                                                               00007880
      II=1                                                              00007890
230   CONTINUE                                                          00007900
      XMIDO   = .5*(XOFF(I) + XOFF(I+1))                                00007910
      YMIDO   = .5*(YOFF(I) + YOFF(I+1))                                00007920
      DLSH=0.5*(DLSP(I)+DLSP(I+1))                                      00001020
      DO 210  J = 1,NODTOT                                              00007960
      FLOG=0.0                                                          00007970
      FTAN=PI                                                           00007980
      IF(J.EQ.I .AND. DLSH.LT.0.0001)GO TO 209          
      DXJ     = XMIDO - X(J)                                            00008000
      DXJP    = XMIDO - X(J+1)                                          00008010
      DYJ     = YMIDO - Y(J)                                            00008020
      DYJP    = YMIDO - Y(J+1)                                          00008030
C     FLOG IS LN(R(I,J+1)/R(I,J)), SEE EQ. (5.3.12)
      FLOG    = .5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))        00008040
C     FTAN IS BETA(I,J), SEE EQ. (5.3.12)
      FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)              00008050
C     CTIMTJ IS COS(THETA(I)-THETA(J))
 209  CTIMTJ  = CTE(II)*COSTHE(J) + STE(II)*SINTHE(J)                   00008060
C     STIMTJ IS SIN(THETA(I)-THETA(J))
      STIMTJ  = STE(II)*COSTHE(J) - CTE(II)*SINTHE(J)                   00008070
C     SEE EQ. (5.5.10B)
      AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)                      00008080
C     SEE EQ. (5.5.10A)
      B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)                      00008090
C     ELEMENTS OF THE COEFFICIENT MATRIX, A(N+1,J), SEE EQ. (5.5.11B)
      A(KUTTA,J) = A(KUTTA,J) - B                                       00008140
C     ELEMENTS OF THE COEFFICIENT MATRIX, A(N+1,N+1), SEE EQ. (5.5.11B)
      A(KUTTA,KUTTA) = A(KUTTA,KUTTA) +AA                               00008150
210   CONTINUE
      if (ICYCLE .ne. 1 ) then      
        DO 211 J=1,NW                                                   
        DXJ     = XMIDO - XW(J)                                         
        DXJP    = XMIDO - XW(J+1)                                       
        DYJ     = YMIDO - YW(J)                                        
        DYJP    = YMIDO - YW(J+1)                                      
        FLOG    = .5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))     
        FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)           
        CTIMTJ  = CTE(II)*COSW(J) + STE(II)*SINW(J)                   
        STIMTJ  = STE(II)*COSW(J) - CTE(II)*SINW(J)                   
        B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)                 
        WAQ     = WAQ-B*QWM(J)                                      
211     CONTINUE
      endif                                                          
      IF(I.EQ.NODTOT) GO TO 220                                         00008290
      I=NODTOT                                                          00008300
      II=2                                                              00008310
      GO TO 230                                                         00008320
220   CONTINUE                                                          00008330
C     ELEMENTS OF VECTOR B FOR I=N+1, SEE EQ. (5.5.11B)        
      BV(KUTTA) = - (CTE(1) + CTE(2))*COSALF                            00008350
     +            - (STE(1) + STE(2))*SINALF - WAQ                      00008360
      RETURN                                                            00008380
      END                                                               00008390
      SUBROUTINE CLCM(SINALF,COSALF)                                    00017920
      COMMON /BOD/ NODTOT,X(201),Y(201),                  
     +             XMID(200),YMID(200),COSTHE(200),SINTHE(200)                  
      COMMON /CPD/ UE(250),CP(250),UEW(50,2),CPW(50,2)                  0001622
      COMMON /COMP / FMACH,C2,C6,C7,C8                                  00000310
      CFX     = 0.0                                                     
      CFY     = 0.0                                                     
      CM      = 0.0                                                     
C                                                                       00008450
C     COMPRESSIBILITY CORRECTION BASED ON KARMAN-TSIEN FORMULA
C
      BETA    = SQRT(1.-FMACH**2)                                        
      DO  I =1,NODTOT                                                  
      CP(I)   = CP(I)/(BETA+(FMACH**2/(1.+BETA))*CP(I)*0.5)              
      V2      = AMAX1(0.0,(1.0-(1.0+C8*CP(I))**(1./C7))/C6+1.0)          
      UE(I)   = SQRT(V2)                                                 
      END DO                                                           
C                                                                        
      DO 100  I = 1,NODTOT                                              
      DX      = X(I+1) - X(I)                                           
      DY      = Y(I+1) - Y(I)                                  
      CFX     = CFX + CP(I)*DY                                          
      CFY     = CFY - CP(I)*DX                                          
      CM      = CM + CP(I)*(DX*XMID(I) + DY*YMID(I))                          
 100  CONTINUE                                                          
      CL      = CFY*COSALF - CFX*SINALF                                 
      WRITE (6,1000) CL,CM                                           
 1000 FORMAT(//,'    CL =',F10.5,'    CM =',F10.5)     
      RETURN                                                            
      END                                                               
      SUBROUTINE GAUSS(M)                                        
      COMMON /COF/ A(201,201),B(201,1),N                                        
      DO 100  K = 1,N-1                                               
      KP      = K + 1                                                   
      DO 100  I = KP,N                                               
      R       = A(I,K)/A(K,K)                                        
      DO 200  J = KP,N                                                   
 200  A(I,J)  = A(I,J) - R*A(K,J)                                      
      DO 100  J = 1,M
 100  B(I,J)    = B(I,J) - R*B(K,J)
      DO 300  K = 1,M
      B(N,K) = B(N,K)/A(N,N)                            
      DO 300  I = N-1,1,-1                                          
      IP      = I + 1                                                   
      DO 400  J = IP,N                                              
 400  B(I,K)    = B(I,K) - A(I,J)*B(J,K)                                 
 300  B(I,K)    = B(I,K)/A(I,I)                                           
      RETURN                                                            
      END                                                               
      SUBROUTINE OBKUTA(SINALF,COSALF,INDEX)                            00000690
      COMMON /BOD/ NODTOT,X(201),Y(201),                  
     +             XMID(200),YMID(200),COSTHE(200),SINTHE(200)                  
      COMMON /COF/ A(201,201),BV(201),KUTTA                                    
      COMMON /NUM/ PI,PI2INV, ICYCLE                                            
      COMMON/COUPL/DLSP(201),VNP(201),QW(51),DELW(51,2)                 00002280
      COMMON /OFFBOD/ XOFF(201),YOFF(201),SOFF(200),COFF(200)           00002330
      COMMON /OBKU/ STE(2),CTE(2),VNTE(2),VTTE(2),VTOT(2)               00000790
      COMMON /WWAKE/ NW,XW(51),YW(51),SINW(51),COSW(51),QWM(51)         00007230
      DIMENSION Q(200)
      DO 50   I = 1,NODTOT                                              00000860
 50   Q(I)    = BV(I)                                                   00000870
      GAMMA   = BV(KUTTA)                                               00000880
      I=1                                                               00000920
      II=1                                                              00000930
   55 CONTINUE                                                          00000950
      XMIDO   = .5*(XOFF(I) + XOFF(I+1))                                00000960
      YMIDO   = .5*(YOFF(I) + YOFF(I+1))                                00000970
C     CONTRIBUTION TO VT(I) FROM FREESREAM VELOCITY, SEE EQ. (5.5.4B)
      VTANG   = COSALF*CTE(II) + SINALF*STE(II)                         00000980
C     CONTRIBUTION TO VN(I) FROM FREESREAM VELOCITY, SEE EQ. (5.5.4A)
      VNOFF   = SINALF*CTE(II) - COSALF*STE(II)                         00000990
      DLSH=0.5*(DLSP(I)+DLSP(I+1))                                      00001020
      DO 120  J = 1,NODTOT                                              00001120
      FLOG    = 0.0                                                     00001140
      FTAN    = PI                                                      00001150
      IF(J.EQ.I.AND.DLSH.LT.0.0001)GO TO 100  
 300  DXJ     = XMIDO - X(J)                                            00001170
      DXJP    = XMIDO - X(J+1)                                          00001180
      DYJ     = YMIDO - Y(J)                                            00001190
      DYJP    = YMIDO - Y(J+1)                                          00001200
C     FLOG IS LN(R(I,J+1)/R(I,J)), SEE EQ. (5.3.12)
      FLOG    = .5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))        00001210
C     FTAN IS BETA(I,J), SEE EQ. (5.3.12)
      FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)              00001220
C     CTIMTJ IS COS(THETA(I)-THETA(J))
 100  CTIMTJ  = CTE(II)*COSTHE(J) + STE(II)*SINTHE(J)                   00001230
C     STIMTJ IS SIN(THETA(I)-THETA(J))
      STIMTJ  = STE(II)*COSTHE(J) - CTE(II)*SINTHE(J)                   00001240
C     AA IS BT(I,J)=AN(I,J), SEE EQ. (5.5.7B)
      AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)                      00001250
C     B IS -AT(I,J)=BN(I,J), SEE EQ. (5.5.7C)
      B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)                      00001260
C     CONTRIBUTION TO VT(I) FROM SINGULARITIES, SEE EQ. (5.5.4B)
      VTANG   = VTANG - B*Q(J) + GAMMA*AA                               00001270
C     CONTRIBUTION TO VN(I) FROM SINGULARITIES, SEE EQ. (5.5.4A)
      VNOFF   = VNOFF + GAMMA*B + AA*Q(J)                               00001280
 120  CONTINUE                                                          00001290
C                                                                       00009880
      if (ICYCLE .ne. 1 ) then  
	  DO J=1,NW                                                     
        DXJ     = XMIDO - XW(J)                                           
        DXJP    = XMIDO - XW(J+1)                                         
        DYJ     = YMIDO - YW(J)                                           
        DYJP    = YMIDO - YW(J+1)                                         
        FLOG    = .5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))        
        FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)              
        CTIMTJ  = CTE(II)*COSW(J) + STE(II)*SINW(J)                       
        STIMTJ  = STE(II)*COSW(J) - CTE(II)*SINW(J)                       
        AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)                      
        B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)                      
        VTANG   = VTANG-B*QWM(J)                                          
        VNOFF   = VNOFF+AA*QWM(J)                                         
       END DO 
      endif                                                        
      VTTE(II)= VTANG                                                   00001460
      VNTE(II)= VNOFF                                                   00001450
      VTOT(II)= SQRT(VTANG*VTANG+VNOFF*VNOFF)
      IF(II.EQ.2) GO TO 130                                             00001470
      I=NODTOT                                                          00001480
      II=2                                                              00001490
      GO TO 55                                                          00001500
130   INDEX=0                                                           00001510
      IF(ABS(VTOT(1)-VTOT(2)).LT.0.00001)INDEX=1                        00001520
      RETURN                                                            00001540
      END                                                               00001550
      SUBROUTINE VPDIS(SINALF,COSALF)                                   00016120
      COMMON /BOD/ NODTOT,X(201),Y(201),                  
     +             XMID(200),YMID(200),COSTHE(200),SINTHE(200)                  
      COMMON /COF/ A(201,201),BV(201),KUTTA                                    
      COMMON /CPD/ UE(250),CP(250),UEW(50,2),CPW(50,2)                  00016220
      COMMON /NUM/ PI,PI2INV, ICYCLE         
      COMMON/COUPL/DLSP(201),VNP(201),QW(51),DELW(51,2)                 00002280
      COMMON /OFFBOD/ XOFF(201),YOFF(201),SOFF(200),COFF(200)           00002330
      COMMON /WWAKE/ NW,XW(51),YW(51),SINW(51),COSW(51),QWM(51)
	COMMON /COMP / FMACH,C2,C6,C7,C8              
      DIMENSION Q(200)
      WRITE (6,1000) 
      DO 50   I = 1,NODTOT                                              00016390
 50   Q(I)    = BV(I)                                                   00016400
      GAMMA   = BV(KUTTA)                                               00016410
      DO 130  I = 1,NODTOT                                              00016450
      XMIDO   = .5*(XOFF(I) + XOFF(I+1))                                00016480
      YMIDO   = .5*(YOFF(I) + YOFF(I+1))                                00016490
C     CONTRIBUTION TO VT(I) FROM FREESREAM VELOCITY, SEE EQ. (5.5.4B)
      VTANG   = COSALF*COFF(I) + SINALF*SOFF(I)                         00016500
C     CONTRIBUTION TO VN(I) FROM FREESREAM VELOCITY, SEE EQ. (5.5.4A)
      VNOFF   = SINALF*COFF(I) - COSALF*SOFF(I)                         00016510
      DLSH=0.5*(DLSP(I)+DLSP(I+1))                                      00016550
      DO 120  J = 1,NODTOT                                              00016650
      FLOG    = 0.0                                                     00016670
      FTAN    = PI                                                      00016680
      IF (J .EQ. I .AND. DLSH. LT. 0.0001) GO TO 100    
 300  DXJ     = XMIDO - X(J)                                            00016700
      DXJP    = XMIDO - X(J+1)                                          00016710
      DYJ     = YMIDO - Y(J)                                            00016720
      DYJP    = YMIDO - Y(J+1)                                          00016730
C     FLOG IS LN(R(I,J+1)/R(I,J)), SEE EQ. (5.3.12)
      FLOG    = .5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))        00016740
C     FTAN IS BETA(I,J), SEE EQ. (5.3.12)
      FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)              00016750
C     CTIMTJ IS COS(THETA(I)-THETA(J))
 100  CTIMTJ  = COFF(I)*COSTHE(J) + SOFF(I)*SINTHE(J)                   00016760
C     STIMTJ IS SIN(THETA(I)-THETA(J))
      STIMTJ  = SOFF(I)*COSTHE(J) - COFF(I)*SINTHE(J)                   00016770
C     AA IS BT(I,J)=AN(I,J), SEE EQ. (5.5.7B)
      AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)                      00016780
C     B IS -AT(I,J)=BN(I,J), SEE EQ. (5.5.7C)
      B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)                      00016790
C     CONTRIBUTION TO VT(I) FROM SINGULARITIES, SEE EQ. (5.5.4B)
      VTANG   = VTANG - B*Q(J) + GAMMA*AA                               00016800
C     CONTRIBUTION TO VN(I) FROM SINGULARITIES, SEE EQ. (5.5.4A)
      VNOFF   = VNOFF + GAMMA*B + AA*Q(J)                               00016810
 120  CONTINUE                                                          00016820
C                                                                       00012710
      
	DO 121 J=1,NW                                                     
        DXJ     = XMIDO - XW(J)                                        
        DXJP    = XMIDO - XW(J+1)                                       
        DYJ     = YMIDO - YW(J)                                       
        DYJP    = YMIDO - YW(J+1)                                     
        FLOG    =.5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))     
        FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)         
        CTIMTJ  = COFF(I)*COSW(J) + SOFF(I)*SINW(J)                   
        STIMTJ  = SOFF(I)*COSW(J) - COFF(I)*SINW(J)                   
        AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)                  
        B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)                 
        VTANG   = VTANG-B*QWM(J)                                     
        VNOFF   = VNOFF+AA*QWM(J)                                      
 121  CONTINUE
                                             
C
      UE(I)=SQRT(VNOFF*VNOFF+VTANG*VTANG)                               00016990
      IF(VTANG.LT.0.0) UE(I)=-UE(I)                                     00017000
      CP(I)   = 1.0 - UE(I)*UE(I)
      BETA    = SQRT(1.-FMACH**2)                                                                                       
      CPP  = CP(I)/(BETA+(FMACH**2/(1.+BETA))*CP(I)*0.5)              
      V2      = AMAX1(0.0,(1.0-(1.0+C8*CPP)**(1./C7))/C6+1.0)          
      UEE  = SQRT(V2)                                                 
      IF ( UE(I) .LT. 0 ) UEE= -UEE
	                     	                                      
      WRITE (6,1050) I,XMID(I),YMID(I),Q(I),GAMMA,CP(I),UEE                   
 130  CONTINUE                                                          00017090
 1000 FORMAT(/,4X,'J',4X,'X(J)',6X,'Y(J)',6X,'Q(J)',5X,'GAMMA',5X,              
     + 'CP(J)',6X,'V(J)',/)                                             
 1050 FORMAT(I5,6F10.5)                                                 
      RETURN                                                            
      END                                                               
      SUBROUTINE VPDWK(SINALF,COSALF)                                   00014170
      COMMON /BOD/ NODTOT,X(201),Y(201),                  
     +             XMID(200),YMID(200),COSTHE(200),SINTHE(200)                  
      COMMON /COF/ A(201,201),BV(201),KUTTA                                    
      COMMON /CPD/ UE(250),CP(250),UEW(50,2),CPW(50,2)                  
      COMMON /NUM/ PI,PI2INV, ICYCLE                                            
      COMMON/COUPL/DLSP(201),VNP(201),QW(51),DELW(51,2)                 
      COMMON /OFFBOD/ XOFF(201),YOFF(201),SOFF(200),COFF(200)           
      COMMON /WWAKE/ NW,XW(51),YW(51),SINW(51),COSW(51),QWM(51)         00007230
      COMMON /COMP / FMACH,C2,C6,C7,C8                                  00014350
      DIMENSION Q(200), XWMID(51), YWMID(51)
      DO 50   I = 1,NODTOT                                              00014390
 50   Q(I)    = BV(I)                                                   00016400
      GAMMA   = BV(KUTTA)
	WRITE(6, 1000)                                              
      DO 130  I = 1,NW                                                  00014450
         DO 125  K = 1,2                                                00014460
            XMIDO   = .5*(XW(I) + XW(I+1))                              00014470
            YMIDO   = .5*(YW(I) + YW(I+1))                              00014480
            VTANG   = COSALF*COSW(I) + SINALF*SINW(I)                   00014490
            VNOFF   = SINALF*COSW(I) - COSALF*SINW(I)                   00014500
            DLSH=0.5*(DELW(I,K)+DELW(I+1,K))                            00014510
            IF (K.EQ.1) THEN                                            00014520
               XMIDO   = XMIDO+SINW(I)*DLSH                             00014530
               YMIDO   = YMIDO-COSW(I)*DLSH                             00014540
            ELSE                                                        00014550
               XMIDO   = XMIDO-SINW(I)*DLSH                             00014560
               YMIDO   = YMIDO+COSW(I)*DLSH                             00014570
            ENDIF                                                       00014580
            DO 120  J = 1,NODTOT                                        00014620
               DXJ     = XMIDO - X(J)                                   00014630
               DXJP    = XMIDO - X(J+1)                                 00014640
               DYJ     = YMIDO - Y(J)                                   00014650
               DYJP    = YMIDO - Y(J+1)                                 00014660
               FLOG    =.5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))00014670
               FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)     00014680
               CTIMTJ  = COSW(I)*COSTHE(J) + SINW(I)*SINTHE(J)          00014690
               STIMTJ  = SINW(I)*COSTHE(J) - COSW(I)*SINTHE(J)          00014700
               AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)             00014710
               B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)             00014720
               VTANG   = VTANG-B*Q(J)+GAMMA*AA                          00014730
               VNOFF   = VNOFF + GAMMA*B + AA*Q(J)                      00014740
 120        CONTINUE                                                    00014750
C
               DO 121 J=1,NW                                            00014780
                  FLOG    = 0.0                                         00014790
                  FTAN    = PI                                          00014800
                  IF (J .EQ. I .AND. DLSH. LT. 0.0001) GO TO 100        
                  DXJ     = XMIDO - XW(J)                               00014820
                  DXJP    = XMIDO - XW(J+1)                             00014830
                  DYJ     = YMIDO - YW(J)                               00014840
                  DYJP    = YMIDO - YW(J+1)                             00014850
                  FLOG    = .5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+     00014860
     1                      DYJ*DYJ))                                   00014870
                  FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)  00014880
 100              CTIMTJ  = COSW(I)*COSW(J) + SINW(I)*SINW(J)           00014890
                  STIMTJ  = SINW(I)*COSW(J) - COSW(I)*SINW(J)           00014900
                  AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)          00014910
                  B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)          00014920
                  VTANG   = VTANG-B*QWM(J)                              00014930
                  VNOFF   = VNOFF+AA*QWM(J)                             00014940
 121           CONTINUE                                                 00014950
C
            UEW(I,K)   = SQRT(VNOFF*VNOFF+VTANG*VTANG)                  00014980
            CPW(I,K)   = 1.0 - UEW(I,K)*UEW(I,K)                        00014990
 125     CONTINUE                                                       00015000
 130  CONTINUE
                                                             
C                                                                       00015020
C  COMPRESSIBILITY CORRECTION BASED ON KARMAN-TSIEN FORMULA
C                                                                       00015040
         BETA    = SQRT(1.-FMACH**2)                                    00015275
         DO 157 K = 1,2                                                 00015280
            DO 157 I = 1,NW                                             00015290
               CPW(I,K) = CPW(I,K)/(BETA+(FMACH**2/(1.+BETA))*CPW(I,K)  00015300
     1                 *0.5)                                            00015310
               V2   = AMAX1(0.0,(1.0-(1.0+C8*CPW(I,K))**(1./C7))/C6+1.0)00008713
               UEW(I,K) = SQRT(V2)
		     
  157    CONTINUE                                                       00015330
C
      
      DO I =1, NW
	   XWMID(I) = 0.5*(XW(I)+XW(I+1))
	   YWMID(I) = 0.5*(YW(I)+YW(I+1))    
	   WRITE (6,1050) I,XWMID(I),YWMID(I),CPW(I,1),CPW(I,2),
     +    UEW(I,1), UEW(I,2)
	ENDDO

      If ( NW .GT. 0 ) Then
	write (6, *)
	write (6, *)
      write (6, *) "  I      XW        YW"
	do I = 1, NW+1
	   write(6, 1060) I,XW(I), YW(I) 
	enddo
	write(6,*)
      write(6,*)
	end if

 1050 FORMAT(I5,6F10.5)
 1060 FORMAT(I5,2F10.5)
 1000  FORMAT(/,4X,'J',3X,'XWMID(J)',2X,'YWMID(J)',2X,'CPW(J,1)'
     +  ,2X,'CPW(J,2)',2X,              
     + 'UEW(J,1)',2X,'UEW(J,2)',/)  
      RETURN                                                            00015360
      END                                                               00015370

C-------------------------------------------------------------------------------------------
      SUBROUTINE WAKEDS(SINALF,COSALF)
      PARAMETER (NB=5,NXB=201,NXBT=NB*NXB,NXW=51,NSOLV=15)
      COMMON /BOD/ NODTOT, X(201), Y(201),
     +      XMID(200),YMID(200),COSTHE(200),SINTHE(200)
      COMMON /NUM  / PI,PI2INV
      COMMON/COUPL/DLSP(201),VNP(201),QW(51),DELW(51,2)
      COMMON /WWAKE/ NW,XW(51),YW(51),SINW(51),COSW(51),
     +               QWM(51)
	COMMON /COF/ A(201,201),BV(201),KUTTA

      DIMENSION XWN(51),YWN(51)
C
      
      GAMMA=BV(NODTOT+1)
      DO 100 I=1,NW
         QWM(I)=0.5*(QW(I)+QW(I+1))
100   CONTINUE
C
      
      XWN(1)= (X(1)+X(NODTOT+1)) / 2.0
      YWN(1)= (Y(1)+Y(NODTOT+1)) / 2.0
      DSTE  = SQRT((X(NODTOT+1)-X(NODTOT))**2+(Y(NODTOT+1)-
     1        Y(NODTOT))**2)
      DL    = AMAX1(0.005,DSTE)
      THW=0.5*(ATAN2(Y(1)-Y(2),X(1)-X(2))+ATAN2(Y(NODTOT+1)-
     1         Y(NODTOT),X(NODTOT+1)-X(NODTOT)))
C
      DO 500 I=1,NW
         IT=0
 120     COSWK=COS(THW)
         SINWK=SIN(THW)
         XWN(I+1)=XWN(I)+DL*COSWK
         YWN(I+1)=YWN(I)+DL*SINWK
         XMIDO=0.5*(XWN(I)+XWN(I+1))
         YMIDO=0.5*(YWN(I)+YWN(I+1))
         VTANG=COSALF*COSWK+SINALF*SINWK
         VN=SINALF*COSWK-COSALF*SINWK
         J = 0
         DO 150 JJ=1,NODTOT
         J = J+1
            DXJ     = XMIDO - X(JJ)
            DXJP    = XMIDO - X(JJ+1)
            DYJ     = YMIDO - Y(JJ)
            DYJP    = YMIDO - Y(JJ+1)
            FLOG    = .5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))
            FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)
            CTIMTJ  = COSWK*COSTHE(J) + SINWK*SINTHE(J)
            STIMTJ  = SINWK*COSTHE(J) - COSWK*SINTHE(J)
            AA      = PI2INV*(FTAN*CTIMTJ + FLOG*STIMTJ)
            B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)
            VTANG   = VTANG-B*BV(J)+GAMMA*AA
            VN      = VN+GAMMA*B+BV(J)*AA
150      CONTINUE
C
         IF (NW.NE.0) THEN
            DO 160 J=1,NW
               FLOG=0.0
               FTAN=PI
               IF(J.EQ.I)GO TO 170
               DXJ     = XMIDO - XW(J)
               DXJP    = XMIDO - XW(J+1)
               DYJ     = YMIDO - YW(J)
               DYJP    = YMIDO - YW(J+1)
               FLOG    =.5*ALOG((DXJP*DXJP+DYJP*DYJP)/(DXJ*DXJ+DYJ*DYJ))
               FTAN    = ATAN2(DYJP*DXJ-DXJP*DYJ,DXJP*DXJ+DYJP*DYJ)
170            CTIMTJ  = COSWK*COSW(J) + SINWK*SINW(J)
               STIMTJ  = SINWK*COSW(J) - COSWK*SINW(J)
               AA      = PI2INV*(FTAN*CTIMTJ+FLOG*STIMTJ)
               B       = PI2INV*(FLOG*CTIMTJ - FTAN*STIMTJ)
               VTANG   = VTANG-B*QWM(J)
               VN      = VN+AA*QWM(J)
160         CONTINUE
         ENDIF
C
         IF(ABS(VN).GE.0.0001) THEN
            IF (IT.GT.40) GOTO 501
            IT=IT+1
            DTHE=ATAN2(VN,VTANG)
            THW=THW+DTHE
            GO TO 120
         ENDIF
  501    CONTINUE
         DL=1.15*DL
  500 CONTINUE
C
      DO 60 I=1,NW+1
         XW(I)=XWN(I)
         YW(I)=YWN(I)
   60 CONTINUE
      DO 70 I=1,NW
         DX=XW(I+1)-XW(I)
         DY=YW(I+1)-YW(I)
         DIST=SQRT(DX*DX+DY*DY)
         SINW(I)=DY/DIST
         COSW(I)=DX/DIST
   70 CONTINUE
      RETURN
 1000 FORMAT(3X,'WARNING: WAKE NOT EXACTLY A DIVIDING STREAMLINE ',
     &       'FOR ELEMENT',I3)
      END
