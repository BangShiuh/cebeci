      PROGRAM pinterface
C***********************************************************************
C  This program finds the upper surface and lower surface and generates 
C  input data for inverse boundary layer program (iblp2d) and stability 
C  (stp2d) programs
C  
C
C  Input:  the output file of the panel3 (HSPM with visicous effects) code
C  OUTPUT: thwaitesInp.txt       (input file for Thwaites method)
C          headInp.txt           (input file for Head's method)
C          iblp2dInp.txt         (input file for iblp2d program)
C          stp2dInp.txt          (input file for stp2d program)
C
      real*8  X(201), Y(201), Q(201),GAMMA(201),CP(201),UE(201),alpha
      real*8  XC(201), YC(201), UEC(201)
	real*8  XCU(201), YCU(201), UECU(201)
      real*8  XCL(201), YCL(201), UECL(201)
      real*8 S(302)
	real*8  XW(101), YW(101), UE1W(101), UE2W(101)
	CHARACTER*160 RLINE
      CHARACTER*21 junk
      character*80 PANEL
	character*80 PanInp, lowername, uppername
	Character   SCYCLE
	CHARACTER*2 SCYCLE2

	IFLOW     = 0
      method    = 0
      IGEOMETRY = 0
	ICYCLE = 1

      WRITE(6,*) "Flow parameters"
      WRITE(6,*) "ENTER Reference length BIGL = (1. for example)"
      READ(5, *) BIGL
      WRITE(6,*) "Enter reference velocity"
      READ(5, *) UREF
      WRITE(6,*) "Enter  kinematic viscosity ",
     1           "CNU = (1.6e-4 for example)"
     	READ(5,*) CNU
	WRITE(6,*) "Enter 1st x-station (NX0) for stability calculation"
	READ(5,*) NX0 
	Write(6,*) "Enter number of instability lines to be calcualted"
      READ(5,*) IXT

C********* MODIFICATION *********

C Ask for the initial estimates for alphar and omegar
	Write(6,*) "ENTER estimation for real(alpha) - Example : alphar =
     f		 0.1"
      READ(5,*) alfr
	Write(6,*) "ENTER estimation for real(omega) - Example : omegar =
     f		 0.035"
      READ(5,*) omegr

C***** End of Modification *****


C       CNU=1.6E-4
C       UREF = 400.0*1.6
C       BIGL  =1
C	 NX0 = 40
C	 IXT = 5  

C      UREF = RL*CNU/BIGL
      RL = UREF * BIGL/CNU 
      
	I = 0
C Generate input data for IBLP2D program
      ETAE = 8.0
      VGP = 1.05
      IF (VGP .LT. 1.0 ) VGP = 1.0
	IF (IFLOW .EQ. 1 ) THEN
         XCTR = 0.037
	ELSE
	   XCTR = 0.0
	ENDIF
      IF ( XCTR .LE. 0.0 ) THEN
C Laminar flow computation
         XCTR = BIGL+1
         DETA = 0.1
      ELSE
C turbulent flow computation with smaller initial step deltaEta
         DETA = 0.01
      ENDIF
	DETA = 0.01
      P2 = 1

C===============================================
C for both surfaces
	WRITE(6,*)"Enter input file name (its format must be ",
     *         "the same as the panel output file eg. panelOut.txt)"
      READ (5, *) PANEL
C      PANEL = "N0012a0.out"

      NUMUP = 0
      OPEN(55, file=PANEL,status='old')
      READ (55,*)
      READ (55,*)
      READ (55,444) junk,alpha
      READ (55,*)
      READ (55,*)
      READ (55,*)
      READ (55,*)
 444  FORMAT(A21,F10.5)

C  Read outputs along the geometry
      I = 0
 201  READ(55,'(A160)',ERR=2020,END=2020) RLINE
      IF(RLINE(1:20) .EQ. '                    ')THEN
	   IF ( J .EQ. I ) THEN
	      NODTOT = I
	      GOTO 202
	   ELSE
	       WRITE(6,*) "Data Format Errors in Input file: ",PANEL
	       GOTO 2020
	   ENDIF
	ELSE
	   I = I + 1
	   READ(RLINE,*) J, X(I),Y(I), Q(I),GAMMA(I),CP(I),UE(I)
         IF ( UE(I) .LE. 0.0 ) NUMUP = I
	   GOTO 201
	ENDIF

C Read the output along the wake    
 202  READ (55,*)
	read (55,*)
	II = 0
 203  READ(55,'(A160)',ERR=2020,END=2020) RLINE
      IF (RLINE(1:20) .EQ. '                    ')THEN
	   IF ( JJ .EQ. II ) THEN
	        NWAKE = II
	        GOTO 2020
	   ELSE 
	        WRITE(6,*) "Data Format Errors in Input file: ",PANEL
	        GOTO 2020
	   ENDIF	  
	ELSE
	   II = II + 1
	   READ(RLINE,*) JJ, XW(II),YW(II),CPW1,CPW2,UE1W(II),UE2W(II) 
	   GOTO 203
	ENDIF  
 2020 CONTINUE
        
 401	READ(55,'(A160)',ERR=490,END=490) RLINE
	if ( RLINE(1:8) .EQ. '    CL =') then
	   read(55,*) 
	   read(55,*) PanInp
	   goto 490
	endif
      goto 401
 490  close(55)
      
C read panel code input data          
      open(unit=55, file=PanInp, status="OLD")
      READ (55, *)
      READ (55,*) NODTOT_pan, NW_pan

	  
C===================================================
C  Generate the output data
C	     
      if (ABS(UE(NUMUP)) .LE. 0.0001 ) then
         XCL(1) = X(NUMUP)
	   YCL(1) = Y(NUMUP)
	   UECL(1) = 0.0
	   Do I = 2, NUMUP
	     J = NUMUP+1-I
	     XCL(I)  = X(J)
		 YCL(I)  = Y(J) 
           UECL(I) = -UE(J)
	   END DO
	   NLower = NUMUP

         XCU(1) = X(NUMUP)
	   YCU(1) = Y(NUMUP)
	   UECU(1) = 0.0
	   Do I = NUMUP+1, NODTOT 
	     J = I +1 - NUMUP 
	     XCU(J)  = X(I)
		 YCU(J)  = Y(I) 
           UECU(J) = UE(I)
	   END DO
	   NUpper = NODTOT + 1 - NUMUP

	elseif( ABS(UE(NUMUP+1)) .le. 0.0001 ) then
         XCL(1) = X(NUMUP+1)
	   YCL(1) = Y(NUMUP+1)
	   UECL(1) = 0.0
	   Do I = 2, NUMUP+1
	     J = NUMUP+2-I
	     XCL(I)  = X(J)
		 YCL(I)  = Y(J) 
           UECL(I) = -UE(J)
	   END DO
	   NLower = NUMUP+1

         XCU(1) = X(NUMUP+1)
	   YCU(1) = Y(NUMUP+1)
	   UECU(1) = 0.0
	   Do I = NUMUP+2, NODTOT 
	     J = I - NUMUP 
	     XCU(J)  = X(I)
		 YCU(J)  = Y(I) 
           UECU(J) = UE(I)
	   END DO
	   NUpper = NODTOT - NUMUP
	else
	   X_stag = X(NUMUP) - UE(NUMUP)*(X(NUMUP+1)-X(NUMUP))
     *              /(UE(NUMUP+1)-UE(NUMUP))
 	   Y_stag = Y(NUMUP) - UE(NUMUP)*(Y(NUMUP+1)-Y(NUMUP))
     *              /(UE(NUMUP+1)-UE(NUMUP))
         XCL(1) = X_stag
	   YCL(1) = Y_stag
	   UECL(1) = 0.0
	   Do I = 2, NUMUP+1
	     J = NUMUP+2-I
	     XCL(I)  = X(J)
		 YCL(I)  = Y(J) 
           UECL(I) = -UE(J)
	   END DO
	   NLower = NUMUP+1

         XCU(1) = X_stag
	   YCU(1) = Y_stag
	   UECU(1) = 0.0
	   Do I = NUMUP+1, NODTOT 
	     J = I +1 - NUMUP 
	     XCU(J)  = X(I)
		 YCU(J)  = Y(I) 
           UECU(J) = UE(I)
	   END DO
	   NUpper = NODTOT + 1 - NUMUP  	   
	endif

      KASE = 1
      KDIS = 0
	ISWPT = 20

C generate lower surface
	WRITE(6,*) " "
	WRITE(6,*) "============================================"
      WRITE(6,*) "Generate iblp2d input file for ",
     +        "lower surface of airfoil  with "
	WRITE(6,*) "NXT,   NXTE,  NXS,   NPT,  ISWPT"
      WRITE(6,*) NLower+NWAKE, NLower, 14, 81, ISWPT
	WRITE(6,*) " RL            XCTR     ETAE     VGP     DETA     P2"
	WRITE(6,8811) RL, XCTR, ETAE,VGP,DETA,P2

	WRITE(6,*) " "

	lowername(1:9) = "lower.inp"

      write(6,*) lowername
	OPEN(13,file=lowername)
	WRITE(13,*)" NXT   NXTE    NXS    NPT   ISWPT  isep  istab  icont"
	WRITE(13,8888)  NLower+NWAKE, NLower, 14, 81, ISWPT, 1, 1, 1
	WRITE(13,*) " RL            XCTR     ETAE     VGP     DETA     P2"

	WRITE(13,8811) RL, XCTR, ETAE,VGP,DETA,P2

	write(13,*) "xbl(i),i=1,nxt"
      WRITE(13,88) ( XCL(I), I=1, NLower)
	If (NWAKE .GT. 0 ) WRITE(13,88) ( XW(I), I=1,NWAKE )
	write(13,*) "ybl(i),i=1,nxt"
      WRITE(13,88) ( YCL(I), I=1, NLower)
	If (NWAKE .GT. 0 ) write(13,88) ( YW(I), I=1, NWAKE)
	write(13,*) "ubl(i),i=1,nxt"
      WRITE(13,88) (UECL(I), I=1, NLower)
	If (NWAKE .GT. 0 ) write(13,88) (UE1W(I), I=1, NWAKE)  
	WRITE(13,*) " nx0  ixt    uinf     bigl     alfar    omegar"

C********* MODIFICATION *********

C Write the initial estimates for alphar and omegar
	WRITE(13,7777)   nx0, ixt,  UREF,   BIGL,  alfr,     omegr 

C***** End of Modification *****

	WRITE(13,*)
      WRITE(13,*)
	

	write(13,*) PANEL, "Angle of attack =", alpha

      WRITE(13,*)
	WRITE(13,*)
      WRITE(13,*) "-----------Data Explanation-------------"
      WRITE(13,*) "NXT NXTE NPT NXS ISWPT ISEP ISTAB ICONT"
      WRITE(13,*) " RL,XCTR,ETAE,VGP,DETA(1),P2(1)"
      WRITE(13,*) " X(1)  ...  X(NXT)"
      WRITE(13,*) " Y(1)  ...  Y(NXT)"
      WRITE(13,*) " UE(1) ...  UE(NXT)"
      WRITE(13,*) " NX0, IXT, UINF, BIGL, ALFAR, OMEGAR"

      WRITE(13,*) " NXT    ---  The number of x-stations,including wake"
     1     ," points if there are" 
      WRITE(13,*) "             if there are no wake points, nxt = nxte"
      WRITE(13,*) " NXTe   ---  The number of x-stations on the body" 
      WRITE(13,*) " NPT   ---  The number of normal grid points"
      WRITE(13,*) " nxs   ---  x - station where cals switch to "
     1            ,"inverse method"
      WRITE(13,*) " iswpt ---  total number of inv/vis"  
      WRITE(13,*) " isep  ---  index:"
      WRITE(13,*) "            = 1, turn on transition if laminar flow "
     1            ,"occurs"
      WRITE(13,*) "            = 0, flow remains laminar even if "
     1            ,"laminar flow occurs"
      WRITE(13,*) " istab ---  index for stability calcualtion"
      WRITE(13,*) "            = 1, yes. two input files for stability "
     1            ,"calculation"   
      WRITE(13,*) "                      will be generated, but user "
     1            ,"has to input"
      WRITE(13,*) "                      nx0,ixt,uinf,bigl,alfar, omega"                     
      WRITE(13,*) "                      nx0  : nx station stability "
     1            ,"calculation starts"
      WRITE(13,*) "                      ixt  : no of lines stability "
      WRITE(13,*) "                             amplification will be "
     1            ,"calculated"
      WRITE(13,*) "                      unif : freestream velocity" 
      WRITE(13,*) "                      bigl : ref. length on which "
     1            ,"Rey is calculated"
      WRITE(13,*) "                      alfar; real value of stability"
     1            ," wave number"
      WRITE(13,*) "                      omega: real value of "
     1            ,"stability frequency" 
      WRITE(13,*) "                        both alfar & omega are "
     1            ,"normalized values"
      WRITE(13,*) "            = 0, no.  no file will be generated. "
     1            ,"user does not"
      WRITE(13,*) "                      have to input nx0,ixt,uinf,"
     1            ,"bigl,alfar, omega" 
      WRITE(13,*) " icont ---  index for using continuation to guess "
     1            ,"alfar, omegar"
      WRITE(13,*) "            = 1, yes. continuation method will be "
     1           ,          "used to determine"
      WRITE(13,*) "                 alfar, omegar. user can input "
     1            ,                  "arbitrary alfar, omegar" 
      WRITE(13,*) "            = 0, no. continuation method will not be"
     1            ,                 " used to determine"
      WRITE(13,*) "                 alfar, omegar. user has to guess "
     1            ,                  "alfar, omegar as"
      WRITE(13,*) "                 accurate as possible at nx0"
      WRITE(13,*) " X(I)  ---  The x-coordinate of geometry"
      WRITE(13,*) " Y(I)  ---  The y-coordinate of geometry"
      WRITE(13,*) " UE(I) ---  The external velocity"
      WRITE(13,*) " RL    ---  Reynolds number"
      WRITE(13,*) " XCTR  ---  Transition location"
      WRITE(13,*) " ETAE  ---  The boundary layer thickness"
      WRITE(13,*) " VGP   ---  The variable grid paramemter"
      WRITE(13,*) " DETA(1)--  The size of first normal grid"
      WRITE(13,*) " P2(1)  --  Initial pressure at first x-station"
      WRITE(13,*) "c"
      WRITE(13,*) "c the following variables are used for stability "
     1            ,"cal. only"
      WRITE(13,*) "c"
      WRITE(13,*) "NX0, IXT, UINF, BIGL, ALFAR, OMEGAR"
      WRITE(13,*)" NX0    --  1st NX station stability to be calculated"
      WRITE(13,*) " IXT    --  no of instability lines to be calcualted"
      WRITE(13,*) "            to determine the transition location"
      WRITE(13,*) " UINF   --  Refernece velocity used to normalize the"
     1            ," vel. scale" 
      WRITE(13,*) " BIGL   --  Refernece length used to normalize the "
     1            ,"length scale"
      WRITE(13,*) " ALFAR  --  Normalized instability wave number"
      WRITE(13,*) " OMEGAR --  Normalized instability wave frequency"
 
      CLOSE(13)

C generate the upper surface
	WRITE(6,*) " "
	WRITE(6,*) "============================================"
      WRITE(6,*) "Generate iblp2d input file for ",
     +       "upper surface airfoil  with "
	WRITE(6,*) "NXT,   NXTE,  NXS,   NPT,  ISWPT"
      WRITE(6,*) NUpper+NWAKE, NUpper, 14, 81, ISWPT
	WRITE(6,*) " RL            XCTR     ETAE     VGP     DETA     P2"
	WRITE(6,8811) RL, XCTR, ETAE,VGP,DETA,P2

	WRITE(6,*) " "

	
	uppername(1:9) = "upper.inp" 
	
      write(6,*) uppername 
      OPEN(13,file=uppername) 
	WRITE(13,*)" NXT   NXTE    NXS    NPT   ISWPT  isep  istab  icont"
	WRITE(13,8888)  NUpper+NWAKE, NUpper, 14, 81, ISWPT,  1,  1,   1
	WRITE(13,*) " RL            XCTR     ETAE     VGP     DETA     P2"
	WRITE(13,8811) RL, XCTR, ETAE,VGP,DETA,P2

	write(13,*) "xbl(i),i=1,nxt"
      WRITE(13,88) ( XCU(I), I=1, NUpper)
	If (NWAKE .GT. 0 ) WRITE(13,88) ( XW(I), I=1,NWAKE )
	write(13,*) "ybl(i),i=1,nxt"
      WRITE(13,88) ( YCU(I), I=1, NUpper)
	If (NWAKE .GT. 0 ) write(13,88) ( YW(I), I=1, NWAKE)
	write(13,*) "ubl(i),i=1,nxt"
      WRITE(13,88) (UECU(I), I=1, NUpper)
      If (NWAKE .GT. 0 ) write(13,88) (UE2W(I), I=1, NWAKE)
	WRITE(13,*) " nx0  ixt    uinf     bigl     alfar    omegar"

C********* MODIFICATION *********

C Write the initial estimates for alphar and omegar
	WRITE(13,7777)   nx0, ixt,  UREF,   BIGL,  alfr,     omegr 
	
C***** End of Modification *****
	
	WRITE(13,*)
      WRITE(13,*)
	

	write(13,*) PANEL, "Angle of attack =", alpha
	
      WRITE(13,*)
	WRITE(13,*)
      WRITE(13,*) "-----------Data Explanation-------------"
      WRITE(13,*) "NXT NXTE NPT NXS ISWPT ISEP ISTAB ICONT"
      WRITE(13,*) " RL,XCTR,ETAE,VGP,DETA(1),P2(1)"
      WRITE(13,*) " X(1)  ...  X(NXT)"
      WRITE(13,*) " Y(1)  ...  Y(NXT)"
      WRITE(13,*) " UE(1) ...  UE(NXT)"
      WRITE(13,*) " NX0, IXT, UINF, BIGL, ALFAR, OMEGAR"

      WRITE(13,*) " NXT    ---  The number of x-stations,including wake"
     1     ," points if there are" 
      WRITE(13,*) "             if there are no wake points, nxt = nxte"
      WRITE(13,*) " NXTe   ---  The number of x-stations on the body" 
      WRITE(13,*) " NPT   ---  The number of normal grid points"
      WRITE(13,*) " nxs   ---  x - station where cals switch to "
     1            ,"inverse method"
      WRITE(13,*) " iswpt ---  total number of inv/vis"  
      WRITE(13,*) " isep  ---  index:"
      WRITE(13,*) "            = 1, turn on transition if laminar flow "
     1            ,"occurs"
      WRITE(13,*) "            = 0, flow remains laminar even if "
     1            ,"laminar flow occurs"
      WRITE(13,*) " istab ---  index for stability calcualtion"
      WRITE(13,*) "            = 1, yes. two input files for stability "
     1            ,"calculation"   
      WRITE(13,*) "                      will be generated, but user "
     1            ,"has to input"
      WRITE(13,*) "                      nx0,ixt,uinf,bigl,alfar, omega"                     
      WRITE(13,*) "                      nx0  : nx station stability "
     1            ,"calculation starts"
      WRITE(13,*) "                      ixt  : no of lines stability "
      WRITE(13,*) "                             amplification will be "
     1            ,"calculated"
      WRITE(13,*) "                      unif : freestream velocity" 
      WRITE(13,*) "                      bigl : ref. length on which "
     1            ,"Rey is calculated"
      WRITE(13,*) "                      alfar; real value of stability"
     1            ," wave number"
      WRITE(13,*) "                      omega: real value of "
     1            ,"stability frequency" 
      WRITE(13,*) "                        both alfar & omega are "
     1            ,"normalized values"
      WRITE(13,*) "            = 0, no.  no file will be generated. "
     1            ,"user does not"
      WRITE(13,*) "                      have to input nx0,ixt,uinf,"
     1            ,"bigl,alfar, omega" 
      WRITE(13,*) " icont ---  index for using continuation to guess "
     1            ,"alfar, omegar"
      WRITE(13,*) "            = 1, yes. continuation method will be "
     1           ,          "used to determine"
      WRITE(13,*) "                 alfar, omegar. user can input "
     1            ,                  "arbitrary alfar, omegar" 
      WRITE(13,*) "            = 0, no. continuation method will not be"
     1            ,                 " used to determine"
      WRITE(13,*) "                 alfar, omegar. user has to guess "
     1            ,                  "alfar, omegar as"
      WRITE(13,*) "                 accurate as possible at nx0"
      WRITE(13,*) " X(I)  ---  The x-coordinate of geometry"
      WRITE(13,*) " Y(I)  ---  The y-coordinate of geometry"
      WRITE(13,*) " UE(I) ---  The external velocity"
      WRITE(13,*) " RL    ---  Reynolds number"
      WRITE(13,*) " XCTR  ---  Transition location"
      WRITE(13,*) " ETAE  ---  The boundary layer thickness"
      WRITE(13,*) " VGP   ---  The variable grid paramemter"
      WRITE(13,*) " DETA(1)--  The size of first normal grid"
      WRITE(13,*) " P2(1)  --  Initial pressure at first x-station"
      WRITE(13,*) "c"
      WRITE(13,*) "c the following variables are used for stability "
     1            ,"cal. only"
      WRITE(13,*) "c"
      WRITE(13,*) "NX0, IXT, UINF, BIGL, ALFAR, OMEGAR"
      WRITE(13,*)" NX0    --  1st NX station stability to be calculated"
      WRITE(13,*) " IXT    --  no of instability lines to be calcualted"
      WRITE(13,*) "            to determine the transition location"
      WRITE(13,*) " UINF   --  Refernece velocity used to normalize the"
     1            ," vel. scale" 
      WRITE(13,*) " BIGL   --  Refernece length used to normalize the "
     1            ,"length scale"
      WRITE(13,*) " ALFAR  --  Normalized instability wave number"
      WRITE(13,*) " OMEGAR --  Normalized instability wave frequency"


 
      CLOSE(13)

      
 443  FORMAT(3I3,2F10.5,F10.7)
 555  FORMAT(3F10.5)
 501  FORMAT(I5,6F10.5)
 88   FORMAT(5(F16.6, 1X))



      CLOSE(6)
	CLOSE(5)
	PRINT*," "
	PRINT*,"Calculations are successfully completed."
C       PRINT*,"The output is saved in ", OUTPUT_NAME
	PRINT*," "
 	PRINT*,"Hit any key to close this DOS-window."
 	READ(5,*) 	 
  
 8000 FORMAT(20I3)
 8101 FORMAT(7F10.5)
 8010 FORMAT(2F10.5, F10.1)
 8888 FORMAT(8(I5, 2X))
 8811 FORMAT(F12.2, 1X, 5(F8.4,1X))
 7777 FORMAT(I4, 2X, I3, 3X, 4(F8.4, 1X))
      STOP
      END
