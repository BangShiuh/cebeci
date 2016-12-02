C Panel and Inverse Boundary-Layer Method
C     program PANEL_IBL
 	parameter (nbt=201)
      parameter (nangle=25)
	Character*80 input_name,  PanelOut_name, iblOut_name,
     +    summary_name 
	real alphas(nangle), x(nbt),y(nbt)
	Integer nodtot,jangle
	Dimension Cl_Res(25), Cd_Res(25), Cm_Res(25)
	
c
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
      IF (IDOT .GT. 7) IDOT = IDOT-3
C     print*, "IDOT=", IDOT 
      DO I = 1, IDOT
	   PanelOut_name(I:I) = input_name(I:I) 
         iblOut_name(I:I) = input_name(I:I)
	   summary_name(I:I) = input_name(I:I)
	ENDDO
	if ( IDOT .EQ. Iend ) THEN
	   PanelOut_name(IDOT+1:IDOT+9) = "Panel.out"
         iblOut_name(IDOT+1:IDOT+7) = "ibl.out"
	   summary_name(IDOT+1:IDOT+11) = "Summary.out"
	ELSE
	   PanelOut_name(IDOT:IDOT+8) = "Panel.out"
         iblOut_name(IDOT:IDOT+6) = "ibl.out"
	   summary_name(IDOT:IDOT+10) = "Summary.out"
	ENDIF
	print*, "Your input file is ", input_name
      print*, "This program will generate three output files:"
	print*, "1. Panel output file name ",   PanelOut_name
      print*, "2. IBL output file name",      iblOut_name 
	print*, "3. Summary output file name",  summary_name 
	print*, "It is runing. It may take several minutes"




      open(unit=8, file=input_name,status='old')
      open(unit=7, file=PanelOut_name)
      open(unit=9, file=summary_name)
      open(unit=12,file=iblOut_name)

      rewind 8
      rewind 7
      rewind 9
      rewind 12

	
	call FGeoInp(nodtot, jangle, x, y, 
     +           fmach, rl, iterative,alphas)
     
	numPanel = nodtot
      jjangle = jangle
	Lxy = 201
	LAA = 25
c       
 	jjangle = min0(15,jangle)
      call PANEL_IBL(numPanel,jjangle, 
     +         x, y, Lxy, alphas, LAA, fmach, rl,
     +         PanelOut_name, iblOut_name, summary_name,
     +         Cl_Res, Cd_Res, Cm_Res,iterative)
     
               close(6)
C              close(5)
               close(7)
               close(8)
               close(9)
               close(12)

	PRINT*," "
	PRINT*,"Calculations are successfully completed."
C      PRINT*,"The output is saved in ", OUTPUT_NAME
	PRINT*," "
 	PRINT*,"Hit any key to close this DOS-window."
 	READ(5,*) 	 

      Stop 
	end     

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine FGeoInp(nodtot, jangle, x, y, 
     +                   fmach, rl, iterative,alphas)
      parameter (nbt=201)
      parameter (nangle=25)
      real alphas(nangle), x(nbt),y(nbt)
      Integer nodtot,jangle
 
      read (8,*)
      read (8,*) nodtot
      read (8,*)
      read (8,*) (x(i),i=1,nodtot+1)
      read (8,*)
      read (8,*) (y(i),i=1,nodtot+1)
      read (8,*)
      read (8,*) jangle
      read (8,*)
      read (8,*) (alphas(j),j=1,jangle)
      read (8,*)
      read (8,*) fmach
      read (8,*)
      read (8,*) rl
C     read (8,*)
C     read (8,*) iterative
      iterative = 1
c      close(8)

	Return
	End


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE PANEL_IBL(numPanel,jjangle, 
     +         xgeom, ygeom, Lxy, attackangles, LAA, S_mach, Reynolds,
     +         PanelOut_name, iblOut_name, summary_name,
     +         Cl_Res, Cd_Res, Cm_Res,iterative)
c
c   main program
c
c   logical control pan-ibl interaction
c
c  -------------------------------------------------------------
c  note : nbt - max total grids on the body;
c         nwt - max total grids along the wake (including t. e.)
c         nbwt- max total grids on the body & wake
c         ncyc- max iteraction cycles per angle of attack
c         nangle - max angles of attack allowed per calculation
c --------------------------------------------------------------
c
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      parameter (nbwtm1=250,ncyc=35,nangle=25)
      common /bod/ nodtot,x(nbt),y(nbt),xmid(nbtm1),ymid(nbtm1),
     +             costhe(nbtm1),sinthe(nbtm1)
      common /cpd/ ue(nbwtm1),cp(nbwtm1),uew(nwtm1,2),
     +             cpw(nwtm1,2)
      common /num/ pi,pi2inv, irestart
      common /coupl/dlsp(nbt),vnp(nbt),qw(nwt),delw(nwt,2)
     1              ,uedlspp(nbt),uedelw(nwt,2)
      common /comp / fmach,c2,c6,c7,c8
      common /wwake/ nw,xw(nwt),yw(nwt),sinw(nwt),cosw(nwt),qwm(nwt)
     1               ,xwmid(nwt),ywmid(nwt)
      common /blinp/ nt(2),nl(2),cdf(2),rl,
     1               xcbl(nbwt,2),ycbl(nbwt,2),uebl(nbwt,2),
     1               dlspbl(nbwt,2),vnpbl(nbwt,2),uedlspbl(nbwt,2),
     1               uedlsp(nbwt),cfs(nbwt),relax,cfsbl(nbwt,2)
      common /pmisc1/ clsave(ncyc,nangle),cmsave(ncyc,nangle),dcl,
     1                cdsave(ncyc,nangle),
     1                alphas(nangle),icyclets(nangle),dcll(ncyc,nangle)
      logical clconv

	Dimension xgeom(LXY), ygeom(LXY), attackangles(LAA) 
	CHARACTER*(*)  PanelOut_name,  
     1              iblOut_name, summary_name
      real*4 Cl_Res(jjangle), Cd_Res(jjangle), Cm_Res(jjangle) 
C
C clean everything
c 
	do i = 1, nbt
	  dlsp(i)=0.0
	  vnp(i)=0.0
	end do
	do i = 1, nwt
	   qw(i)=0.
	   delw(i,1)=0
	   delw(i,2)=0
	end do 
        do i = 1,nbwt
           uedlsp(i) = 0.0
        enddo
c
      nodtot = numPanel
	fmach = S_mach
      rl = Reynolds
	do i= 1, numPanel+1
	  x(i) = xgeom(i)
	  y(i) = ygeom(i)
	enddo
      jangle = jjangle
	do i = 1, jangle
	  alphas(i) = attackangles(i)
	end do

c  descriptions of the input data
c
c  nodtot  - total number of panels. the airfoil-defined grid points 
c            should be equal to nodtot + 1
c  jangle  - total # of angles of attack interested.   
c  x,y     - normalized airfoil cartesian coordinates (with chord length) 
c            ; the data are input clockwise from the lower trailing edge to 
c            upper trailing edge
c  alphas  - values of "jangle" angles of attack; they should be input in
c            ascending order, i.e. from lower to higher values
c  fmach   - free-stream mach number; since the viscous flow calculations  
c            are based on incompressible flows, high mach number combined 
c            with high angles of attack are not recommended
c  rl      - Reynolds number based on the freestream velocity & chord 
c            length; for inviscid flow calculations only, set rl = 0.0
c  iterative - 0  there is no pan - ibl interation
c              1  cal. is done by performing pan - ibl interation
c
      if(iterative .eq. 0) then
         icyclet = 1
         iswpt   = 25
      else
         icyclet = 30
         iswpt   = 8
      endif
c

      if(rl .gt. 0.0) then 
         if(icyclet .gt. 1) then
            write(7,3001)
3001        format(30x,'  panel_ibl iterative results '/
     1             30x,'  ***************** ')
         else
            write(7,3102)
3102        format(30x,'  panel_ibl non_iterative results '/
     1             30x,'  ***************** ')
         endif


         write(12,3002)
3002     format(//30x,'  inverse boundary layer results '/
     1          30x,'  ****************************** ')
      else
         write(7,3003)
3003     format(//30x,'  panel results '/
     1          30x,'  ************* ')
      endif

c
c  regrid panel distribution
 
      ntot = nodtot+1
      call regrid(ntot,x,y)
      nodtot = ntot - 1
c
c  set parameters
c
      interp  = 0
      dcl     = 0.0005    ! converged criterior for cl
      nw      = 40        ! # of wake grids including t. e.
      irestart= 0
c
c  loop for pan_ibl calculations
c
      iang   = 0
      iobkk  = 0
101   iang   = iang + 1
      alpha  = alphas(iang)
      do icycle = 1,icyclet
c
c  calculate inviscid flow
c
         icyc = icycle
         call panel(alpha,iang,icycle,iobkk,cl,cm)
c
         do i = 1,nw
            xwmid(i) = 0.5*(xw(i) + xw(i+1))
            ywmid(i) = 0.5*(yw(i) + yw(i+1))
         enddo
c
         clsave(icycle,iang) = cl
         cmsave(icycle,iang) = cm
         if(rl .le. 0.0) goto 99 
c
c  check for convergence of the currecnt angle attack iteration
c
         clconv = .false.
         if(icycle .ge. 2) then
            dcll(icyc,iang)=clsave(icyc,iang)-clsave(icyc-1,iang)
            if(abs(dcll(icyc,iang)) .lt. dcl) clconv = .true.
         else
            dcll(1,iang) = 0.0
         endif
c
c  determine relaxation parameter for updating solutions
c
         if(clconv) then
            relax = 1.0
         else
            if(icycle .le. 15) then
               relax = amax1(1.0,1.10-0.05*icycle)
            else
     1         relax = amax1(0.7,1.0-0.05*(icycle-15))
            endif
            if(icycle .ge. 3) then
               dclprod = dcll(icyc,iang) * dcll(icyc-1,iang)
               if(dclprod .lt. 0.0 ) then
                  if(abs(dcll(icyc,iang)) .ge. 
     1                   abs(dcll(icyc-1,iang)) ) then
                      relax = amin1(0.7,relax)
                  endif
               else
                  if(abs(dcll(icyc,iang)) .le.
     1                   abs(dcll(icyc-1,iang)) ) then
                      relax = amax1(1.0,relax * 1.05) 
                  endif
               endif 
            endif 
         endif       
         if(icycle .ge. icyclet) clconv = .true.

c
c  separate the airfoil into upper & lower surfaces for ibl calculations
c
         call uplows(interp,ic,nbwt,nodtot,xmid,ymid,ue,dlsp,uedlsp,nl,
     1               xcbl,ycbl,uebl,dlspbl,uedlspbl)
c
c  bl calculation using ibl code
c
         do k = 1,2
            nt(k) = nl(k)
            if(nw .gt. 1) then
               nt(k) = nl(k) + nw
               j     = nl(k)
               do i = 1,nw
                  j = j + 1
                  xcbl(j,k) = xwmid(i)
                  ycbl(j,k) = ywmid(i)
                  uebl(j,k) = uew(i,k)
                  dlspbl(j,k) = delw(i,k)
                  uedlspbl(j,k)= uedelw(i,k)
               enddo
            endif
            call blmain(istop,irestart,k,icycle,alpha,clconv,nt(k),
     1           nl(k),rl,cdf(k),xcbl(1,k),ycbl(1,k),uebl(1,k),
     1           dlspbl(1,k),vnpbl(1,k),uedlspbl(1,k),cfsbl(1,k),
     1           fmach,iswpt,iang)
            if(istop .eq. 1) then
               write(7,7000) alpha,icycle,k
7000           format('sols. fail at alpha =',f8.3,' cyc=',i3,
     1                'surface = ',i2/
     1                'cals. stop') 
               jangle = iang-1 
               goto 199
            endif
         enddo
         cdsave(icycle,iang) = cdf(1) + cdf(2)
         if(icycle .gt. 2) then
            relaxcd = 0.75
            cdsave(icycle,iang) = relaxcd*cdsave(icycle,iang) +
     1                       (1.-relaxcd)*cdsave(icycle-1,iang) 
         endif
c
c
c  option uses the neumann output control points as the b. l. grid point
c  dlspbl & vnpbl are calcualted at the grid points and no interpolation
c  is needed
c  merge upper & lower blowing & displacement thickness for pan code
c
         i   = 0
c
c   lower surface
         do j = nl(1),2,-1
            i = i + 1
            dlsp(i) = relax*dlspbl(j,1)+(1.-relax)*dlsp(i)
            vnp (i) = relax*vnpbl(j,1) +(1.-relax)*vnp(i)
            uedlsp(i)= relax*uedlspbl(j,1)+(1.-relax)*uedlsp(i)
            cfs(i)   = cfsbl(j,1)
         enddo
c
         if(ic .le. 2) then
            i = i + 1
            dlspi = 0.5 * (dlspbl(2,1) + dlspbl(2,2))
            vnpi  = 0.5 * (vnpbl (2,1) + vnpbl (2,2))
            uedlspi = 0.5 * (uedlspbl (2,1) + uedlspbl (2,2))
            dlsp(i) = relax*dlspi + (1.-relax)*dlsp(i)
            vnp(i)  = relax*vnpi + (1.-relax)*vnp(i)
            uedlsp(i) = relax*uedlspi + (1.-relax)*uedlsp(i)
            cfs(i)  = 0.5 * (cfsbl(2,1) + cfsbl(2,2))
         endif
c
c  upper surface
         do j = 2,nl(2)
            i = i  + 1
            dlsp(i) = relax*dlspbl(j,2)+(1.-relax)*dlsp(i)
            uedlsp(i)= relax*uedlspbl(j,2)+(1.-relax)*uedlsp(i)
            vnp (i) = relax*vnpbl (j,2)+(1.-relax)*vnp(i)
            cfs (i) = cfsbl (j,2)
         enddo
         if( i .ne. nodtot ) then
             irror = 1
             write(6,*) 'b.l. & pan grids not compatible; cal. stops'
	       close(6)
	       close(5)
	       close(7)
	       close(8)
	       close(9)
	       close(12)
	       return
C             stop
         endif
c
         if(nw .ne. 0) then
c
c  wake data
            do k = 1,2
               nwbl = nt(k) - nl(k)
               if(nw .ne. nwbl) then
                  ierror = 2
                  write(6,*) 'bl & pan wake grids  not compatible;',
     1                          ' cal. stops'
               endif
               i   = 0
               j   = nl(k)
               do i = 1,nw
                  j = j + 1
                  delw(i,k) = relax*dlspbl(j,k)+(1.-relax)*delw(i,k)
                  uedelw(i,k)=relax*uedlspbl(j,k)+(1.-relax)*
     1                           uedelw(i,k)
               enddo
            enddo
            i    = 0
            j1   = nl(1)
            j2   = nl(2)
            do i = 1,nw
               j1 = j1 + 1
               j2 = j2 + 1
               qw(i)=relax*(vnpbl(j1,1)+vnpbl(j2,2))+(1.-relax)*qw(i)
            enddo
         endif
c
         irestart = 1
c
c  converge for currecnt angle attack iteration
c
         if(clconv) goto 99
      enddo
c
99    icyclets(iang) = icyc
c
c ----------------------------------------------------------------------
c
c  write the results for each angle of attack in the unit 7
c
      if(rl .le. 0.0) then
c
c  inviscid flow only
c
         write(7,1001) alpha,fmach,rl,icyc 
         write(7,1002) cl,cm
         write(7,4003)
         do i = 1,nodtot
            write(7,4006) i,xmid(i),ymid(i),cp(i),ue(i)
         enddo
         if(nw .gt.0) then
            write(7,4005)
            do i = 1,nw
               write(7,4006) i,xwmid(i),ywmid(i),uew(i,1),uew(i,2) 
            enddo
         endif

      else
         write(7,1001) alpha,fmach,rl,icyc 
         write(7,1002) cl,cm
         write(7,1003)
         do i = 1,nodtot
            write(7,1006) i,xmid(i),ymid(i),cp(i),ue(i),dlsp(i),vnp(i)
     1                    ,cfs(i)
         enddo
         if(nw .gt.0) then
            write(7,1005)
            do i = 1,nw
               write(7,1006) i,xwmid(i),ywmid(i),uew(i,1),uew(i,2),
     1                    delw(i,1),delw(i,2),qw(i)
            enddo
         endif
      endif
c
      if(iterative .eq. 0) then
         irestart = 0
c
c  reset b. l. data = 0.0
c
         do i = 1,nbt
            dlsp(i) = 0.0
            vnp(i)  = 0.0
            cfs(i)  = 0.0
         enddo
         do i = 1,nwt
            delw(i,1) = 0.0
            delw(i,2) = 0.0
            qw(i)     = 0.0
c           uedelw(i,1) = 0.0 
c           uedelw(i,2) = 0.0 
         enddo
c        do i = 1,nbwt
c           uedlsp(i) = 0.0
c        enddo
      endif
c
      if(iang .lt. jangle ) goto 101
c
c  write lifthistory data
c
199   continue
       write(9,299) rl,fmach 
       if(rl .le. 0.0) goto 305
299   format( ' alpha vs cl, cd & cm for Rc=',e13.4,' Mach=',f10.5,/)     
       write(9,*) " cycle    cl       dcl      cd        cm      ",
     1            " alpha"
       do j = 1,jangle
          alpha = alphas(j)
          do i = 1,icyclets(j)
             write(9,300) i,clsave(i,j),dcll(i,j),cdsave(i,j),
     1                   cmsave(i,j),alpha
300         format(i5,5f10.5)
          enddo
         write(9,301)
301      format( 5x," -----------------------------------------------")
      enddo
c
c  write the summary of alpha ~ cl etc. data
c
305   continue
      write(9,*)
      write(9,*) "   j    alpha      cl        cd        cm"
      write(9,*) "   -    -----      --        --        --"
      do j = 1,jangle
         alpha = alphas(j)
         cl    = clsave(icyclets(j),j)
	   cl_res(j) = cl
         if(rl .le. 0.0) then
            cd = 0.0
         else
            cd    = cdsave(icyclets(j),j)
         endif
	   cd_res(j) = cd
         cm    = cmsave(icyclets(j),j)
	   cm_res(j) = cm
         write(9,302) j,alpha,cl,cd,cm 
302      format(i5,4f10.5)
      enddo
c
1001  format(/'alpha =',f8.4,' fmach = ',f6.3,' rl = ',e13.3, 
     1        ' cyc  =',i5) 
1002  format(' cl = ',f10.5,'  cm = ',f10.5,/)
1003  format(' results on the airfoil'/
     1  '   i     xmid      ymid        cp        ue      dlsp',    
     1  '      vw         cfs ')
1005  format(' results along the wake'/
     1  '   i    xwmid     ywmid      uewl      uewu     delwl',
     1  '     delwu       qw ')
1006  format(i4,7f10.6)
4003  format(' results on the airfoil'/
     1  '   i     xmid      ymid        cp        ue ')   
4005  format(' results along the wake'/
     1  '   i    xwmid     ywmid      uewl      uewu ')

4006  format(i4,4f10.6)

C	Print*," "
C	print*,"The computation is successfully finished",
C     +       " Enter any key to close the program window"
C	read(5,*)
c	       close(6)
cc	       close(5)
c   	       close(7)
c	       close(8)
c 	       close(9)
c	       close(12)
	       return
C     stop
      end
 
      subroutine uplows(intp,ic,nbwt,nodtot,x,y,ue,dlsp,uedlsp,nl,xcbl,
     1                  ycbl,uebl,dlspbl,uedlspbl)
c
c  this subroutine performs the function of separating the airfoil into
c  the upper & lower surface for b.l. calculations
c
      dimension x(nodtot),y(nodtot),ue(nodtot),dlsp(nodtot),
     1       uedlsp(nodtot),xcbl(nbwt,2),ycbl(nbwt,2),uebl(nbwt,2),
     1       dlspbl(nbwt,2),uedlspbl(nbwt,2),nl(2)
c     dimension x(*),y(*),ue(*),dlsp(*),xcbl(*),ycbl(8),uebl(*),
c    1          dislbl(*),nl(*)
c
      do i = 1,nodtot
         if(ue(i) .le. 0.0) numup = i
      enddo
c
      istag = 0
      if (abs(ue(numup)) .le. 0.0001 ) then
c
c  lower surface
c
         xcbl(1,1) = x(numup)
         ycbl(1,1) = y(numup)
         uebl(1,1) = 0.0
         dlspbl(1,1) = dlsp(numup)
         uedlspbl(1,1)= uedlsp(numup)
         i         = 1
         do j = numup-1,1,-1
            i = i+1
            xcbl(i,1)  = x(j)
            ycbl(i,1)  = y(j)
            uebl(i,1) = -ue(j)
            dlspbl(i,1) = dlsp(j)
            uedlspbl(i,1)= uedlsp(j)
         enddo
         nl(1) = i
c
         xcbl(1,2) = x(numup)
         ycbl(1,2) = y(numup)
         uebl(1,2) = 0.0
         dlspbl(1,2)= dlsp(numup)
         uedlspbl(1,2)= uedlsp(numup)
         j          = 1
         do i = numup+1, nodtot
            j = j +1
            xcbl(j,2)  = x(i)
            ycbl(j,2)  = y(i)
            uebl(j,2)  = ue(i)
            dlspbl(j,2)= dlsp(i)
            uedlspbl(j,2)= uedlsp(i)
         enddo
         nl(2) = j
         ic    = 1
cc       write(6,*) 'ic,numup,nodtot,nl1,nl2 =',ic,numup,nodtot,nl(1),nl
 
      elseif( abs(ue(numup+1)) .le. 0.0001 ) then
         i      = 1
         xcbl(1,1) = x(numup+1)
         ycbl(1,1) = y(numup+1)
         uebl(1,1) = 0.0
         dlspbl(1,1)= dlsp(numup+1)
         uedlspbl(1,1)= uedlsp(numup+1)
         do j = numup,1,-1
            i = i + 1
            xcbl(i,1)  = x(j)
            ycbl(i,1)  = y(j)
            uebl(i,1) = -ue(j)
            dlspbl(i,1) = dlsp(j)
            uedlspbl(i,1) = uedlsp(j)
         end do
         nl(1) = i
 
         j      = 1
         xcbl(1,2) = x(numup+1)
         ycbl(1,2) = y(numup+1)
         uebl(1,2) = 0.0
         dlspbl(1,2)= dlsp(numup+1)
         uedlspbl(1,2)= uedlsp(numup+1)
         do i = numup+2, nodtot
            j = j + 1
            xcbl(j,2)  = x(i)
            ycbl(j,2)  = y(i)
            uebl(j,2) = ue(i)
            dlspbl(j,2)= dlsp(i)
            uedlspbl(j,2)= uedlsp(i)
         end do
         nl(2) = j
         ic    = 2
      else
         x_stag = x(numup) - ue(numup)*(x(numup+1)-x(numup))
     *              /(ue(numup+1)-ue(numup))
         y_stag = y(numup) - ue(numup)*(y(numup+1)-y(numup))
     *              /(ue(numup+1)-ue(numup))
         d_stag = dlsp(numup) - ue(numup)*(dlsp(numup+1)-dlsp(numup))
     *              /(ue(numup+1)-ue(numup))
         ued_stag=uedlsp(numup)-ue(numup)*(uedlsp(numup+1)-
     *            uedlsp(numup))  /(ue(numup+1)-ue(numup))
 
         xcbl(1,1) = x_stag
         ycbl(1,1) = y_stag
         uebl(1,1) = 0.0
         dlspbl(1,1) = d_stag
         uedlspbl(1,1)= ued_stag
         i    = 1
         do j = numup,1,-1
            i = i + 1
            xcbl(i,1)  = x(j)
            ycbl(i,1)  = y(j)
            uebl(i,1) = -ue(j)
            dlspbl(i,1)= dlsp(j)
            uedlspbl(i,1)= uedlsp(j)
         end do
         nl(1) = i
 
         xcbl(1,2) = x_stag
         ycbl(1,2) = y_stag
         dlspbl(1,2) = d_stag
         uedlspbl(1,2)= ued_stag
         uebl(1,2) = 0.0
         j         = 1
         do i = numup+1, nodtot
            j = j +1
            xcbl(j,2)  = x(i)
            ycbl(j,2)  = y(i)
            uebl(j,2) = abs(ue(i))
            dlspbl(j,2)= dlsp(i)
            uedlspbl(j,2)= uedlsp(i)
         end do
         nl(2) = j
         ic    = 3
      endif
      return
      end
 
      subroutine panel(alpha,iangle,icycle,iobkk,cl,cm)
c
c  subroutine controls the inviscid flow calculation
c
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      parameter (nbwtm1=250,ncyc=30,nangle=15)
 
      common /bod/ nodtot,x(nbt),y(nbt),xmid(nbtm1),ymid(nbtm1),
     +             costhe(nbtm1),sinthe(nbtm1)
      common /num/ pi,pi2inv, irestart
      common/coupl/dlsp(nbt),vnp(nbt),qw(nwt),delw(nwt,2)
     1              ,uedlsp(nbt),uedelw(nwt,2)
 
      common /offbod/ xoff(nbt),yoff(nbt),soff(nbtm1),coff(nbtm1)
      common /obku/ ste(2),cte(2),vnte(2),vtte(2),vtot(2)
     +              ,vntes(2),vttes(2),vtots(2)
     +              ,stes(2),ctes(2)
 
      common /comp / fmach,c2,c6,c7,c8
      common /wwake/ nw,xw(nwt),yw(nwt),sinw(nwt),cosw(nwt),qwm(nwt)
     1               ,xwmid(nwt),ywmid(nwt)
 
ccc   dimension title(20),dbody(nbtm1),sbgeo(nbt),cbgeo(nbt)
      dimension dbody(nbtm1),sbgeo(nbt),cbgeo(nbt)
      dimension dlspgrid(nbt)
c
      pi      = 3.1415926585
      pi2inv  = .5/pi
c
      cosalf  = cos(alpha*pi/180.)
      sinalf  = sin(alpha*pi/180.)
      gam    = 1.4
      c2     = 0.5 * (gam-1.0)
      c6     = c2 * fmach**2
      c7     = gam/(gam-1.0)
      c8     = 0.5 * gam * fmach**2
      do 100 i=1,nodtot
c     xmi and ymi, see eq. (5.3.12)
         xmid(i) = .5*(x(i) + x(i+1))
         ymid(i) = .5*(y(i) + y(i+1))
         dx=x(i+1)-x(i)
         dy=y(i+1)-y(i)
         dist=sqrt(dx*dx+dy*dy)
c     see eq. (5.3.2)
         sinthe(i)=dy/dist
         costhe(i)=dx/dist
         dbody(i)=dist
100   continue
      if (irestart .eq. 1 ) then
         do i=1,nw
            dx=xw(i+1)-xw(i)
            dy=yw(i+1)-yw(i)
            dist=sqrt(dx*dx+dy*dy)
            sinw(i)=dy/dist
            cosw(i)=dx/dist
cc          qwm(i)=0.5*(qw(i)+qw(i+1))
            qwm(i)=qw(i)
         enddo
      else
         do i = 1,nw
            qw(i)  = 0.0
            qwm(i) = 0.0
         enddo
      endif
c
      do 200 i=2,nodtot
      sbgeo(i)=(sinthe(i)*dbody(i-1)+sinthe(i-1)*dbody(i))/(dbody(i-1)+
     +          dbody(i))
200   cbgeo(i)=(costhe(i)*dbody(i-1)+costhe(i-1)*dbody(i))/(dbody(i-1)+
     +          dbody(i))
      sbgeo(1)=2.0*sinthe(1)-sbgeo(2)
      sbgeo(nodtot+1)=2.0*sinthe(nodtot)-sbgeo(nodtot)
      cbgeo(1)=2.0*costhe(1)-cbgeo(2)
      cbgeo(nodtot+1)=2.0*costhe(nodtot)-cbgeo(nodtot)
c
c  kcc - note that dlsp are stored at mid points
c
      if(irestart .eq. 1 ) then
         do j = 2,nodtot
            dlspgrid(j) = 0.5 * (dlsp(j-1) + dlsp(j))
         enddo
         dlspgrid(1) = dlspgrid(2) + (dlspgrid(2)-dlspgrid(3))
     1              /(x(2)-x(3))*(x(1)-x(2))
         dlspgrid(nodtot+1) = dlspgrid(nodtot)+ (dlspgrid(nodtot)
     1       -dlspgrid(nodtot-1))/(x(nodtot)-x(nodtot-1))*
     1       (x(nodtot+1)-x(nodtot))
      else
         do j = 1,nodtot+1
            dlspgrid(j) = 0.0
            dlsp(j)     = 0.0
         enddo
      endif
c
cckcc
      do 220 j=1,nodtot+1
c     xoff(j)=x(j)-sbgeo(j)*dlsp(j)
      xoff(j)=x(j)-sbgeo(j)*dlspgrid(j)
c     yoff(j)=y(j)+cbgeo(j)*dlsp(j)
      yoff(j)=y(j)+cbgeo(j)*dlspgrid(j)
220   continue
      cfac    = amin1(1.0,0.5 + 0.5*(iangle-1+icycle))
cc    cfac    = 0.5
c     write(44,*) 'iangle,icycle = ',iangle,icycle
c     write(44,*) 'it,ste1,ste2,cte1,cte2,vt1,vt2,vn1,vn2,dvtot'
 
      do 210 i=1,nodtot
      dx=xoff(i+1)-xoff(i)
      dy=yoff(i+1)-yoff(i)
      dist=sqrt(dx*dx+dy*dy)
c     see eq. (5.3.2)
      soff(i)=dy/dist
      coff(i)=dx/dist
210   continue
c
c  added by kcc
ccc      ste(1)=soff(1)
ccc      ste(2)=soff(nodtot)
ccc      cte(1)=coff(1)
ccc      cte(2)=coff(nodtot)
      if(iobkk.eq.0) then
         ste(1)=soff(1)
         ste(2)=soff(nodtot)
         cte(1)=coff(1)
         cte(2)=coff(nodtot)
      else
         ste(1)=stes(1)*cfac + (1.-cfac)*soff(1)
         ste(2)=stes(2)*cfac + (1.-cfac)*soff(nodtot)
         cte(2)=ctes(2)*cfac + (1.-cfac)*coff(nodtot)
         cte(1)=ctes(1)*cfac + (1.-cfac)*coff(1)
      endif
      stes(1) = ste(1)
      stes(2) = ste(2)
      ctes(2) = cte(2)
      ctes(1) = cte(1)
c
c  kcc - end
c
      it=1
 240  continue
      call coef(sinalf,cosalf)
      call gauss(1)
      call obkuta(dvtot,sinalf,cosalf,index)
      call wakeds(sinalf,cosalf)
c     write(44,1111) it,ste(1),ste(2),cte(1),cte(2)
c    *        ,vtte(1),vtte(2),vnte(1),vnte(2),dvtot
1111  format(i3,9f10.6)
 
c
c  added by kcc
      if(it.eq.1) then
         do ii = 1,2
            vttes(ii) = vtte(ii)
            vntes(ii) = vnte(ii)
            vtots(ii) = vtot(ii)
         enddo
      endif
c  kcc - end
c
      if(index.eq.1) go to 230
c
c  added by kcc
ccc   if(it.gt.10) stop
      if(it.gt.10) then
         dvtmax = -999.
         do ii = 1,2
            vtm = amax1(1.0,vtot(ii))
            dvt = abs(vtots(ii) - vtot(ii))
            if(dvt.gt. dvtmax) dvtmax = dvt
         enddo
         if(dvtmax .lt. 0.0001) goto 230
      endif
c
c  kcc - end
      if(it.gt.20) stop 10
      it=it+1
      do 250 ii=1,2
      steo=ste(ii)
      cteo=cte(ii)
      ste(ii)=(steo*vtte(ii)+cteo*vnte(ii))/vtot(ii)
      cte(ii)=(cteo*vtte(ii)-steo*vnte(ii))/vtot(ii)
250   continue
      ste(1)=-ste(1)
      cte(1)=-cte(1)
c
c  added kcc
      if(it .gt. 5) then
         rexte = 0.75
         do ii = 1,2
            ste(ii)=rexte * ste(ii)+(1.-rexte)*stes(ii)
            cte(ii)=rexte * cte(ii)+(1.-rexte)*ctes(ii)
         enddo
      endif
c
c  kcc  - end
 
c
c  save data for relaxation
c
      stes(1) = ste(1)
      stes(2) = ste(2)
      ctes(2) = cte(2)
      ctes(1) = cte(1)
      do ii = 1,2
         vttes(ii) = vtte(ii)
         vntes(ii) = vnte(ii)
         vtots(ii) = vtot(ii)
      enddo
 
      go to 240
230   continue
c
c  save data for the new calculations
      stes(1) = ste(1)
      stes(2) = ste(2)
      ctes(2) = cte(2)
      ctes(1) = cte(1)
      iobkk = 1
c
      call vpdis(sinalf,cosalf)
      call vpdwk(sinalf,cosalf)
      call clcm(sinalf,cosalf,cl,cm)
c
c                 stop
      return
      end
c
      subroutine coef(sinalf,cosalf)
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      parameter (nbwtm1=250,ncyc=30,nangle=15)
 
      common /bod/ nodtot,x(nbt),y(nbt),xmid(nbtm1),ymid(nbtm1),
     +             costhe(nbtm1),sinthe(nbtm1)
      common /cof/ a(nbt,nbt),bv(nbt),kutta
      common /num/ pi,pi2inv, irestart
      common/coupl/dlsp(nbt),vnp(nbt),qw(nwt),delw(nwt,2)
     1              ,uedlsp(nbt),uedelw(nwt,2)
 
      common /offbod/ xoff(nbt),yoff(nbt),soff(nbtm1),coff(nbtm1)
      common /obku/ ste(2),cte(2),vnte(2),vtte(2),vtot(2)
     +              ,vntes(2),vttes(2),vtots(2)
     +              ,stes(2),ctes(2)
 
      common /wwake/ nw,xw(nwt),yw(nwt),sinw(nwt),cosw(nwt),qwm(nwt)
     1               ,xwmid(nwt),ywmid(nwt)
 
      kutta   = nodtot + 1
      do 120  i = 1,nodtot
      a(i,kutta) = 0.0
      do 110  j = 1,nodtot
      flog    = 0.0
      ftan    = pi
      if (j .eq. i)     go to 100
      dxj     = xmid(i) - x(j)
      dxjp    = xmid(i) - x(j+1)
      dyj     = ymid(i) - y(j)
      dyjp    = ymid(i) - y(j+1)
c     flog is ln(r(i,j+1)/r(i,j)), see eq. (5.3.12)
      flog    = .5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
c     ftan is beta(i,j), see eq. (5.3.12)
      ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
c     ctimtj is cos(theta(i)-theta(j))
 100  ctimtj  = costhe(i)*costhe(j) + sinthe(i)*sinthe(j)
c     stimtj is sin(theta(i)-theta(j))
      stimtj  = sinthe(i)*costhe(j) - costhe(i)*sinthe(j)
c     elements of the coefficient matrix, a(i,j), see eq. (5.4.1a)
      a(i,j)  = pi2inv*(ftan*ctimtj + flog*stimtj)
      b       = pi2inv*(flog*ctimtj - ftan*stimtj)
c     elements of the coefficient matrix, a(i,n+1), see eq. (5.4.1b)
      a(i,kutta) = a(i,kutta) + b
 110  continue
      waq=0.0
	if (irestart .ne. 1 ) then
        do 111 j=1,nw
        dxj     = xmid(i) - xw(j)
        dxjp    = xmid(i) - xw(j+1)
        dyj     = ymid(i) - yw(j)
        dyjp    = ymid(i) - yw(j+1)
        flog    =.5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
        ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
        ctimtj  = costhe(i)*cosw(j) + sinthe(i)*sinw(j)
        stimtj  = sinthe(i)*cosw(j) - costhe(i)*sinw(j)
        aa      = pi2inv*(ftan*ctimtj + flog*stimtj)
 111    waq     = waq+aa*qwm(j)
      endif
c     elements of vector b for i=1,...,n, see eq. (5.5.1)
      bv(i) = sinthe(i)*cosalf - costhe(i)*sinalf
ckcc    +        +0.5*(vnp(i)+vnp(i+1))-waq
     +        +vnp(i)-waq
 120  continue
200   do 90   j = 1,kutta
 90   a(kutta,j)   = 0.0
      waq=0.0
      i=1
      ii=1
230   continue
      xmido   = .5*(xoff(i) + xoff(i+1))
      ymido   = .5*(yoff(i) + yoff(i+1))
ckcc     dlsh=0.5*(dlsp(i)+dlsp(i+1))
      dlsh=dlsp(i)
      do 210  j = 1,nodtot
      flog=0.0
      ftan=pi
      if(j.eq.i .and. dlsh.lt.0.0001)go to 209
      dxj     = xmido - x(j)
      dxjp    = xmido - x(j+1)
      dyj     = ymido - y(j)
      dyjp    = ymido - y(j+1)
c     flog is ln(r(i,j+1)/r(i,j)), see eq. (5.3.12)
      flog    = .5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
c     ftan is beta(i,j), see eq. (5.3.12)
      ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
c     ctimtj is cos(theta(i)-theta(j))
 209  ctimtj  = cte(ii)*costhe(j) + ste(ii)*sinthe(j)
c     stimtj is sin(theta(i)-theta(j))
      stimtj  = ste(ii)*costhe(j) - cte(ii)*sinthe(j)
c     see eq. (5.5.10b)
      aa      = pi2inv*(ftan*ctimtj + flog*stimtj)
c     see eq. (5.5.10a)
      b       = pi2inv*(flog*ctimtj - ftan*stimtj)
c     elements of the coefficient matrix, a(n+1,j), see eq. (5.5.11b)
      a(kutta,j) = a(kutta,j) - b
c     elements of the coefficient matrix, a(n+1,n+1), see eq. (5.5.11b)
      a(kutta,kutta) = a(kutta,kutta) +aa
210   continue
      if (irestart .ne. 1 ) then
        do 211 j=1,nw
        dxj     = xmido - xw(j)
        dxjp    = xmido - xw(j+1)
        dyj     = ymido - yw(j)
        dyjp    = ymido - yw(j+1)
        flog    = .5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
        ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
        ctimtj  = cte(ii)*cosw(j) + ste(ii)*sinw(j)
        stimtj  = ste(ii)*cosw(j) - cte(ii)*sinw(j)
        b       = pi2inv*(flog*ctimtj - ftan*stimtj)
        waq     = waq-b*qwm(j)
211     continue
      endif
      if(i.eq.nodtot) go to 220
      i=nodtot
      ii=2
      go to 230
220   continue
c     elements of vector b for i=n+1, see eq. (5.5.11b)
      bv(kutta) = - (cte(1) + cte(2))*cosalf
     +            - (ste(1) + ste(2))*sinalf - waq
      return
      end
      subroutine clcm(sinalf,cosalf,cl,cm)
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      parameter (nbwtm1=250,ncyc=30,nangle=15)
 
      common /bod/ nodtot,x(nbt),y(nbt),xmid(nbtm1),ymid(nbtm1),
     +             costhe(nbtm1),sinthe(nbtm1)
      common /cpd/ ue(nbwtm1),cp(nbwtm1),uew(nwtm1,2),
     +             cpw(nwtm1,2)
      common /comp / fmach,c2,c6,c7,c8
      cfx     = 0.0
      cfy     = 0.0
      cm      = 0.0
c
c     note that correction has been done in sub. vpdis & vpwdis
c     compressibility correction based on karman-tsien formula
c
c     beta    = sqrt(1.-fmach**2)
c     do  i =1,nodtot
c     uei     = ue(i)
c     cp(i)   = cp(i)/(beta+(fmach**2/(1.+beta))*cp(i)*0.5)
c     v2      = amax1(0.0,(1.0-(1.0+c8*cp(i))**(1./c7))/c6+1.0)
c     ue(i)   = sqrt(v2)
c     if(uei .lt. 0) ue(i) =- ue(i)         !kcc
c     end do
c
      do 100  i = 1,nodtot
      dx      = x(i+1) - x(i)
      dy      = y(i+1) - y(i)
      cfx     = cfx + cp(i)*dy
      cfy     = cfy - cp(i)*dx
      cm      = cm + cp(i)*(dx*xmid(i) + dy*ymid(i))
 100  continue
      cl      = cfy*cosalf - cfx*sinalf
      return
      end
      subroutine gauss(m)
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      common /cof/ a(nbt,nbt),b(nbt,1),n
      do 100  k = 1,n-1
      kp      = k + 1
      do 100  i = kp,n
      r       = a(i,k)/a(k,k)
      do 200  j = kp,n
 200  a(i,j)  = a(i,j) - r*a(k,j)
      do 100  j = 1,m
 100  b(i,j)    = b(i,j) - r*b(k,j)
      do 300  k = 1,m
      b(n,k) = b(n,k)/a(n,n)
      do 300  i = n-1,1,-1
      ip      = i + 1
      do 400  j = ip,n
 400  b(i,k)    = b(i,k) - a(i,j)*b(j,k)
 300  b(i,k)    = b(i,k)/a(i,i)
      return
      end
      subroutine obkuta(dvtot,sinalf,cosalf,index)
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      parameter (nbwtm1=250,ncyc=30,nangle=15)
 
      common /bod/ nodtot,x(nbt),y(nbt),xmid(nbtm1),ymid(nbtm1),
     +             costhe(nbtm1),sinthe(nbtm1)
      common /cof/ a(nbt,201),bv(nbt),kutta
      common /num/ pi,pi2inv, irestart
      common/coupl/dlsp(nbt),vnp(nbt),qw(nwt),delw(nwt,2)
     1              ,uedlsp(nbt),uedelw(nwt,2)
 
      common /offbod/ xoff(nbt),yoff(nbt),soff(nbtm1),coff(nbtm1)
      common /obku/ ste(2),cte(2),vnte(2),vtte(2),vtot(2)
     +              ,vntes(2),vttes(2),vtots(2)
     +              ,stes(2),ctes(2)
 
      common /wwake/ nw,xw(nwt),yw(nwt),sinw(nwt),cosw(nwt),qwm(nwt)
     1               ,xwmid(nwt),ywmid(nwt)
 
      dimension q(nbtm1)
      do 50   i = 1,nodtot
 50   q(i)    = bv(i)
      gamma   = bv(kutta)
      i=1
      ii=1
   55 continue
      xmido   = .5*(xoff(i) + xoff(i+1))
      ymido   = .5*(yoff(i) + yoff(i+1))
c     contribution to vt(i) from freesream velocity, see eq. (5.5.4b)
      vtang   = cosalf*cte(ii) + sinalf*ste(ii)
c     contribution to vn(i) from freesream velocity, see eq. (5.5.4a)
      vnoff   = sinalf*cte(ii) - cosalf*ste(ii)
ckcc       dlsh=0.5*(dlsp(i)+dlsp(i+1))
      dlsh=dlsp(i)
      do 120  j = 1,nodtot
      flog    = 0.0
      ftan    = pi
      if(j.eq.i.and.dlsh.lt.0.0001)go to 100
 300  dxj     = xmido - x(j)
      dxjp    = xmido - x(j+1)
      dyj     = ymido - y(j)
      dyjp    = ymido - y(j+1)
c     flog is ln(r(i,j+1)/r(i,j)), see eq. (5.3.12)
      flog    = .5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
c     ftan is beta(i,j), see eq. (5.3.12)
      ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
c     ctimtj is cos(theta(i)-theta(j))
 100  ctimtj  = cte(ii)*costhe(j) + ste(ii)*sinthe(j)
c     stimtj is sin(theta(i)-theta(j))
      stimtj  = ste(ii)*costhe(j) - cte(ii)*sinthe(j)
c     aa is bt(i,j)=an(i,j), see eq. (5.5.7b)
      aa      = pi2inv*(ftan*ctimtj + flog*stimtj)
c     b is -at(i,j)=bn(i,j), see eq. (5.5.7c)
      b       = pi2inv*(flog*ctimtj - ftan*stimtj)
c     contribution to vt(i) from singularities, see eq. (5.5.4b)
      vtang   = vtang - b*q(j) + gamma*aa
c     contribution to vn(i) from singularities, see eq. (5.5.4a)
      vnoff   = vnoff + gamma*b + aa*q(j)
 120  continue
c
      if (irestart .ne. 1 ) then
	  do j=1,nw
             dxj   = xmido - xw(j)
             dxjp  = xmido - xw(j+1)
             dyj   = ymido - yw(j)
             dyjp  = ymido - yw(j+1)
             flog  = .5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
             ftan  = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
             ctimtj= cte(ii)*cosw(j) + ste(ii)*sinw(j)
             stimtj= ste(ii)*cosw(j) - cte(ii)*sinw(j)
             aa    = pi2inv*(ftan*ctimtj + flog*stimtj)
             b     = pi2inv*(flog*ctimtj - ftan*stimtj)
             vtang = vtang-b*qwm(j)
             vnoff = vnoff+aa*qwm(j)
          end do
      endif
c
      vtte(ii)= vtang
      vnte(ii)= vnoff
      vtot(ii)= sqrt(vtang*vtang+vnoff*vnoff)
      if(ii.eq.2) go to 130
      i=nodtot
      ii=2
      go to 55
130   index=0
      dvtot = abs(vtot(1)-vtot(2))
      if(abs(vtot(1)-vtot(2)).lt.0.000025) index=1
      return
      end
      subroutine vpdis(sinalf,cosalf)
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      parameter (nbwtm1=250,ncyc=30,nangle=15)
 
      common /bod/ nodtot,x(nbt),y(nbt),xmid(nbtm1),ymid(nbtm1),
     +             costhe(nbtm1),sinthe(nbtm1)
      common /cof/ a(nbt,nbt),bv(nbt),kutta
      common /cpd/ ue(nbwtm1),cp(nbwtm1),uew(nwtm1,2),
     +             cpw(nwtm1,2)
      common /num/ pi,pi2inv, irestart
      common/coupl/dlsp(nbt),vnp(nbt),qw(nwt),delw(nwt,2)
     1              ,uedlsp(nbt),uedelw(nwt,2)
      common /offbod/ xoff(nbt),yoff(nbt),soff(nbtm1),coff(nbtm1)
      common /wwake/ nw,xw(nwt),yw(nwt),sinw(nwt),cosw(nwt),qwm(nwt)
     1               ,xwmid(nwt),ywmid(nwt)
 
      common /comp / fmach,c2,c6,c7,c8
      dimension q(nbtm1)
      do 50   i = 1,nodtot
 50   q(i)    = bv(i)
      gamma   = bv(kutta)
      do 130  i = 1,nodtot
      xmido   = .5*(xoff(i) + xoff(i+1))
      ymido   = .5*(yoff(i) + yoff(i+1))
c     contribution to vt(i) from freesream velocity, see eq. (5.5.4b)
      vtang   = cosalf*coff(i) + sinalf*soff(i)
c     contribution to vn(i) from freesream velocity, see eq. (5.5.4a)
      vnoff   = sinalf*coff(i) - cosalf*soff(i)
ckcc         dlsh=0.5*(dlsp(i)+dlsp(i+1))
      dlsh=dlsp(i)
      do 120  j = 1,nodtot
      flog    = 0.0
      ftan    = pi
      if (j .eq. i .and. dlsh. lt. 0.0001) go to 100
 300  dxj     = xmido - x(j)
      dxjp    = xmido - x(j+1)
      dyj     = ymido - y(j)
      dyjp    = ymido - y(j+1)
c     flog is ln(r(i,j+1)/r(i,j)), see eq. (5.3.12)
      flog    = .5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
c     ftan is beta(i,j), see eq. (5.3.12)
      ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
c     ctimtj is cos(theta(i)-theta(j))
 100  ctimtj  = coff(i)*costhe(j) + soff(i)*sinthe(j)
c     stimtj is sin(theta(i)-theta(j))
      stimtj  = soff(i)*costhe(j) - coff(i)*sinthe(j)
c     aa is bt(i,j)=an(i,j), see eq. (5.5.7b)
      aa      = pi2inv*(ftan*ctimtj + flog*stimtj)
c     b is -at(i,j)=bn(i,j), see eq. (5.5.7c)
      b       = pi2inv*(flog*ctimtj - ftan*stimtj)
c     contribution to vt(i) from singularities, see eq. (5.5.4b)
      vtang   = vtang - b*q(j) + gamma*aa
c     contribution to vn(i) from singularities, see eq. (5.5.4a)
      vnoff   = vnoff + gamma*b + aa*q(j)
 120  continue
c
 
	do 121 j=1,nw
        dxj     = xmido - xw(j)
        dxjp    = xmido - xw(j+1)
        dyj     = ymido - yw(j)
        dyjp    = ymido - yw(j+1)
        flog    =.5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
        ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
        ctimtj  = coff(i)*cosw(j) + soff(i)*sinw(j)
        stimtj  = soff(i)*cosw(j) - coff(i)*sinw(j)
        aa      = pi2inv*(ftan*ctimtj + flog*stimtj)
        b       = pi2inv*(flog*ctimtj - ftan*stimtj)
        vtang   = vtang-b*qwm(j)
        vnoff   = vnoff+aa*qwm(j)
 121  continue
 
c
      ue(i)=sqrt(vnoff*vnoff+vtang*vtang)
      vmax    = sqrt(1.+1./c6) -0.0025
      cpmin   = -1./c8 + 0.001
c
      if(vtang.lt.0.0) ue(i)=-ue(i)
c
c  mach number correction using the Karman-Tsien formula
      cp(i)   = 1.0 - ue(i)*ue(i)
      cpi     = cp(i)
      beta    = sqrt(1.-fmach**2)
      clamda  = fmach**2 /(1.+beta**2)
      cp(i)   = cpi/(beta + clamda*cpi*0.5)
      if(cp(i) .lt. cpmin) cp(i) = cpmin
cc    cpp  = cp(i)/(beta+(fmach**2/(1.+beta)*cp(i)*0.5)
      if(cp(i) .gt. cpmin) then
         v2 = amax1(0.0,(1.0-(1.0+c8*cp(i))**(1./c7))/c6+1.0)
      else
         v2 = vmax ** 2
      endif
c
      uee  = sqrt(v2)
      if ( ue(i) .lt. 0 ) uee= -uee
      ue(i)= uee
 130  continue
      return
      end
c
      subroutine vpdwk(sinalf,cosalf)
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      parameter (nbwtm1=250,ncyc=30,nangle=15)
 
      common /bod/ nodtot,x(nbt),y(nbt),xmid(nbtm1),ymid(nbtm1),
     +             costhe(nbtm1),sinthe(nbtm1)
      common /cof/ a(nbt,201),bv(nbt),kutta
      common /cpd/ ue(nbwtm1),cp(nbwtm1),uew(nwtm1,2),
     +             cpw(nwtm1,2)
      common /num/ pi,pi2inv, irestart
      common/coupl/dlsp(nbt),vnp(nbt),qw(nwt),delw(nwt,2)
     1              ,uedlsp(nbt),uedelw(nwt,2)
      common /offbod/ xoff(nbt),yoff(nbt),soff(nbtm1),coff(nbtm1)
      common /wwake/ nw,xw(nwt),yw(nwt),sinw(nwt),cosw(nwt),qwm(nwt)
     1               ,xwmid(nwt),ywmid(nwt)
 
      common /comp / fmach,c2,c6,c7,c8
      dimension q(nbtm1)
      do 50   i = 1,nodtot
 50   q(i)    = bv(i)
      gamma   = bv(kutta)
      do 130  i = 1,nw
         do 125  k = 1,2
            xmido   = .5*(xw(i) + xw(i+1))
            ymido   = .5*(yw(i) + yw(i+1))
            vtang   = cosalf*cosw(i) + sinalf*sinw(i)
            vnoff   = sinalf*cosw(i) - cosalf*sinw(i)
ckcc             dlsh=0.5*(delw(i,k)+delw(i+1,k))
            dlsh=delw(i,k)
            if (k.eq.1) then
               xmido   = xmido+sinw(i)*dlsh
               ymido   = ymido-cosw(i)*dlsh
            else
               xmido   = xmido-sinw(i)*dlsh
               ymido   = ymido+cosw(i)*dlsh
            endif
            do 120  j = 1,nodtot
               dxj     = xmido - x(j)
               dxjp    = xmido - x(j+1)
               dyj     = ymido - y(j)
               dyjp    = ymido - y(j+1)
               flog    =.5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
               ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
               ctimtj  = cosw(i)*costhe(j) + sinw(i)*sinthe(j)
               stimtj  = sinw(i)*costhe(j) - cosw(i)*sinthe(j)
               aa      = pi2inv*(ftan*ctimtj + flog*stimtj)
               b       = pi2inv*(flog*ctimtj - ftan*stimtj)
               vtang   = vtang-b*q(j)+gamma*aa
               vnoff   = vnoff + gamma*b + aa*q(j)
 120        continue
c
               do 121 j=1,nw
                  flog    = 0.0
                  ftan    = pi
                  if (j .eq. i .and. dlsh. lt. 0.0001) go to 100
                  dxj     = xmido - xw(j)
                  dxjp    = xmido - xw(j+1)
                  dyj     = ymido - yw(j)
                  dyjp    = ymido - yw(j+1)
                  flog    = .5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+
     1                      dyj*dyj))
                  ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
 100              ctimtj  = cosw(i)*cosw(j) + sinw(i)*sinw(j)
                  stimtj  = sinw(i)*cosw(j) - cosw(i)*sinw(j)
                  aa      = pi2inv*(ftan*ctimtj + flog*stimtj)
                  b       = pi2inv*(flog*ctimtj - ftan*stimtj)
                  vtang   = vtang-b*qwm(j)
                  vnoff   = vnoff+aa*qwm(j)
 121           continue
c
            uew(i,k)   = sqrt(vnoff*vnoff+vtang*vtang)
            cpw(i,k)   = 1.0 - uew(i,k)*uew(i,k)
 125     continue
 130  continue
 
c
c  compressibility correction based on karman-tsien formula
c
         beta    = sqrt(1.-fmach**2)
         do 157 k = 1,2
            do 157 i = 1,nw
               cpw(i,k) = cpw(i,k)/(beta+(fmach**2/(1.+beta))*cpw(i,k)
     1                 *0.5)
               v2   = amax1(0.0,(1.0-(1.0+c8*cpw(i,k))**(1./c7))/c6+1.0)
               uew(i,k) = sqrt(v2)
  157    continue
      return
      end
      subroutine wakeds(sinalf,cosalf)
      parameter (nbt=201,nwt=51,nbwt=251,nbtm1=200,nwtm1=50)
      parameter (nbwtm1=250,ncyc=30,nangle=15)
 
      common /bod/ nodtot,x(nbt),y(nbt),xmid(nbtm1),ymid(nbtm1),
     +             costhe(nbtm1),sinthe(nbtm1)
      common /num  / pi,pi2inv,irestart
      common/coupl/dlsp(nbt),vnp(nbt),qw(nwt),delw(nwt,2)
     1              ,uedlsp(nbt),uedelw(nwt,2)
      common /wwake/ nw,xw(nwt),yw(nwt),sinw(nwt),cosw(nwt),qwm(nwt)
     1               ,xwmid(nwt),ywmid(nwt)
      common /cof/ a(nbt,nbt),bv(nbt),kutta
 
      dimension xwn(nwt),ywn(nwt)
c
 
      gamma=bv(nodtot+1)
      do 100 i=1,nw
ckcc        qwm(i)=0.5*(qw(i)+qw(i+1))
         qwm(i) = qw(i)
100   continue
c
 
      xwn(1)= (x(1)+x(nodtot+1)) / 2.0
      ywn(1)= (y(1)+y(nodtot+1)) / 2.0
c     dste  = sqrt((x(nodtot+1)-x(nodtot))**2+(y(nodtot+1)-
c    1        y(nodtot))**2)
      dste  = sqrt((x(nodtot)-x(nodtot-1))**2+(y(nodtot)-
     1        y(nodtot-1))**2)
cc    dl    = amax1(0.005, dste)
      dl    = amax1(0.003, dste)
c
      thw=0.5*(atan2(y(1)-y(2),x(1)-x(2))+atan2(y(nodtot+1)-
     1         y(nodtot),x(nodtot+1)-x(nodtot)))
c
      do 500 i=1,nw
         it=0
 120     coswk=cos(thw)
         sinwk=sin(thw)
         xwn(i+1)=xwn(i)+dl*coswk
         ywn(i+1)=ywn(i)+dl*sinwk
         xmido=0.5*(xwn(i)+xwn(i+1))
         ymido=0.5*(ywn(i)+ywn(i+1))
         vtang=cosalf*coswk+sinalf*sinwk
         vn=sinalf*coswk-cosalf*sinwk
         j = 0
         do 150 jj=1,nodtot
         j = j+1
            dxj     = xmido - x(jj)
            dxjp    = xmido - x(jj+1)
            dyj     = ymido - y(jj)
            dyjp    = ymido - y(jj+1)
            flog    = .5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
            ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
            ctimtj  = coswk*costhe(j) + sinwk*sinthe(j)
            stimtj  = sinwk*costhe(j) - coswk*sinthe(j)
            aa      = pi2inv*(ftan*ctimtj + flog*stimtj)
            b       = pi2inv*(flog*ctimtj - ftan*stimtj)
            vtang   = vtang-b*bv(j)+gamma*aa
            vn      = vn+gamma*b+bv(j)*aa
150      continue
c
         if (nw.ne.0) then
            do 160 j=1,nw
               flog=0.0
               ftan=pi
               if(j.eq.i)go to 170
               dxj     = xmido - xw(j)
               dxjp    = xmido - xw(j+1)
               dyj     = ymido - yw(j)
               dyjp    = ymido - yw(j+1)
               flog    =.5*alog((dxjp*dxjp+dyjp*dyjp)/(dxj*dxj+dyj*dyj))
               ftan    = atan2(dyjp*dxj-dxjp*dyj,dxjp*dxj+dyjp*dyj)
170            ctimtj  = coswk*cosw(j) + sinwk*sinw(j)
               stimtj  = sinwk*cosw(j) - coswk*sinw(j)
               aa      = pi2inv*(ftan*ctimtj+flog*stimtj)
               b       = pi2inv*(flog*ctimtj - ftan*stimtj)
               vtang   = vtang-b*qwm(j)
               vn      = vn+aa*qwm(j)
160         continue
         endif
c
         if(abs(vn).ge.0.0001) then
            if (it.gt.40) goto 501
            it=it+1
            dthe=atan2(vn,vtang)
            thw=thw+dthe
            go to 120
         endif
  501    continue
         dl=1.15*dl
  500 continue
c
      do 60 i=1,nw+1
         xw(i)=xwn(i)
         yw(i)=ywn(i)
   60 continue
      do 70 i=1,nw
         dx=xw(i+1)-xw(i)
         dy=yw(i+1)-yw(i)
         dist=sqrt(dx*dx+dy*dy)
         sinw(i)=dy/dist
         cosw(i)=dx/dist
   70 continue
      return
 1000 format(3x,'warning: wake not exactly a dividing streamline ',
     &       'for element',i3)
      end
c
c  regrid panel for better pan_ibl solutions
c
      subroutine regrid(ntot,x,y)
      parameter (nb = 201)
      dimension x(*),y(*)
      dimension xout(nb),yout(nb),sint(nb),sout(nb),
     1          s(nb),d1(nb),d2(nb),d3(nb),d4(nb),
     1          xtemp(nb),ytemp(nb)
c
c  search for x min
c
      xmin  = 9999.
      do i = 1,ntot
        if(x(i).lt. xmin) then
           xmin = x(i)
           imin = i
        endif
      enddo
c
      nt  = ntot
      nh  = max0(76,nt/2 + 1)
      ang = 7.0
      call modcos(nh,ang,s)
c
c  lower surface
      sint(1) = 0.0
      do i = 1,imin
         xtemp(i) = x(imin-i+1)
         ytemp(i) = y(imin-i+1)
      enddo
      do i = 2,imin
         sint(i) = sint(i-1)+sqrt((xtemp(i)-xtemp(i-1))**2+
     1                  (ytemp(i)-ytemp(i-1))**2)
      enddo
      do i = 1,nh
         sout(i)=sint(imin)*s(i)/s(nh)
      enddo
      call diff3 (imin,sint,xtemp,d1,d2,d3,0)
      call intrp3(imin,sint,xtemp,d1,d2,d3,nh,sout,xout,d4)
      call diff3 (imin,sint,ytemp,d1,d2,d3,0)
      call intrp3(imin,sint,ytemp,d1,d2,d3,nh,sout,yout,d4)
      do i = 1,nh
         xtemp(i) = xout(nh-i+1)
         ytemp(i) = yout(nh-i+1)
      enddo
      do i = 1,nh
         xout(i) = xtemp(i)
         yout(i) = ytemp(i)
      enddo
c
c  upper surface
c
      sint(1) = 0.0
      ii     = 1
      do i = imin+1,ntot
         ii     = ii + 1
         sint(ii) = sint(ii-1)+sqrt((x(i)-x(i-1))**2+(y(i)-
     1                  y(i-1))**2)
      enddo
      in  = ii
      do i = 1,nh
         sout(i)=sint(in)*s(i)/s(nh)
      enddo
      call diff3 (in,sint,x(imin),d1,d2,d3,0)
      call intrp3(in,sint,x(imin),d1,d2,d3,nh,sout,xout(nh),d4)
      call diff3 (in,sint,y(imin),d1,d2,d3,0)
      call intrp3(in,sint,y(imin),d1,d2,d3,nh,sout,yout(nh),d4)
      ntot = 2*nh - 1
      do i = 1,ntot
         x(i) = xout(i)
         y(i) = yout(i)
      enddo
      return
      end
c
      subroutine modcos(n, ang, x)
c
c  ********************************************************************
c  *                                                                  *
c  *    generate b.l. x-wise grid using modified cosine distribution  *
c  *                                                                  *
c  ********************************************************************
c
      dimension x(n)
      data crad/57.2957795/, bpi/3.14159265/
      save
c
c   ------------------------------------------------------------
c
      nn     = 2 * n - 1
      en     = float((nn-1)/2)
      tho    = ang/crad
      cto1   = 1. + cos(tho)
      dth    = (bpi - tho) / en
      fi     = float(n - 2)
      do 10 i=n,nn
         fi     = 1.0 + fi
         ii     = i - n + 1
         xii    = tho + fi * dth
         x(ii)  = (1.0 + cos(xii))/cto1
10    continue
      x1     = x(1)
      xn     = x(n)
      ch     = xn -x1
      fn1    = float(n-1)
      do 20 i=1,n
         x(i)   = (x(i)-x1)/ch
20    continue
      if(x(2).lt.0.35 * x(3)) x(2) = 0.35 * x(3)
c
      return
      end
 
      subroutine diff3 (n,x,f,fp,fpp,fppp,iend)
c
c  ********************************************************************
c  *                                                                  *
c  *    calculate 1st, 2nd & 3rd derivatives of the function f        *
c  *                                                                  *
c  ********************************************************************
c
      dimension x(n),f(n),fp(n),fpp(n),fppp(n)
c
c  first derivatives using weighted angles
c
      do 5 i = 1,n-1
         dx     = x(i+1)-x(i)
         df     = f(i+1)-f(i)
         fpp(i) = atan2(df,dx)
   5  continue
c
      do 10 i = 2,n-1
         dx1    = x(i) - x(i-1)
         dx2    = x(i+1) - x(i)
         ang    = (dx2*fpp(i-1)+dx1*fpp(i))/(dx1+dx2)
         fp(i)  = tan(ang)
   10 continue
      if(iend.eq.0) then
c
c  iend    = 0 , extrapolate for end values
c
         fp(1)   = 2.*tan(fpp(1))   - fp(2)
         fp(n)   = 2.*tan(fpp(n-1)) - fp(n-1)
      else
c
c  iend    = 1, derivatives are continuous across ends
c
         dx1     = x(2) - x(1)
         dxn     = x(n) - x(n-1)
         ang     = (dx1*fpp(n-1) + dxn*fpp(1)) / (dx1 + dxn)
         fp(1)   = tan(ang)
         fp(n)   = fp(1)
      endif
c
c  second & third derivatives using cubic fits
c
      do 30 i = 2,n
         i1      = i - 1
         dx      = x(i) - x(i1)
         dx3     = - dx ** 3 / 6.0
         df      = f(i) - f(i1)
         dfp     = fp(i) - fp(i1)
         pfp     = fp(i) + fp(i1)
         fppp(i1)= (2.*df - dx*pfp) / dx3
         fpp(i1) = dfp/dx - dx*fppp(i1)/2.
   30 continue
      fppp(n) = fppp(n-1)
      fpp (n) = fpp (n-1) + fppp(n-1) * (x(n) - x(n-1))
c
      return
      end
      subroutine intrp3 (n1,x1,f1,fp1,fpp1,fppp1,n2,x2,f2,fp2)
c
c  ********************************************************************
c  *                                                                  *
c  *    this subroutine performs cubic interpolation                  *
c  *                                                                  *
c  ********************************************************************
c
c  function description :
c
c     given the values of a function (f1) and its derivatives
c     at n1 values of the independent variable (x1)
c
c     find the values of the function (f2) & derivative fp2
c     at n2 values of the independent variable (x2)
c
c     x1 and x2 must be in ascending order
c
      dimension x1(n1),f1(n1),fp1(n1),fpp1(n1),fppp1(n1),
     1          x2(n2),f2(n2),fp2(n2)
c
      dx1min  = x1(n1) - x1(1)
      tl      = dx1min
      do 10 i = 2,n1
         dx   = x1(i) - x1(i-1)
         dx1min = amin1(dx1min, dx)
   10 continue
      eps     = amin1(1.e-7*tl, 1.e-4 * dx1min )
      jt      = 2
      do 50 i = 1,n2
         do 20 j = jt,n1
            if (x1(j).ge.x2(i)) go to 30
   20    continue
         j       = n1
   30    if (abs(x1(j)-x2(i)) .lt. eps) then
            f2(i) = f1(j)
            fp2(i)= fp1(j)
         else
            jt      = j
            j1      = j - 1
            dxx     = x2(i)-x1(j1)
            dxx2    = dxx*dxx/2.
            dxx3    = dxx2*dxx/3.
            f2(i)   = f1(j1)+dxx*fp1(j1)+dxx2*fpp1(j1)+dxx3*fppp1(j1)
            fp2(i)  = fp1(j1)+dxx*fpp1(j1)+dxx2*fppp1(j1)
         endif
   50 continue
c
      return
      end
 
C********************************************************************
C  Inverse Boundary Layer Program
C*******************************************************************
c***************************************************************************
c  control subroutine for ibl calculation   
c  the subroutine is called from 'main' program
c
      subroutine blmain(istop,irestart,isurf,icyc,alpha,clconv,nxtin,
     1     nxtein,rlin,cdf,xcbl,ycbl,uebl,dlspbl,vnpbl,uedispbl,cfsbl, 
     1     fmach,iswptin,iang) 
c
c  explanation of the arguments
c     irestart   = 1, if the caculation starts using the viscous inf.
c                from the previsous calculation, otherwise ie 0
c     isurf      = 1, lower surface; = 2, upper surface 
c     icyc       cycle counter
c     alpha      angle of attack
c     clconv     logical for cl convergence in pbn_ibl interaction 
c     nxtin      total number of b.l. streamwise grid points including wake
c     nxtein     total number of b.l. streamwise grid points on airfoil only
c     rlin       Reynolds number based on airfoil chord length & freestream
c                velocity
c     xcbl,ycbl,dispbl,vnpbl    b.l. x-, y- coordinates, displacement
c                thickness , blowing velocity & skin friction coeff. 
c
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blprc2/ ff(my),uu(my),ww(my),vv(my)
      common /blc2so/ delf(my),delu(my),delv(my),delw(my)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blgtyy/ p1(mx),p2(mx)
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),dp(mx),gi,ddeps,om
      common /bledy1/ xtr,gamtr(mx),ffs(mx),alfas(mx),edv(my)
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx),uedls(mx)
      common /blmisc/ wcal,conv,omeg,iexceed 
      dimension       d1(mx2),vnp(mx),work(mx) 
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c the variables passed from the calling subroutine (main)
c
      dimension xcbl(nxtin),ycbl(nxtin),uebl(nxtin),dlspbl(nxtin),
     1          vnpbl(nxtin),uedispbl(nxtin),cfsbl(nxtin)
c +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cc
      character*2  st
      logical wcal,clconv,conv,lgc1,lgc2,lgc3 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  lagrange interpolation or extrapolation
c
      grang(x1,x2,x3,y1,y2,y3,x0)= (x0-x2)*(x0-x3)/(x1-x2)/(x1-x3)*y1
     1      +(x0-x1)*(x0-x3)/(x2-x1)/(x2-x3)*y2+(x0-x1)*(x0-x2)
     2      /(x3-x1)/(x3-x2)*y3
c  lagrange differentiation
c
C     diff3p_2nd(y1,y2,y3,x1,x2,x3,x4) = 
C    1                Y1*(2.0*X4-X2-X3)/(X1-X2)/(X1-X3) +
C    1                Y2*(2.0*X4-X1-X3)/(X2-X1)/(X2-X3) +
C    2                Y3*(2.0*X4-X1-X2)/(X3-X1)/(X3-X2)

c
c  convert input variables (arguments) to variables used in the ibl cal.
c
      nxt  = nxtin
      nxte = nxtein
      npt  = 101
      rl   = rlin
      iswpt = iswptin
      etae = 8.0
      vgp  = 1.12
      deta(1) = 0.015
      xtr  = 0.995
      xtr  = amin1(xtr,xbl(nxte-5))
      icycle = icyc
      rmi  = fmach
      do i = 1,nxt
         xbl(i) = xcbl(i)
         ybl(i) = ycbl(i)
         ubl(i) = uebl(i)
         dls(i) = dlspbl(i)
         uedls(i) = uedispbl(i)
      enddo
c
      call input(irestart,iang) 
      call ivpl
      call hic(nxt,rl,x(1))
c 
      kqwik = 0 
      ktecp = 0
      ddeps = 0.0025 
      iswp  = 0
      itv   = 1
      nx    = 0
      ntr   = nxte + 1
      conv  = .false.
      omeg  = 1.0
   50 iswp  = iswp + 1
      if(iswp .ge. iswpt) conv = .true.
      wcal  = .false.
      iexceed = 0
  100 nx    = nx + 1
      rx    = uevs(nx)*x(nx)*rl
      sqrx  = sqrt(rx)
      if(nx .gt. 1) then
         do j=1,npt
            ff(j)  = f(j,1)
            uu(j)  = u(j,1)
            ww(j)  = w(j,1)
            vv(j)  = v(j,1)
            f(j,1) = f(j,2)
            u(j,1) = u(j,2)
            v(j,1) = v(j,2)
            w(j,1) = w(j,2)
            b(j,1) = b(j,2)
         enddo
      endif
c
      if (itv .eq. 0) then
         sum2     = 0.
         do k = nx+1,nxt+1
            sum2     = sum2 + cc(k,nx)*(dd(k)-dp(k))
         enddo
         sum1     = 0.
         do k = 3,nx-1
           sum1     = sum1 + cc(k,nx)*(dd(k)-dp(k)) 
         enddo
         gi       = ubl(nx) + sum1 + sum2 - cc(nx,nx)*dp(nx) 
      end if
c
      it    = 0
      igrow  = 0
  120 it    = it+1
      if(it .ge. 20) then
         write(12,9999) alpha,icycle,isurf,iswp,nx
         iexceed = iexceed 
         if(iexceed .ge. 5) then
            istop = 1 
            return
         else
            go to 999
         endif
      endif
c
      if(nx .ge. ntr) then
         call eddy 
      else
         if(nx .gt. 3) then
            isep = 0
            if(v(1,2) .lt. 0.0) isep = 1
            call trns(isep,itran)
            if(itran .eq.1) then
               if(isep .eq. 1) goto 100
               igrow = 0
               it = 0
               goto 120
            endif
         endif
      endif
c
      iback    = 0
      if(itv .eq. 0 ) then
         iback = 1
         if(iswp .gt. 2) then
            iback = 2
            lgc1  = .false.
            if(nx .le. nxs + 3) lgc1 = .true.
            lgc2  = .false.
            if(wcal) then
               if(nx .lt. nxte+3) lgc2  = .true.
               imin = ismin(np,u(1,2),1)
               umin = u(imin,2)
c              if(umin .lt. -0.05) lgc2 = .true.
c              if(umin .lt. -0.0 ) lgc2 = .true.
            endif
            lgc3  = .false.
            if(nx .eq. ntr .or. nx.eq.ntr+1 ) lgc3 = .true.
            if(lgc1 .or. lgc2 .or. lgc3) iback = 1
         endif
      endif
c  
      itmin   = max0(iback,1)
      if(wcal) then
         nadd  = 4
         itmin = max0(itmin,(nxte+nadd-nx))
      endif

      call coef4(iback) 
      call solv4(wcal,istop) 
      if(istop .eq. 1) return
c
      if (v (1,2) .le. 0.0 .and. itv .eq. 1) then
         write(12,201) nx,x(nx),ubl(nx),p1(nx),p2(nx),v(1,2)
         write(12,202) (i,x(i),xbl(i),ybl(i),p1(i),p2(i),i=1,nxt)
202      format('i,x,xbl,ybl,p1,p2'/(i4,5e13.5))
         istop = 1
         return
      endif
c
      if(it .le. itmin) goto 120
      if(nx .le. nxte) then
c  convergence of the solutions - b. l. 
         if(nx .lt. ntr ) then
            vweps = 0.0001
            if(v(1,2) .lt. 0.0) vweps = 0.001
            if(abs(delv(1)).gt. vweps) goto 120
         else
            v1  = amax1(1.0, v(1,2))
            if(abs(delv(1)/v1) .gt. 0.01) goto 120
         end if
      else
c  convergence of the solutions - wake
         uweps  = 0.0025
cc       if(u(1,2).lt.-0.005) uweps = amin1(0.005,0.0025-0.02*u(1,2))
         if(u(1,2).lt.-0.005) uweps = amax1(0.0025,-0.025*u(1,2))
         if(abs(delu(1)).gt.uweps) goto 120
      endif
      if(np .lt. npt .and. igrow .lt. 3) then
         if( abs(v(np,2)) .gt. 0.0010) then
            igrow = igrow + 1
c  boundary layer growth
            call growth(1)
            it    = 1
            go to 120
         endif
      endif
  999 continue
      call output
      if(nx .lt. nxs ) go to 100
      if (nx .eq. nxs) then
        itv     = 0
        call swtch(0)
        go to 100
      end if
c
c  logic to speed up the convergence
c
      ddb   = 0.020 * ( dd(nx) - db(nx))
      do i = nx+1,nxt
         ddfact  = 1.0
         if(i.gt.nxte) ddfact = x(nxte) / x(i)
         dd(i) = db(i) + ddb * ddfact
      enddo
c
      if(nx .lt. nxt  ) then
         if(nx.eq.nxte) then
             call wakepr
             wcal = .true.
         endif
         goto 100
      endif
c
      call amean(nxs,nxt,x,uevs,d1,1)
      if (iswp .lt. iswpt) then
c
c  check for convergence
c
         dmax = dd(nxs)
         ddmax= abs(dd(nxs) - db(nxs))
         do i = nxs+1,nxt
            dmax = amax1(dmax,dd(i))
            ddmax= amax1(ddmax,abs(dd(i)-db(i)))
         enddo
         if(ddmax/dmax.lt.ddeps) conv = .true.
         if(conv) goto 150 
c
         if(iswp .gt. 1) then
c
c  over-relax d values to speed-up convergence
c
            do i= nxs,nxt
               dd(i)=dd(i)*(1.0+omeg*(uep(i)/uevs(i)-1.0))
               uep(i) = uevs(i)
               db(i)= dd(i)
            enddo
c
         else 
            do i= nxs,nxt
               dd(i)=dd(i)*(1.0+omeg*(uevs(i)/ubl(i)-1.0))
               uep(i)= uevs(i)
               db(i)= dd(i)
            enddo
         endif
c
         dd(nxt+1)=grang(x(nxt-2),x(nxt-1),x(nxt),
     1            dd(nxt-2), dd(nxt-1),dd(nxt),x(nxt+1))
         nx     = nxs
         call swtch(1)
         goto  50
      end if
      
c
150   call amean(nxs,nxte,x,cfs,d1,1)
c
c  calculate blowing velocity, vnp
c
      do i=1, NXT
         uedls(i) = uevs(i)*dls(i)
      end do
c
c  calculate form drag
c
c     dlsm = 0.5 * (dls(nxt)+dls(nxt-1))
c     them = 0.5 * (theta(nxt)+theta(nxt-1))
c     hh   = dlsm / them
c     uevsm= 0.5* (uevs(nxt) + uevs(nxt-1))
c     cdf1 = 2.0 * them*uevsm**((hh+5.)*0.5) 
c
c     i    = nxt-1
c     dlsm = 0.5 * (dls(i)+dls(i-1))
c     them = 0.5 * (theta(i)+theta(i-1))
c     hh   = dlsm / them
c     uevsm= 0.5* (uevs(i) + uevs(i-1))
c     cdf2 = 2.0 * them*uevsm**((hh+5.)*0.5)
c
c     cdf  = (cdf1 + cdf2 ) / 2.
c
c
c  calculate form drag
c
      i    = nxte
      dlsm = dls(i)
      them = theta(i)
      hh   = dlsm / them
      uevsm= uevs(i)
      cdf1 = 2.0 * them*uevsm**((hh+5.0+0.4*rmi**2)*0.5)
c
      i    = i-1
      dlsm = dls(i)
      them = theta(i)
      hh   = dlsm / them
      uevsm= uevs(i)
      cdf2 = 2.0 * them*uevsm**((hh+5.+4.*rmi**2)*0.5)
      dx1  = x(i+1) - x(i)
c
      i    = i-1
      dlsm = dls(i)
      them = theta(i)
      hh   = dlsm / them
      uevsm= uevs(i)
      cdf3 = 2.0 * them*uevsm**((hh+5.+4.*rmi**2)*0.5)
      dx2  = x(i+1) - x(i)
      cdf12= 0.5 * ( cdf1 + cdf2)
      cdf23= 0.5 * ( cdf3 + cdf2)
      cdf  = (cdf12 * dx2 + cdf23*dx1)/(dx1 + dx2)

c
c  option to calculate vnp
c
      call amean(1,nxt,x,uedls,work,1)
      call diff1(nxt, x, uedls, vnp, work)
      call amean(1,nxte,x,vnp,work,1)
c
c  write summary of b.l. solutions
c
      if(clconv) then
         write(12,9140) alpha,icycle,isurf,xtr,ntr,nxs,iswp  
         write(12,9130) (i,x(i),xbl(i),ubl(i),uevs(i),cfs(i),uc(i),
     1      theta(i),dls(i),dd(i),dp(i),VNP(i),npsav(i),itsav(i),
     1      i=2,nxt)
      endif
c
      do i = 1,nxt
         dlspbl(i) = dls(i)
         vnpbl(i)  = vnp(i)
         uedispbl(i) = uedls(i)
         cfsbl(i)    = cfs(i)
      enddo
c
c  write data into qwikplot,tecpplot  
c
      if(kqwik .eq. 0) return
      iqwik = 88
      itecp = 87
      open(88,file='qwik_file',form='unformatted')
      open(87,file='tecp_file',form='formatted')
      rewind iqwik
      rewind itecp
      write(st,'(i2)') isurf
      write(iqwik) 'surf= ',st,0.0 
      write(iqwik) 'x       ',nxt,(x(i),     i=1, nxt)
      write(iqwik) 'xbl     ',nxt,(xbl(i),   i=1, nxt)
      write(iqwik) 'ubl     ',nxt,(ubl(i),   i=1, nxt)
      write(iqwik) 'theta   ',nxt,(theta(i), i=1, nxt)
      write(iqwik) 'dls     ',nxt,(dls(i),   i=1, nxt)
      write(iqwik) 'uevs    ',nxt,(uevs(i),  i=1, nxt)
c
      write(iqwik) 'xbl_body',nxte,(xbl(i),  i=1, nxte)
      write(iqwik) 'cfs     ',nxte,(cfs(i),  i=1, nxte)
      write(iqwik) 'vw      ',nxte,(vw(i),   i=1, nxte)
      if(nxt .gt. nxte) then
         ii = nxt - nxte
         write(iqwik) 'xbl_wake',ii,(xbl(i), i=nxte+1, nxt)
         write(iqwik) 'uc      ',ii,(uc(i),  i=nxte+1, nxt)
      endif
	close(iqwik)
c
c  tecplot file
c
      write(itecp,99)
99    format('title= "2d inverse b. l. summary"')
      write(itecp,110)
110   format('variables="s","x/c","uei","delst","theta","vw"',
     1       '"cf","uev" ')
      write(itecp,121) isurf, nxte 
121   format('zone','isurf =',i3,' i=',i3,' F = point')
      do i = 1,nxte
         write(itecp,122) x(i),xbl(i),ubl(i),dls(i),theta(i),
     1            vw(i),cfs(i),uevs(i)
122      format(8f10.6)
      enddo
      close(itecp)
c
      return  
c
  201 format (' flow separates in the standard mode; print data',
     1        ' for program check '/1h ,2x,'nx = ',i3,2x,'x = ',
     2        f10.5,2x,'ubl = ',f10.5,2x,'p1 = ',f10.5,2x,'p2 = ',
     3        f10.5,2x,'vw = ',f10.5,/)
 9130 format (3h nx,6x,2h x,7x,2hxc,5x,4huei ,4x,4huev ,8x,
     &    3h cf,8x,3huc ,6x,5htheta,8x,3hdls,8x,3h dd,8x,3h dp
     &    8x,3h vw,4x,3h np,1x,2hit/(i3,2f10.6,2f8.4,7e11.3,i5,i3))
 9140 format (/,'b. l. solutions for a cyc surf xtr ntr nxs iswp = ',
     &        f8.3,2i3,f8.5,3i5/)
 9100 format(//'** summary of boundary layer solutions '/
     &        'nxs =',i3,4x,'** transition location : ntr =  ',i3)
 9999 format('iter. exceeds max at a,cyc,surf,swp,nx =',f10.5,4i5/)
      end

      subroutine input(irestart,iang) 
      parameter(mx=201,my=151,mx2=mx*2)
      common /blgtyy/ p1(mx),p2(mx)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),dp(mx),gi,ddeps,om
      common /blinp3/ rl, rx, sqrl, sqrx
      common /bledy1/ xtr,gamtr(mx),ffs(mx),alfas(mx),edv(my)
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx),uedls(mx)
      dimension d1(mx)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      grang(x1,x2,x3,y1,y2,y3,x0)= (x0-x2)*(x0-x3)/(x1-x2)/(x1-x3)*y1
     1      +(x0-x1)*(x0-x3)/(x2-x1)/(x2-x3)*y2+(x0-x1)*(x0-x2)
     2      /(x3-x1)/(x3-x2)*y3
c
c     read (5,*  )
c     read (5,*  ) nxt,nxte,nxs,npt,iswpt 
c     read (5,*  )
c     read (5,*  ) rl,xtr,etae,vgp,deta(1), p2(1)
c     read (5,*  )
c     read (5,*  ) (xbl(i),i=1,nxt)
c     read (5,*  )
c     read (5,*  ) (ybl(i),i=1,nxt)
c     read (5,*  )
c     read (5,*  ) (ubl(i),i=1,nxt)
c
c  generate b.l. grid in the normal direction
c
      if ((vgp-1.0) .gt. 0.001) then    
        np = alog((etae/deta(1))*(vgp-1.0)+1.0)/alog(vgp)+1.0001
      else    
        np = etae/deta(1)+1.0001
      end if
      eta(1) = 0.0
      do j=2,npt
        deta(j)  = vgp*deta(j-1)
        eta(j)   = eta(j-1)+deta(j-1)
        a(j)     = 0.5*deta(j-1)
      enddo
c
      x(1) = 0.0
      do i = 2,nxt
         x(i) = x(i-1)+sqrt((xbl(i)-xbl(i-1))**2+(ybl(i)-ybl(i-1))**2)
      enddo
c
      umx  = -999.0
      do i = 2,nxte
         if(ubl(i) .gt. umx) then
            imx = i
            umx = ubl(i)
         endif
      enddo
c     call amean(imx,nxte,x,ubl,d1,1)
      call amean(1  ,nxte,x,ubl,d1,1)
      if(nxt .gt. nxte) then
         call amean(nxte-3,nxte+3,x,ubl,d1,2)
      endif
c    
      do i=1,nxt
         uevs(i) = ubl(i)
         uep(i)  = ubl(i)
      enddo
c  initial deltastar & dd distribution
      sqrl = sqrt(rl)
      if(irestart.eq.1) then
c
c  dls & dd are initially set to the values of the previous cycle
c
         do i = 1,nxt
            dp(i) = uedls(i)*sqrl 
            db(i) = dp(i)
            dd(i) = dp(i)
         enddo
      else
c
c  initial guess of dd & dls
c
         dd(1) = 0.0
         dls(1)= 0.0
         db(1) = 0.0
         dp(1) = 0.0
cc	 if(iang .gt. 1) then 
cc	    do i = 1,nxt
cc  	       db(i) = uedls(i)*sqrl
cc 	       dd(i) = db(i)
c	       dp(i) = 0.0
c           enddo
c        else
            dd(1) = 0.0
	    dls(1)= 0.0
	    db(1) = 0.0
	    dp(1) = 0.0
            do i = 2,nxte
               rx     = rl*ubl(i)*x(i)
               dls(i) = 0.046875*x(i)/rx**0.2
               dd(i)  = ubl(i)*dls(i)*sqrl
               db(i)  = dd(i)
               dp(i)  = 0.0
            enddo
c        endif
c
         do i=nxte+1,nxt
            dd(i) = dd(nxte)/exp( (x(i)-x(nxte))/(x(nxt)-x(nxte)))
            db(i) = dd(i)
            dp(i) = 0.0
         enddo
      endif
c
      x(nxt+1)  = 3.0 * x(nxt) -3.0 * x(nxt-1) + x(nxt-2)
      dd(nxt+1) =grang(x(nxt-2),x(nxt-1),x(nxt),
     1              dd(nxt-2),dd(nxt-1),dd(nxt),x(nxt+1))
      dp(nxt+1) =grang(x(nxt-2),x(nxt-1),x(nxt),
     1              dp(nxt-2),dp(nxt-1),dp(nxt),x(nxt+1))

      call diff1(nxt,x,ubl,p2,p1) 
      p2(1) = 1.0
      p1(1) = 0.5*(p2(1) + 1.0)
      do i=2,nxt 
         p2(i) = x(i)/ubl(i)*p2(i)
         p1(i) = 0.5*(1.0+p2(i))
      enddo
c
c  search for nxs, location where calculation switch from the standard to
c  inverse mode
c
      do i = 6,nxte
         ii  = i
         if(ubl(i+1) .lt. ubl(i)) goto 10 
      enddo
10    nxs  = min0(12,ii-4) 
       
      return                                                            
c5000 format(a)
c6040 format(/' etae=',f8.3,3x,'vgp=',f7.3,3x,'deta(1)=',f7.3/)
      end                                                               

      subroutine ivpl
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      etanpq= 0.25*eta(np)
      etau15= 1.5/eta(np)
      do 30 j=1,np
         etab  = eta(j)/eta(np)
         etab2 = etab**2
         f(j,2)= etanpq*etab2*(3.0-0.5*etab2)
         u(j,2)= 0.5*etab*(3.0-etab2)
         w(j,2)= 1.0
         v(j,2)= etau15*(1.0-etab2)
         b(j,2) = 1.0
   30 continue
      return
      end

      subroutine hic(nxt,rl,x)
      parameter(mx=201,my=151,mx2=mx*2)
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),dp(mx),gi,ddeps,om
      dimension       x(*),e(mx)
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      pi      = 4. * atan(1.0 )
      pisqrl  = pi * sqrt( rl )
      do 65 i = 2, nxt
         e (1)   = 0.
         do 60 j= 2, i-1
            e(j)= alog(abs((x(i)-x(j-1))/(x(i)-x(j))))
     1            /(x(j)-x(j-1))
   60    continue
         j      = i
         r1     = (x(j+1) -x(j)) /(x(j+1) - x(j-1))
         e(j)   = (r1*alog(abs((x(i)-x(j-1))/(x(i)-x(j+1))))+2.0)
     1            /(x(j) - x(j-1))
         j      = i + 1
         r1     = (x(j-1) - x(j-2)) / (x(j) - x(j-2))
         e(j  ) = (r1*alog(abs((x(i)-x(j-2))/(x(i)-x(j))))-2.0)
     1            /(x(j) - x(j-1))
         do 61 j= i+2,nxt+1
            e(j)= alog(abs((x(i)-x(j-1))/(x(i)-x(j))))
     1            /(x(j) - x(j-1))
   61    continue
         do 62 j =2,nxt+1
            cc(j-1,i) = ( e(j-1) - e(j) ) / pisqrl
   62    continue
         cc(nxt+1,i) = e(nxt+1) / pisqrl
   65 continue
      return
      end
c
      subroutine eddy
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /bledy1/ xtr,gamtr(mx),ffs(mx),alfas(mx),edv(my)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blcwak/ edsave(my),ytedge(my),wscale,ydelwk
      dimension       fint(my),yy(my),edintp(my),du(my)
c   --------------------------------------------------------
      if(itv .ne. 0) then
         sqrey = sqrx
      else
         sqxc  = sqrt(x(nx))
         sqrey = sqxc * sqrl
      endif
c
      iwk = 1
      if(nx .le. nxte) iwk = 0
c
c  calculate fint
      gamtrnx = gamtr(nx)
      call gamcal(iwk,np,gamtrnx,eta,u(1,2),fint,du)
c
      if(iwk .eq. 1) goto 1000
c
c  calculate alfa
c
      if(nx .eq. ntr) then
         alfas(nx) = 0.0168
         ffs(nx)   = 1.0
      else
         if(it .eq. 1) then
            alfas(nx) = alfas(nx-1)
            ffs(nx)   = ffs(nx-1)
         else
            call calfa
         endif
      endif
c
c  calculate equivalent 'v(1,2)', um, in the damping term
c  this is necessary because the original pplus in the CS model is
c  not applicable to the flow with seapration
c
      jm       = 1
      swall    = abs(v(1,2))
      ttm      = -9999.0
      if(nx.gt.ntr .or. it.gt.1) then
         wsrex= 0.75
         do j = 1,np
            tt = v(j,2)*(wsrex * edv(j) +(1.-wsrex))
            if(tt.gt. ttm) then
               ttm = tt
               jm  = j    
            endif
         enddo
      endif
      smax     = amax1(swall, ttm)
      um       = sqrt(smax)
c
      rk    = 0.40
      cya   = sqrt(sqrey)/26.0
      edviic= rk**2 * sqrey * gamtr(nx)
      cedv  = sqrey*gamtr(nx)*u(np,2)
      edvoo = alfas(nx)*(eta(np)-f(np,2)/u(np,2))*cedv
      do j = 1,np
         edv(j)= edvoo
      enddo
      do j = 1,np
         yoa   = amin1(50.0,cya*um*eta(j))
         el    = 1.0 - exp(-yoa)
         edvin = edviic*abs(v(j,2))*(eta(j)*el)**2
         if(edvin .gt. edvoo) goto 190
         edv(j)= edvin
      enddo
c
  190 do j=1,np
         edv(j) = edv(j) * fint(j)
         b(j,2) = 1.0+edv(j)
      enddo
      return
c
c  wake region
c
 1000 continue
c
      edvoo  = 0.064 * (eta(np)-f(np,2)/u(np,2)) * sqrey * u(np,2)
      cexpt  = amin1(20.,(x(nx)-x(nxte))/(wscale*ydelwk))
      expt   = 1.0 - exp(-cexpt)
      yy(1) = 0.0
      yoeta = x(nx) / sqrey
      do j = 2, np
         yy(j)  = yoeta * eta(j)
         if(yy(j) .lt. ytedge(npt)) nn = np
      enddo
      call lntp(npt,ytedge,edsave,nn,yy,edintp)
      if(nn .lt. np) then
         do j = nn+1,np
            edintp(j) = edintp(nn)
         enddo
      endif
c
      do j = 1,np
         edvwj = edvoo * fint(j)
         edv(j)= amax1(0.0,edintp(j)+(edvwj-edintp(j)) * expt)
         b(j,2)= 1.0+edv(j)
      enddo
      return
      end

      subroutine coef4(iback) 
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blprc2/ ff(my),uu(my),ww(my),vv(my)
      common /blcs30/ s1(my),s2(my),s3(my),s4(my),s5(my),s6(my),
     1                s7(my),s8(my),r1(my),r2(my),r3(my),r4(my)
      common /blcs31/ gamma1,gamma2
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blgtyy/ p1(mx),p2(mx)
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),dp(mx),gi,ddeps,om
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(itv .eq. 0) goto 200
c
c  standard method - use central difference method
c
      cel = 0.0                                                         
      if(nx .gt. 1) cel = 0.5*(x(nx)+x(nx-1))/(x(nx)-x(nx-1))           
      l1  = itv
      l2  = 1 - itv
      p1p = l1*p1(nx) + cel + 0.5*l2                                       
c
      do 100 j=2,np
c  present station
        usb   = 0.5*(u(j,2)**2+u(j-1,2)**2)
        wsb   = 0.5*(w(j,2)**2+w(j-1,2)**2)
        ub    = 0.5*(u(j,2)+u(j-1,2))
        vb    = 0.5*(v(j,2)+v(j-1,2))
        wb    = 0.5*(w(j,2)+w(j-1,2))
        fb    = 0.5*(f(j,2)+f(j-1,2))
        fvb   = 0.5*(f(j,2)*v(j,2)+f(j-1,2)*v(j-1,2))
        derbv = (b(j,2)*v(j,2)-b(j-1,2)*v(j-1,2))/deta(j-1)
        flare = 1.0
        if(ub .lt. 0.0) flare = 0.0
        p2p = l1*p2(nx) + cel*flare
        if(nx.eq.1) then
          crb = -p2(nx)
          cfb = 0.0
          cvb = 0.0
        else
c   previous station
          cfb   = 0.5*(f(j,1)+f(j-1,1))
          cvb   = 0.5*(v(j,1)+v(j-1,1))
          cusb  = 0.5*(u(j,1)**2+u(j-1,1)**2)
          cwsb  = 0.5*(w(j,1)**2+w(j-1,1)**2)
          cfvb  = 0.5*(f(j,1)*v(j,1)+f(j-1,1)*v(j-1,1))
          cderbv= (b(j,1)*v(j,1)-b(j-1,1)*v(j-1,1))/deta(j-1)
          clb = cderbv + l1*(p1(nx-1)*cfvb + p2(nx-1)*(1.0-cusb)) +
     &         0.5*l2*cfvb 
          crb = -clb - l1*p2(nx) - cel*flare*cusb + cel*cfvb
     &         + cel*l2*cwsb
        endif
c  definitions of sj
        s1(j) = b(j,2)/deta(j-1) + 0.5*p1p*f(j,2) -0.5*cel*cfb
        s2(j) = -b(j-1,2)/deta(j-1)+0.5*p1p*f(j-1,2)-0.5*cel*cfb      
        s3(j) = 0.5*(p1p*v(j,2) + cel*cvb)                                
        s4(j) = 0.5*(p1p*v(j-1,2) + cel*cvb)                              
        s5(j) = -p2p*u(j,2)
        s6(j) = -p2p*u(j-1,2)
        s7(j) = l2*cel*wb
        s8(j) = l2*cel*wb
c  definitions of rj
        r1(j)  = f(j-1,2) - f(j,2) + ub*deta(j-1)
        r2(j) = crb - (derbv + p1p*fvb - p2p*usb) 
     &         -cel*(fb*cvb-vb*cfb + l2*wsb)
        r3(j-1)= u(j-1,2) - u(j,2) + vb*deta(j-1)
        r4(j-1)= w(j-1,2) - w(j,2) 
  100 continue
      goto 300
  200 continue
c
c  inverse method
c
      if(iback .eq. 2) then
c
c  3-p backward difference
c
         cx1  = x(nx)*(x(nx)-x(nx-1))/(x(nx-2)-x(nx-1))/
     1          (x(nx-2)-x(nx))
         cx2  = x(nx)*(x(nx)-x(nx-2))/(x(nx-1)-x(nx-2))/
     1          (x(nx-1)-x(nx))
         cx3  = x(nx)*(2.*x(nx)-x(nx-1)-x(nx-2))/(x(nx)-x(nx-2))/
     1          (x(nx)-x(nx-1))
      else
c
c  2-p backward difference
         cx1  = 0.0
         cx2  = - x(nx)/(x(nx)-x(nx-1))
         cx3  = -cx2
      endif
c
      p1(nx) = 0.5 
      p1h    = 0.5 * p1(nx)
      do j=2,np
c  present station
        usb   = 0.5*(u(j,2)**2+u(j-1,2)**2)
        wsb   = 0.5*(w(j,2)**2+w(j-1,2)**2)
        ub    = 0.5*(u(j,2)+u(j-1,2))
        vb    = 0.5*(v(j,2)+v(j-1,2))
        wb    = 0.5*(w(j,2)+w(j-1,2))
        fb    = 0.5*(f(j,2)+f(j-1,2))
        fvb   = 0.5*(f(j,2)*v(j,2)+f(j-1,2)*v(j-1,2))
        derbv = (b(j,2)*v(j,2)-b(j-1,2)*v(j-1,2))/deta(j-1)
        flare = 1.0
        if(ub .lt. 0.0) flare = 0.0
c
c   previous station
c
        cfb   = 0.5*(f(j,1)+f(j-1,1))
        cusb  = 0.5*(u(j,1)**2+u(j-1,1)**2)
        cwsb  = 0.5*(w(j,1)**2+w(j-1,1)**2)
c
c  at nx - 2
c
        ffb   = 0.5*(ff(j)  + ff(j-1))
        uusb  = 0.5*(uu(j)**2  + uu(j-1)**2)
        wwsb  = 0.5*(ww(j)**2  + ww(j-1)**2)
c
c  defination of derivatives
c
        xdfdx = cx1*ffb + cx2*cfb + cx3 * fb
        vxdfdx= 0.5*cx3
        xdusdx= cx1 * uusb + cx2 * cusb + cx3 * usb
        vxdusdxj= cx3*u(j,2)
        vxdusdxj1= cx3*u(j-1,2)
        xdwsdx= cx1 * wwsb + cx2 * cwsb + cx3 * wsb
        vxdwsdxj= cx3*w(j,2)
        vxdwsdxj1= cx3*w(j-1,2)
c
c  definitions of sj
        s1(j) = b(j,2)/deta(j-1) + p1h*f(j,2) +0.5*xdfdx
        s2(j) =-b(j-1,2)/deta(j-1)+p1h*f(j-1,2)+0.5*xdfdx
        s3(j) = p1h*v(j,2)   + vb * vxdfdx
        s4(j) = p1h*v(j-1,2) + vb * vxdfdx
        s5(j) = -0.5 * vxdusdxj * flare 
        s6(j) = -0.5 * vxdusdxj1* flare
        s7(j) =  0.5 * vxdwsdxj
        s8(j) =  0.5 * vxdwsdxj1 
c  definitions of rj
        r1(j)  = f(j-1,2) - f(j,2) + ub*deta(j-1)
        r2(j)  = -(derbv+p1(nx)*fvb+(0.5*xdwsdx-0.5*xdusdx*
     +             flare+vb*xdfdx))
        r3(j-1)= u(j-1,2) - u(j,2) + vb*deta(j-1)
        r4(j-1)= w(j-1,2) - w(j,2)
      enddo  
c
  300 continue
c  boundary conditions
      r1(1)  = 0.0
      r2(1)  = 0.0
      r4(np) = w(np,2)-u(np,2)
      if (itv .eq. 1 ) then
c  standard mode
         gamma1 = 0.0
         gamma2 = 1.0
         r3(np) = 0.0
      else
c  inverse mode
         gamma1 = cc(nx,nx) * sqrt(x(nx))
         gamma2 = 1.0 - eta(np) * gamma1
         r3(np) = gi - (gamma2*w(np,2) +gamma1*f(np,2))
c
c  residue relaxtion at the first wake station
cc
         if(nx.eq.nxte+1) then
            ism  = ismin(1,u(1,2),np)
            umin = u(ism,2)
            if(umin .lt. -0.005) then 
               reit = amin1(1.0,0.25*it) 
            else
               reit = amin1(1.0,0.50*it)
            endif
c
            do j = 1,np
               r1(j) = reit*r1(j)
               r2(j) = reit*r2(j)
               r3(j) = reit*r3(j)
               r4(j) = reit*r4(j)
            enddo
         endif
      endif
      return
      end
c
      subroutine swtch(index)
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx),uedls(mx)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blgtyy/ p1(mx),p2(mx)
      common /blsave/ fsave(my),usave(my),vsave(my),bsave(my),
     1                wsave(my),npsave
c    -----------------------------------------------------------
      if(index.gt.0) goto 40
      wnp   = ubl(nx)
      ak    = sqrt(1.0/wnp)
      do j=1,npt 
         u(j,2) = wnp*u(j,2)
         v(j,2) = wnp*v(j,2)/ak
         f(j,2) = wnp*f(j,2)*ak
         fsave(j)  = f(j,2)
         usave(j)  = u(j,2)
         vsave(j)  = v(j,2)
         bsave(j)  = b(j,2)
      enddo
      do j = 1,npt 
         w(j,2)   = u(np,2)
         wsave(j) = w(j,2)
      enddo
      do 30 j= 1,npt
         eta(j) = ak*eta(j)
   30 continue
      do 35 j=2,npt
         deta(j-1)= eta(j)-eta(j-1)
         a(j)     = 0.5*deta(j-1)
   35 continue
      npsave = np
      do i = nxs+1,nxt
         p2(i) = 0.0
         p1(i) = 0.5
      enddo
      return
   40 np     = npsave
      do j=1,npt 
         f(j,2) = fsave(j)
         u(j,2) = usave(j)
         w(j,2) = wsave(j)
         v(j,2) = vsave(j)
         b(j,2) = bsave(j)
      enddo
      return
      end

      subroutine wakepr
      parameter(mx=201,my=151,mx2=mx*2,my2=my*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /bledy1/ xtr,gamtr(mx),ffs(mx),alfas(mx),edv(my)
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blinp8/ d1(mx2),d2(mx2),d3(mx2),d4(mx2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blcwak/ edsave(my),ytedge(my),wscale,ydelwk
      dimension       y(my2),uu(my2),vv(my2)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (nx.ge.ntr) then
        do 20 j=1,np
          edv(j) = edv(j) /gamtr(nx)
   20   continue
c     else
c       gamtr(nx) = 1.0
c       ntr       = nx
c       it        = 1
c       call eddy
c       do 30 j =1,np
c         edsave(j) = edv(j)
c         b(j,2) = 1.0  
c  30   continue
c       gamtr(nx) = 0.0
c       ntr       = nx + 1
      endif
      do 40 j=1,np
         jj    = np-1+j
         uu(jj)= u (j,2)
         y(jj) = eta(j)
   40 continue
      do 50 j = 2,np
         jj     = np-j+1
         uu(jj) = u (j,2)
         y(jj)  = -eta(j)
   50 continue
      jn     = 2*np-1
      call amean(np-3,np+3,y,uu,d1,2)
      call amean(np-5,np+5,y,uu,d1,1)
      call diff1(jn,y,uu,vv,d1)
      do 60 jj= np,jn
         j       = jj - np + 1
         u (j ,2)= uu(jj)
         v (j ,2)= vv(jj)
         edsave(j)= edv(j)
   60 continue
      edsanp  = edsave(np)
      do 65 j = np,npt
         edsave(j)= edsanp
   65 continue
      v(1,2)  = 0.0
      call integ(np ,eta, u(1,2), v(1,2), f(1,2))
      ytedge(1 )  = 0.0
      cy     = (1.0/sqrl) *  sqrt(x(nx) )
      udel   = 0.995 * u(np,2)
      do 230 j = 2,np
         ytedge(j)  = cy * eta(j)
         if ( u(j,2) .lt. udel ) jj  = j
  230 continue
      do 240 j = np+1,npt
         ytedge(j) = cy * eta(j)
  240 continue
      ydelwk = ytedge(jj)+(ytedge(jj+1)-ytedge(jj))/(u(jj+1,2)-u(jj,2))*
     1         (udel-u(jj,2))
      wscale  = 20.0
      return
      end

      subroutine integ (n,x,f,fp,h)
      dimension x(n),f(n),fp(n),h(n)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      h(1)=0.0
      do 10 i=2,n
         i1  = i-1
         dx  = x(i)-x(i1)
         df  = f(i)-f(i1)
         dfp = fp(i)-fp(i1)
         fpf =fp(i)+fp(i1)
         dx2 =dx*dx/2.0
         dx3 =dx2*dx/3.0
         dx4 =dx3*dx/4.0
         fppp=(-2.*df+dx*fpf)/dx3
         fpp =dfp/dx - dx*fppp/2.0
         h(i)=h(i1) + f(i1)*dx + fp(i1)*dx2 + fpp*dx3 + fppp*dx4
   10 continue
      return
      end

      subroutine diff1 (n,x,f,fp,fpp)
      dimension x(n),f(n),fp(n),fpp(n)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do 5 i = 1,n-1
         dx    = x(i+1)-x(i)
         df    = f(i+1)-f(i)
         fpp(i)= atan2(df,dx)
   5  continue
      do 10 i  = 2,n-1
         dx1   = x(i) - x(i-1)
         dx2   = x(i+1) - x(i)
         ang   = (dx2*fpp(i-1)+dx1*fpp(i))/(dx1+ dx2)
         fp(i)  = tan(ang)
   10 continue
      fp(1)   = 2.*tan(fpp(1)) - fp(2)
      fp(n)   = 2.*tan(fpp(n-1))- fp(n-1)
      return
      end

      function diff3p_2nd (Y1,Y2,Y3,X1,X2,X3,X4)
C  diff3p_2nd is a three-point (x1, x2, x3) and 2nd order
C  difference formula  used to calculate the derivative of y with respect 
C  to x at x=x4.
         diff3p_2nd = Y1*(2.0*X4-X2-X3)/(X1-X2)/(X1-X3) +
     1                Y2*(2.0*X4-X1-X3)/(X2-X1)/(X2-X3) +
     2                Y3*(2.0*X4-X1-X2)/(X3-X1)/(X3-X2)
	return
	end

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

      subroutine amean(ns,nd,x,y,d,it)
      dimension x(nd),y(nd),d(nd)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((nd-ns).gt.2 .and. it.gt.0) then
         do 20 k=1,it
            do 10 i = ns,nd
               d(i) = y(i)
   10       continue
            do 15 i=ns+1,nd-1
               f1    = (x(i + 1) -x(i))/(x(i+1) - x(i-1))
               y(i)  = 0.5*(f1*d(i-1)+d(i)+(1.0-f1)*d(i+1))
   15       continue
   20    continue
      endif
      return
      end

      subroutine solv4(wcal,istop)                                                  
      parameter(mx=201,my=151,mx2=mx*2)
      common /blcs31/ gamma1,gamma2
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blc2so/ delf(my),delu(my),delv(my),delw(my)
      common /blcs30/ s1(my),s2(my),s3(my),s4(my),s5(my),s6(my),
     1                s7(my),s8(my),r1(my),r2(my),r3(my),r4(my)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx),uedls(mx)
      common /blmiss/ delvw(25)
      dimension     a11(my),a12(my),a13(my),a14(my),a21(my),            
     1              a22(my),a23(my),a24(my),g11(my),g12(my),            
     1              g13(my),g14(my),g21(my),g22(my),g23(my),            
     1              g24(my),w1(my),w2(my),w3(my),w4(my)                 
      logical wcal
c 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c 
      a11(1)   = 1.0                                                    
      a12(1)   = 0.0                                                    
      a13(1)   = 0.0                                                    
      a14(1)   = 0.0                                                    
      if(nx .gt.nxte) then
c
c  wake flow: v = 0.0
         a21(1)   = 0.0
         a22(1)   = 0.0
         a23(1)   = 1.0
         a24(1)   = 0.0
      else
c
c b .l. - u = 0.0
         a21(1)   = 0.0                                                    
         a22(1)   = 1.0                                                    
         a23(1)   = 0.0                                                    
         a24(1)   = 0.0                                                    
      endif
c
      w1(1)    = r1(1)                                                  
      w2(1)    = r2(1)                                                  
      w3(1)    = r3(1)                                                  
      w4(1)    = r4(1)                                                  
      do 10 j = 2,np                                                    
      aa1      = a13(j-1)-a(j)*a12(j-1)                                 
      aa2      = a23(j-1)-a(j)*a22(j-1)                                 
      aa3      = s2(j)-a(j)*s6(j)                                       
      det      = aa2*a11(j-1)-aa1*a21(j-1)                              
      ajs      = a(j)**2                                                
      g11(j)   = -(aa2+a21(j-1)*ajs)/det                                
      g12(j)   = (a11(j-1)*ajs+aa1)/det                                 
      g13(j)   = a12(j-1)*g11(j)+a22(j-1)*g12(j)+a(j)                   
      g14(j)   = a14(j-1)*g11(j)+a24(j-1)*g12(j)                        
      g21(j)   = (s4(j)*aa2-a21(j-1)*aa3)/det                           
      g22(j)   = (a11(j-1)*aa3-s4(j)*aa1)/det                           
      g23(j)   = a12(j-1)*g21(j)+a22(j-1)*g22(j)-s6(j)                  
      g24(j)   = a14(j-1)*g21(j)+a24(j-1)*g22(j)-s8(j)                  
      a11(j)   = 1.0                                                    
      a12(j)   = -a(j)-g13(j)                                           
      a13(j)   = a(j)*g13(j)                                            
      a14(j)   = -g14(j)                                                
      a21(j)   = s3(j)                                                  
      a22(j)   = s5(j)-g23(j)                                           
      a23(j)   = s1(j)+a(j)*g23(j)                                      
      a24(j)   = s7(j)-g24(j)                                           
      w1(j)    = r1(j) -g11(j)*w1(j-1)-g12(j)*w2(j-1)-w3(j-1)*g13(j)    
     1           -g14(j)*w4(j-1)                                        
      w2(j)    = r2(j) -g21(j)*w1(j-1)-g22(j)*w2(j-1)-w3(j-1)*g23(j)    
     1           -g24(j)*w4(j-1)                                        
      w3(j)    = r3(j)                                                  
      w4(j)    = r4(j)                                                  
   10 continue                                                          
      d = gamma1*(a13(np)*a24(np)-a14(np)*a23(np)-a12(np)*a23(np)+      
     1      a13(np)*a22(np)) + gamma2*(a11(np)*a23(np)-a13(np)*a21(np)) 
      df = w3(np)*(a13(np)*a24(np)-a14(np)*a23(np)-a12(np)*a23(np)+     
     1      a13(np)*a22(np)) - gamma2*(w4(np)*(a12(np)*a23(np)-a13(np)* 
     2      a22(np)) - w1(np)*a23(np)+w2(np)*a13(np))                   
      du = w4(np)*(gamma1*(a24(np)*a13(np)-a23(np)*a14(np)) + gamma2*   
     1      (a11(np)*a23(np)-a21(np)*a13(np)))+ gamma1*(w2(np)*a13(np)- 
     2      w1(np)*a23(np)) - w3(np)*(a21(np)*a13(np)-a23(np)*a11(np))  
      dv = gamma1*(w1(np)*a24(np)-w2(np)*a14(np)) + w3(np)*(a21(np)*    
     1      a14(np)-a11(np)*a24(np)) - gamma2*(w1(np)*a21(np)-w2(np)*   
     2      a11(np)) +w4(np)*(gamma1*(a22(np)*a14(np)-a24(np)*a12(np))  
     3      +gamma2*(a21(np)*a12(np)-a11(np)*a22(np)))+ gamma1*(w1(np)* 
     4      a22(np)-w2(np)*a12(np)) + w3(np)*(a21(np)*a12(np)-a11(np)*  
     5      a22(np))                                                    
      dw = gamma1*(w2(np)*a13(np)-w1(np)*a23(np)+w4(np)*(a12(np)*a23(np)
     1      -a13(np)*a22(np))) + w3(np)*(a11(np)*a23(np)-a13(np)*       
     2      a21(np))                                                    
c 
      det = d                                                           
      delf(np) = df/det                                                 
      delu(np) = du/det                                                 
      delv(np) = dv/det                                                 
      delw(np) = dw/det                                                 
c
      do j = np-1,1,-1
         cc1      = delu(j+1)-w3(j)-a(j+1)*delv(j+1)                       
         cc2      = delw(j+1)-w4(j)                                        
         cc3      = a13(j)-a(j+1)*a12(j)                                   
         cc4      = w1(j)-a12(j)*cc1-a14(j)*cc2                            
         cc5      = a23(j)-a(j+1)*a22(j)                                   
         cc6      = w2(j)-a22(j)*cc1-a24(j)*cc2                            
         deno     = a11(j)*cc5-a21(j)*cc3                                  
         delf(j)  = (cc4*cc5-cc3*cc6)/deno                                 
         delv(j)  = (a11(j)*cc6-a21(j)*cc4)/deno                           
         delw(j)  = cc2                                                    
         delu(j)  = cc1-a(j+1)*delv(j)                                     
      enddo
c
c  under-relaxation for numerical stability 
c
      rdel   = 1.0
      if(nx .gt. nxte) then
         delvw(it) = delu(1)
      else
         delvw(it) = delv(1)
      endif 
c
      imin = ismin(np,u(1,2),1)
      umin = u(imin,2)

      if(wcal) then
        if(nx .eq. nxte+1) then
           rdel = amin1(1.0,0.35+0.25*(it-1))
        else
           if(nx .le. nxte+2 .or. umin .lt. -0.005)
     1        rdel = amin1(1.0,0.5+0.25*(it-1))
        endif
      endif
c
      if(nx.eq.ntr .or. nx.eq.ntr+1 .or. nx.eq.ntr+2)
     1    rdel = amin1(1.0,0.5+0.25*(it-1))
c
      if(it.gt.5) then
         prod    = delvw(it) * delvw(it-1)
         if(prod.lt.0.0 ) rdel=amin1(rdel,0.75)
      endif
c
c     if(umin.lt.-0.005) rdel = amin1(rdel,0.85)
      do 30 j = 1,np                                                    
         f(j,2) = f(j,2)+rdel*delf(j)                                       
         u(j,2) = u(j,2)+rdel*delu(j)                                       
         v(j,2) = v(j,2)+rdel*delv(j)                                       
         w(j,2) = w(j,2)+rdel*delw(j)                                       
   30 continue                                                          
      if(nx .le.nxte) then
c
c  b. l.
         u(1,2) = 0.0                                                    
      else
c
c  wake
         v(1,2) = 0.0
      endif
      if(itv .eq. 0) uevs(nx) = u(np,2)
      istop = 0
      if(nx .ge. ntr) call edgchk(umin,istop) 
c
      return                                                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      end                                                               
c
      subroutine edgchk(umin,istop) 
c
c  check for edge balance & smooth the velocity profiles if necessary
c
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blgtyy/ p1(mx),p2(mx)

      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx),
     1                uc(mx),uep(mx),uedls(mx)
      dimension d1(my),d2(my)
c
c    check for maximum v
c
      iope = 1
      iops = 0
      istop= 0
      vmax = v(1,2)
      jp   = 1 
      do j = 1,np
         if(v(j,2) .gt. vmax) then
            vmax = v(j,2)
            jp = j
         endif
      enddo
c
c   check for smoothness
c
      do j = jp+1,np-1
         jj   = j
         if(u(j,2) .gt. u(np,2) .or. v(j,2) .lt. 0.0) go to 40
      enddo
      return
c
   40 j1     = jj-1
c
c  smooth velocity profiles near the edge of b. l. to eliminate
c  unnecessary growth because of the balance there 
c
      if(j1 .le. 1) then
         write(12,*) 'print in sub. edgchk'
         write(12,*) 'solutions diverge. cal. stop at nx=',nx 
         write(12,*) 'j,eta,f(j,2),u(j,2),v(j,2),w(j,2),b(j,2),j=1,np'
         do j = 1,np
            write(12,10) j,eta(j),f(j,2),u(j,2),v(j,2),w(j,2),b(j,2)
10          format(i5,6e13.5)
         enddo
         write(12,202) (i,x(i),xbl(i),ybl(i),ubl(i),uevs(i),
     1                  p1(i),p2(i),i=1,nxt)
202      format('i,x,xbl,ybl,ubl,uev,p1,p2'/(i4,7e13.5))
         istop = 1
         return
      endif
c
      call amean(j1,np,eta,f(1,2),d1,1)
      call amean(j1,np,eta,u(1,2),d1,1)
      if(iope .eq. 1) then
c        do j = 1,np
c           d2(j) = v(j,2)
c        enddo
c        call amean(j1,np,eta,d2,d1,1)
c 
         djp    = deta(j1-1)
         vjp    = (u(j1,2)-u(j1-1,2))/djp
         do j=j1+1,np
            djm      = deta(j-1)
            vjm      = (u(j,2)-u(j-1,2))/djm
            v(j-1,2) = (vjp*djm + vjm*djp)/(djm+djp)
            vjp      = vjm
            djp      = djm
         enddo
         v(np,2)= -v(np-1,2)+2.*(u(np,2)-u(np-1,2))/deta(np-1)
c        do j = 1,np
c           v(j,2) = (d2(j) + v(j,2))*0.5
c        enddo
      else 
         call amean(j1,np,eta,u(1,2),d1,1)
      endif
cc
c  check for flow separation
c  smooths velocity profiles in the flow reversal region
c  to enhence numerical stability. this logic is well
c  justified because a flare approximation has been
c  used in the difference formulation where there is
c  a region of flow separation
c
      if(umin .ge. -0.005) return
c
      do j = np,2,-1
         jj  = j
         if(u(j,2) .lt. 0.0) goto 20
      enddo
      return
   20 j1    = min0(np,jj + 2 )
      call amean(1,j1 ,eta,f(1,2),d1,1)
      call amean(1,j1 ,eta,u(1,2),d1,1)
      if(iops .eq. 1) then
         do j = 1,np
            d2(j) = v(j,2)
         enddo
         call amean(1,j1 ,eta,d2,d1,1)

         do j= j1-1,2,-1
            djp      = deta(J  )
            vjp      = (u(j+1,2)-u(j,2))/djp
            djm      = deta(j-1)
            vjm      = (u(j,2)-u(j-1,2))/djm
            v(j  ,2) = (vjp*djm + vjm*djp)/(djm+djp)
         enddo
         v(1,2)= -v(2,2)+2.0*(u(2,2)-u(1,2))/deta(1)  
         do j = 1,np
            v(j,2) = 0.5 * (d2(j) + v(j,2))
         enddo
      else
         call amean(1,j1 ,eta,v(1,2),d1,1)
      endif
      return
      end
c
      subroutine output
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blprc2/ ff(my),uu(my),ww(my),vv(my)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),dp(mx),gi,ddeps,om
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx),uedls(mx)
      dimension       d1(my)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      npsav(nx)  = np
      itsav(nx)  = it
      if(nx .eq. 1 ) then
         cfs(nx) = 0.0
         call growth(2) 
         return
      endif
      sum  = 0.0
      do 5 j =1,np
         d1(j) = u(j,2) * (1.0-u(j,2)/u(np,2))
    5 continue
      do 10 j=2,np
         sum  = sum  + a(j)*(d1(j)  + d1(j-1))
   10 continue
      if(itv .ne. 0) then
         rth      = sqrx * sum
         dls(nx)= x(nx)/sqrx*(eta(np)- f(np,2))
         theta(nx)= x(nx)/sqrx*sum
         rds      = sqrx * (eta(np) - (f(np,2)-f(1,2)))
         cfs(nx)= 2.0 * v (1,2) / sqrx
         vw(nx) = v(1,2)
      else
         sqxc     = sqrt(x(nx))
         dls(nx)  = sqxc/sqrl*(eta(np)- f(np,2)/u(np,2))
         theta(nx)= sqrt(x(nx)/rl)*sum/u(np,2)
         rth      = sqrl * sum * sqxc
         rds      = rx * (dls(nx)/x(nx))
         if(nx.le.nxte) then
            vw(nx) = v(1,2)
            cfs(nx)= 2.0*v(1,2)/(sqrl*sqxc *u(np,2) ** 2)
            uc(nx) = 0.0
         else
            uc(nx) = u(1,2)
            cfs(nx)= 0.0
            vw(nx) = 0.0
         endif
         uevs(nx) =  u(np,2)
      endif
      dd(nx)   = dls(nx) * uevs(nx) * sqrl
      if(itv .ne. 0) db(nx) = dd(nx)
      vw(nx)   =  v(1,2)
      if(nx .gt. ntr) call calfa
      call growth(2)
      return
      end
c
c
      subroutine calfa
c
c  function: calculate alpa based on simpson'argument and cebeci &
c            kcc correlation
c  subroutine is called from subroutine eddy
c
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /bledy1/ xtr,gamtr(mx),ffs(mx),alfas(mx),edv(my)
      common /blgtyy/ p1(mx),p2(mx)
      common /blinp3/ rl, rx, sqrl, sqrx
c
c   modify outer eddy based on simpson suggestion
c
c   step 1 calculate (du/dx)/(du/dy) at the location where shear is
c   maximum
c
      if(itv .ne. 0) then
         sqrey = sqrx
      else
         sqxc  = sqrt(x(nx))
         sqrey = sqxc * sqrl
      endif
c
      tm  = -9999.
      do j=1,np
         tt = edv(j)*v(j,2)
         if(tt .gt. tm) then
            jm = j
            tm = tt
         endif
      enddo
      vm       = v(jm,2)
      if(vm.le.0.0 .or. tm.le.0.0) return
c
      du    = x(nx)*(u(jm,2)-u(jm,1))/(x(nx) - x(nx-1))
      dudx  = du+p2(nx)*u(jm,2)+0.5*eta(jm)*v(jm,2)*(p2(nx)-1.0) 
      rr    = dudx/v(jm,2)/sqrey
c
c   step 2 : calculate (uu - vv)/uv
c
      rt        = amax1(0.0,v(1,2)/tm)
      if(rt.gt.1.0) then
         cr     = (1.0 + rt) /rt
      else
         cr     = 6.0 /(1.0+2.0*rt*(2.0-rt))
      endif
c
c   step 3 : calculate ff
c
      fr   = cr * rr
c
      fr   = amin1(fr,0.50)
      fr   = amax1(fr,-1.50)
      rex  = 0.50
      ffs(nx)   = (1.- rex)*ffs(nx) + (1.0 -fr)*rex
      alfas(nx) = 0.0168 / ffs(nx)**1.5
c
c  limit on the change of alfas from the previous station
c
      nxx      = nx - ntr
      factor   = amin1(0.50,0.15+0.10*nxx)
ccc   factor   = amin1(0.75,0.15+0.10*nxx)
      alfasu   = alfas(nx-1) * (1.+factor)
      alfasl   = alfas(nx-1) * (1.-factor)
      alfas(nx)= amin1(alfasu, alfas(nx))
      alfas(nx)= amax1(alfasl, alfas(nx))
      ffs(nx)  = (0.0168/alfas(nx))**(1./1.5)
c
      return
      end

      subroutine gamcal(iwk,np,gamtr,eta,u,fint,du)
c
c  calculate intermittency
c
      dimension eta(np),u(np),fint(np),du(np)
c
      dimension hh(13),ydel(13),sgdel(13)
      data hh/1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4/
      data ydel/0.50,0.665,0.8,0.843,0.864,0.88,0.89,0.898,0.905,
     1          0.913,0.921,0.928,0.936/
      data sgdel/0.180,0.156,0.142,0.136,0.133,0.131,0.1285,0.1279,
     1           0.1276,0.1273,0.127,0.1267,0.1264/

c
c  --------------------------------------------------------------------
c
      udel    = 0.995 * (u(np)-u(1))
      do j = 2 , np
         jj    = j
         du(j)    = u(j) - u(1)
         if (u(j) .gt.udel) goto 20
      enddo
   20 edel   = eta(jj-1)+(eta(jj)-eta(jj-1))/(du(jj)-du(jj-1))*
     1         (udel-du(jj-1))
c
      if(iwk .eq. 0) go to 50
c
c  cal fint for wake flow
c
      do j = 1,np
         etadel  = amin1(1.25,eta(j)/edel) 
         fint(j) = 1. / (1.+5.5*etadel**6)
      enddo 
      return
c
c   calculate fint ( gamma in c. s. model) using head's correlation
c
c        calc. of h = dels/theta (based on x-value for 3d case)
c
50    continue
      uj    = u(1)/u(np)
      termp = (1.-uj)
      termpt= uj*(1.-uj)
      sumu  = 0.0
      sumt  = 0.0
      do j  = 2,np
         aj    = 0.5*(eta(j)-eta(j-1))
         uj    = amin1(1.0,u(j) / u(np))
         term  = 1.- uj
         sumu  = sumu + aj*(term+termp)
         termp = term
         termt = uj * (1.- uj)
         sumt  = sumt + aj*(termt+termpt)
         termpt= termt
      enddo
      hcal  = sumu / sumt
      heq   = 1.4
      sqgmr = gamtr**2
      hcal  = sqgmr * hcal + (1.-sqgmr)* heq
      if(hcal.le. hh(1)) then
         ydelz  = ydel(1)
         sgdelz = sgdel(1)
      elseif(hcal.ge. hh(13)) then
         ydelz  = ydel(13)
         sgdelz = sgdel(13)
      else
         do j  = 2,13
            jj = j-1
            if(hcal .lt. hh(j) ) goto 60
         enddo
60       ydelz  = ydel(jj)+(ydel(jj+1)-ydel(jj))*(hcal-hh(jj))
     1              /(hh(jj+1)-hh(jj))
         sgdelz = sgdel(jj)+(sgdel(jj+1)-sgdel(jj))*(hcal-hh(jj))
     1              /(hh(jj+1)-hh(jj))
      endif
c
      sq2   = sqrt(2.0)
      do j  = 1,np
         etadel = amin1(1.25,eta(j)/edel)
         z      = (etadel-ydelz)/(sq2*sgdelz)
         z1     = abs(z)
         t      = 1.0/(1.0 + 0.3275911*z1)
         erfz   = 1.0 + t*(-0.254829592 +t*(0.284496736 +t*
     1          (-1.421413741+ t*(1.453152027 - 1.061405429*t))))*
     1          exp(-z1**2)
         if(z .lt. 0.0) erfz=-erfz
         fint(j)= 0.5*(1.0 - erfz)
      enddo
c
      return
      end
c
      subroutine trns(isep,itran)
c
c  determine the transition location
c
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blprc2/ ff(my),uu(my),ww(my),vv(my)
      common /bledy1/ xtr,gamtr(mx),ffs(mx),alfas(mx),edv(my)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx),
     1                uc(mx),uep(mx),uedls(mx)

      dimension d1(my)
c
c  ------------------------------------------------------------------
c
      itran = 0
      if(isep .eq. 1) then
c
c  laminar flow separation, move back one station (for numerical
c  stability) and set transtion location there
c
c        write(6,*) 'lam. separation at nx = ',nx
         nx    = nx-1
         np    = npsav(nx)
         do j  = 1,npt 
            f(j,2) = f(j,1)
            u(j,2) = u(j,1)
            v(j,2) = v(j,1)
            w(j,2) = w(j,1)
            f(j,1) = ff(j)
            u(j,1) = uu(j)
            v(j,1) = vv(j)
            w(j,1) = ww(j)
         enddo
         itran = 1
cc       write(6,*) 'lam. separation at nx = ',nx           
      else
c
c  search for ntr corresponding to the input xtr
c
         if ( xtr .gt. 0.0 ) then
            if (xbl(nx).ge. xtr) itran = 1
         else
c
            if (xbl(nx) .lt. abs(xtr)) itran = 1
         endif
c
         if(itran .eq. 0 .and. xtr .le. xbl(nxte)) then
c
c  check transition location based on michel method 
c
            do j =1,np
               d1(j) = u(j,2) * (1.0-u(j,2)/u(np,2))
            enddo
            sum = 0.0
            do j= 2,np
               sum  = sum  + a(j)*(d1(j)  + d1(j-1))
            enddo
            if(itv .ne. 0) then
               rx       = x(nx)* ubl(nx) * rl
               rth      = sqrt(rx) * sum
            else
               sqx      = sqrt(x(nx))
               rth      = sqrl * sum * sqx 
               rx       = x(nx)* u(np,2) * rl
            endif
            rm  = 1.174*(1.0 + 22400.0/rx)*rx**(0.46)
            if(rm .lt. rth) itran = 1
c           if(itran .eq. 1) 
c    1         write(6,*) 'nx,iswp,rth,rm',nx,iswp,rth,rm
         endif
      endif
c
      if(itran .eq. 0) return
c
c  calcualte gamtr
c
      ntr     = nx
      xtr     = xbl(ntr)
      rxntr   = x(ntr-1)* ubl(ntr-1) * rl
c
c  modified by kcc on 9/5/2002
c
      if(rxntr .lt. 0.6e06) then
         csqr = 213.0*(alog10(rxntr) - 4.7323)
         csqr = amax1(csqr,360.0)
      else
         csqr = 222.7
         crl  = (rxntr / 0.6e+06)**1.5
         csqr = amin1(3600.,crl * csqr)
         csqr = amax1(csqr,360.)
      endif
c
      cgg    = 3.0/csqr
c
c     cgg    = 1. /1200
      ggft    = cgg *rl**2/rxntr**1.34*ubl(ntr-1)**3
      ueintg = 0.0
      u1     = 1.0/ubl(ntr-1)
      nocal  = 0
      do i = ntr,nxte
         if(nocal .eq. 1) then
            gamtr(i) = 1.0
         else
            u2      = 1.0/ubl(i)
            ueintg  = ueintg+0.5 *(u1+u2)*(x(i)-x(i-1))
            u1      = u2
            gg      = ggft*ueintg*(x(i)-x(ntr-1))
            if(gg .gt. 10.0) then
               gamtr(i)= 1.0
               nocal   = 1
            else
               gamtr(i)= 1.0-exp(-gg)
            endif
            alfas(i) = 0.0168
            ffs(i)   = 1.0
         endif
      enddo
c
c     nocal = 1
      if(nocal .eq. 0 ) then
c
c  adjust gamtr so that at the trailing edge, gamtr = 1.0
c
         dgamtr = (1.0 - gamtr(nxte))/(x(nxte) - x(ntr-1))
         do i = ntr,nxte
            gamtr(i) = gamtr(i) + dgamtr * (x(i) - x(ntr-1))
            gamtr(i) = amin1(1.0,gamtr(i))
         enddo
      endif
c
      return
      end
c
      subroutine growth(index) 
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my),etae,vgp

c
      npo   = np
      if(index .eq. 1) then
         np    = min0(npt,np+2)
         npend = np
      else
         npend = npt
      endif
c
      do j=npo+1,npend
         v(j,2)= 0.0
         u(j,2)= u(npo,2)
         f(j,2)= f(j-1,2)+u(j,2)*deta(j-1)
         b(j,2)= b(npo,2)
         w(j,2)= w(npo,2)
      enddo
      return
      end
c
      function ismin(n,a,m)

      dimension a(m,n)
      amin  = a(m,1)
      ismin = 1
      do 10 i=1,n
         if (a(m,i) .lt. amin) then
            amin  = a(m,i)
            ismin = i
         end if
   10 continue
      return
      end

