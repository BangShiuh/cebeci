c***************************************************************************
c     main
C  Inverse Boundary Layer Program
C
      parameter(mx=201,my=151,mx2=mx*2)
      complex alfa,omega
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blprc2/ ff(my),uu(my),vv(my),ww(my)
      common /blc2so/ delf(my),delu(my),delv(my),delw(my)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blgtyy/ p1(mx),p2(mx)
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),gi,ddeps,om
      common /bledy1/ gamtr(mx),ffs(mx),alfas(mx),edv(my),isep 
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx)
      common /blstab/ istab,icont,nx0,ixt,uinf,bigl,alfar, omegar,
     1                ustab(my,mx),vstab(my,mx),etastab(my)
      dimension       rdelst(mx),h(mx),d1(mx2)
      dimension       rdelst_cri(14),h_cri(14)
      data rdelst_cri/12940.0,10920.0,8890.0,7680.0,6230.0,4550.0, 
     1                2830.0,1380.0,865.0,520.0,318.0,199.0,138.0,
     1                67.0/
      data h_cri/2.216,2.240,2.274,2.297,2.325,2.362,2.411,2.481,     
     1          2.529,2.591,2.676,2.801,2.963,4.029/
      character*80 dd_file
      logical conv,lgc1,lgc2,lgc3 
      CHARACTER*80 input_name, output_name, StabInput
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  lagrange interpolation or extrapolation
      grang(x1,x2,x3,y1,y2,y3,x0)= (x0-x2)*(x0-x3)/(x1-x2)/(x1-x3)*y1
     1      +(x0-x1)*(x0-x3)/(x2-x1)/(x2-x3)*y2+(x0-x1)*(x0-x2)
     2      /(x3-x1)/(x3-x2)*y3
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
      if ( IDOT .GE. 5 ) IDOT = IDOT-3 
      DO I = 1, IDOT
         output_name(I:I) = input_name(I:I)
	   StabInput(I:I) = input_name(I:I)
	ENDDO
	if ( IDOT .EQ. Iend ) THEN
         OUTPUT_NAME(IDOT+1:IDOT+7)   = "Out.txt"
	   StabInput(IDOT+1:IDOT+10) = "StbInp.txt"
	ELSE
         OUTPUT_NAME(IDOT:IDOT+6)   = "Out.txt"
	   StabInput(IDOT:IDOT+9) = "StbInp.txt"
	ENDIF
	print*, "Your input file is ", input_name
      print*,"This program bl.exe will generate several output files:"
	print*, "1. output file name ",   OUTPUT_NAME
      print*, "2. Generated stability input:",  StabInput 
	print*, "3. Velocity Profile in file fort.1"

      open(unit=5, file=input_name, status="OLD")
      open(unit=6, file=output_name)

      call input 
      call ivpl
      call hic(nxt,rl,x(1))
      ddeps = 0.0025
      iswp  = 0
      itv   = 1
      nx    = 0
      conv  = .false.
      om    = 1.5
   50 iswp  = iswp + 1
  100 nx    = nx + 1
      rx    = uevs(nx)*x(nx)*rl
      sqrx  = sqrt(rx)
      write(6,6000) nx,xbl(nx),x(nx)
      if(nx.gt.1) then
         do j=1,npt 
            ff(j)  = f(j,1)
            uu(j)  = u(j,1)
            vv(j)  = v(j,1)
            ww(j)  = w(j,1)
            f(j,1) = f(j,2)
            u(j,1) = u(j,2)
            v(j,1) = v(j,2)
            w(j,1) = w(j,2)
            b(j,1) = b(j,2)
         enddo
      end if
c
      if (itv .eq. 0) then
         sum2     = 0.
         do k = nx+1,nxt+1
            sum2     = sum2 + cc(k,nx)*dd(k)
         enddo
         sum1     = 0.
         do k = 3,nx-1
           sum1     = sum1 + cc(k,nx)*dd(k)
         enddo
         gi       = ubl(nx) + sum1 + sum2
      end if
c
      it    = 0
      igrow  = 0
  120 it    = it+1
      if(it .ge. 15) then
         write(6,9999) nx
         go to 999
      endif
      if(nx .ge. ntr) call eddy
      iback    = 0
      if(itv .eq. 0 ) then
         iback = 1
         if(iswp .gt. 3) then
            iback = 2
            lgc1  = .false.
            if(nx .le. nxs + 3) lgc1 = .true.
            lgc2  = .false.
            if(nx .eq. nxte+1 .or. nx .eq. nxte+2) lgc2  = .true.
            lgc3  = .false.
            if(nx .eq. ntr .or. nx.eq.ntr+1 ) lgc3 = .true.
            if(lgc1 .or. lgc2 .or. lgc3) iback = 1
         endif
      endif
c  
      call coef(iback) 
      call solv4
      if(nx .le. nxte) then
         write(6,*) 'iswp,it,v1,delv1 =',iswp, it, v(1,2), delv(1) 
      else
         write(6,*) 'iswp,it,u1,delu1 =',iswp, it, u(1,2), delu(1)
      endif
c
c  flow separation. check for transtion.
c
      if (v(1,2) .le. 0.0) then
         if(itv .eq. 1) then
            write(6,201) nx,x(nx),ubl(nx),p1(nx),p2(nx),v(1,2)
            stop
         else
            if(nx .lt. ntr ) then
               if(isep .eq. 1) then
c
c  turn on transition. and calculations continue
c
                  call trns(1) 
                  it   = 1
                  igrow= 0
                  goto 100
               else
c
c  calculations continue as a laminar flow until the specified
c  transition location is met
               endif
            endif
         endif
      endif
c
      if(nx .le. nxte) then
c  convergence of the solutions - b. l. 
         if(nx .lt. ntr ) then
            if(abs(delv(1)).gt. 0.0001) goto 120
         else
            v1  = amax1(1.0, v(1,2))
            if(abs(delv(1)/v1) .gt. 0.01) goto 120
         end if
      else
c  convergence of the solutions - wake
         if(abs(delu(1)).gt. 0.0025) goto 120
      endif
      if(np .lt. npt .and. igrow .lt. 3) then
         if( abs(v(np,2)) .gt. 0.0010) then
            call growth(1) 
            igrow = igrow + 1
            it    = 1
            go to 120
         endif
      endif
  999 continue
      call output
      write (6,9000) (j,eta(j),f(j,2),u(j,2),v(j,2),b(j,2),j=1,np-1,4) 
      j  = np
      write (6,9001) j,eta(j),f(j,2),u(j,2),v(j,2),b(j,2) 
      if(nx .lt. nxs ) go to 100
      if (nx .eq. nxs) then
        itv     = 0
        call swtch(0)
        go to 100
      end if
c
c  logic to speed up the convergence
c
c     if(iswp .lt. 5) then
         ddb   = 0.05 * ( dd(nx) - db(nx))
         do i = nx+1,nxt
            ddfact  = 1.0
            if(i.gt.nxte) ddfact = x(nxte) / x(i)
            dd(i) = db(i) + ddb * ddfact
         enddo
c     endif
c
      if(nx .lt. nxt  ) then
         if(nx.eq.nxte) then
             call wakepr
         endif
         goto 100
      endif
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
c
c        write (8,9110) iswp, ddmax/dmax 
c        write (8,9120) (i,xbl(i),uep(i),uevs(i),cfs(i),uc(i),
c    1                theta(i),dls(i),dd(i),db(i),i=2,nxt)

         if(conv) then
            write(6,111) iswp,ddmax/dmax 
 111        format(' inv/vis calculations converge at sweep = ',
     1             i3,' ddmax/dmax =',f10.5,' cal. stop')
            goto 150
         endif
c
         if(iswp .gt. 1) then
c
c  over-relax d values to speed-up convergence
c
            do i= nxs,nxt
               dd(i)=dd(i)*(1.0+om*(uep(i)/uevs(i)-1.0))
               uep(i) = uevs(i)
               db(i)= dd(i)
            enddo
c
         else 
            do i= nxs,nxt
               dd(i)=dd(i)*(1.0+om*(uevs(i)/ubl(i)-1.0))
               uep(i)= uevs(i)
               db(i)= dd(i)
            enddo
         endif
c
c        call amean(nxs,nxt,x,dd,d1,1)
c
         dd(nxt+1)=grang(x(nxt-2),x(nxt-1),x(nxt),
     1            dd(nxt-2), dd(nxt-1),dd(nxt),x(nxt+1))
         nx     = nxs
         call swtch(1)
         goto  50
      end if
c
150   call amean(nxs,nxte,x,cfs,d1,1)
      write(6,9100)  nxs,ntr
c        
      write(6,9110) iswp, ddmax/dmax
9110  format( 'sweep = ',i3, ' ddmax/dmax =',f10.5)
      write (6,9130) (i,x(i),xbl(i),ubl(i),uevs(i),cfs(i),uc(i),
     1                theta(i),dls(i),dd(i),db(i),i=2,nxt)
      write(6,*)
      write(6,*) 'normal termination in main '
      
      iqwik = 88
      itecp = 87
c     iddfl = 86
      open(88,file='qwik_file',form='unformatted')
      open(87,file='tecp_file',form='formatted')
c     open(86,file='iblp_out.ddfile',form='formatted')

c
      rewind iqwik
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
c  write unit 1 file for stability calculations
c
      if(istab .eq. 1) then
c
         if(nx0 .gt. ntr-1) then
            write(6,*) '1st station for stabiity calculation ',
     1        ' falls in turbulent region. no files are',
     1        ' generated for stability analysis'
            goto 1100
         endif
c
c  write unit 2 file as standard input file for stability
c  calculations
c
         do i = 2,ntr-1
            h(i) = dls(i) / theta(i)
            rdelst(i) = rl * uevs(i) * dls(i)
         enddo
c
990      do i = 1,14
            ii  = i
            if(h(nx0) .lt. h_cri(i)) goto 1000 
         enddo
1000     rdel_c = rdelst_cri(ii-1)+(rdelst_cri(ii)-
     1      rdelst_cri(ii-1))/(h_cri(ii)-h_cri(ii-1))*
     1      (h(nx0)-h_cri(ii-1))
         if(rdelst(nx0) .lt. rdel_c) then
            nx0 = nx0 + 1    
            if(nx0.gt.ntr-1) then
               write(6,*) 'no stability cal. can be performed',
     1           ' no files generated for stability analysis'
               goto 1100
            endif
            goto 990
         endif
c
         ixt    = min0(ixt,(ntr-1) - nx0)
         alfai  = 0.0
         omegai = 0.0
         alfa   = cmplx (alfar,alfai)
         omega  = cmplx (omegar,omegai)
         open(2,file=StabInput)        
         write(2,*) ' nxt   nx0   ixt  icont'
         write(2,910) ntr-1,nx0,ixt,icont 
         write(2,*) '   uinf      bigl       rl'
         write(2,920) uinf,bigl,rl
         write(2,*) '  alfa                 omega'
         write(2,930) alfa ,omega
         write(2,*) 'x/bigl(i) = 1,nxt'
         write(2,940) (xbl(i),i=1,ntr-1)
         write(2,*) 's/bigl(i), i = 1,nxt'
         write(2,940) (x(i),i=1,ntr-1)
         write(2,*) 'ue/uinf(i),i=1,nxt)' 
         write(2,940) (uevs(i),i=1,ntr-1)
910      format(4i5)
920      format(2f10.3,f10.1) 
930      format(4f10.5)
940      format(7f10.6)
c
         do i = 1,ntr-1
            np  = npsav(i)
            write(1,666) i,np
            write(1,888) (etastab(j),j=1,np)
            write(1,888) (ustab(j,i),j=1,np)
            write(1,888) (vstab(j,i),j=1,np)
         enddo
 666     FORMAT(2I5)
 888     FORMAT(6E14.6)

      endif
c
        
c
c  tecplot file
c
1100  continue
      write(itecp,99)
99    format('title= "2d inverse b. l. summary"')
      write(itecp,110)
110   format('variables="s","x/c","uei","delst","theta","vw"',
     1       '"cf","uev" ')
      write(itecp,121) nxte 
121   format('zone',' i=',i5,' F = point')
      do i = 1,nxte
         write(itecp,122) x(i),xbl(i),ubl(i),dls(i),theta(i),
     1            vw(i),cfs(i),uevs(i)
122      format(8f10.6)
      enddo
      close(itecp)

      close(5)
	close(6)
	close(2)

	PRINT*," "
	PRINT*,"Calculations are successfully completed."
C      PRINT*,"The output is saved in ", OUTPUT_NAME
	PRINT*," "
 	PRINT*,"Hit any key to close this DOS-window."
 	READ(5,*) 	
c
      stop
  201 format (' flow separates in the standard mode; print data',
     1        ' for program check '/1h ,2x,'nx = ',i3,2x,'x = ',
     2        f10.5,2x,'ubl = ',f10.5,2x,'p1 = ',f10.5,2x,'p2 = ',
     3        f10.5,2x,'vw = ',f10.5,/)
 6000 format(/'nx=',i3,2x,'x/c=',f10.5,3x,'s =',f10.5) 
 9000 format('  j ',4x,'eta',9x,'f',13x,'u',13x,'v',13x,'b'/
     &       (i3,f8.3,4e14.5))
 9001 format(i3,f8.3,4e14.5)
 9120 format (3h nx,5x,2hxc,5x,4huep ,4x,4huvis,8x,3h cf,8x,3huc ,6x,
     &    5htheta,8x,3hdls,8x,3h dd,8x,3h db/(i3,f9.6,2f8.4,6e11.3)) 
 9130 format (3h nx,6x,2h x,7x,2hxc,5x,4huei ,4x,4huev ,8x,
     &       3h cf,8x,3huc ,6x,5htheta,8x,3hdls,8x,3h dd,8x,3h db
     &       /(i3,2f10.6,2f8.4,6e11.3))

 9100 format(//'** summary of boundary layer solutions '/
     &        'nxs =',i3,4x,'** transition location : ntr =  ',i3)
 9999 format('iteration exceeds max  at  nx =',i3/)
      end

      subroutine input 
      parameter(mx=201,my=151,mx2=mx*2)
      common /blgtyy/ p1(mx),p2(mx)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),gi,ddeps,om
      common /blinp3/ rl, rx, sqrl, sqrx
      common /bledy1/ gamtr(mx),ffs(mx),alfas(mx),edv(my),isep 
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx)
      common /blstab/ istab,icont,nx0,ixt,uinf,bigl,alfar, omegar, 
     1                ustab(my,mx),vstab(my,mx),etastab(my)
      dimension ddin(mx),dlsin(mx),xin(mx)
      character*80 ddp_file 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      grang(x1,x2,x3,y1,y2,y3,x0)= (x0-x2)*(x0-x3)/(x1-x2)/(x1-x3)*y1
     1      +(x0-x1)*(x0-x3)/(x2-x1)/(x2-x3)*y2+(x0-x1)*(x0-x2)
     2      /(x3-x1)/(x3-x2)*y3
c
      read (5,*  )
      read (5,*  ) nxt,nxte,nxs,npt,iswpt,isep,istab,icont 
      read (5,*  )
      read (5,*  ) rl,xtr,etae,vgp,deta(1), p2(1)
      read (5,*  )
      read (5,*  ) (xbl(i),i=1,nxt)
      read (5,*  )
      read (5,*  ) (ybl(i),i=1,nxt)
      read (5,*  )
      read (5,*  ) (ubl(i),i=1,nxt)
c
      if(istab .eq. 1) then
c
c  input additional data for generating input file for bl stability
c  calcualtion
c
         read(5,*)
         read(5,*) nx0,ixt,uinf,bigl,alfar,omegar
      endif
c
      if ((vgp-1.0) .gt. 0.001) then    
        np = alog((etae/deta(1))*(vgp-1.0)+1.0)/alog(vgp)+1.0001
      else    
        np = etae/deta(1)+1.0001
      end if
	write (6, *) "Eta direction dimension np = ", np
      if (np .gt. npt) then
        write (6,9000)                                                    
        stop 
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
c     determine transition location                                     
      xmin = 9999.0
      do i = 1,nxte
         if(xbl(i) .lt. xmin) then
            imin = i
            xmin = xbl(i)
         endif
      enddo
c
      if(xtr .lt. xbl(nxte)) then
         do 30 i=imin,nxte                                                    
            if (xbl(i) .ge. xtr) go to 40                                  
   30    continue                                                          
   40    ntr = i                                                           
c     calculate gamtr                                                   
         rxntr   = x(ntr-1)* ubl(ntr-1) * rl                                
         ggft    = rl**2/rxntr**1.34*ubl(ntr-1)**3                          
         ueintg  = 0.0                                                     
         u1      = 0.5/ubl(ntr-1)/ 1200.0                                   
         do i = ntr,nxt                                                
            u2      = 0.5/ubl(i)/1200.0                                     
            ueintg  = ueintg+(u1+u2)*(x(i)-x(i-1))                         
            u1      = u2                                                   
            gg      = ggft*ueintg*(x(i)-x(ntr-1))                          
            if(gg .gt. 10.0) then                                          
               gamtr(i)=1.0                                                 
            else                                                           
               gamtr(i)= 1.0-exp(-gg)                                       
            endif                                                          
            alfas(i) = 0.0168
         enddo
      else
         ntr = nxt + 1
      endif
c
      write(6,6000) rl,ntr,xtr
      write (6,6040) etae,vgp,deta(1)
c
      do i=1,nxt
         uevs(i) = ubl(i)
         uep(i)  = ubl(i)
      enddo
c  initial deltastar distribution
      sqrl = sqrt(rl)
      dd(1) = 0.0
      dls(i)= 0.0
      do i = 2,nxte
         rx     = rl*x(i)
         dls(i) = 0.046875*x(i)/rx**0.2
         dd(i)  = ubl(i)*dls(i)*sqrt(rl)
         db(i)  = dd(i)
      enddo
      do i=nxte+1,nxt
         dd(i) = dd(nxte)/exp( (x(i)-x(nxte))/(x(nxt)-x(nxte)))
         db(i) = dd(i)
      enddo
c
      x(nxt+1)  = 3.0 * x(nxt) -3.0 * x(nxt-1) + x(nxt-2)
      dd(nxt+1) =grang(x(nxt-2),x(nxt-1),x(nxt),
     1              dd(nxt-2),dd(nxt-1),dd(nxt),x(nxt+1))
      p21   = p2(1)
      call diff1(nxt,x,ubl,p2,p1) 
      p2(1) = p21
      p1(1) = 0.5*(p2(1) + 1.0)
      do i=2,nxt 
         p2(i) = x(i)/ubl(i)*p2(i)
         p1(i) = 0.5*(1.0+p2(i))
      enddo
      return                                                            
 5000 format(a)
 6000 format(/' rl =',e10.3,5x,'ntr =',i5,5x,'xctr =',f7.3/)                     
 6040 format(/' etae=',f8.3,3x,'vgp=',f7.3,3x,'deta(1)=',f7.3/)
 9000 format('np exceeded npt - program terminated')               
      end                                                               

      subroutine ivpl
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
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
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),gi,ddeps,om
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
      common /bledy1/ gamtr(mx),ffs(mx),alfas(mx),edv(my),isep 
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blcwak/ edsave(my),ytedge(my),wscale,ydelwk
      dimension       fint(my),yy(my),edintp(my)
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
      call gamcal(iwk,np,gamtrnx,eta,u(1,2),fint)
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
            tt = v(j,2)*(wsrex +(1.-wsrex) * edv(j))
            if(tt.gt. ttm) then
               ttm = tt
               jm  = j    
            endif
         enddo
      endif
      smax     = amax1(swall, ttm)
      um       = sqrt(smax)
c
      rk    = 0.41
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

      subroutine coef(iback) 
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blprc2/ ff(my),uu(my),vv(my),ww(my)
      common /blcs30/ s1(my),s2(my),s3(my),s4(my),s5(my),s6(my),
     1                s7(my),s8(my),r1(my),r2(my),r3(my),r4(my)
      common /blcs31/ gamma1,gamma2
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blgtyy/ p1(mx),p2(mx)
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),gi,ddeps,om
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
        s1(j) = b(j,2)/deta(j-1) + 0.5*p1p*fb -0.5*cel*cfb
        s2(j) = -b(j-1,2)/deta(j-1) + 0.5*p1p*fb - 0.5*cel*cfb      
        s3(j) = 0.5*(p1p*vb + cel*cvb)                                
        s4(j) = 0.5*(p1p*vb + cel*cvb)                              
        s5(j) = -p2p*ub
        s6(j) = -p2p*ub
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
     1                uc(mx),uep(mx)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blgtyy/ p1(mx),p2(mx)
      common /blsave/ fsave(my),usave(my),vsave(my),bsave(my),
     1                wsave(my),npsave
c    -----------------------------------------------------------
      if(index.gt.0) goto 40
      wnp   = ubl(nx)
      ak    = sqrt(1.0/wnp)
      do j=1,np
         u(j,2) = wnp*u(j,2)
         v(j,2) = wnp*v(j,2)/ak
         f(j,2) = wnp*f(j,2)*ak
         fsave(j)  = f(j,2)
         usave(j)  = u(j,2)
         vsave(j)  = v(j,2)
         bsave(j)  = b(j,2)
      enddo
      do j = 1,np
         w(j,2) = u(np,2)
         wsave(j)  = w(j,2)
      enddo
      do 30 j= 1,npt
         eta(j) = ak*eta(j)
   30 continue
      do 35 j=2,npt
         deta(j-1)= eta(j)-eta(j-1)
         a(j)  = 0.5*deta(j-1)
   35 continue
      npsave = np
      do i = nxs+1,nxt
         p2(i) = 0.0
         p1(i) = 0.5
      enddo
      return
   40 np     = npsave
      do j=1,np
         f(j,2) = fsave(j)
         u(j,2) = usave(j)
         w(j,2) = wsave(j)
         v(j,2) = vsave(j)
         b(j,2) = bsave(j)
      enddo
      return
      end

      subroutine wakepr
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /bledy1/ gamtr(mx),ffs(mx),alfas(mx),edv(my),isep 
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blinp8/ d1(mx2),d2(mx2),d3(mx2),d4(mx2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blcwak/ edsave(my),ytedge(my),wscale,ydelwk
      dimension       y(mx),uu(mx),vv(mx)
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

      subroutine solv4                                                  
      parameter(mx=201,my=151,mx2=mx*2)
      common /blcs31/ gamma1,gamma2
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blc2so/ delf(my),delu(my),delv(my),delw(my)
      common /blcs30/ s1(my),s2(my),s3(my),s4(my),s5(my),s6(my),
     1                s7(my),s8(my),r1(my),r2(my),r3(my),r4(my)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx)
      dimension     a11(my),a12(my),a13(my),a14(my),a21(my),            
     1              a22(my),a23(my),a24(my),g11(my),g12(my),            
     1              g13(my),g14(my),g21(my),g22(my),g23(my),            
     1              g24(my),w1(my),w2(my),w3(my),w4(my)                 
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
      do 30 j = 1,np                                                    
         f(j,2) = f(j,2)+delf(j)                                       
         u(j,2) = u(j,2)+delu(j)                                       
         v(j,2) = v(j,2)+delv(j)                                       
         w(j,2) = w(j,2)+delw(j)                                       
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
      call edgchk
c
      return                                                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      end                                                               
c
      subroutine edgchk
c
c  check for edge balance & smooth the velocity profiles if necessary
c
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      dimension d1(my)

c    check for maximum v
c
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
   40 j1     = jj-1
c
c  smooth velocity profiles near the edge of b. l. to eliminate
c  unnecessary growth because of the balance there 
c
      call amean(j1,np,eta,f(1,2),d1,1)
      call amean(j1,np,eta,u(1,2),d1,1)
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
      return
      end
c
      subroutine output
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),gi,ddeps,om
      common /blinp3/ rl, rx, sqrl, sqrx
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx), 
     1                uc(mx),uep(mx)
      common /blstab/ istab,icont,nx0,ixt,uinf,bigl,alfar, omegar,
     1                ustab(my,mx),vstab(my,mx),etastab(my)
      dimension       d1(my)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      npsav(nx)  = np
      itsav(nx)  = it
      if(nx .eq. 1 ) then
         cfs(nx) = 0.0
         do j = 1,npt
            etastab(j) = eta(j)
         enddo
         goto 1000
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
c
1000  continue
      call growth (2)
c
c  save data for stability calculation
c
      if(nx .lt. ntr) then
         if(itv .eq. 1) then
            do j = 1,np
               ustab(j,nx) = u(j,2)
               vstab(j,nx) = v(j,2)
            enddo
         else
            wnp   = u(np,2)
            ak    = sqrt(1.0/wnp)
            do j = 1,np
               ustab(j,nx) = u(j,2) / wnp
               vstab(j,nx) = v(j,2) / wnp * ak
            enddo
         endif
      endif
      return
      end
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
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /bledy1/ gamtr(mx),ffs(mx),alfas(mx),edv(my),isep 
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
      factor   = amin1(0.75,0.15+0.10*nxx)
      alfasu   = alfas(nx-1) * (1.+factor)
      alfasl   = alfas(nx-1) * (1.-factor)
      alfas(nx)= amin1(alfasu, alfas(nx))
      alfas(nx)= amax1(alfasl, alfas(nx))
      ffs(nx)  = (0.0168/alfas(nx))**(1./1.5)
c
      return
      end

      subroutine gamcal(iwk,np,gamtr,eta,u,fint)
c
c  calculate intermittency
c
      dimension eta(np),u(np),fint(np)
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
      udel    = 0.995 * u(np)
      do j = 2 , np
         jj    = j
         if (u(j) .gt.udel) goto 20
      enddo
   20 edel   = eta(jj-1)+(eta(jj)-eta(jj-1))/(u(jj)-u(jj-1))*
     1         (udel-(u(jj-1)-u(1)))
c
      if(iwk .eq. 0) go to 50
c
c  cal fint for wake flow
c
      do j = 1,np
         fint(j) = 1. / (1.+5.5*(eta(j)/edel)**6)
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
      subroutine trns(jsep)         
      parameter(mx=201,my=151,mx2=mx*2)
      common /blc001/ nx,np,nxte,nxt,npt,ntr,nxs,it,itv,iswp,iswpt
      common /blprc1/ f(my,2),u(my,2),v(my,2),w(my,2),b(my,2)
      common /blprc2/ ff(my),uu(my),vv(my),ww(my)
      common /blgtyy/ p1(mx),p2(mx)
      common /blgrd1/ x(mx),eta(my),a(my),deta(my)
      common /blhlbs/ cc(mx,mx),dd(mx),db(mx),gi,ddeps,om
      common /blinp3/ rl, rx, sqrl, sqrx
      common /bledy1/ gamtr(mx),ffs(mx),alfas(mx),edv(my),isep 
      common /blrc01/ xbl(mx),ubl(mx),ybl(mx),cfs(mx),theta(mx),
     1                dls(mx),uevs(mx),npsav(mx),itsav(mx),vw(mx),
     1                uc(mx),uep(mx)
c  -------------------------------------------------------------------  
c
c  determine the transtion location
c
      itran = 0
      if(jsep .eq. 1) then
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
      endif
c
      if(itran .eq. 1) then
         ntr = nx
c     calculate gamtr
         rxntr   = x(ntr-1)* ubl(ntr-1) * rl
         ggft    = rl**2/rxntr**1.34*ubl(ntr-1)**3
         ueintg  = 0.0
         u1      = 0.5/ubl(ntr-1)/ 1200.0
         do i = ntr,nxt
            u2      = 0.5/ubl(i)/1200.0
            ueintg  = ueintg+(u1+u2)*(x(i)-x(i-1))
            u1      = u2
            gg      = ggft*ueintg*(x(i)-x(ntr-1))
            if(gg .gt. 10.0) then
               gamtr(i)=1.0
            else
               gamtr(i)= 1.0-exp(-gg)
            endif
            alfas(i) = 0.0168
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



