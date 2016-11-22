c	Xavier Dumont

c	NACA airfoil generator

c	the goal of this program is to generate the shape of a
c	NACA airfoil from the NACA four-digit series given the
c	four numbers of its definition


	Program main

	integer nmax, nsur, i, nwake, column, counta, count
	parameter (pi = 4.0*atan(1.0))
	real chord,m,p,t,h,x,yc,yt,dyc,theta,xu,yu,xl,yl,mc,pm,aoa,mach
	dimension x(1000), yc(1000), yt(1000), dyc(1000), theta(1000),
     f		xu(1000), yu(1000), xl(1000), yl(1000), xa(1000), ya(1000)

	write(6,*) "NACA airfoil generator"

C	Getting the information about the airfoil

	chord = 1

	write(6,*) "Enter NACA airfoil 1st digit"
	read(5,*) m
	write(6,*) "Enter NACA airfoil 2nd digit"
	read(5,*) p
	write(6,*) "Enter NACA airfoil 3rd digit and 4th digit"
	read(5,*) t

	write(6,*) "Enter number of points to be calculated (even number)"
	read(5,*) nmax

	nsur = nmax/2
	nwake = 0

	write(6,*) "Enter the angle of attack"
	read(5,*) aoa
	write(6,*) "Enter the mach number"
	read(5,*) mach

C	Generate values of x from 0 to the maximum chord

	h = chord/nsur
	x(1) = 0
	x(nsur+1) = chord

	do i = 2, nsur
        x(i) = x(i-1)+(x(nsur+1)-x(1))/nsur
	enddo

c       create concentration at the edge

c	do i = 2, nsur
c        x(i) = (1.0-cos(PI*x(i))) / 2.0
c	enddo


C	Compute the mean camber line coordinates

	sep = p*chord/10
	mc = m/100
	pm = p/10

	do i = 1, nsur+1
		if (x(i).lt.sep) then
			yc(i) = mc*(2*pm*x(i)-x(i)**2)/(pm**2)
            dyc(i) = 2.0*mc*(pm-x(i))/(pm**2)
		else
			yc(i) = mc*((1-2*pm)+2*pm*x(i)-x(i)**2)/((1-pm)**2)
            dyc(i) = 2.0*mc*(pm-x(i))/((1-pm)**2)
		endif
	enddo

C	Calculate the thickness distribution

	do i = 1, nsur+1
		yt(i) = 0.01*t*(0.2969*sqrt(x(i))-0.1260*x(i)-0.3516*x(i)**2+
     f				0.2843*x(i)**3-0.1015*x(i)**4)/0.2
	enddo

C	Determine the final coordinates

	do i = 1, nsur+1
		theta(i) = atan(dyc(i))
	enddo

c		upper surface

	do i = 1, nsur+1                               
		xu(i) = x(i) - yt(i)*sin(theta(i))
		yu(i) = -6.77928*xu(i)**6 + 22.8241*xu(i)**5 - 30.2556*xu(i)**4 + 
     f          19.8253 * xu(i)**3 - 6.88211*xu(i)**2 + 1.26759*xu(i)
	enddo

c		lower surface

	do i = 1, nsur+1
		xl(i) = x(i) + yt(i)*sin(theta(i)) 
		yl(i) = 3.08687*xl(i)**6 - 12.0141*xl(i)**5 + 17.0528*xl(i)**4 - 
     f          11.4338 * xl(i)**3 + 3.9684*xl(i)**2 - 0.66017*xl(i)

	enddo

C	Generate the output file

	do i = 1, nsur+1
		xa(i) = xl(nsur+2-i)
		ya(i) = yl(nsur+2-i)
		xa(i+nsur) = xu(i)
		ya(i+nsur) = yu(i)
	enddo

	open (unit=8, file='naca.txt')

	write(8,1000)
	write(8,1050) nmax, nwake

	write(8,1100) nmax+1

	count = 0
	column = 6
	counta = nmax/column
	
	do i = 0, counta
		write(8,1150) (xa(column*(i)+j), j = 1,column)
	    count = count + column
	enddo

	count = 0

	write(8,1200) nmax+1
	
	do i = 0, counta
		write(8,1150) (ya(column*(i)+j), j = 1,column)
		count = count + column
	enddo
c	nb = nmax - count
c	write(8,1150) (ya(count+1+j), j = 1,nb)

	write(8,1250)
	write(8,1150) aoa
	write(8,1300)
	write(8,1150) mach

	close(unit=8)

	stop

 1000 format ('  NODTOT   NWAKE')
 1050 format ('   ',i4,'   ',i4)
 1100 format (' X(I), I=1...   ',i4,'            Airfoil x-coordinate')
 1150 format (6f10.6)
 1200 format (' Y(I), I=1...   ',i4,'            Airfoil y-coordinate')
 1250 format (' angle of attack')
 1300 format (' Mach number')

	end
