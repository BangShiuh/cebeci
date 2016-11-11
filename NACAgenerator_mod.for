c	the goal of this program is to generate the shape of a
c	NACA airfoil from the NACA four-digit series given the
c	four numbers of its definition

c   Note: the input data should have the order of 
c   lower trailing edge to leading edge
c   upper leading edge to trailing edge 
c   source /opt/intel/bin/compilervars.sh intel64


	Program main
    
	integer nmax, nsur, i, nwake, column, counta, count
	real chord,m,p,t,h,x,yc,yt,dyc,theta,xu,yu,xl,yl,mc,pm,aoa,mach
	dimension x(1000), yc(1000), yt(1000), dyc(1000), theta(1000),
     f		xu(1000), yu(1000), xl(1000), yl(1000), xa(1000), ya(1000)
      character*80 input_name

c     enter the file's name
      write(*,*) 'enter the files name'
      read(5,*) input_name

c     enter how many points of data
      write(*,*) 'enter how many points of data'
      read(*,*) nmax

	  write(6,*) "Enter the angle of attack"
	  read(5,*) aoa
	  write(6,*) "Enter the mach number"
	  read(5,*) mach

c     not sure what is nwake 
      nwake = 0

c     read the file from here
      open(unit=7, file=input_name,status='OLD')
	      do i = 1, nmax
		      read(7,2000) xa(i), ya(i)
	      enddo
      CLOSE(unit=7)
 2000 format(2f10.5)

C	Generate the output file
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

 1000 format ('  NOTDOT   NWAKE')
 1050 format ('   ',i4,'   ',i4)
 1100 format (' X(I), I=1...   ',i4,'            Airfoil x-coordinate')
 1150 format (' ',6f10.6)
 1200 format (' Y(I), I=1...   ',i4,'            Airfoil y-coordinate')
 1250 format (' angle of attack')
 1300 format (' Mach number')

	end
