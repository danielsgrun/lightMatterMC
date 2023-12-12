! PROGRAMA PARA TRABALHO FINAL DE METODOS COMPUTACIONAIS DA FISICA B
! BY: DANIEL S GRUN

! RUDIMENTOS DE SIMULACOES BIDIMENSIONAIS DE N PARTICULAS

real function u(a,b)  ! FUNCAO ACELERACAO
real a,b
real, parameter :: g = 1e-3  ! PARAMETRO PARA EVITAR SINGULARIDADE
u = -4*3.14**2/((a**2 + b**2 + g)**(3./2.))
end function u

real function rand(p) ! GERADOR DE CONGRUENCIA PARA COND. INICIAIS
implicit none
integer, parameter :: a = 1029, m = 1048576, c = 221591
integer :: p
!integer :: seed
integer, save :: seed
if (p .ne. 0) seed = p
seed = mod((a*seed+c),m)
rand = float(seed)/m

end function rand


program grav ! PROGRAMA PRINCIPAL
  implicit none
  real,dimension(:),allocatable :: x,y,vx,vy,ax,ay,mass,a
  real,dimension(:),allocatable :: x1,y1,vx1,vy1,ax1,ay1
  real,dimension(:),allocatable :: x2,y2,vx2,vy2,ax2,ay2
  real,dimension(:),allocatable :: x3,y3,vx3,vy3,ax3,ay3
  real,dimension(:),allocatable :: x4,y4,vx4,vy4,ax4,ay4
  real,external :: u,rand 
  real w1,w2,w,x_aux,y_aux,dt,t,tmax,rmax,r,r1,pp,v
  integer i,j,k,n,p,contador,m,c
  open(unit=1,file="inic.dat")
  open(unit=2,file="evol.dat")
  open(unit=3,file="final.dat")
  
  read*, n,tmax,p

  contador = 0
  dt = 1e-04    ! MODIFICAR CASO O NUMERO DE PARTICULAS SEJA > 5
  t = 0
  m = 1048576
  c = 221591

  allocate(x(n))
  allocate(x1(n))
  allocate(x2(n))
  allocate(x3(n))
  allocate(x4(n))
                   ! DIMENSIONANDO OS ARRAYS DE POSICAO,
  allocate(y(n))   ! VELOCIDADE, ACELERACAO E MASSA
  allocate(y1(n))
  allocate(y2(n))
  allocate(y3(n))
  allocate(y4(n))
    
  allocate(vx(n))
  allocate(vx1(n))
  allocate(vx2(n))
  allocate(vx3(n))
  allocate(vx4(n))

  allocate(vy(n))
  allocate(vy1(n))
  allocate(vy2(n))
  allocate(vy3(n))
  allocate(vy4(n))

  allocate(ax(n))
  allocate(ax1(n))
  allocate(ax2(n))
  allocate(ax3(n))
  allocate(ax4(n))

  allocate(ay(n))
  allocate(ay1(n))
  allocate(ay2(n))
  allocate(ay3(n))
  allocate(ay4(n))

  allocate(a(n))

  allocate(mass(n))

  call system('rm *.png')
  call system('rm *.jpg')
  call system('rm *.gif')  

  rmax = 5*n
  pp = rand(p)

  do i=1,n  

    w1 = rand(0) + float(1)/float(m)
    r1 = rmax*w1

    w1 = rand(0)
    w = 2*int(2*w1) - 1

    w1 = rand(0) + float(1)/float(m)
    x(i) = w*r1*(w1)
    
    w1 = rand(0)
    w = 2*int(2*w1) - 1    

    y(i) = w*sqrt(r1**2 - x(i)**2)

    write(1,*) x(i),y(i)
    
    w1 = rand(0) 
    mass(i) = 1e+2*(w1 + 0.01)

    ax(i) = 0
    ay(i) = 0
    a(i) = 0

    do j=1,n

      x_aux = x(i)-x(j)
      y_aux = y(i)-y(j) 
      ax(i) = ax(i) + u(x_aux,y_aux)*x_aux*mass(j)
      ay(i) = ay(i) + u(x_aux,y_aux)*y_aux*mass(j)      

    enddo

    a(i) = sqrt(ax(i)**2 + ay(i)**2)    
    v = sqrt(r1*a(i))
    vx(i) = - v*y(i)/r1   
    vy(i) = v*x(i)/r1
    
    write(1,*) x(i),y(i)

  enddo
  
  do while (t .le. tmax)            
    
    open(unit=4,file='teste.dat')

    do i=1,n 

      ax1(i) = 0
      ay1(i) = 0
      ax2(i) = 0
      ay2(i) = 0
      ax3(i) = 0
      ay3(i) = 0
      ax4(i) = 0
      ay4(i) = 0    

      vx1(i) = vx(i)
      vy1(i) = vy(i)
      x1(i) = x(i)
      y1(i) = y(i) 

      do j=1,n
        if (j .ne. i) then
          x_aux = x1(i)-x1(j)
          y_aux = y1(i)-y1(j)
          ax1(i) = ax1(i) + u(x_aux,y_aux)*x_aux*mass(j)
          ay1(i) = ay1(i) + u(x_aux,y_aux)*y_aux*mass(j)
        endif
      enddo

      x2(i) = x(i) + vx1(i)*dt/2.
      y2(i) = y(i) + vy1(i)*dt/2.
    
      vx2(i) = vx(i) + ax1(i)*dt/2
      vy2(i) = vy(i) + ay1(i)*dt/2
    
      do j=1,n
        if (j .ne. i) then
          x_aux = x2(i)-x2(j)
          y_aux = y2(i)-y2(j)
          ax2(i) = ax2(i) + u(x_aux,y_aux)*x_aux*mass(j)
          ay2(i) = ay2(i) + u(x_aux,y_aux)*y_aux*mass(j)
        endif
      enddo

      x3(i) = x(i) + vx2(i)*dt/2.
      y3(i) = y(i) + vy2(i)*dt/2.
    
      vx3(i) = vx(i) + ax2(i)*dt/2
      vy3(i) = vy(i) + ay2(i)*dt/2
    
      do j=1,n
        if (j .ne. i) then
          x_aux = x3(i)-x3(j)
          y_aux = y3(i)-y3(j)
          ax3(i) = ax3(i) + u(x_aux,y_aux)*x_aux*mass(j)
          ay3(i) = ay3(i) + u(x_aux,y_aux)*y_aux*mass(j)
        endif
      enddo
      
      x4(i) = x(i) + vx3(i)*dt/2.
      y4(i) = y(i) + vy3(i)*dt/2.
    
      vx4(i) = vx(i) + ax3(i)*dt/2.
      vy4(i) = vy(i) + ay3(i)*dt/2.
    
      do j=1,n
        if (j .ne. i) then
          x_aux = x4(i)-x4(j)
          y_aux = y4(i)-y4(j)
          ax4(i) = ax4(i) + u(x_aux,y_aux)*x_aux*mass(j)
          ay4(i) = ay4(i) + u(x_aux,y_aux)*y_aux*mass(j)
        endif
      enddo      
    
      x(i) = x(i) + 1./6.*(vx1(i) + 2*(vx2(i)+vx3(i)) + vx4(i))*dt
      y(i) = y(i) + 1./6.*(vy1(i) + 2*(vy2(i)+vy3(i)) + vy4(i))*dt
      vx(i) = vx(i) + 1./6.*(ax1(i) + 2*(ax2(i) + ax3(i)) + ax4(i))*dt
      vy(i) = vy(i) + 1./6.*(ay1(i) + 2*(ay2(i) + ay3(i)) + ay4(i))*dt    

      write(4,*) x(i),y(i)         
      if (t .gt. tmax - dt) then
        write(3,*) x(i),y(i)     
      endif
    enddo

    t = t + dt

    write(2,*) x(1),y(1),x(2),y(2),x(3),y(3), &
&     x(4),y(4),x(5),y(5),x(6),y(6), &
&     x(7),y(7),x(8),y(8),x(9),y(9),x(10),y(10)

    contador = contador + 1
    if (mod(contador,7) == 0) then   
      call system ('python plot.py') ! CHAMA O PLOT.PY PARA REALIZAR O
    end if                           ! PLOT DAS PARTICULAS A CADA
                                     ! INSTANTE
    close(unit=4)

  enddo

!  call system ('convert -delay 3 -loop 0 *.jpg simulation.gif') ! DESCOMENTAR
!  call system ('rm *.jpg')                                 ! PARA FAZER UM 
!  call system ('eog simulation.gif')                       ! GIF COM AS
                                                            ! IMAGENS PLOTADAS
end program grav



    



    

    
  
  
    
  
