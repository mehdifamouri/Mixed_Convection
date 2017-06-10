	IMPLICIT NONE
	INTEGER :: i,j,iteration,ii,jj,l,nx,ny,Iterationpp,qq
      REAL :: u,v,u1,v1,u2,v2,p,p1,pp1,pp2,x,xx,yy,y,t,t1,sy,hf
	/,re,pr,Ra,Pe,gr,
	/ae,aw,an,as,ap,fe,fw,fn,fs,de,dw,dn,ds,su,sv,st,spp,di,dj,
     /qu,qv,qp,qt,qpp,dx0,dy0,dx,dy,dxi,dxii,dyj,dyjj,
     /bRes,aCf,aNu,al,
	/ux,uy,vx,vy,tx,ty,FFI,HTI,Sg,BE,aOmg,Br,
	/SumFFI,SumHTI,SumSG,SumBe,
	/b1,b2,udif,vdif,pdif,tdif,difpp,gu1,gv1,gpp1,gp1,gt1

	dimension u(150,150),u1(150,150),u2(150,150)
	/,v(150,150),v1(150,150),v2(150,150)
	/,p(150,150),p1(150,150)
     /,pp1(150,150),pp2(150,150)
     /,di(150,150),dj(150,150)
     /,t(150,150),t1(150,150)
	/,sy(150,150),hf(150,150)
	/,FFI(150,150),HTI(150,150),Sg(150,150),BE(150,150)
	/,x(150),y(150),xx(150),yy(150)


c========== Input Data Ra Re Nx Ny =======================	
	write(*,*)"RE,Pr,Ra,GAMA="
	read(*,*)re,pr,Ra,al
	write(*,*)"Nx,ny="
	read(*,*)b1,b2
c      re=100.
	Ra=10.**(Ra)
c	Pr=.7
	al=0.
	Pe=re*pr
      Gr=ra/pr
      al=al*3.14/180.
	aOmg=.1
	Br=.1
	nx=101*b1
	ny=75*b2
      gu1=.3
	gv1=.3
	gpp1=1.
	gp1=.3
	gt1=.1
	qu=5.
	qv=5.
	qp=5.
	qt=5.
	qpp=6.
c=======     X     Y     dx     dy  ======================      
	dx0=1./(nx-1)
      dy0=dx0*1.

	dx=dx0
 	dy=dy0

	write(*,93)Re,Pr,Ra,Gr
93    format("  Re=",f6.0,"    Pr=",f7.3,"     Ra=",f7.0,"    Gr=",f8.0)
	x(0)=0.
	x(1)=dx0/2.
	do i=2,nx-1
     	x(i)=x(i-1)+dx0
	end do
	x(nx)=x(nx-1)+dx0/2.

	y(0)=0.
	y(1)=dy0/2.
	do j=2,ny-1
     	y(j)=y(j-1)+dy0
c	write(*,*)j,y(j)
	end do
	y(ny)=y(ny-1)+dy0/2.

	dx=dx0
	dy=dy0
	nx=nx
	ny=ny
	xx=x
	yy=y
	write(*,94)nx,ny,dx,dy
94    format("  Nx=",i6,"    Ny=",i6,"      dX=",f7.4,"    dY=",f7.4)
c========  Initial Guess    ===============================
	u=0.7
	u1=0.7
	v=.1
	v1=.1

	p1=0.
	p=0.
	t=0.0
	t1=0.0	
      do i=1,nx
	t(i,j)=.5
	t1(i,j)=.5
	end do
	i =1
      do j=2,ny
	u(i,j)=1.
	u1(i,j)=1.
	end do
	bRes=1
c===============================================
c=========   START    START  ===================
c===============================================
222	iteration=iteration+1
	if (1.*iteration/1000.eq.int(1.*iteration/1000)) then
	write(*,92)iteration/1000,ii,jj,
     /udif,vdif,pdif,tdif,aNu,aCf,abs(bRes)
92	format(" >> ",i5,i3,i3,7e9.3)
	write(*,*)
	if (iteration.lt.30)goto 300
c==============================================================================================================
c=================== Post Proccesing ==========================================================================
c=================== Post Proccesing ==========================================================================
c==========================================================================================================
 	open(1,file="1 Velocity.plt")
	open(2,file="2 Pressure.plt")
	open(3,file="3 Tempearture.plt")
	open(4,file="4 Nusselt.txt")
	open(5,file="5 Firection Coefficient.txt")
	open(6,file="6 Stream Function.plt")
	open(7,file="7 Heat Function.plt")
	open(8,file="8 FFI.plt")
	open(9,file="9 HTI.plt")
	open(10,file="10 Sg.plt")
	open(11,file="11 Be.plt")
	open(12,file="12 Details.txt")


	write(1,*)'Zone  i=',ny-1,'    j=',nx-1
	write(2,*)'Zone  i=',ny-1,'    j=',nx-1
	write(3,*)'Zone  i=',ny-1,'    j=',nx-1
	write(6,*)'Zone  i=',ny-1,'    j=',nx-1
	write(7,*)'Zone  i=',ny-1,'    j=',nx-1
	write(8,*)'Zone  i=',ny-1,'    j=',nx-1
	write(9,*)'Zone  i=',ny-1,'    j=',nx-1
	write(10,*)'Zone  i=',ny-1,'    j=',nx-1
	write(11,*)'Zone  i=',ny-1,'    j=',nx-1

	do i=1,nx-1
	do j=1,ny-1
	write(1,*)x(i),y(j),u(i,j),v(i,j)
	write(3,*)x(i),y(j),t(i,j)
	write(6,*)x(i),y(j),sy(i,j)
	write(7,*)x(i),y(j),Pe*u(i,j)*t(i,j)-(t(i+1,j)-t(i-1,j))/(2*dx)
	/,Pe*v(i,j)*t(i,j)-(t(i,j+1)-t(i,j-1))/(2*dx)
	write(8,*)x(i),y(j),FFI(i,j)
	write(9,*)x(i),y(j),HTI(i,j)
	write(10,*)x(i),y(j),Sg(i,j)
	write(11,*)x(i),y(j),Be(i,j)
	if((i.eq.0).or.(i.eq.nx).or.(j.eq.0).or.(j.eq.ny)) goto 4
	write(2,*)x(i),y(j),p(i,j)
4	end do
	end do

	close(1)
	close(2)
	close(3)
	close(6)
	close(7)
	close(8)
	close(9)
	close(10)
	close(11)


c======== Nusselt======================================
	aNu=0.
      do i=2,nx
	b1=1./t(i,0)
	write(4,*)x(i),b1
      aNu=aNu+b1*dx
	end do			  
	write(4,*)
	write(4,*)"Nusselt Average=",aNu

c========== Firection Factor===========================
	aCf=0.
      do i=1,nx
	dyj=dy/2
	dyjj=dy+dy/2.
	b1=(-u(i,2)*dyj**2+u(i,1)*dyjj**2-u(i,0)*(dyjj**2-dyj**2))
	//(dyj*dyjj*(dyjj-dyj))
      aCf=aCf+b1*dx/(Re/2)
	write(5,*)x(i),b1/(Re/2)
	end do
	write(5,*)
	write(5,*)"Cf Average=",acf
	write(5,*)
	write(5,*)"====== Bondry Layer REsults =============="
      do ii=2,nx
     	write(5,*)x(ii),.664*(re*x(ii))**(-.5)
      end do
      close(4)
	close(5)
c========= Stream Function================================
c========= Stream Function================================
	do ii=1,100
      do i=1,nx-1
	do j=1,ny-1
	
      b1=(v(i+1,j)-v(i,j))/(dx)
	b1=b1+(u(i,j+1)-u(i,j))/(dy)
	sy(i,j)=(sy(i,j+1)+sy(i,j-1)+sy(i+1,j)+sy(i-1,j)-b1*dy**2)/4

	end do
	end do
c=============Botom wall
	j=0
	do i=1,nx
	sy(i,j)=0.
      end do
c=============Top Boundary
	j=ny
	do i=1,nx
	sy(i,j)=sy(i,j-1)*0
      end do
c=============Inlet 
	i=1
	do j=1,ny
c	sy(i,j)=(j-1.)*dy
	sy(i,j)=sy(i+1,j)
      end do
c=============Outlet 
	i=nx
	do j=1,ny
	sy(i,j)=sy(i-1,j)
      end do

	end do 
c========= Heat Function================================
c========= Heat Function================================
c	do ii=1,10000
c      do i=1,nx-1
c	do j=1,ny-1
	
      b1=(v(i+1,j)*t(i+1,j)-u(i-1,j)*t(i-1,j))/(2.*dx)
	b1=b1+(u(i,j+1)*t(i,j+1)-u(i,j-1)*t(i,j+1))/(2.*dy)
	b1=b1*Pe
	hf(i,j)=(hf(i,j+1)+hf(i,j-1)+hf(i+1,j)+hf(i-1,j)-b1*dy**2)/4

c	end do
c	end do
c=============Botom wall
c	j=0
c	do i=1,nx
c	hf(i,j)=0.
c      end do
c=============Top Boundary
	j=ny
c	do i=1,nx
c	hf(i,j)=hf(i,j-1)
c      end do
c=============Inlet 
	i=1
c	do j=1,ny
c	sy(i,j)=(j-1.)*dy
	hf(i,j)=hf(i+1,j)
c      end do
c=============Outlet 
	i=nx
c	do j=1,ny
	hf(i,j)=hf(i-1,j)
c      end do

c	end do 
c=============== ENTROPY GENERATION ==============
c=============== ENTROPY GENERATION ==============
      do i = 1,nx-1
 	do j = 1,ny-1

      tx=(t(i+1,j)-t(i-1,j))/(2*dx)
	ty=(t(i,j+1)-t(i,j-1))/(2*dy)

	ux=(u(i+1,j)-u(i,j))/dx
	uy=(u(i,j+1)-u(i,j))/dy

	vx=(v(i+1,j)-u(i,j))/dx
	vy=(v(i,j+1)-v(i,j))/dy
            
	HTI(i,j)=tx**2+ty**2
      FFI(i,j)=2*((ux)**2+(vy)**2)+(uy+vx)**2
 	end do
	end do
c=============Botom wall
	j=1
	do i=1,nx
	FFI(i,j)=FFI(i,j+1)
	HTI(i,j)=HTI(i,j+1)
      end do
c=============Top Boundary
	j=ny
	do i=1,nx
	FFI(i,j)=FFI(i,j-1)
	HTI(i,j)=HTI(i,j-1)
      end do
c=============Inlet 
	i=1
	do j=1,ny
	FFI(i,j)=FFI(i+1,j)
	HTI(i,j)=HTI(i+1,j)
      end do
c=============Outlet 
	i=nx
	do j=1,ny+1
	FFI(i,j)=FFI(i-1,j)
	HTI(i,j)=HTI(i-1,j)
      end do
c================Dimensionless Parameter
      do i=1,nx-1
	do j=1,ny-1
	FFI(i,j)=Br*FFI(i,j)/(1.+aOmg*t(i,j))
	HTI(i,j)=aOmg**2*HTI(i,j)/(1.+aOmg*t(i,j))**2
	Sg(i,j)=FFI(i,j)+HTI(i,j)
	Be(i,j)=HTI(i,j)/(FFI(i,j)+HTI(i,j))
      end do							
	end do
c================Sum FFI HTI ...
	SumFFI=0.
 	SumHTI=0.
 	SumBe=0.
      do i=1,nx-1
	do j=1,ny-1
	dx=x(i+1)-x(i)
	dy=y(j+1)-y(j)
	b1=(FFI(i,j)+FFI(i+1,j)+FFI(i,j+1)+FFI(i+1,j+1))
	b2=(HTI(i,j)+HTI(i+1,j)+HTI(i,j+1)+HTI(i+1,j+1))
	SumFFI=SumFFI+.25*b1*dx*dy
	SumHTI=SumHTI+.25*b2*dx*dy
	b1=(Be(i,j)+Be(i+1,j)+Be(i,j+1)+Be(i+1,j+1))
	SumBe=SumBe+.25*b1*dx*dy
     	end do
	end do
	SumSG=SumFFI+SumHTI


c================Details=================
	write(12,94)nx,ny,dx,dy
	write(12,95)qu,qv,qt,qp
	write(12,*)" Iteration=",iteration
      write(12,*)" Resiual b=",bRes
	write(12,*)"======== Iso Flux ========================="
	write(12,*)"Agel=",al*180/3.14 
      write(12,*)"  Pr=",pr
	write(12,*)"  Pe=",Pe
	write(12,*)"  Gr=",Gr
 	write(12,*)"  Ra=",Ra
	write(12,*)" RaX=",Ra*sin(al)
	write(12,*)" RaY=",Ra*cos(al)
	write(12,*)"===================================="
	write(12,*)"  Cf=",acf
	write(12,*)"  Nu=",aNu
      write(12,*)"==================================="
	write(12,*)"Omega=",aOmg
	write(12,*)"   Br=",Br
	write(12,*)"  FFI=",SumFFI
	write(12,*)"  HTI=",SumHTI
	write(12,*)"   Sg=",SumSg
	write(12,*)"   Be=",SumBe


     
      close(12)
95    format("   u-Er= ",f3.0,"   v-Er= ",
     /f3.0,"   T-Er= ",f3.0,"  P-Er= ",f3.0)
c=============================================================================================
	end if
c=============================================================================================
c=============================================================================================
300   l=l	
	if(bRes.lt.1e-3)then
	gu1=.1
	gv1=.1
	gp1=.1
	gt1=.1
	end if
	dx=dx0
	dy=dy0
c===============   UUUUU   =======================
c===============   UUUUU   =======================
	do i=2,nx-1
	do j=1,ny-1
      Fw=(u1(i,j)+u1(i-1,j))/2.*dy
	Fe=(u1(i+1,j)+u1(i,j))/2.*dy
	Fs=(v1(i,j)+v1(i-1,j))/2.*dx
	Fn=(v1(i,j+1)+v1(i-1,j+1))/2.*dx

	Dw=(1./Re)/(dx)*dy
	De=(1./Re)/(dx)*dy
	Ds=(1./Re)/(dy)*dx
	Dn=(1./Re)/(dy)*dx
	
	if (j.eq.1)Ds=Ds*2.
 	if (j.eq.ny-1)Dn=Dn*2.

	su=-(p1(i,j)-p1(i-1,j))*dy
	if(iteration.gt.7000)su=su+Ra*sin(al)/(Pe*Re)*(t(i,j)-t(i-1,j))*dx

      ae=-Fe/2.+De
	aw=+Fw/2.+Dw
	an=-Fn/2.+Dn
	as=+Fs/2.+Ds
	ap=-Fw/2.+Fe/2.-Fs/2.+Fn/2.+De+Dw+Dn+Ds
	b1=ae*u1(i+1,j)+aw*u1(i-1,j)+an*u1(i,j+1)+as*u1(i,j-1)+su
      u2(i,j)=u1(i,j)+gu1*(-u1(i,j)+(b1/ap))
	di(i,j)=dy/(ap)

	end do
	end do
c=============Botom wall
	j=0
	do i=1,nx
	u2(i,j)=0.
      end do
c=============Top Boundary
	j=ny
	do i=1,nx
	u2(i,j)=u2(i,j-1)
      end do
c=============Inlet 
	i=1
	do j=1,ny
	u2(i,j)=1.
      end do
c=============Outlet 
	i=nx
	do j=1,ny
	u2(i,j)=u2(i-1,j)
      end do
c===============   VVVVV   ===================================
c===============   VVVVV   ===================================
	do i=1,nx-1
	do j=2,ny-1

      Fw=(u1(i,j)+u1(i,j-1))/2.*dy
	Fe=(u1(i+1,j)+u1(i+1,j-1))/2.*dy
	Fs=(v1(i,j-1)+v1(i,j))/2.*dx
	Fn=(v1(i,j)+v1(i,j+1))/2.*dx

	Dw=(1./Re)/(dx)*dy
	De=(1./Re)/(dx)*dy
	Ds=(1./Re)/(dy)*dx
	Dn=(1./Re)/(dy)*dx

	if (i.eq.1)Dw=Dw*2.
	if (i.eq.nx-1)De=De*2.
      
	sv=-(p1(i,j)-p1(i,j-1))*dy
	if(iteration.gt.7000)sv=sv+Ra*cos(al)/(Pe*Re)*(t(i,j)-t(i,j-1))*dy

      ae=-Fe/2.+De
	aw=+Fw/2.+Dw
	an=-Fn/2.+Dn
	as=+Fs/2.+Ds
	ap=-Fw/2.+Fe/2.-Fs/2.+Fn/2.+De+Dw+Dn+Ds
	b1=ae*v1(i+1,j)+aw*v1(i-1,j)+an*v1(i,j+1)+as*v1(i,j-1)+sv
	v2(i,j)=v1(i,j)+gv1*(-v1(i,j)+(b1/ap))

	dj(i,j)=dx/(ap)

	end do
	end do
c=============Botom wall
	j=1
	do i=1,nx
	v2(i,j)=0.
      end do
c=============Top Boundary
	j=ny
	do i=1,nx
	v2(i,j)=v2(i,j-1)
      end do
c=============Inlet 
	i=0
	do j=1,ny
	v2(i,j)=0.
      end do
c=============Outlet 
	i=nx
	do j=2,ny
	v2(i,j)=v2(i-1,j)
      end do
c===============   di  di   =======================
c===============   di  di   =======================
	do i=1,nx-1
	do j=1,ny-1

      Fw=(u2(i,j)+u2(i-1,j))/2.*dy
	Fe=(u2(i+1,j)+u2(i,j))/2.*dy
	Fs=(v2(i,j)+v2(i-1,j))/2.*dx
	Fn=(v2(i,j+1)+v2(i-1,j+1))/2.*dx

	Dw=(1./Re)/(dx)*dy
	De=(1./Re)/(dx)*dy
	Ds=(1./Re)/(dy)*dx
	Dn=(1./Re)/(dy)*dx
 
	if (j.eq.1)Ds=Ds*2.
 	if (j.eq.ny-1)Dn=Dn*2.

      ae=-Fe/2.+De
	aw=+Fw/2.+Dw
	an=-Fn/2.+Dn
	as=+Fs/2.+Ds
	ap=-Fw/2.+Fe/2.-Fs/2.+Fn/2.+De+Dw+Dn+Ds
    
   	di(i,j)=dy/(ap)

	end do
	end do
c===============  dj  dj   ===================================
c===============  dj  dj   ===================================
	do i=1,nx-1
	do j=1,ny-1

      Fw=(u2(i,j)+u2(i,j-1))/2.*dy
	Fe=(u2(i+1,j)+u2(i+1,j-1))/2.*dy
	Fs=(v2(i,j-1)+v2(i,j))/2.*dx
	Fn=(v2(i,j)+v2(i,j+1))/2.*dx

	Dw=(1./Re)/(dx)*dy
	De=(1./Re)/(dx)*dy
	Ds=(1./Re)/(dy)*dx
	Dn=(1./Re)/(dy)*dx

	if (i.eq.1)Dw=Dw*2.
	if (i.eq.nx-1)De=De*2.
	

      ae=-Fe/2.+De
	aw=+Fw/2.+Dw
	an=-Fn/2.+Dn
	as=+Fs/2.+Ds
	ap=-Fw/2.+Fe/2.-Fs/2.+Fn/2.+De+Dw+Dn+Ds

	dj(i,j)=dx/(ap)

	end do
	end do

c========   Perusser Coorection   =============
c========   Perusser Coorection   =============
	
      pp1=0.
	pp2=0.
	bRes=0.
	Iterationpp=0
111	Iterationpp=Iterationpp+1
      do i=1,nx-1
	do j=1,ny-1
	ae=di(i+1,j)*dy
	aw=di(i,j)*dy
	an=dj(i,j+1)*dx
	as=dj(i,j)*dx
	if (i.eq.1) aw=0.
	if (i.eq.nx-1) ae=0.
	if (j.eq.1) as=0.
	if (j.eq.ny-1) an=0.
     	ap=ae+aw+an+as
	spp=+u2(i,j)*dy-u2(i+1,j)*dy+v2(i,j)*dx-v2(i,j+1)*dx
	b1=ae*pp1(i+1,j)+aw*pp1(i-1,j)+an*pp1(i,j+1)+as*pp1(i,j-1)+spp
	if (bRes.lt.abs(spp)) bRes=abs(spp)
      pp2(i,j)=pp1(i,j)+gpp1*(-pp1(i,j)+b1/(ap))

      end do
	end do
c	pp2(1,1)=pp2(1,2)
c	pp2(2,1)=pp2(2,2)

	if (Iterationpp.gt.15)goto 11
c=============Error 
	do i=2,nx-1
	do j=2,ny-1
      if (abs(pp2(i,j)).ge.1e-9) then
      difpp=Abs((pp2(i,j)-pp1(i,j))/pp2(i,j))/gpp1
	else
	difpp=0.
	end if
      If (difpp-10.**(-qpp) .gt.0.0)goto 1
      end do
	end do
	goto 11
c=============Repelace 
1	do i=1,nx
	do j=1,ny
	pp1(i,j)=pp2(i,j)
	end do
	end do
	goto 111
11     L=L
c========   Correct Velocity&Pressure   =============
c========   Correct Velocity&Pressure   =============
	do i=1,nx-1
	do j=1,ny-1
	p(i,j)=p1(i,j)+gp1*pp2(i,j)
      end do
	end do

	do i=1,nx
	do j=1,ny
	u2(i,j)=u2(i,j)+di(i,j)*(pp2(i-1,j)-pp2(i,j))
	v2(i,j)=v2(i,j)+dj(i,j)*(pp2(i,j-1)-pp2(i,j))
	
	u(i,j)=u1(i,j)+1*(-u1(i,j)+u2(i,j))
	v(i,j)=v1(i,j)+1*(-v1(i,j)+v2(i,j))
      end do		
	end do
c=============Botom wall
	j=1
	do i=1,nx
c	p(i,j)=p(i,j+1)
	u(i,j-1)=0.
	v(i,j)=0.
	v(i,j-1)=0.
      end do
c=============Top Boundary
	j=ny
	do i=1,nx
c	p(i,j)=1.
	u(i,j)=u(i,j-1)
	v(i,j)=v(i,j-1)
      end do
c=============Inlet 
	i=1
	do j=0,ny
c	p(i,j)=1.
	u(i,j)=1.
	u(i-1,j)=1.
	v(i-1,j)=0.
       end do
c=============Outlet 
	i=nx
	do j=1,ny
c	p(i,j)=p(i-1,j)
	u(i,j)=u(i-1,j)
	v(i,j)=v(i-1,j)
       end do
c==================================
c	p(1,1)=1.
c	p(2,1)=p(2,2)
c===============   TTTTTTTT   ===================================
c===============   TTTTTTTT   ===================================
	Pe=re*pr
 	if(iteration.lt.10000) Pe=re*pr/1.2  
      if(iteration.lt.9000) Pe=re*pr/1.5 
      if(iteration.lt.5000) Pe=re*pr/2
      if(iteration.lt.4000) Pe=re*pr/3 
      if(iteration.lt.3000) Pe=re*pr/6 
      if(iteration.lt.2000) Pe=re*pr/10 
	if(Pr.lt.3)Pe=re*pr

	do i=1,nx-1
	do j=1,ny-1

      Fw=u(i,j)*dy
	Fe=u(i+1,j)*dy
	Fs=v(i,j)*dx
	Fn=v(i,j+1)*dx

	Dw=(1./Pe)/(dx)*dy
	De=(1./Pe)/(dx)*dy
	Ds=(1./Pe)/(dy)*dx
	Dn=(1./Pe)/(dy)*dx

	if (j.eq.1)Ds=Ds*2.
	if (j.eq.ny-1)Dn=Dn*2.
	if (i.eq.1)Dw=Dw*2.
	if (i.eq.nx-1)De=De*2.
	
      st=0.
	if (j.eq.1)st=dx/Pe
	if (j.eq.1)Ds=0.

      ae=-Fe/2.+De
	aw=+Fw/2.+Dw
	an=-Fn/2.+Dn
	as=+Fs/2.+Ds
	ap=-Fw/2.+Fe/2.-Fs/2.+Fn/2.+De+Dw+Dn+Ds
	b1=ae*t1(i+1,j)+aw*t1(i-1,j)+an*t1(i,j+1)+as*t1(i,j-1)+st
	t(i,j)=t1(i,j)+gt1*(-t1(i,j)+(b1/ap))
	end do
	end do
c=============Botom wall
	j=0
	do i=1,nx
	t(i,j)=t(i,j+1)+dy/2.
      end do
c=============Top Boundary
	j=ny
	do i=1,nx
	t(i,j)=0.
      end do
c=============Inlet 
	i=0
	do j=1,ny
	t(i,j)=0.
      end do
c=============Outlet 
	i=nx
	do j=0,ny
	t(i,j)=t(i-1,j)
      end do
     
c==========  Restart     Reastart   ==================
c== ========  Restart     Reastart   ==================
32    l=l
      do i = 1,nx
     	do j = 1,ny

	if (abs(u(i,j)).ge.1e-9) then
      udif=Abs((u(i,j)-u1(i,j))/u(i,j))/gu1
	else
	udif=0.
	end if
	if (abs(v(i,j)).ge..1e-9)then
	vdif=Abs((v(i,j)-v1(i,j))/v1(i,j))/gv1
	else
	vdif=0.
	end if
	if (abs(p(i,j)).ge..1e-9)then
	pdif=Abs((p(i,j)-p1(i,j))/p(i,j))/gp1
	else
	pdif=0.
	end if
	if (abs(t1(i,j)).ge..1e-9)then
	tdif=Abs((t(i,j)-t1(i,j))/t(i,j))/gt1
	else
	tdif=0.
	end if

	ii=i
	jj=j
	
	If (udif-10.**(-qu) .gt.0.0 )goto 2
	If (vdif-10.**(-qv) .gt.0.0)goto 2
	If (pdif-10.**(-qp) .gt.0.0)goto 2
	If (tdif-10.**(-qt) .gt.0.0)goto 2


	end do
	end do
	goto 22
c=============Repelace 
2     l=l
      do i=0,nx
	do j=0,ny
	p1(i,j)=p(i,j)
	u1(i,j)=u(i,j)
	v1(i,j)=v(i,j)
	t1(i,j)=t(i,j)
	end do
	end do
	goto 222
22    l=l
c=====================================================
c 	open(1,file="1 Velocity.plt")
c	open(2,file="2 Pressure.plt")
c	write(1,*)'Zone  i=',ny,'    j=',nx
c	write(2,*)'Zone  i=',ny-1,'    j=',nx-1
c	do i=1,nx
c	do j=1,ny
c	write(1,*)x(i),y(j),u(i,j),v(i,j)
c	end do
c	end do
c
c	do i=1,nx-1
c 	do j=1,ny-1
c 	write(2,*)x(i),y(j),p(i,j)
c	end do
c	end do

c	close(2)
91    format(5f12.4)          
      stop

	end