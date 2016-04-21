SUBROUTINE F(k,beita,flag)
USE GLOBAL_VARIATION
IMPLICIT NONE
!-------------------------variation introduction--------------------------------------
!flag=0: stand for F(rou)
!flag=1: stands for F(u)
!flag=2: stands for F(v)
!flag=3: stands for F(p)
!k(m,n): is the intermediate term for LDDRK term and the subscript changes from 0~4
!and k is taken as both the input variations and output variations
!-------------------------------------------------------------------------------------
 integer::i,j,m,n,flag,stat_x,stat_y,j_temp
 real(kind=8)::k(m,n),k1(m,n)
 real(kind=8)::E1(m,n),F1(m,n),E2(14,n),F2(14,n),E3(m,14),F3(m,14)
 real(kind=8)::u0(7),v0(7),cos_thita,sin_thita,v_thita
 real(kind=8)::x_coor(m),y_coor(n)
 real(kind=8)::beita
 real(kind=8),external::DRP7
!
!
 do i=1,m
	 x_coor(i)=-x_range/2.0d0+(i-1)*x_range/(m-1)
 end do
!
 do j=1,n
	y_coor(j)=-y_range/2.0d0+(j-1)*y_range/(n-1)
 end do
!
!
!interior region
 if(flag==0)then
	 E1=Ma*(rou+beita*k)+u
 	 F1=v
 elseif(flag==1)then
	 E1=Ma*(u+beita*k)+p
 	 F1=0.0d0
 elseif(flag==2)then
	 E1=Ma*(v+beita*k)
 	 F1=p
 elseif(flag==3)then
 	 E1=Ma*(p+beita*k)+u
 	 F1=v
 end if
!inflow and outflow region
!inflow
 do i=1,7
	do j=1,n-6
		cos_thita=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2)
       	sin_thita=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2)
		v_thita=(Ma*cos_thita+dsqrt(1-Ma**2*sin_thita**2))
 		if(flag==0)then
			E2(i,j-3)=v_thita*cos_thita*(rou(i,j)+beita*k(i,j))
			F2(i,j-3)=v_thita*sin_thita*(rou(i,j)+beita*k(i,j))
		elseif(flag==1)then
			E2(i,j-3)=v_thita*cos_thita*(u(i,j)+beita*k(i,j))
			F2(i,j-3)=v_thita*sin_thita*(u(i,j)+beita*k(i,j))
 		elseif(flag==2)then
			E2(i,j-3)=v_thita*cos_thita*(v(i,j)+beita*k(i,j))
			F2(i,j-3)=v_thita*sin_thita*(v(i,j)+beita*k(i,j))
		else
			E2(i,j-3)=v_thita*cos_thita*(p(i,j)+beita*k(i,j))
			F2(i,j-3)=v_thita*sin_thita*(p(i,j)+beita*k(i,j))
		end if
	end do
 end do
!outflow
 do i=n-6,n
	do j=1,m-6
		cos_thita=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2)
        sin_thita=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2)
		v_thita=(Ma*cos_thita+dsqrt(1-Ma**2*sin_thita**2))
		if(flag==0)then
			E2(i-n+14,j-3)=Ma*(rou(i,j)+beita*k(i,j))+v_thita*cos_thita*p(i,j)-Ma*p(i,j)
			F2(i-n+14,j-3)=v_thita*sin_thita*p(i,j)
		elseif(flag==1)then
			E2(i-n+14,j-3)=Ma*(u(i,j)+beita*k(i,j))+p(i,j)
			F2(i-n+14,j-3)=0.0d0
		elseif(flag==2)then
			E2(i-n+14,j-3)=Ma*(v(i,j)+beita*k(i,j))+p(i,j)
			F2(i-n+14,j-3)=p(i,j)
		else
			E2(i-n+14,j-3)=v_thita*cos_thita*(p(i,j)+beita*k(i,j))
			F2(i-n+14,j-3)=v_thita*sin_thita*(p(i,j)+beita*k(i,j))
 		end if
 	end do
 end do
!upper and lower region
 do i=1,m
 	j=1
	do while(j<=7 .or. j>=n-6)
		cos_thita=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2)
        sin_thita=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2)
		v_thita=(Ma*cos_thita+dsqrt(1-Ma**2*sin_thita**2))
		if(j<7)then
			if(flag==0)then
				E3(i,j)=v_thita*cos_thita*(rou(i,j)+beita*k(i,j))
				F3(i,j)=v_thita*sin_thita*(rou(i,j)+beita*k(i,j))
			elseif(flag==1)then
				E3(i,j)=v_thita*cos_thita*(u(i,j)+beita*k(i,j))
				F3(i,j)=v_thita*sin_thita*(u(i,j)+beita*k(i,j))
			elseif(flag==2)then
				E3(i,j)=v_thita*cos_thita*(v(i,j)+beita*k(i,j))
				F3(i,j)=v_thita*sin_thita*(v(i,j)+beita*k(i,j))
			else
				E3(i,j)=v_thita*cos_thita*(p(i,j)+beita*k(i,j))
				F3(i,j)=v_thita*sin_thita*(p(i,j)+beita*k(i,j))
			end if
			j=j+1
		elseif(j==7)then
			if(flag==0)then
				E3(i,j)=v_thita*cos_thita*(rou(i,j)+beita*k(i,j))
				F3(i,j)=v_thita*sin_thita*(rou(i,j)+beita*k(i,j))
			elseif(flag==1)then
				E3(i,j)=v_thita*cos_thita*(u(i,j)+beita*k(i,j))
				F3(i,j)=v_thita*sin_thita*(u(i,j)+beita*k(i,j))
			elseif(flag==2)then
				E3(i,j)=v_thita*cos_thita*(v(i,j)+beita*k(i,j))
				F3(i,j)=v_thita*sin_thita*(v(i,j)+beita*k(i,j))
			else
				E3(i,j)=v_thita*cos_thita*(p(i,j)+beita*k(i,j))
				F3(i,j)=v_thita*sin_thita*(p(i,j)+beita*k(i,j))
			end if
			j=j+n-13
		else
			if(flag==0)then
				E3(i,j-n+14)=v_thita*cos_thita*(rou(i,j)+beita*k(i,j))
				F3(i,j-n+14)=v_thita*sin_thita*(rou(i,j)+beita*k(i,j))
			elseif(flag==1)then
				E3(i,j-n+14)=v_thita*cos_thita*(u(i,j)+beita*k(i,j))
				F3(i,j-n+14)=v_thita*sin_thita*(u(i,j)+beita*k(i,j))
			elseif(flag==2)then
				E3(i,j-n+14)=v_thita*cos_thita*(v(i,j)+beita*k(i,j))
				F3(i,j-n+14)=v_thita*sin_thita*(v(i,j)+beita*k(i,j))
			else
				E3(i,j-n+14)=v_thita*cos_thita*(p(i,j)+beita*k(i,j))
				F3(i,j-n+14)=v_thita*sin_thita*(p(i,j)+beita*k(i,j))
			end if
			j=j+1
		end if
	end do
 end do
!
!interior region: get k(m,n)
 do i=4,m-3
 	do j=4,n-3
		u0=E1((i-3):(i+3),j)
		v0=F1(i,(j-3):(j+3))
		stat_x=-3
		stat_y=-3
		k1(i,j)=-DRP7(u0,delta_x,stat_x)-DRP7(v0,delta_x,stat_y)
	end do
 end do
!
!inflow region
 do i=1,3
	do j=4,n-3
		u0=E2(1:7,j)
		stat_x=1-i
		v0=F2(i,(j-3):(j+3))
		stat_y=-3
		k1(i,j)=-DRP7(u0,delta_x,stat_x)-DRP7(v0,delta_x,stat_y)
	end do
 end do
!
!outflow region
 do i=1,3
	do j=4,n-3
		u0=E2(8:14,j)
		stat_x=-3-i
		v0=F2(i+11,(j-3):(j+3))
		stat_y=-3
		k1(i+m-3,j)=-DRP7(u0,delta_x,stat_x)-DRP7(v0,delta_x,stat_y)
		cos_thita=x_coor(i+11)/dsqrt(x_coor(i+11)**2+y_coor(j)**2)
        sin_thita=y_coor(j)/dsqrt(x_coor(i+11)**2+y_coor(j)**2)
		v_thita=(Ma*cos_thita+dsqrt(1-Ma**2*sin_thita**2))
		if(flag==3)then
			k1(i+m-3,j)=k1(i+m-3,j)-p(i+11,j)*v_thita/(2.0d0*sqrt(x_coor(i+11)**2+y_coor(j)**2))
		end if
	end do
 end do
!
!upper and lower region
 do i=1,m
	j=1
	do while(j<=3 .or. ((j>=n-2) .and. (j<=n)))
		if(j<=3)then
			v0=F3(i,7:1:-1)
			stat_y=j-7
			j_temp=j
		else
			v0=F3(i,8:14)
			stat_y=n-6-j
			j_temp=j-n+14
		end if
		if(i<=3)then
			u0=E3(1:7,j_temp)
			stat_x=1-i
		elseif(i>=m-2)then
			u0=E3((m-6):m,j_temp)
			stat_x=m-6-i
		else
			u0=E3((i-3):(i+3),j_temp)
			stat_x=-3
		end if
		k1(i,j)=-DRP7(u0,delta_x,stat_x)-DRP7(v0,delta_x,stat_y)
		if(j<3)then
			j=j+1
		elseif(j=3)then
			j=j+n-5
		else
			j=j+1
		end if
	end do
 end do
 k=k1
END SUBROUTINE F


!F(rou)
!interior region
	E=Ma*(rou+beita*k)+u
	F=v
	do i=4,m-3
 		do j=4,n-3
			u0=E((i-3):(i+3),j)
			v0=F(i,(j-3):(j+3))
			stat_x=-3
			stat_y=-3
			k(i,j)=-DRP7(u0,delta_x,stat_x)-DRP7(v0,delta_x,stat_y)
		end do
	end do
!inflow region: radiation boundary condition
 	do i=1,7
		do j=4,n-3
			cos_thita=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2)
        	sin_thita=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2)
			v_thita=(Ma*cos_thita+dsqrt(1-Ma**2*sin_thita**2))
			E(i,j)=v_thita*cos_thita*(rou(i,j)+beita*k(i,j))
			F(i,j)=v_thita*sin_thita*(rou(i,j)+beita*k(i,j))
		end do
 	end do
!
	do i=1,3
		do j=4,n-3
			u0=E(1:7,j)
			stat_x=1-i
			v0=F(i,(j-3):(j+3))
			stat_y=-3
			k(i,j)=-DRP7(u0,delta_x,stat_x)-DRP7(v0,delta_x,stat_y)
		end do
	end do
!upper and lower region: radiation boundary condition
 do i=1,m
	j=1
	do while(j<=7 .or. j>=n-6)
		cos_thita=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2)
        sin_thita=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2)
		v_thita=(Ma*cos_thita+dsqrt(1-Ma**2*sin_thita**2))
		E(i,j)=v_thita*cos_thita*(rou(i,j)+beita*k(i,j))
		F(i,j)=v_thita*sin_thita*(rou(i,j)+beita*k(i,j))
		j=j+1
	end do
 end do
!
 do i=1,m
	j=1
	do while(j<=3 .or. j>=n-2)
		if(j<=3)then
			v0=F(i,7:1:-1)
			stat_y=j-7
		else
			v0=F(i,(n-6):n)
			stat_y=n-6-j
		end if
		if(i<=3)then
			u0=E(1:7,j)
			stat_x=1-i
		elseif(i>=m-2)then
			u0=E((m-6):m,j)
			stat_x=m-6-i
		else
			u0=E((i-3):(i+3),j)
			stat_x=-3
		end if
		k(i,j)=-DRP7(u0,delta_x,stat_x)-DRP7(v0,delta_x,stat_y)
		j=j+1
	end do
 end do
!outer region: outflow boundary condition
 do i=n-6,n
	do j=4,m-3
		cos_thita=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2)
        sin_thita=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2)
		v_thita=(Ma*cos_thita+dsqrt(1-Ma**2*sin_thita**2))
		E(i,j)=Ma*(rou(i,j)+beita*k(i,j))+v_thita*cos_thita*p(i,j)-Ma*p(i,j)
		F(i,j)=v_thita*sin_thita*p(i,j)
 	end do
 end do
!
 do i=1,3
	do j=4,n-3
		u0=E((m-6):m,j)
		stat_x=m-6-i
		v0=F(i,(j-3):(j+3))
		stat_y=-3
		k(i,j)=-DRP7(u0,delta_x,stat_x)-DRP7(v0,delta_x,stat_y)
		cos_thita=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2)
        sin_thita=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2)
		v_thita=(Ma*cos_thita+dsqrt(1-Ma**2*sin_thita**2))
		k(i,j)=-p(i,j)*v_thita/(2.0d0*sqrt(x_coor(i)**2+y_coor(j)**2))
	end do
 end do
