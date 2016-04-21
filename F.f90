SUBROUTINE F(k_rou,k_u,k_v,k_p,rou,u,v,p,m,n,Ma,delta_x,x_range,y_range)
IMPLICIT NONE
 integer::i,j,m,n,flag
 real(kind=8)::k_rou(m,n),k_u(m,n),k_v(m,n),k_p(m,n)
 real(kind=8)::rou(m,n),u(m,n),v(m,n),p(m,n),x_coor(m),y_coor(n)
 real(kind=8)::Deriv_E1(m,n),Deriv_E2(m,n),Deriv_E3(m,n),Deriv_E4(m,n)
 real(kind=8)::Deriv_F1(m,n),Deriv_F2(m,n),Deriv_F3(m,n),Deriv_F4(m,n)  
 real(kind=8)::Radi_coeff_cos(m,n),Radi_coeff_sin(m,n)
 real(kind=8)::Ma,delta_x,x_range,y_range,cos_thita(m,n),sin_thita(m,n)
!
 do i=1,m
	 x_coor(i)=-x_range/2.0d0+(i-1)*x_range/(m-1)
 end do
!
 do j=1,n
	y_coor(j)=-y_range/2.0d0+(j-1)*y_range/(n-1)
 end do
!
 k_rou=0.0d0
 k_u=0.0d0
 k_v=0.0d0
 k_p=0.0d0
!
 flag=0
 Call Partial_Deriv(Deriv_E1,-rou,m,n,delta_x,flag)
 Call Partial_Deriv(Deriv_E2,-u,m,n,delta_x,flag)
 Call Partial_Deriv(Deriv_E3,-v,m,n,delta_x,flag)
 Call Partial_Deriv(Deriv_E4,-p,m,n,delta_x,flag)
!
 flag=1
 Call Partial_Deriv(Deriv_F1,-rou,m,n,delta_x,flag)
 Call Partial_Deriv(Deriv_F2,-u,m,n,delta_x,flag)
 Call Partial_Deriv(Deriv_F3,-v,m,n,delta_x,flag)
 Call Partial_Deriv(Deriv_F4,-p,m,n,delta_x,flag)
!
 do i=4,m-3
 	do j=4,n-3
		k_rou(i,j)=Ma*Deriv_E1(i,j)+Deriv_E2(i,j)+Deriv_F3(i,j);
    	k_u(i,j)=Ma*Deriv_E2(i,j)+Deriv_E4(i,j);
        k_v(i,j)=Ma*Deriv_E3(i,j)+Deriv_F4(i,j);
        k_p(i,j)=Ma*Deriv_E4(i,j)+Deriv_E2(i,j)+Deriv_F3(i,j);
	end do
 end do
!
!%inlet region(three rows)
do i=1,m
    do j=1,n
        cos_thita(i,j)=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2);
        sin_thita(i,j)=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2);
        Radi_coeff_cos(i,j)=(Ma*cos_thita(i,j)+dsqrt(1-Ma**2*sin_thita(i,j)**2))*cos_thita(i,j);
        Radi_coeff_sin(i,j)=(Ma*cos_thita(i,j)+dsqrt(1-Ma**2*sin_thita(i,j)**2))*sin_thita(i,j);
        !Radi_coeff_cos(i,j)=(Ma*cos_thita(i,j)+dsqrt(1-Ma**2*sin_thita(i,j)**2));
        !Radi_coeff_sin(i,j)=(Ma*cos_thita(i,j)+dsqrt(1-Ma**2*sin_thita(i,j)**2));
    end do
end do
 do i=1,3
    do j=1,n
        k_rou(i,j)=Radi_coeff_cos(i,j)*Deriv_E1(i,j)+Radi_coeff_sin(i,j)*Deriv_F1(i,j);
        k_u(i,j)=Radi_coeff_cos(i,j)*Deriv_E2(i,j)+Radi_coeff_sin(i,j)*Deriv_F2(i,j);
        k_v(i,j)=Radi_coeff_cos(i,j)*Deriv_E3(i,j)+Radi_coeff_sin(i,j)*Deriv_F3(i,j);
        k_p(i,j)=Radi_coeff_cos(i,j)*Deriv_E4(i,j)+Radi_coeff_sin(i,j)*Deriv_F4(i,j);
        !k_rou(i,j)=Radi_coeff_cos(i,j)*(Deriv_E1(i,j)/cos_thita(i,j)+rou(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_u(i,j)=Radi_coeff_cos(i,j)*(Deriv_E2(i,j)/cos_thita(i,j)+u(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_v(i,j)=Radi_coeff_cos(i,j)*(Deriv_E3(i,j)/cos_thita(i,j)+v(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_p(i,j)=Radi_coeff_cos(i,j)*(Deriv_E4(i,j)/cos_thita(i,j)+p(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
    end do
end do
!
!%upper and lower region
 do i=4,(m-3)
    do j=1,3
        k_rou(i,j)=Radi_coeff_cos(i,j)*Deriv_E1(i,j)+Radi_coeff_sin(i,j)*Deriv_F1(i,j);
        k_u(i,j)=Radi_coeff_cos(i,j)*Deriv_E2(i,j)+Radi_coeff_sin(i,j)*Deriv_F2(i,j);
        k_v(i,j)=Radi_coeff_cos(i,j)*Deriv_E3(i,j)+Radi_coeff_sin(i,j)*Deriv_F3(i,j);
        k_p(i,j)=Radi_coeff_cos(i,j)*Deriv_E4(i,j)+Radi_coeff_sin(i,j)*Deriv_F4(i,j);
        !k_rou(i,j)=Radi_coeff_cos(i,j)*(Deriv_E1(i,j)/cos_thita(i,j)+rou(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_u(i,j)=Radi_coeff_cos(i,j)*(Deriv_E2(i,j)/cos_thita(i,j)+u(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_v(i,j)=Radi_coeff_cos(i,j)*(Deriv_E3(i,j)/cos_thita(i,j)+v(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_p(i,j)=Radi_coeff_cos(i,j)*(Deriv_E4(i,j)/cos_thita(i,j)+p(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
    end do
    do j=(n-2),n
        k_rou(i,j)=Radi_coeff_cos(i,j)*Deriv_E1(i,j)+Radi_coeff_sin(i,j)*Deriv_F1(i,j);
        k_u(i,j)=Radi_coeff_cos(i,j)*Deriv_E2(i,j)+Radi_coeff_sin(i,j)*Deriv_F2(i,j);
        k_v(i,j)=Radi_coeff_cos(i,j)*Deriv_E3(i,j)+Radi_coeff_sin(i,j)*Deriv_F3(i,j);
        k_p(i,j)=Radi_coeff_cos(i,j)*Deriv_E4(i,j)+Radi_coeff_sin(i,j)*Deriv_F4(i,j);
!
        !k_rou(i,j)=Radi_coeff_cos(i,j)*(Deriv_E1(i,j)/cos_thita(i,j)+rou(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_u(i,j)=Radi_coeff_cos(i,j)*(Deriv_E2(i,j)/cos_thita(i,j)+u(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_v(i,j)=Radi_coeff_cos(i,j)*(Deriv_E3(i,j)/cos_thita(i,j)+v(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
        !k_p(i,j)=Radi_coeff_cos(i,j)*(Deriv_E4(i,j)/cos_thita(i,j)+p(i,j)/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2)));
    end do
end do
!%%
!%outlet region
do i=(m-2),m
    do j=1,n
        !cos_thita=x_coor(i)/dsqrt(x_coor(i)**2+y_coor(j)**2);
        !sin_thita=y_coor(j)/dsqrt(x_coor(i)**2+y_coor(j)**2);
        k_p(i,j)=Radi_coeff_cos(i,j)*cos_thita(i,j)*Deriv_E4(i,j)+Radi_coeff_sin(i,j)*sin_thita(i,j)*Deriv_F4(i,j)&
            &-p(i,j)*(Ma*cos_thita(i,j)+dsqrt(1-Ma**2*sin_thita(i,j)**2))/(2.0d0*dsqrt(x_coor(i)**2+y_coor(j)**2));
        k_rou(i,j)=k_p(i,j)+Ma*Deriv_E1(i,j)-Ma*Deriv_E4(i,j);
        k_u(i,j)=Ma*Deriv_E2(i,j)+Deriv_E4(i,j);
        k_v(i,j)=Ma*Deriv_E3(i,j)+Deriv_F4(i,j);
    end do
end do
END SUBROUTINE F
