SUBROUTINE LDDRK(rou_next,u_next,v_next,p_next,rou,u,v,p,m,n,Ma,delta_x,delta_t,x_range,y_range)
implicit none
!!%the LDDRK scheme to get the next value for variables rou,u,v,p
!%
 integer::i,j,m,n
 real(kind=8)::rou_next(m,n),u_next(m,n),v_next(m,n),p_next(m,n)
 real(kind=8)::k_rou1(m,n),k_u1(m,n),k_v1(m,n),k_p1(m,n)
 real(kind=8)::k_rou2(m,n),k_u2(m,n),k_v2(m,n),k_p2(m,n)
 real(kind=8)::k_rou3(m,n),k_u3(m,n),k_v3(m,n),k_p3(m,n)
 real(kind=8)::k_rou4(m,n),k_u4(m,n),k_v4(m,n),k_p4(m,n)
 real(kind=8)::rou(m,n),u(m,n),v(m,n),p(m,n)
 real(kind=8)::rou1(m,n),u1(m,n),v1(m,n),p1(m,n)
 real(kind=8)::rou2(m,n),u2(m,n),v2(m,n),p2(m,n)
 real(kind=8)::rou3(m,n),u3(m,n),v3(m,n),p3(m,n)
 real(kind=8)::Ma,delta_x,delta_t,x_range,y_range
 real(kind=8)::w(4)=(/0.1630296d0,0.348012d0,0.3259288d0,0.1630296d0/)
 real(kind=8)::beita(3)=(/0.5d0,0.5d0,1.0d0/)
!%get k1
 Call F(k_rou1,k_u1,k_v1,k_p1,rou,u,v,p,m,n,Ma,delta_x,x_range,y_range);
 k_rou1=k_rou1*delta_t;
 k_u1=k_u1*delta_t;
 k_v1=k_v1*delta_t;
 k_p1=k_p1*delta_t;
 rou1=rou+beita(1)*k_rou1;
 u1=u+beita(1)*k_u1;
 v1=v+beita(1)*k_v1;
 p1=p+beita(1)*k_p1;
!%get k2 
 Call F(k_rou2,k_u2,k_v2,k_p2,rou1,u1,v1,p1,m,n,Ma,delta_x,x_range,y_range);
 k_rou2=k_rou2*delta_t;
 k_u2=k_u2*delta_t;
 k_v2=k_v2*delta_t;
 k_p2=k_p2*delta_t;
 rou2=rou+beita(2)*k_rou2;
 u2=u+beita(2)*k_u2;
 v2=v+beita(2)*k_v2;
 p2=p+beita(2)*k_p2;
!%get k3
 Call F(k_rou3,k_u3,k_v3,k_p3,rou2,u2,v2,p2,m,n,Ma,delta_x,x_range,y_range);
 k_rou3=k_rou3*delta_t;
 k_u3=k_u3*delta_t;
 k_v3=k_v3*delta_t;
 k_p3=k_p3*delta_t;
 rou3=rou+beita(3)*k_rou3;
 u3=u+beita(3)*k_u3;
 v3=v+beita(3)*k_v3;
 p3=p+beita(3)*k_p3;
!%get k4
 Call F(k_rou4,k_u4,k_v4,k_p4,rou3,u3,v3,p3,m,n,Ma,delta_x,x_range,y_range);
 k_rou4=k_rou4*delta_t;
 k_u4=k_u4*delta_t;
 k_v4=k_v4*delta_t;
 k_p4=k_p4*delta_t;
!%get the next flow field
 rou_next=rou+w(1)*k_rou1+w(2)*k_rou2+w(3)*k_rou3+w(4)*k_rou4;
 u_next=u+w(1)*k_u1+w(2)*k_u2+w(3)*k_u3+w(4)*k_u4;
 v_next=v+w(1)*k_v1+w(2)*k_v2+w(3)*k_v3+w(4)*k_v4;
 p_next=p+w(1)*k_p1+w(2)*k_p2+w(3)*k_p3+w(4)*k_p4;
end SUBROUTINE LDDRK
