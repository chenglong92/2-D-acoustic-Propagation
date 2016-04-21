Module BASIC_VARIATION
IMPLICIT NONE
 real(kind=8)::Ma,x_range,y_range,delta_x,delta_t
 parameter(Ma=0.5,x_range=200,y_range=200,delta_x=1d0,delta_t=0.0569d0)
 integer::m,n,t_step
 parameter(m=int(x_range/delta_x)+1,n=int(y_range/delta_x)+1,t_step=5000)
END Module BASIC_VARIATION
!
MODULE GLOBAL_VARIATION
 USE BASIC_VARIATION
IMPLICIT NONE
 real(kind=8)::rou(m,n),u(m,n),v(m,n),p(m,n)
 COMMON /origin/ rou,u,v,p
END MODULE GLOBAL_VARIATION
!
MODULE GLOBAL_DERIV
 USE BASIC_VARIATION
IMPLICIT NONE
 real(kind=8)::k_rou(m,n),k_u(m,n),k_v(m,n),k_p(m,n)
 COMMON /origin/ k_rou,k_u,k_v,k_p
MODULE GLOBAL_DERIV
!
MODULE LDDRK_K1234
USE BASIC_VARIATION
IMPLICIT NONE
 real(kind=8)::k_rou1(m,n),k_u1(m,n),k_v1(m,n),k_p1(m,n)
 real(kind=8)::k_rou2(m,n),k_u2(m,n),k_v2(m,n),k_p2(m,n)
 real(kind=8)::k_rou3(m,n),k_u3(m,n),k_v3(m,n),k_p3(m,n)
 real(kind=8)::k_rou4(m,n),k_u4(m,n),k_v4(m,n),k_p4(m,n)
 COMMON /k_rou1234/ k_rou1,k_rou2,k_rou3,k_rou4
 COMMON /k_u1234/ k_u1,k_u2,k_u3,k_u4
 COMMON /k_v1234/ k_v1,k_v2,k_v3,k_v4
 COMMON /k_p1234/ k_p1,k_p2,k_p3,k_p4
END MODULE LDDRK_K1234
