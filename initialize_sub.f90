SUBROUTINE initialize(rou,u,v,p,m,n,x_range,y_range)
IMPLICIT NONE
!%intialize the perturbation field
!%The initial perturbations include acoustic, entropy and vorticity wave
 integer::i,j,m,n
 real(kind=8)::x,y,x_range,y_range
 real(kind=8)::rou(m,n),u(m,n),v(m,n),p(m,n)
 do i=1,m
    do j=1,n
        x=-x_range/2.0+(i-1)*x_range/(1.0*(m-1))
        y=-y_range/2.0+(j-1)*y_range/(1.0*(n-1))
        rou(i,j)=exp(-log(2.0d0)*(x**2+y**2)/9.0d0)+0.1*exp(-log(2.0d0)*((x-67)**2+y**2)/25.0d0)
        u(i,j)=0.04*y*exp(-log(2.0d0)*((x-67)**2+y**2)/25.0d0)
        v(i,j)=-0.04*(x-67)*exp(-log(2.0d0)*((x-67)**2+y**2)/25.0d0)
        p(i,j)=exp(-log(2.0d0)*(x**2+y**2)/9.0d0)
    end do
end do
end SUBROUTINE initialize
