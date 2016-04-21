SUBROUTINE Partial_Deriv(Deriv_u,u,m,n,delta_x,flag)
implicit none
!%this function can reture derivatives for all points
!%using 7-points DRP scheme and the optimized backward or forward
!%difference scheme
!%flag=0: stands for the x direction derivative for variable u
!%flag=1: stands for the y direction derivative for variable u
!%delta_x: stands for the distance between two grids
 integer::i,j,m,n,flag
 real(kind=8)::Deriv_u(m,n),u(m,n)
 real(kind=8)::delta_x
 real(kind=8),external::DRP7
 if(flag==0)then
!%get the x-direction derivatives for u
    do i=1,m
        do j=1,n
            if(i<=3)then
                 Deriv_u(i,j)=DRP7(u(1:7,j),delta_x,1-i);
            elseif(i>=m-2)then
                Deriv_u(i,j)=DRP7(u((m-6):m,j),delta_x,m-6-i);
            else
                Deriv_u(i,j)=DRP7(u((i-3):(i+3),j),delta_x,-3);
            end if
        end do
    end do
 elseif(flag==1)then
!%get the y-direction derivatives for u
    do i=1,m
        do j=1,n
            !if(j<=3)then
            !    Deriv_u(i,j)=DRP7(u(i,7:-1:1),delta_x,j-7);
            !elseif(j>=n-2)then
            !    Deriv_u(i,j)=DRP7(u(i,(n-6):n),delta_x,n-6-j);
            !else
            !    Deriv_u(i,j)=DRP7(u(i,(j-3):(j+3)),delta_x,-3);
            !end if
	    if(j<=3)then
                Deriv_u(i,j)=DRP7(u(i,1:7),delta_x,1-j);
            elseif(j>=n-2)then
                Deriv_u(i,j)=DRP7(u(i,(n-6):n),delta_x,n-6-j);
            else
                Deriv_u(i,j)=DRP7(u(i,(j-3):(j+3)),delta_x,-3);
            end if
        end do
    end do
end if
end SUBROUTINE Partial_Deriv
