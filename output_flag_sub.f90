subroutine output_flag(i,delta_t,flag)
implicit none
integer::i,j,m,flag
parameter(m=8)
integer::outs(m)=(/1,3,5,10,20,40,60,80/)
real(kind=8)::delta_t
!
flag=0
do j=1,m
	if(i==int(outs(j)/delta_t))then
		flag=1
	end if
end do
end subroutine output_flag
