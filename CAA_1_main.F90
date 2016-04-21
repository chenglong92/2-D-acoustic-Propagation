PROGRAM CAA_1
IMPLICIT NONE
 integer::i,j,k,kk
 real(kind=8)::Ma,x_range,y_range,delta_x,delta_t
 parameter(Ma=0.5,x_range=200,y_range=200,delta_x=1d0,delta_t=0.0569d0)
 integer::m,n,t_step
 parameter(m=int(x_range/delta_x)+1,n=int(y_range/delta_x)+1,t_step=5000)
 real(kind=8)::rou(m,n),u(m,n),v(m,n),p(m,n)
 real(kind=8)::rou_next(m,n),u_next(m,n),v_next(m,n),p_next(m,n)
 character(len=80)::filename,str1,g,gg,ggg,gggg
 real(kind=8)::start,finish
!
 call cpu_time(start) 
 Call initialize(rou,u,v,p,m,n,x_range,y_range)
 g='.dat'
 gg='rou_'
 ggg='u_'
 gggg='v_'
 kk=0
 do i=1,t_step
    Call LDDRK(rou_next,u_next,v_next,p_next,rou,u,v,p,m,n,Ma,delta_x,delta_t,x_range,y_range)
    rou=rou_next
    u=u_next
    v=v_next
    p=p_next
 	if(i==1 .or. mod(i,100)==0)then
    write(*,*)'i=',i
		kk=kk+1
        if(kk/10>=1)then
     	 	write(str1,'(I2)')kk
        else
            write(str1,'(I1)')kk
        end if
 		write(filename,*)(trim(gg)//trim(str1)//trim(g))
		open(unit=10,file=trim(filename))
		write(filename,*)(trim(ggg)//trim(str1)//trim(g))
		open(unit=11,file=trim(filename))
		write(filename,*)(trim(gggg)//trim(str1)//trim(g))
		open(unit=12,file=trim(filename))
		do j=1,m
			do k=1,n
				if(k==n)then
					write(unit=10,fmt='(f18.8)',advance='yes')rou(j,k)
                    write(unit=11,fmt='(f18.8)',advance='yes')u(j,k)
                    write(unit=12,fmt='(f18.8)',advance='yes')v(j,k)
				else
 					write(unit=10,fmt='(f18.8)',advance='no')rou(j,k)
                    write(unit=11,fmt='(f18.8)',advance='no')u(j,k)
                    write(unit=12,fmt='(f18.8)',advance='no')v(j,k)
				end if
			end do
		end do
		close(unit=10)
        close(unit=11)
        close(unit=12)
	end if
end do
call cpu_time(finish)
write(*,*)'time:',finish-start
!%matlabpool close;
END PROGRAM CAA_1
