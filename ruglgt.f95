!forrit til að rugla arfgerðargreiningar

program ruglgt
implicit none
integer :: i,j,k,l,m,n,o
integer :: nogt,noall, noall_1
integer :: norugl
integer :: ierr
integer :: einn=1

integer,allocatable :: animal(:)

real :: ruglpct
real :: af(10),wgt(10)
real,allocatable :: gt(:,:)

character*30 :: gtfile
character*30 :: form2

open(10,file='gtrugl_control.txt')

read(10,*)gtfile
read(10,*)ruglpct
read(10,*)noall
noall_1 = noall-1
do i=1, noall_1
  read(10,*)af(i)
enddo
close(10)

do j=2,noall_1
  af(j)=af(j)+af(j-1)
enddo

open(11,file=gtfile)
ierr=0
i=0
do while(ierr.eq.0)
  read(11,*,iostat=ierr)
  i=i+1
enddo
nogt=i-1
allocate(animal(nogt),gt(noall_1,nogt))
rewind(11)
gt=0.0
do i=1,nogt
  read(11,*)animal(i),einn,gt(1:noall_1,i)
enddo
close(11)

norugl=nint(ruglpct*nogt)
i=0
do while(i.le.norugl)
  wgt=0.0
  call sample(nogt,j)
  call base(k,noall,af(1:noall))
  wgt(k)=wgt(k)+1.0
  call base(k,noall,af(1:noall))
  wgt(k)=wgt(k)+1.0
  if(sum(abs(wgt-gt(:,j))).gt.0.01)then
    i=i+1
    gt(:,j)=wgt
  endif
enddo

write(form2,'(a10,i1,a10)')'(i8,1x,i2,',noall_1,'(1x,f4.2))'

open(21,file='phenotypes_rugl.txt')
do i=1,nogt
  write(21,form2)animal(i),einn,gt(1:noall_1,i)
enddo
close(21)
end program ruglgt

subroutine sample(n,m)
!finding animals at random 
  integer :: n
  real :: rn
  call random_number(rn)
  rn=rn*real(n)+0.5
  m=nint(rn)

end subroutine sample

subroutine base(k,l,af)             !out value, number of alleles, vector of allele frequencies. 
  ! Sampling alleles for base animals without known parents
  integer :: l,k
  real :: rn
  real :: af(l)
  call random_number(rn)
  do k=1,l
    if(rn.lt.af(k))exit
  enddo
end subroutine base