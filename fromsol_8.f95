! Code for reading SOL file and pairing with the test phenotypes
! Version 6 prints out specific groups based on year read from 
! a seperate file

program fromsol
implicit none
integer :: i,j,k,l,m,n
integer :: noped,notr
integer :: ierr
integer :: readsol(7)
integer :: id
integer :: dummy(4)
integer, allocatable :: ar(:)

real :: readebv
real :: intercept(3,7)
real, allocatable :: tbv(:,:), ebv(:,:,:)

logical,allocatable :: gt(:)

character*30 :: form1

notr=4


!Read the true "phenotypes" from file
open(11,file='../phenotypes_all.txt')
ierr=0
i=0
do while (ierr.eq.0)
  read(11,*,iostat=ierr)
  i=i+1
enddo
noped = i-1
allocate(tbv(notr,noped),ebv(3,notr,noped),gt(noped),ar(noped))
tbv=0
ebv=0
rewind(11)
open(13,file='../pedout_kyn.txt')
do i=1,noped
  read(11,*)id,(tbv(j,i),j=1,notr)
  read(13,*)dummy(1:3),ar(i)
enddo

close(13)
close(11)
!Read the SOL file

open(12,file='rktest90.SOL')
ierr=0
gt=.FALSE.
do while(ierr.eq.0)
  read(12,*,iostat=ierr)readsol(:),readebv
  if(readsol(1).eq.2)intercept(1,readsol(2))=readebv
  if(readsol(1).eq.4)ebv(1,readsol(2),readsol(5))=readebv
  if(readsol(1).eq.4.and.readsol(6).gt.0)gt(readsol(5))=.TRUE.
enddo
close(12)

open(12,file='rktest95.SOL')
ierr=0

do while(ierr.eq.0)
  read(12,*,iostat=ierr)readsol(:),readebv
  if(readsol(1).eq.2)intercept(2,readsol(2))=readebv
  if(readsol(1).eq.4)ebv(2,readsol(2),readsol(5))=readebv
enddo
close(12)

open(12,file='rktest99.SOL')
ierr=0

do while(ierr.eq.0)
  read(12,*,iostat=ierr)readsol(:),readebv
  if(readsol(1).eq.2)intercept(3,readsol(2))=readebv
  if(readsol(1).eq.4)ebv(3,readsol(2),readsol(5))=readebv
enddo
close(12)

write(form1,'(a4,i1,a10)')'(i8,',4*notr,'(1x,f6.3))'
open(21,file='forcor_2021.txt')
open(22,file='forcor_2016.txt')
do i=1,noped
  do j=1,notr
    ebv(1,j,i)=ebv(1,j,i)+intercept(1,j)
	ebv(2,j,i)=ebv(2,j,i)+intercept(2,j)
	ebv(3,j,i)=ebv(3,j,i)+intercept(3,j)
  enddo
  if(.not.gt(i).and.ar(i).gt.2020)then
    write(21,form1)i,ebv(:,1:notr,i),tbv(1:notr,i)
  elseif(.not.gt(i).and.ar(i).gt.2015)then
    write(22,form1)i,ebv(:,1:notr,i),tbv(1:notr,i)
  endif
enddo
close(21)
close(22)

stop
end program fromsol