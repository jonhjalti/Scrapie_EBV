! Forrit til að rugla ættskrá gáfulega
! Á að lesa upplýsingar um bu og ar 
! Setur foreldra af öðrum grip frá sama ári.
! Annað hvort bara vitlaus faðir eða báðir foreldrar

program wrong_pedigree
implicit none
integer :: i,j,k,l,m,n,o
integer :: noped
integer :: norugl
integer :: buar(500)
integer :: arvis(3,100)
integer :: tmppar(2)
integer,allocatable :: ped(:,:),ar(:),bu(:)

real :: ruglpct,vixlpct

character*30 :: pedfile
character*30 :: bufile

open(10,file='rugl_control.txt')

read(10,*)pedfile
read(10,*)bufile
read(10,*)noped
read(10,*)ruglpct
read(10,*)vixlpct

close(10)

allocate(ped(3,noped),ar(noped),bu(noped))

open(11,file=pedfile)
open(12,file=bufile)
do i=1,noped
  read(11,*)ped(:,i),ar(i)
  read(12,*)bu(i)
enddo
close(11)
close(12)

!Find first and last animal from each year
j=1
l=0
do i=4831, noped
  if(ar(i).ne.l)then
    arvis(1,j)=ar(i)
    arvis(2,j)=i
    if(j.gt.1)arvis(3,j-1)=i-1
    j=j+1
    l=ar(i)
  endif
  if(j.eq.100)then
    write(*,*)'ERROR'
	exit
  endif
enddo
write(*,*)j, ' ar'

open(22,file='wrongsire.txt')

norugl=nint(ruglpct*noped)
write(*,*)'A maximum of ',norugl,' sires to change'
o=0
do i=1,norugl
  call sample(noped,k)
  if(ped(2,k).eq.0)cycle
  l=1
  do while(arvis(1,l).lt.ar(k))
    l=l+1
    if(l.eq.100)then
      write(*,*)'ERROR'
	  exit
    endif
  enddo
  m=0
  do j=arvis(2,l), arvis(3,l)
    if(bu(k).eq.bu(j).and.ped(2,j).gt.0)then
	  m=m+1
	  buar(m)=j
	endif
	if(m.eq.499)then
      write(*,*)'stort bu-ar',ar(k),bu(k),k
	  exit
    endif
  enddo
  if(m.eq.0)cycle
  do j=1,m 
    call sample(m,n)
	if(ped(2,k).ne.ped(2,buar(n)))then
	  write(22,*)k,buar(n),ped(2,k),ped(2,buar(n))
	  ped(2,k)=ped(2,buar(n))
	  o=o+1
	  exit
	endif
  enddo
enddo
write(*,*)o,' sires changed'
close(22)

open(23,file='parswap.txt')
norugl=nint(vixlpct*noped)
write(*,*)'A maximum of ',norugl,' parents to swap'
o=0
do i=1,norugl
  call sample(noped,k)
  if(ped(2,k).eq.0.and.ped(3,k).eq.0)cycle
  l=1
  do while(arvis(1,l).lt.ar(k))
    l=l+1
    if(l.eq.100)then
      write(*,*)'ERROR'
	  exit
    endif
  enddo
  m=0
  do j=arvis(2,l), arvis(3,l)
    if(bu(k).eq.bu(j))then
	  m=m+1
	  buar(m)=j
	endif
  enddo
  if(m.eq.0)cycle
  do j=1,m 
    call sample(m,n)
	if(k.ne.buar(n))then
	  write(23,*)k,n,buar(n),ped(2:3,k),ped(2:3,buar(n))
	  tmppar=ped(2:3,k)
	  ped(2:3,k)=ped(2:3,buar(n))
	  ped(2:3,buar(n))=tmppar
	  o=o+1
	  exit
	endif
  enddo
enddo
write(*,*)o,' parents swaps done'
close(23)

open(21,file='ped_rug.txt')
do i=1,noped
  write(21,*)ped(:,i),ar(i)
enddo
close(21)

end program wrong_pedigree

subroutine sample(n,m)
!finding animals at random 
  integer :: n
  real :: rn
  call random_number(rn)
  rn=rn*real(n)+0.5
  m=nint(rn)

end subroutine sample