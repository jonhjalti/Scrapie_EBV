!Program for simulating genotypes for a single locus based on provided pedigree
!The pedigree should be with running numbers, 1 to number of animals, ordered with oldest animals first
!Missing parents should be coded with integers <1
!Costom made version for making seperate training files
!Can have two extra colomns in the pedigree file, one for sex and on for year
!Animals with known genotype can be chosen based on sex and year. 

program simgt
implicit none
integer :: i,j,k,l,m,n,o
integer :: ierr
integer :: noped
integer :: ngr,noout
integer :: ar
integer,allocatable :: miss(:)
integer,allocatable :: ped(:,:),gt(:,:),fourth(:),year(:)

real :: proptest(10)
real,allocatable :: af(:,:)
real,allocatable :: rgt(:,:)

character*40 :: pedfile
character*40 :: outfile1,outfile2
character*30 :: form1, form2
character*1  :: rletter,yletter
character*10 :: scencode(10)

logical :: test,testcontrol(10),testyear(10)

!Control file
!Includes the name of the pedigree file
!and allele frequencies to use in the base generations
open(10,file='control.txt')

read(10,'(a)')pedfile             !Name of pedigree file
read(10,*)ngr                  !Number of different allele frequencies
read(10,*)l                   !Number of alleles
allocate(af(ngr,l),miss(ngr))
do i=1,ngr
  read(10,*)miss(i)
  do j=1,l-1
    read(10,*)af(i,j)
  enddo
enddo
read(10,*)noout
if(noout.gt.10)then
  noout=10
  write(*,*)'Max 10 outputs'
endif 
 
testcontrol=.FALSE.
testyear=.FALSE.

do i=1,noout
  read(10,*)scencode(i)     !Name of scenario to write to file name
  read(10,'(a1)')rletter
  read(10,'(a1)')yletter
  if(rletter.eq."M")testcontrol(i)=.TRUE.
  if(yletter.eq."Y")testyear(i)=.TRUE.
  read(10,*)proptest(i)
  if(testyear(i))read(10,*)ar
enddo
close(10)

do i=1,ngr
  do j=2,l-1
    af(i,j)=af(i,j)+af(i,j-1)
  enddo
enddo
!Reading pedigree file to count the lines
ierr=0
i=0
open(11,file=pedfile)
do while (ierr.eq.0)
  read(11,*,iostat=ierr)
  i=i+1
enddo
noped=i-1
write(*,*)'Finished reading control file'
write(*,*)ngr,'unknown groups:'
do i=1, ngr
  write(*,*)miss(i), (af(i,j),j=1,l-1)
enddo
write(*,*)noped,'read from pedigree'

!Now reading the pedigree file again, saving the information
rewind(11)
allocate(ped(3,noped),gt(2,noped),rgt(l,noped),fourth(noped),year(noped))
do i=1,noped
  if(any(testcontrol).and.any(testyear))then
    read(11,*)ped(:,i),year(i),fourth(i)
  elseif(any(testcontrol))then
    read(11,*)ped(:,i),fourth(i)
  elseif(any(testyear))then
    read(11,*)ped(:,i),year(i)
  else
    read(11,*)ped(:,i)
  endif
  if(ped(2,i).gt.i)then
    write(*,*)'Invalid pedigree, sire',ped(2,i),'is less than',i
	exit
  endif
  if(ped(3,i).gt.i)then
    write(*,*)'Invalid pedigree, sire',ped(3,i),'is less than',i
	exit
  endif
enddo
close(11)

!Simulating genotypes

do i=1,noped
  do j=1,2      !over paternal and maternal alleles
    if(ped(j+1,i).lt.1)then
	  do k=1,ngr
	    if(ped(j+1,i).eq.miss(k))then
		  call base(m,l-1,af(k,:))
		  gt(j,i)=m
		endif
	  enddo
	else
	  call mendel(m)
	  gt(j,i)=gt(m,ped(j+1,i))
	endif
  enddo
enddo
rgt=0.0

!Format the genotypes as allele count
write(form1,'(a4,i1,a10)')'(i8,',l,'(1x,f4.2))'
write(form2,'(a10,i1,a10)')'(i8,1x,i2,',l-1,'(1x,f4.2))'
do i=1, noped
  rgt(gt(1,i),i)=1.0
  rgt(gt(2,i),i)=rgt(gt(2,i),i)+1.0
enddo
!write out the genotypes

!open(21,file='genotypes.txt')

open(22,file='phenotypes_all.txt')
do i=1,noped
!  write(21,'(i1,1x,i1)')gt(:,i)
  write(22,form1)ped(1,i),rgt(:,i)
enddo
close(22)

!Finding which anmals should be in the training set for each scenario
do o=1, noout
  write(outfile2,'(a4,a)')scencode(o),'/phenotypes_train.txt'
  open(23,file=outfile2)
  do i=1,noped
    test = .TRUE.
    if(testcontrol(o).and.testyear(o))then
      if(fourth(i).eq.1.and.year(i).lt.ar)call randgroup(proptest(o),test)
    elseif(testcontrol(o))then
      if(fourth(i).eq.1)call randgroup(proptest(o),test)
    elseif(testyear(o))then
      if(year(i).lt.ar)call randgroup(proptest(o),test)
    else
      call randgroup(proptest(o),test)
    endif
    if(test)then
!      write(24,form1)ped(1,i),rgt(:,i)
      continue
    else
      write(23,form2)ped(1,i),1,(rgt(j,i),j=1,l-1)
    endif
  enddo
  !close(21)
!  close(22)
  close(23)
  !close(24)
enddo
stop
end program simgt

!Subroutines

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

subroutine mendel(n)
!For selecting random allele from parents
  integer :: n
  real :: rn
  call random_number(rn)
  n=1
  if(rn.ge.0.50000)n=2

end subroutine mendel

subroutine randgroup(pr,n)
!For selecting random allele from parents
  logical :: n
  real :: rn, pr
  call random_number(rn)
  n=.TRUE.
  if(rn.ge.pr)n=.FALSE.

end subroutine randgroup