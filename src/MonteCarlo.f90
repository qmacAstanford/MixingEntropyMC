program main
use mersenne_twister
use precision, only: dp, eps
implicit none
type(random_stat) rand_stat ! for random number generator
integer NT
integer ii, replicate, seed
integer, parameter :: NT_values(8)= [1000, 10000, 50000, 100000, 200000, 300000, 400000, 500000]
character(LEN=30) :: buffer
character(LEN=30) :: fname
logical isfile, verbose
integer widom_attempts, widom_successes
real(dp) phi
logical integer_steps
real(dp) L_big_real

verbose = .True.
Call get_command_Argument(1,buffer)
read(buffer,*) seed
print*, "seed", seed

Call get_command_Argument(2,buffer)
read(buffer,*) fname
print*, "fname", fname

Call get_command_Argument(3,buffer)
read(buffer,*) integer_steps
print*, "integer_steps", integer_steps

Call get_command_Argument(4,buffer)
read(buffer,*) L_big_real
print*, "L_big_real", L_big_real

!call random_setseed(7343, rand_stat) ! random seed for head node
call random_setseed(seed, rand_stat) ! random seed for head node

do replicate=1,100
    do ii = 1,8
        NT = NT_values(ii)
        call calc_P(rand_stat, NT, verbose, widom_attempts, widom_successes, phi, integer_steps, L_big_real)

        inquire(file = fname, exist = isfile)
        if (isfile) then
            open (unit = 5, file = fname, status ='OLD', POSITION = "append")
        else
            open (unit = 5, file = fname, status = 'NEW')
        endif
        write(5,*) NT, widom_attempts, widom_successes, phi
        close(5)

    enddo
enddo

end program


subroutine calc_P(rand_stat, NT, verbose, widom_attempts, widom_successes, phi,&
                  integer_steps, L_big_real)
use mersenne_twister
use binning
use precision, only: dp, eps
implicit none

type(random_stat), intent(inout) :: rand_stat
integer, intent(in) :: NT
logical, intent(in) :: verbose
integer, intent(out) :: widom_attempts
integer, intent(out) :: widom_successes
real(dp), intent(out) :: phi
logical, intent(in) :: integer_steps  ! On lattice?
real(dp), intent(in) :: L_big_real  ! side length of cubes

! --------------
! Binning variables
! --------------
type(binType) bin
double precision setBinSize(3)
integer setBinShape(3)
double precision maxvals(3)
double precision minvals(3)

real(dp), allocatable, dimension(:,:) :: R_real
real(dp) RP(3)
double precision sideLength

integer ix, iy, iz ! integers for selecting position
integer ii,jj
logical success, isCollision  ! collison happend
integer ind_move ! move counter "time"
logical periodic


! -----------------
! Move variables
! -----------------
integer irnd(1)
integer irnd3(3)
integer step_size ! amplitude of steps
real step_vec(3)
real step_scalar(1)
integer N_moves
integer move_success, move_attempts
integer percent_slide

! ---------------
! Measureing Observables
! ---------------
integer t_start_taking_data
integer lowest_bead
integer bin_counts(10)
real(dp) bin_width
logical trackLowest

! ------------------
!  Variable parameters
! ------------------
!NT = 10000
sideLength = 100.0_dp*L_big_real
! if you changes L_big_real or side length, be sure to change int insert move below!
N_moves = MAX(10000000, 2000*NT)
!integer_steps = .True.
bin_width=1.0_dp
step_size=2
trackLowest=.False.
periodic = .True.


allocate( R_real(3,NT) )

t_start_taking_data = N_moves/2
maxvals=[sideLength, sideLength, sideLength]
minvals=[0.0_dp, 0.0_dp, 0.0_dp]
setBinSize = [sideLength, sideLength, sideLength]
setBinShape = [8,8,8]
widom_attempts = 0
widom_successes = 0
bin_counts=0

! ----------------
! initial condigion
! ----------------
ix=0
iy=0
iz=0
call constructBin(bin,setBinShape,minvals,setBinSize)
lowest_bead = 1
do ii=1,NT
    R_real(:,ii) = ([real(ix,dp), real(iy,dp), real(iz,dp)]+0.5_dp)*L_big_real
    call addBead(bin,R_real,NT,ii)

    if (trackLowest .and. R_real(2,ii)  < R_real(2,lowest_bead)) then
        lowest_bead = ii
    endif

    ix=ix+1
    if (ix >= sideLength/L_big_real - eps) then
        iy=iy+1
        ix=0
    endif
    if (iy >= sideLength/L_big_real - eps) then
        iz=iz+1
        iy=0
    endif
enddo

move_success=0
move_attempts=0
do ind_move=1,N_moves
    if (mod(ind_move,100) < 70) then
        ! -----------
        ! Slide move
        ! -----------
        call random_index(NT, irnd, rand_stat)
        ii=irnd(1)
        if (integer_steps) then
            call random_index(step_size*2+1, irnd3, rand_stat)
            irnd3 = irnd3 - step_size - 1
            RP = R_real(:,ii) + [real(irnd3(1),dp), real(irnd3(2),dp), real(irnd3(3),dp)]
        else
            call random_number(step_scalar, rand_stat)
            call random_index(3, irnd, rand_stat)
            RP = R_real(:,ii)
            RP(irnd(1)) = R_real(irnd(1),ii) + step_size*(step_scalar(1)-0.5_dp)
        endif
        if (periodic) then
            RP(1) = modulo(RP(1),sideLength)
            RP(2) = modulo(RP(2),sideLength)
            RP(3) = modulo(RP(3),sideLength)
        endif
    else
        ! -----------
        ! Random insertion move
        ! -----------
        call random_index(NT, irnd, rand_stat)
        ii=irnd(1)
        if (integer_steps) then
            call random_index(int(sideLength+eps), irnd3, rand_stat)
            if (periodic) then
                RP(1) = real(irnd3(1)-1,dp)
                RP(2) = real(irnd3(2)-1,dp)
                RP(3) = real(irnd3(3)-1,dp)
            else
                print*, "Not written yet ..."
                stop
            endif
        else
            call random_number(step_vec, rand_stat)
            if (periodic) then
                RP = step_vec*sideLength
            else
                RP = (step_vec)*(sideLength-L_big_real) + L_big_real/2.0_dp
            endif
        endif
    endif

    !success=.True.
    success = .not. isCollision(bin, NT, R_real, RP, L_big_real, ii, periodic, sideLength)
    if ((.not. periodic) .and.&
        (MINVAL(RP) < L_big_real/2.0_dp - eps .or.&
         MAXVAL(RP) > sideLength - L_big_real/2.0_dp + eps)) then
        success = .False.
    endif
    if (success) then

        if (trackLowest) then
            if (ii == lowest_bead .and. RP(2) > R_real(2,ii)) then
                lowest_bead = 0
            else if (RP(2) < R_real(2,lowest_bead)) then
                lowest_bead=ii
            endif
        endif

        call removeBead(bin,R_real(:,ii),ii)
        R_real(:,ii) = RP
        if (trackLowest .and. lowest_bead==0) then
            lowest_bead=1
            do jj=1,NT
                if (R_real(2,jj)<R_real(2,lowest_bead)) lowest_bead=jj
            enddo
        endif
        call addBead(bin,R_real,NT,ii)
        move_success = move_success + 1
    endif
    move_attempts = move_attempts + 1

    if (ind_move>t_start_taking_data) then
        !--------------------------
        ! widom insertion
        !--------------------------
        if (integer_steps) then
            call random_index(int(sideLength+eps), irnd3, rand_stat)
            if (periodic) then
                RP(1) = real(irnd3(1)-1,dp)
                RP(2) = real(irnd3(2)-1,dp)
                RP(3) = real(irnd3(3)-1,dp)
            else
                print*, "Not written yet ..."
                stop
            endif
        else
            call random_number(step_vec, rand_stat)
            if (periodic) then
                RP = step_vec*sideLength
            else
                RP = (step_vec)*(sideLength-L_big_real) + L_big_real/2.0_dp
            endif
        endif

        success = .not. isCollision(bin, NT, R_real, RP, L_big_real, -1, periodic, sideLength)
        if (success) widom_successes = widom_successes+1
        widom_attempts = widom_attempts + 1
        !-----------------------------
        ! pressure
        !----------------------------
        if (trackLowest .and. R_real(2,lowest_bead) < L_big_real*0.5_dp + 10*bin_width) then
            jj = ceiling( (R_real(2,lowest_bead) - L_big_real*0.5_dp)/bin_width )
            bin_counts(jj) = bin_counts(jj) + 1
            !print*, "R_real",R_real(2,lowest_bead)
            !print*, "jj", jj
        endif
    endif


enddo
if (verbose) then
    print*, "----------------------------------------------------"
    print*, "Done"
    !do ii=1, NT
    !    print*, R_real(:,ii)
    !enddo
    print*, "Average position", [sum(R_real(1,:))/NT, sum(R_real(2,:))/NT, sum(R_real(3,:))/NT]
    print*, "widom success rate", real(widom_successes,dp)/real(widom_attempts,dp)
    print*, "move success rate", real(move_success,dp)/real(move_attempts,dp)
    print*, "Phi ", real(NT,dp)*(L_big_real**3)/(sideLength**3)
    if (trackLowest) print*, "bin_counts", bin_counts
    print*, "No. of succeses", widom_successes
    print*, "of total", widom_attempts
endif
phi = real(NT,dp)*(L_big_real**3)/(sideLength**3)


deallocate( R_real )
end subroutine
function isCollision(bin, NT, R_real, RP, L_big_real, exclude, periodic, sideLength)
    use binning
    use precision, only: dp, eps
    implicit none
    logical isCollision
    type(binType), intent(in) :: bin
    integer, intent(in) :: NT
    real(dp), intent(in) :: R_real(3,NT)
    real(dp), intent(in) :: RP(3)
    real(dp), intent(in) :: L_big_real
    integer, intent(in) :: exclude  ! don't count collsions with this bead
    logical, intent(in) :: periodic
    real(dp), intent(in) :: sideLength
    logical collision
    integer ix, iy, iz, jj
    real(dp) R_test(3)

    ! --------------
    ! Binning variables
    ! --------------
    integer, parameter :: maxNeighbors=27
    integer nNeighbors !number of neighbors found (so far)
    integer neighbors(maxNeighbors) ! list of bead ID's
    double precision distances(maxNeighbors) ! list of |r-r| values

    nNeighbors = 0 ! Clear list of neighbors
    if (periodic) then
        ! Loop over current an nearest periods to check for neighbors
        do ix=-1,1
            if (ix==-1 .and. RP(1) > L_big_real) cycle
            if (ix==1 .and. RP(1) < sideLength-L_big_real) cycle
            do iy=-1,1
                if (iy==-1 .and. RP(2) > L_big_real) cycle
                if (iy==1 .and. RP(2) < sideLength-L_big_real) cycle
                do iz=-1,1
                    if (iz==-1 .and. RP(3) > L_big_real) cycle
                    if (iz==1 .and. RP(3) < sideLength-L_big_real) cycle
                    R_test(1)=RP(1) + (real(ix,dp)*sideLength)
                    R_test(2)=RP(2) + (real(iy,dp)*sideLength)
                    R_test(3)=RP(3) + (real(iz,dp)*sideLength)
                    call findNeighbors(bin,R_test,L_big_real + eps,R_real,&
                                       NT,maxNeighbors,neighbors,distances,nNeighbors,-1)

                enddo
            enddo
        enddo
    else
        call findNeighbors(bin,RP,L_big_real + eps,R_real,NT,&
                           maxNeighbors,neighbors,distances,nNeighbors,-1)
    endif
    !print*, "number of neighbors", nNeighbors
    do jj = 1,nNeighbors
        if (neighbors(jj) == exclude) cycle
        if (distances(jj) < L_big_real - eps) then ! allow to overlap by eps for rounding error
            isCollision = .True.
            return
        endif
    enddo
    isCollision = .False.
end function

