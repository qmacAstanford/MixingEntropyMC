program main
use mersenne_twister
use binning
use precision, only: dp, eps
implicit none

! --------------
! Binning variables
! --------------
type(binType) bin
double precision setBinSize(3)
integer setBinShape(3)
double precision maxvals(3)
double precision minvals(3)
integer, parameter :: maxNeighbors=27
integer nNeighbors !number of neighbors found (so far)
integer neighbors(maxNeighbors) ! list of bead ID's
double precision distances(maxNeighbors) ! list of |r-r| values

integer NT
real(dp), allocatable, dimension(:,:) :: R_real
real(dp) RP(3)
double precision sideLength
real(dp) L_big_real  ! side length of cubes
real(dp) radius

integer ix, iy, iz ! integers for selecting position
integer ii,jj
logical success, collision  ! collison happend
type(random_stat) rand_stat ! for random number generator
integer ind_move ! move counter "time"
logical integer_steps  ! On lattice?


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
integer widom_attempts, widom_successes
integer t_start_taking_data
integer lowest_bead
integer bin_counts(10)
real(dp) bin_width

! ------------------
!  Variable parameters
! ------------------
NT = 5000
L_big_real=10.0_dp
sideLength = 1000.0_dp
N_moves = 90000000
integer_steps = .False.
bin_width=1.0_dp
step_size=2


allocate( R_real(3,NT) )
call random_setseed(7343, rand_stat) ! random seed for head node

t_start_taking_data = N_moves/2
maxvals=[sideLength, sideLength, sideLength]
minvals=[0.0_dp, 0.0_dp, 0.0_dp]
setBinSize = [sideLength, sideLength, sideLength]
setBinShape = [8,8,8]
radius = 1.75_dp*L_big_real ! > sqrt(3)
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
    
    if (R_real(2,ii)  < R_real(2,lowest_bead)) then
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
    else
        ! -----------
        ! Random insertion move
        ! -----------
        call random_index(NT, irnd, rand_stat)
        ii=irnd(1)
        if (integer_steps) then
            print*, "Not written yet ..."
            stop
        else
            call random_number(step_vec, rand_stat)
            RP = (step_vec)*(sideLength-L_big_real) + L_big_real/2.0_dp
        endif
    endif

    nNeighbors=0
    call findNeighbors(bin,RP,radius,R_real,NT,&
                          maxNeighbors,neighbors,distances,nNeighbors)
    !print*, "number of neighbors", nNeighbors
    do jj = 1,nNeighbors
        if (neighbors(jj) == ii) cycle
        success = collision(R_real(:,neighbors(jj)), RP, L_big_real)
        if (.not. success) then
            !print*, "Rand into jj", jj, " ii", ii
            exit
        endif
    enddo
    if (MINVAL(RP) < L_big_real/2.0_dp - eps .or.&
        MAXVAL(RP) > sideLength - L_big_real/2.0_dp + eps) then
        success = .False.
    endif
    if (success) then

        if (ii == lowest_bead .and. RP(2) > R_real(2,ii)) then
            lowest_bead = 0
        else if (RP(2) < R_real(2,lowest_bead)) then
            lowest_bead=ii
        endif

        call removeBead(bin,R_real(:,ii),ii)
        R_real(:,ii) = RP
        if (lowest_bead==0) then
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
        call random_number(step_vec, rand_stat)
        RP = step_vec*(sideLength-L_big_real) + L_big_real/2.0_dp
        nNeighbors=0
        call findNeighbors(bin,RP,radius,R_real,NT,&
                          maxNeighbors,neighbors,distances,nNeighbors)
        do jj = 1,nNeighbors
            success = collision(R_real(:,neighbors(jj)), RP, L_big_real)
            if (.not. success) exit
        enddo
        if (success) widom_successes = widom_successes+1
        widom_attempts = widom_attempts + 1
        !-----------------------------
        ! pressure
        !----------------------------
        if (R_real(2,lowest_bead) < L_big_real*0.5_dp + 10*bin_width) then
            jj = ceiling( (R_real(2,lowest_bead) - L_big_real*0.5_dp)/bin_width )
            bin_counts(jj) = bin_counts(jj) + 1
            !print*, "R_real",R_real(2,lowest_bead)
            !print*, "jj", jj
        else
            continue
            print*, "not counted"
        endif
    endif


enddo
print*, "Done"
!do ii=1, NT
!    print*, R_real(:,ii)
!enddo
print*, "Average position", [sum(R_real(1,:))/NT, sum(R_real(2,:))/NT, sum(R_real(3,:))/NT]
print*, "widom success rate", real(widom_successes,dp)/real(widom_attempts,dp)
print*, "move success rate", real(move_success,dp)/real(move_attempts,dp)
print*, "Phi ", real(NT,dp)*(L_big_real**3)/(sideLength**3)
print*, "bin_counts", bin_counts
print*, "of total", widom_attempts

deallocate( R_real )
end program
function collision(R1, R2, width)
    use precision, only: dp
    logical collision
    real(dp) R1(3)
    real(dp) R2(3)
    real(dp) width
    integer dd
    collision = .True.
    do dd=1,3
        if (abs(R1(dd)-R2(dd)) < width) then
            !print*, "Collision dd",dd," R_real",R_real(dd,neighbors(jj))," RP", RP(dd)
            continue
        else
            collision = .False.
            !print*, "pass      dd",dd," R_real",R_real(dd,neighbors(jj))," RP", RP(dd)
            return
        endif
    enddo
end function collision

