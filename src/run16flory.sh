#/bin/bash

rm -f out
gfortran mersenne_twister.f90 precision.f03 binning.f03 MonteCarlo.f90 -o out -O3
#gfortran precision.f03 binning.f03 MonteCarlo.f90 mersenne_twister.f90 
nohup ./out 7172 florydata1 .True. 1.0 > florymessage1 & 
nohup ./out 6173 florydata2 .True. 1.0 > florymessage2 &
nohup ./out 5174 florydata3 .True. 1.0 > florymessage3 &
nohup ./out 4175 florydata4 .True. 1.0 > florymessage4 &
nohup ./out 3176 florydata5 .True. 1.0 > florymessage5 &
nohup ./out 2177 florydata6 .True. 1.0 > florymessage6 &
nohup ./out 1178 florydata7 .True. 1.0 > florymessage7 &
nohup ./out 9179 florydata8 .True. 1.0 > florymessage8 &
nohup ./out 8110 florydata9 .True. 1.0 > florymessage9 &
nohup ./out 7111 florydata10 .True. 1.0 > florymessage10 &
nohup ./out 6112 florydata11 .True. 1.0 > florymessage11 &
nohup ./out 5113 florydata12 .True. 1.0 > florymessage12 &
nohup ./out 4114 florydata13 .True. 1.0 > florymessage13 &
nohup ./out 3115 florydata14 .True. 1.0 > florymessage14 &
nohup ./out 2116 florydata15 .True. 1.0 > florymessage15 &
nohup ./out 1117 florydata16 .True. 1.0 > florymessage16 &
