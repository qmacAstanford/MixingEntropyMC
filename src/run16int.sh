#/bin/bash

rm -f out
gfortran mersenne_twister.f90 precision.f03 binning.f03 MonteCarlo.f90 -o out -O3
#gfortran precision.f03 binning.f03 MonteCarlo.f90 mersenne_twister.f90 
nohup ./out 7172 intdata1 .True. 10.0 > intmessage1 & 
nohup ./out 6173 intdata2 .True. 10.0 > intmessage2 &
nohup ./out 5174 intdata3 .True. 10.0 > intmessage3 &
nohup ./out 4175 intdata4 .True. 10.0 > intmessage4 &
nohup ./out 3176 intdata5 .True. 10.0 > intmessage5 &
nohup ./out 2177 intdata6 .True. 10.0 > intmessage6 &
nohup ./out 1178 intdata7 .True. 10.0 > intmessage7 &
nohup ./out 9179 intdata8 .True. 10.0 > intmessage8 &
nohup ./out 8110 intdata9 .True. 10.0 > intmessage9 &
nohup ./out 7111 intdata10 .True. 10.0 > intmessage10 &
nohup ./out 6112 intdata11 .True. 10.0 > intmessage11 &
nohup ./out 5113 intdata12 .True. 10.0 > intmessage12 &
nohup ./out 4114 intdata13 .True. 10.0 > intmessage13 &
nohup ./out 3115 intdata14 .True. 10.0 > intmessage14 &
nohup ./out 2116 intdata15 .True. 10.0 > intmessage15 &
nohup ./out 1117 intdata16 .True. 10.0 > intmessage16 &
