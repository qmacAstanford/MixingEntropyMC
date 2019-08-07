#/bin/bash

rm -f out
gfortran mersenne_twister.f90 precision.f03 binning.f03 MonteCarlo.f90 -o out -O3
#gfortran precision.f03 binning.f03 MonteCarlo.f90 mersenne_twister.f90
nohup ./out 7172 data1 .False. 10.0 > message1 &
nohup ./out 6173 data2 .False. 10.0 > message2 &
nohup ./out 5174 data3 .False. 10.0 > message3 &
nohup ./out 4175 data4 .False. 10.0 > message4 &
nohup ./out 3176 data5 .False. 10.0 > message5 &
nohup ./out 2177 data6 .False. 10.0 > message6 &
nohup ./out 1178 data7 .False. 10.0 > message7 &
nohup ./out 9179 data8 .False. 10.0 > message8 &
nohup ./out 8110 data9 .False. 10.0 > message9 &
nohup ./out 7111 data10 .False. 10.0 > message10 &
nohup ./out 6112 data11 .False. 10.0 > message11 &
nohup ./out 5113 data12 .False. 10.0 > message12 &
nohup ./out 4114 data13 .False. 10.0 > message13 &
nohup ./out 3115 data14 .False. 10.0 > message14 &
nohup ./out 2116 data15 .False. 10.0 > message15 &
nohup ./out 1117 data16 .False. 10.0 > message16 &
