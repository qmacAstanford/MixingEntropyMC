#/bin/bash

rm -f out
gfortran mersenne_twister.f90 precision.f03 binning.f03 MonteCarlo.f90 -o out -O3
#gfortran precision.f03 binning.f03 MonteCarlo.f90 mersenne_twister.f90 
nohup ./out 7172 data1 .False. > message1 & 
nohup ./out 6173 data2 .False. > message2 &
nohup ./out 5174 data3 .False. > message3 &
nohup ./out 4175 data4 .False. > message4 &
nohup ./out 3176 data5 .False. > message5 &
nohup ./out 2177 data6 .False. > message6 &
nohup ./out 1178 data7 .False. > message7 &
nohup ./out 9179 data8 .False. > message8 &
nohup ./out 8110 data9 .False. > message9 &
nohup ./out 7111 data10 .False. > message10 &
nohup ./out 6112 data11 .False. > message11 &
nohup ./out 5113 data12 .False. > message12 &
nohup ./out 4114 data13 .False. > message13 &
nohup ./out 3115 data14 .False. > message14 &
nohup ./out 2116 data15 .False. > message15 &
nohup ./out 1117 data16 .False. > message16 &
