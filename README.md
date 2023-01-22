g++ -O2 -o finito.exe finitoMain.cpp getSingleMatrix.cpp finito.c d_thMul.c d_ppNormVector.c ~/HiC/straw_may_2022/C++/straw.cpp -I ~/HiC/straw_may_2022/C++ -lz -lcurl -lpthread

g++ -O2 -o finitoGW.exe finitoGW.cpp getGWMatrix.cpp finito.c d_thMul.c d_ppNormVector.c ~/HiC/straw_may_2022/C++/straw.cpp -I ~/HiC/straw_may_2022/C++ -lz -lcurl -lpthread
