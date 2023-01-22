Balancing of either chromosome-specific or Genome Wide (eithe full or inter only) contact matrices from hic file.

You need gcc **version 7**. Please make sure that zlib, curl and pthread. You will also need straw (https://github.com/aidenlab/straw). In my case it is installed in ~/HiC/straw_may_2022. Replace it with your location.

To create an executable which computes norm vector for a particular chromosome do

**g++ -O2 -o finito.exe finitoMain.cpp getSingleMatrix.cpp finito.c d_thMul.c d_ppNormVector.c ~/HiC/straw_may_2022/C++/straw.cpp -I ~/HiC/straw_may_2022/C++ -lz -lcurl -lpthread**

Run **./finito.exe** to see usage.

To create an executable which computes Genome Wide (GW) norm vector do

**g++ -O2 -o finitoGW.exe finitoGW.cpp getGWMatrix.cpp finito.c d_thMul.c d_ppNormVector.c ~/HiC/straw_may_2022/C++/straw.cpp -I ~/HiC/straw_may_2022/C++ -lz -lcurl -lpthread**

Run **./finitoGW.exe** to see usage.
