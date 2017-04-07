gcc bluestein_fft_test.c -DDEBUG=1 -g -o bluestein.test
gcc bluestein_gen_params.zz.c -DDEBUG=1 -g -o bluestein.gen
gcc bluestein_crt_conv.zz.c -DDEBUG=1 -g -o bluestein.conv
./bluestein.gen
./bluestein.conv
./bluestein.test
