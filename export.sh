rm -rf bluestein
mkdir bluestein
cp bluestein_* bluestein/
cp make_bluestein.sh bluestein/
cp -r data.bs bluestein/
rm bluestein/data.bs/param*
rm bluestein/data.bs/poly.crt*
rm bluestein/bluestein_z_test.c
rm bluestein/bluestein_norev.h
rm bluestein/bluestein_zz_test.c
rm bluestein/bluestein_*.test 
