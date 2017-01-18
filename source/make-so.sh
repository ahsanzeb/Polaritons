
# f2py
f2py='/opt/local/bin/f2py-3.4'

rm tcorr.so
$f2py -c -m tcorr fortran/tcorr.f


