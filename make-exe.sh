
# python version you want to use:
pythonpath='#!/usr/local/bin/python3.5';
# executable name:
exe='polaritons'

# make a zip of src, add the pythonpath on top, and make it executable 
cd src
zip -r ../x.zip *
cd ..
echo $pythonpath | cat - x.zip > $exe
chmod +x $exe
rm x.zip

### if want to keep the executable to bin:
bin='/Users/panda/bin'
mv $exe $bin
### if need to add the path of bin to PATH:
# echo "PATH=$bin:"'$PATH': >> ~/.bashrc



