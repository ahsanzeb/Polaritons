
cd src
zip -r ../x.zip *
cd ..


echo '#!/usr/local/bin/python3.5' | cat - x.zip > polaritons
chmod +x polaritons
rm x.zip



