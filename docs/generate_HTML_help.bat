cd source 
sphinx-apidoc -T -f -o . ../../chemdistiller
cd ..
call make html
