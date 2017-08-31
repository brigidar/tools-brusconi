
#!/bin/bash

find * -name "*.fasta" | while read file
do
sed -i.bak "s/^[>.* ]/>${file} /g" $file;
sed -i.bak "s/.fasta//g" $file;
done
rm *bak
