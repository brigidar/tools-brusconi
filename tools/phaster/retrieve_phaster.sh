#!/bin/bash
wget --post-file=$1 "http://phaster.ca/phaster_api" -O 'output_1.txt'
cut -f 1 -d "," output_1.txt > job_id.txt; sed -i.bak 's/{"job_id":"//g' job_id.txt; sed -i.bak 's/"//g' job_id.txt; job1="$(cat job_id.txt)";

#checks if there is already an update from PHASTER for the prediction.
minimumsize=100
file=output_1.txt
actualsize=$(wc -c <"$file")
while [ $minimumsize -ge $actualsize ]
do
sleep 300
wget "http://phaster.ca/phaster_api?acc=$job1" -O output_1.txt
file=output_1.txt
actualsize=$(wc -c <"$file")

done

# once it is updated it takes the zip file and puts into a file that will be recognized by Galaxy
if grep -lq '"status":"Complete"' output_1.txt ; then
    cut -f 4 -d "," output_1.txt > phaster.txt
fi

sed -i.bak 's/"zip":"//g' phaster.txt ; sed -i.bak 's/"//g' phaster.txt;
phaster_1="$(cat phaster.txt)";
mkdir test;
wget "$phaster_1" -P ./test/;
unzip ./test/*.zip;
mv phage_regions.fna output1.fasta;




