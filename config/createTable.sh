# Change the base_dir variable to the base path of your input files, if necessary also change the pattern 
# after -name bellow to match the ending of your files.
base_dir="/home/eli4001/lab_data/sharedData/PPCG"
find $base_dir -name "*R1_001.fastq.gz" | grep -v "Undetermined" |sort > F1.txt
find $base_dir -name "*R2_001.fastq.gz" | grep -v "Undetermined" |sort > F2.txt
sed 's/.*\///g' F1.txt | sed 's/_L00.*//g' > sampleNames.txt
sed 's/.*\///g' F1.txt | sed 's/.*\(L00.\).*/\1/g' > unit.txt
paste -d'\t' <(cut -f1 sampleNames.txt) <(cut -f1 unit.txt) <(cut -f1 F1.txt) <(cut -f1 F2.txt) > samples.tsv
sed -i '1s/^/sample_name\tunit_name\tfq1\tfq2\n/' samples.tsv
rm -f F1.txt F2.txt sampleNames.txt unit.txt
