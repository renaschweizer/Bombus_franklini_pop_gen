# Follow **runs of homozygosity** calculations from Robinson et al. 2022: https://github.com/jarobin/vaquitagenomics2022/blob/v1/runs_of_homozygosity/runs_of_homozygosity_vaquita.sh

cd runs_of_homozygosity
cp ../filt_all/bombus_franklini_combined_PASS_n25_maxMissing75_minDP15_noSing.recode.vcf . 

IN=bombus_franklini_combined_PASS_n25_maxMissing75_minDP15_noSing.recode.vcf

# Identify ROH
module load bcftools/1.19
bcftools roh -G 30 -Orz -o ${IN}_bcftoolsROH.txt.gz ${IN}

# Reformat output
zcat ${IN}_bcftoolsROH.txt.gz | tail -n+4 \
| sed 's/# //g' | sed -E 's/\[[0-9]\]//g' | sed 's/ (bp)//g' \
| sed 's/ (average fwd-bwd phred score)//g' \
| tr ' ' '_'> ${IN}_bcftoolsROH.txt


OUT=bombus_franklini_combined_PASS_n25_maxMissing75_minDP15_noSing_pseudohom.vcf.gz

cat ${IN} | head -n 1000 | grep "^#" > head.tmp
# Manually edit head.tmp to add a sample name to the last line ("pseudohom")
cat ${IN} | grep -v "^#" | sed -e 's/$/\t0\/0/g' | cat head.tmp - | gzip > ${OUT}
bcftools roh -G 30 -Orz -o ${OUT}_bcftoolsROH.txt.gz ${OUT}
zcat ${OUT}_bcftoolsROH.txt.gz \
| awk -v s=pseudohom 'BEGIN{sum=0}{if ($2==s){sum+=$6}}END{printf "%s\t%s\n", s, sum}'
# pseudohom	ROH length: 279879617

cat ${IN} | head -n 1000 | grep "^#" | tail -n 1 | cut -f10- | tr '\t' '\n' > samples.list

# Calculate Froh using max length calculated from pseudo-genome
DATA=bombus_franklini_combined_PASS_n25_maxMissing75_minDP15_noSing.recode.vcf_bcftoolsROH.txt.gz
while read -r SAMPLE ; do 
zcat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/279879617, num, sum}'
done < samples.list >> roh_greater1Mb.txt

# Calculate Froh using max length calculated from pseudo-genome, chromosomes only
DATA=bombus_franklini_combined_PASS_n25_maxMissing75_minDP15_noSing.recode.vcf_bcftoolsROH.txt.gz
while read -r SAMPLE ; do 
zcat ${DATA} \
| awk -v s=${SAMPLE} 'BEGIN{sum=0}{if ($2==s && $3~"NC_" && $6>=1e6){sum+=$6; num+=1}}END{printf "%s\t%s\t%s\t%s\n", s, sum/279879617, num, sum}'
done < samples.list >> roh_greater1Mb_chrOnly.txt



# subset to only 18 chromosomes for detectRUNS

zcat bombus_franklini_combined_PASS_n25_maxMissing75_minDP15_noSing.recode.vcf_bcftoolsROH.txt.gz \
| egrep "NC_066344.1|NC_066345.1|NC_066346.1|NC_066347.1|NC_066348.1|NC_066349.1|NC_066350.1|NC_066351.1|NC_066352.1|NC_066353.1|NC_066354.1|NC_066355.1|NC_066356.1|NC_066357.1|NC_066358.1|NC_066359.1|NC_066360.1|NC_066361.1" \
| cat head2.tmp - | gzip - >  bombus_franklini_combined_PASS_n25_maxMissing75_minDP15_noSing.recode.vcf_bcftoolsROH_chr.txt.gz
