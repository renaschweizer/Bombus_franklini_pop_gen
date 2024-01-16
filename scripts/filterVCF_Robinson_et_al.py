#!usr/bin/env python

'''
Custom filtering of VCF files
based on jarobinson vacquita_filterVCF.py
modified to only filter based on depth of coverage

Input: raw VCF
Output: filtered VCF prints to screen
- Sites failing filters are marked as FAIL_? or WARN_? in the 7th column
- Sites where REF is in [A,C,G,T] and ALT is in [A,C,G,T,.] go on to genotype filtering if AD and DP present in FORMAT
- Filtered out genotypes are changed to './.', all others reported
- Sites with non-reference alleles remaining after genotype filtering also filtered based on values in INFO column

Possible usage:

SCRIPT=filterVCF.py
python ${SCRIPT} myfile.vcf.gz | bgzip > myfile_filtered.vcf.gz
tabix -p vcf myfile_filtered.vcf.gz

'''

import sys
import gzip
import re

vcf_file = sys.argv[1]
VCF = gzip.open(vcf_file, 'rt')


# Min depth (1/3x mean)
minD={
'Bombus_franklini_BLX2616':3,
'Bombus_franklini_BLX2735':9,
'Bombus_franklini_BLX2736':13,
'Bombus_franklini_BLX3212':3,
'Bombus_franklini_BLX3213':8,
'Bombus_franklini_BLX3214':4,
'Bombus_franklini_BLX3215':5,
'Bombus_franklini_BLX3216':8,
'Bombus_franklini_BLX3217':14,
'Bombus_franklini_BLX3218':5,
'Bombus_franklini_BLX3219':3,
'Bombus_franklini_BLX3220':13,
'Bombus_franklini_BLX3221':9,
'Bombus_franklini_BLX3222':6,
'Bombus_franklini_BLX3223':7,
'Bombus_franklini_BLX3224':9,
'Bombus_franklini_BLX3225':5,
'Bombus_franklini_BLX3226':3,
'Bombus_franklini_BLX3227':9,
'Bombus_franklini_BLX3228':3,
'Bombus_franklini_BLX3229':7,
'Bombus_franklini_BLX3230':18,
'Bombus_franklini_BLX3231':22,
'Bombus_franklini_BLX3232':17,
'Bombus_franklini_BLX3233':14,
'Bombus_franklini_BLX3234':13,
'Bombus_franklini_BLX3235':8,
'Bombus_franklini_BLX3236':8,
'Bombus_franklini_BLX3237':6,
'Bombus_franklini_BLX3238':8,
'Bombus_franklini_GNS104':20,
'Baff_256':4,
'Baff_257':6,
'Baff_258':6,
'Baff_ESI3Aug211':6,
'Baff_ESI3Aug212':6,
'Baff_ESI5Aug212':4,
'Baff_MHJ725B':5,
'Baff_MJH706C':2,
'Baff_MJH707A':5,
'Baff_MJH725E':4,
'Baff_MJH725G':2,
'Baff_MJH727A':3,
'Baff_MLBLON121':4,
'Baff_MLBnest1':5,
'Baff_MLBnest2':5,
'Baff_MLBnest3':5,
'Baff_MLBNOVA47':5,
'Baff_MLBNOVA48':5,
'Bombus_affinis_aff019':3,
'Bombus_affinis_aff020':3,
'Bombus_affinis_aff021':2,
'Bombus_affinis_aff022':3,
'Bombus_affinis_aff023':2,
'Bombus_affinis_aff024':1,
'Bombus_affinis_aff025':2,
'Bombus_affinis_aff026':2,
'Bombus_affinis_aff027':2,
'Bombus_affinis_aff028':2,
'Bombus_affinis_aff029':2,
'Bombus_affinis_aff030':3,
'Bombus_affinis_aff031':3,
'Bombus_affinis_aff032':3,
'Bombus_affinis_aff033':1,
'Bombus_affinis_aff034':3,
'Bombus_affinis_aff035':2,
'Bombus_affinis_aff036':3,
'Bombus_affinis_aff037':5,
'Bombus_affinis_aff038':4,
'Bombus_affinis_aff039':3,
'Bombus_affinis_aff040':4,
'Bombus_affinis_aff041':4,
'Bombus_affinis_aff042':4,
'Bombus_affinis_aff043':3,
'Bombus_affinis_aff044':2,
'Bombus_affinis_aff045':3,
'Bombus_affinis_aff046':3,
'Bombus_affinis_aff047':2,
'Bombus_affinis_aff048':5,
'Bombus_affinis_aff049':4,
'Bombus_affinis_aff050':4,
'Bombus_affinis_aff051':4,
'Bombus_affinis_aff052':3,
'Bombus_affinis_aff053':5,
'Bombus_affinis_aff054':5,
'Bombus_affinis_aff055':4,
'Bombus_affinis_aff056':3,
'Bombus_affinis_aff057':4,
'Bombus_affinis_BLX3405':2,
'Bombus_affinis_BLX3406':3,
'Bombus_affinis_BLX3407':3,
'Bombus_affinis_BLX3408':3,
'Bombus_affinis_BLX3409':3,
'Bombus_affinis_BLX3410':3,
'Bombus_affinis_BLX3411':3,
'Bombus_affinis_BLX3412':3,
'Bombus_affinis_BLX3413':3}

# Max depth (2x mean)
maxD={
'Bombus_franklini_BLX2616':17,
'Bombus_franklini_BLX2735':56,
'Bombus_franklini_BLX2736':76,
'Bombus_franklini_BLX3212':15,
'Bombus_franklini_BLX3213':51,
'Bombus_franklini_BLX3214':23,
'Bombus_franklini_BLX3215':32,
'Bombus_franklini_BLX3216':51,
'Bombus_franklini_BLX3217':87,
'Bombus_franklini_BLX3218':32,
'Bombus_franklini_BLX3219':20,
'Bombus_franklini_BLX3220':80,
'Bombus_franklini_BLX3221':55,
'Bombus_franklini_BLX3222':33,
'Bombus_franklini_BLX3223':42,
'Bombus_franklini_BLX3224':55,
'Bombus_franklini_BLX3225':30,
'Bombus_franklini_BLX3226':17,
'Bombus_franklini_BLX3227':51,
'Bombus_franklini_BLX3228':20,
'Bombus_franklini_BLX3229':43,
'Bombus_franklini_BLX3230':109,
'Bombus_franklini_BLX3231':129,
'Bombus_franklini_BLX3232':100,
'Bombus_franklini_BLX3233':85,
'Bombus_franklini_BLX3234':75,
'Bombus_franklini_BLX3235':47,
'Bombus_franklini_BLX3236':50,
'Bombus_franklini_BLX3237':35,
'Bombus_franklini_BLX3238':49,
'Bombus_franklini_GNS104':122,
'Baff_256':25,
'Baff_257':36,
'Baff_258':33,
'Baff_ESI3Aug211':38,
'Baff_ESI3Aug212':34,
'Baff_ESI5Aug212':24,
'Baff_MHJ725B':27,
'Baff_MJH706C':14,
'Baff_MJH707A':33,
'Baff_MJH725E':23,
'Baff_MJH725G':14,
'Baff_MJH727A':15,
'Baff_MLBLON121':22,
'Baff_MLBnest1':28,
'Baff_MLBnest2':27,
'Baff_MLBnest3':31,
'Baff_MLBNOVA47':32,
'Baff_MLBNOVA48':29,
'Bombus_affinis_aff019':18,
'Bombus_affinis_aff020':15,
'Bombus_affinis_aff021':13,
'Bombus_affinis_aff022':18,
'Bombus_affinis_aff023':13,
'Bombus_affinis_aff024':8,
'Bombus_affinis_aff025':15,
'Bombus_affinis_aff026':10,
'Bombus_affinis_aff027':14,
'Bombus_affinis_aff028':14,
'Bombus_affinis_aff029':11,
'Bombus_affinis_aff030':19,
'Bombus_affinis_aff031':17,
'Bombus_affinis_aff032':16,
'Bombus_affinis_aff033':8,
'Bombus_affinis_aff034':17,
'Bombus_affinis_aff035':14,
'Bombus_affinis_aff036':17,
'Bombus_affinis_aff037':28,
'Bombus_affinis_aff038':25,
'Bombus_affinis_aff039':17,
'Bombus_affinis_aff040':24,
'Bombus_affinis_aff041':24,
'Bombus_affinis_aff042':26,
'Bombus_affinis_aff043':17,
'Bombus_affinis_aff044':15,
'Bombus_affinis_aff045':17,
'Bombus_affinis_aff046':17,
'Bombus_affinis_aff047':11,
'Bombus_affinis_aff048':28,
'Bombus_affinis_aff049':22,
'Bombus_affinis_aff050':22,
'Bombus_affinis_aff051':27,
'Bombus_affinis_aff052':19,
'Bombus_affinis_aff053':30,
'Bombus_affinis_aff054':27,
'Bombus_affinis_aff055':23,
'Bombus_affinis_aff056':20,
'Bombus_affinis_aff057':23,
'Bombus_affinis_BLX3405':12,
'Bombus_affinis_BLX3406':17,
'Bombus_affinis_BLX3407':16,
'Bombus_affinis_BLX3408':15,
'Bombus_affinis_BLX3409':16,
'Bombus_affinis_BLX3410':20,
'Bombus_affinis_BLX3411':16,
'Bombus_affinis_BLX3412':17,
'Bombus_affinis_BLX3413':17}


# Individual genotype filtering function
#  - Genotypes failing filters are set to missing (./.)
#  - Applies individual min and max depth filters
#  - Filters heterozygotes if the allele balance (REF/DP) is <20% or >80%
#  - Filters homozygotes if more than 10% of the alleles are different type
#  - 'sample' is the sample name
#  - 'GT_entry' is the entire genotype entry for that individual (typically GT:AD:DP:GQ)
#  - 'ADpos' is the position of the AD field in FORMAT (determined below)
#  - 'DPpos' is the position of the DP field in FORMAT (determined below)

def GTfilter(sample, GT_entry, ADpos, DPpos):
    if GT_entry[:1]=='.' : return GT_entry
    else:
        gt=GT_entry.split(':')
        nocall='./.:' + ':'.join(gt[1:])
        #print(gt)
        #print(len(gt))
        #print(gt[0])
        #print(DPpos)
        #print(gt[DPpos])
        if len(gt)>1:
            if gt[0] in ('0/0','0/1','1/1') and gt[DPpos]!='.':
                DP=int(gt[DPpos])
                if minD[sample]<=DP<=maxD[sample]:
                    if gt[0]=='0/0':
                        return GT_entry
#		REF=float(gt[ADpos].split(',')[0])
#                AB=float(REF/DP)
                    elif gt[0]=='0/1':
                        REF=float(gt[ADpos].split(',')[0])
                        AB=float(REF/DP)
                        if 0.2<=AB<=0.8: return GT_entry
                        else: return nocall
                    elif gt[0]=='1/1':
                        REF=float(gt[ADpos].split(',')[0])
                        AB=float(REF/DP)
                        if AB<=0.1: return GT_entry
                        else: return nocall
                    else: return nocall
                else: return nocall
            else: return nocall
        else: return nocall

# Get list of samples in VCF file
samples=[]
for line in VCF:
    if line.startswith('##'):
        pass
    else:
        for i in line.split()[9:]: samples.append(i)
        break


# Go back to beginning of file
VCF.seek(0)


# Write pre-existing header lines & add new lines describing filters being applied
for line0 in VCF:
    if line0.startswith('#'):
        if line0.startswith('##FORMAT'):
            sys.stdout.write('##FILTER=<ID=FAIL_REF,Description="Reference allele not one of [A,C,G,T].">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_ALT,Description="Alternate allele not one of [A,C,G,T,.].">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noADi,Description="AD not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noDPi,Description="DP not present in FORMAT.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_noGT,Description="No called genotypes remain after filtering.">\n')
            sys.stdout.write('##FILTER=<ID=WARN_missing,Description="Excess missingness (>25% of samples uncalled or set to missing).">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_QD,Description="QD < 4.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_FS,Description="FS > 60.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_MQ,Description="MQ < 40.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_MQRankSum,Description="MQRankSum < -12.5.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_ReadPosRankSum,Description="ReadPosRankSum < -8.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_SOR,Description="SOR > 3.0.">\n')
            sys.stdout.write('##FILTER=<ID=FAIL_excessHet,Description="Excess heterozygosity (>75% of genotypes are 0/1).">\n')
            sys.stdout.write('##INFO=<ID=VariantType,Number=.,Type=String,Description="Whether the Variant is monomorphic or SNP).">\n')
            sys.stdout.write(line0)
            break
        #elif line0.startswith('##INFO'):
         #   sys.stdout.write('##INFO=<ID=VariantType,Number=.,Type=String,Description="Whether the Variant is monomorphic or SNP).">\n')
          #  sys.stdout.write(line0)
           # break
        else: sys.stdout.write(line0)


# Go through VCF file line by line to apply filters
for line0 in VCF:
    if line0.startswith('#'):
        sys.stdout.write(line0)
        continue
    line=line0.strip().split('\t')
    #print(line)

### Site filtering:
### Keep any filters that have already been applied
    filter=[]
    if line[6] not in ('.', 'PASS'):
        filter.append(line[6])

### Check REF allele
    if line[3] not in ['A','C','G','T']:       
        filter.append('FAIL_REF') 

### Check ALT allele
    if line[4] not in ['A','C','G','T','.']:
        filter.append('FAIL_ALT') 

### Access INFO field annotations
    if ';' in line[7]:
        INFO=line[7].split(';')
        d=dict(x.split('=') for x in INFO)
    else:
        INFO=line[7]
        if '=' in INFO:
            d={INFO.split('=')[0]:INFO.split('=')[1]}
        else: d={}

### Get the position of AD, DP in genotype fields
    if 'AD' in line[8]:
        ADpos=line[8].split(':').index('AD')
#    else: filter.append('FAIL_noADi')
    else: ADpos=999 

    #print(ADpos)
	
    if 'DP' in line[8]:
        DPpos=line[8].split(':').index('DP')
    else: filter.append('FAIL_noDPi')

### If any filters failed, write out line and continue
    if filter!=[]:
        sys.stdout.write('%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:])) )
        continue

### Genotype filtering:
    GT_list=[]
    for i in range(0,len(samples)):
        GT=GTfilter(samples[i],line[i+9],ADpos,DPpos)
        GT_list.append(GT)

### Recalculate AC, AN, AF for INFO (after this step, modified INFO values will be output)
    REF=2*[x[:3] for x in GT_list].count('0/0') + [x[:3] for x in GT_list].count('0/1')
    ALT=2*[x[:3] for x in GT_list].count('1/1') + [x[:3] for x in GT_list].count('0/1')
    if REF+ALT==0:
        filter.append('FAIL_noGT')
        sys.stdout.write('%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), '\t'.join(line[7:9]), '\t'.join(GT_list)) )
        continue    
    d['AC']=ALT
    d['AN']=REF+ALT
    d['AF']=round(float(ALT)/(float(REF)+float(ALT)), 4)

### Warn if >25% of genotypes missing
    n_missing=[x[:3] for x in GT_list].count('./.')
    if n_missing>0.25*len(samples):
        filter.append('WARN_missing')

### Fail sites with excess heterozygosity (>75% of genotypes are heterozygous)
    n_het=sum(x[:3]=='0/1' for x in GT_list)
    n_called=sum(x[:3]!='./.' for x in GT_list)
    if float(n_het/n_called)>0.75:
        filter.append('FAIL_excessHet')
        
### Set VariantType, outputting sites with just hom. REF genotypes without further filtering
    if ALT==0:
        d['VariantType']='NO_VARIATION'   
        if filter==[]:
            filter.append('PASS')
        sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), ';'.join('{0}={1}'.format(key, val) for key, val in sorted(d.items())), line[8], '\t'.join(GT_list)) )
        continue
    elif REF==0:
        d['VariantType']='NO_VARIATION'
    else:
        d['VariantType']='SNP'

### Fail sites with poor variant metrics
    if 'QD' in d and float(d['QD']) < 4.0:
        filter.append('FAIL_QD')
    if 'FS' in d and float(d['FS']) > 60.0:
        filter.append('FAIL_FS')
    if 'MQ' in d and float(d['MQ']) < 40.0:
        filter.append('FAIL_MQ')
    if 'MQRankSum' in d and float(d['MQRankSum']) < -12.5:
        filter.append('FAIL_MQRankSum')
    if 'ReadPosRankSum' in d and float(d['ReadPosRankSum']) < -8.0:
        filter.append('FAIL_ReadPosRankSum')
    if 'SOR' in d and float(d['SOR']) > 3.0:
        filter.append('FAIL_SOR')

### Write out new line
    if filter==[]:
        filter.append('PASS')
    sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(filter), ';'.join('{0}={1}'.format(key, val) for key, val in sorted(d.items())), line[8], '\t'.join(GT_list)) )


# Close files and exit
VCF.close()
exit()
