#!/bin/bash

trap "exit 2" 2

#Parse my arguments
base=''
fai=''
indexes=''
it=''
output_name=''
outdir=''
fragments=''

while [ -n "$1" ]; do
    case "$1" in
        -b)
            base="$2"
            shift 1
            ;;
        --fai)
            fai="$2"
            shift 1
            ;;
        --indexes)
            indexes="$2"
            shift 1
            ;;
        -i)
            it="$2"
            shift 1
            ;;
        --outdir)
            outdir="$2"
            shift 1
            ;;
        --output_name)
            output_name="$2"
            shift 1
            ;;
        --source)
            source="$2"
            shift 1
            ;;
        --readsu)
            readsu="$2"
            shift 1
            ;;
        --reads1)
            reads1="$2"
            shift 1
            ;;
        --reads2)
            reads2="$2"
            shift 1
            ;;
        --var_fraction)
            var_fraction="$2"
            shift 1
            ;;
        --var_depth)
            var_depth="$2"
            shift 1
            ;;
        --variant_caller)
            variant_caller="$2"
            shift 1
            ;;
        -t)
            threads="$2"
            shift 1
            ;;
        --input_type)
            input_type="$2"
            shift 1
            ;;
        --scoremin)
            scoremin="$2"
            shift 1
            ;;
    esac
    shift 1
done

#cd "$outdir"
itp=$((it - 1))

mkdir it"$it"

if [ ! "$it" -eq 0 ]; then
    error=$(bcftools consensus -f "$base" -o it"$it"/"$output_name"_it"$it".fasta it"$itp"/variants.bcf 2>&1)
    error=$(echo "$error" | grep -v 'Applied' | grep -v 'Note: the --sample' | grep -v 'overlaps')
    if [ -n "$error" ]; then
         "$error" 1>&2
    fi

    base=it"$it"/"$output_name"_it"$it".fasta
    mkdir it"$it"/indexes
    hisat2-build -q -p "$threads" -f it"$it"/"$output_name"_it"$it".fasta it"$it"/indexes/"$output_name"_it"$it"
    samtools faidx it"$it"/"$output_name"_it"$it".fasta
    indexes=it"$it"/indexes/"$output_name"_it"$it"
    fai=it"$it"/"$output_name"_it"$it".fasta.fai
else
    if [ "$fai" = 'NA' ]; then
        samtools faidx -o original_base.fasta.fai "$base"
        fai=original_base.fasta.fai
    fi
    if [ "$indexes" = 'NA' ]; then
        mkdir original_indexes
        hisat2-build -q -f "$base" original_indexes/original_indexes
        indexes=original_indexes/original_indexes
    fi
fi

if [ "$input_type" == 'fasta' ]; then
    hisat2 -x "$indexes" \
           -f -p "$threads" --no-softclip --no-spliced-alignment \
           --score-min "$scoremin" -k 1 --no-unal --norc \
           -U "$source" -S it"$it"/raw.sam \
           2>it"$it"/MAPPING.stderr
elif [ "$input_type" == 'fastq_u' ]; then
    hisat2 -x "$indexes" \
           -q -p "$threads" --no-softclip \
           --pen-canintronlen G,-8,7.5 --pen-noncanintronlen G,-8,7.5 \
           --score-min "$scoremin" -k 1 --no-unal \
           -U "$readsu" -S it"$it"/raw.sam \
           2>it"$it"/MAPPING.stderr
elif [ "$input_type" == 'fastq_p' ]; then
    hisat2 -x "$indexes" \
           -q -p "$threads" --no-softclip \
           --pen-canintronlen G,-8,7.5 --pen-noncanintronlen G,-8,7.5 \
           --score-min "$scoremin" -k 1 --no-unal \
           -1 "$reads1" -2 "$reads2" -S it"$it"/raw.sam \
           2>it"$it"/MAPPING.stderr
fi



error=$(samtools sort -l 9 -O BAM -@ "$threads" -o it"$it"/sorted.bam it"$it"/raw.sam 2>&1)
error=$(echo "$error" | grep -v 'merging from')
if [ -n "$error" ]; then
    echo "$error" 1>&2
fi

samtools index it"$it"/sorted.bam

rm it"$it"/raw.sam

#Process substitution doesn't work in shell, so I need to create my
#regions file separately
if [ "$variant_caller" -eq 'freebayes' ]; then
    fasta_generate_regions.py "$fai" 100000 > it"$it"/freebayes_regions.tmp
    freebayes-parallel it"$it"/freebayes_regions.tmp "$threads" -f "$base" --min-alternate-fraction "$var_fraction" --min-alternate-count "$var_depth" it"$it"/sorted.bam > it"$it"/variants.vcf
    rm it"$it"/freebayes_regions.tmp
] elif [ "$variant_caller" -eq 'bcftools' ]; then
    bcftools mpileup it"$it"/sorted.bam -a FORMAT/AD,FORMAT/DP --threads "$threads" -f "$base" it"$it"/sorted.bam |
        bcftools call -O v -m -v --threads 8 |
        bcftools filter --threads "$threads" -i "FORMAT/AD[0:1] >= ${var_depth} && (FORMAT/AD[0:1] / (FORMAT/AD[0:0] + FORMAT/AD[0:1])) >= ${var_fraction} > it"$it"/variants.vcf
fi

bcftools view -O bcf -o it"$it"/variants.bcf it"$it"/variants.vcf
bcftools index it"$it"/variants.bcf
