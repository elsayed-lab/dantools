#!/usr/bin/bash

trap "exit" SIGINT

#This script will finish much of the post-processing I have to do
#after the iterations are technically done

while [ -n "$1" ]; do
    case "$1" in
        --ref)
            ref="$2"
            shift 1
            ;;
        --alt)
            alt="$2"
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
        --contigs)
            IFS=',' read -r -a contigs <<< "$2"
            shift 1
            ;;
        --ref_ends)
            IFS=',' read -r -a ref_ends <<< "$2"
            shift 1
            ;;
        --alt_ends)
            IFS=',' read -r -a alt_ends <<< "$2"
            shift 1
            ;;
        --breaks)
            IFS=',' read -r -a breaks <<< "$2"
            shift 1
            ;;
        --aln)
            aln="$2"
            shift 1
            ;;
        --threads)
            threads="$2"
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
        --input_type)
            input_type="$2"
            shift 1
            ;;
    esac
    shift 1;
done
#cd "$outdir"

hisat2-build -q -p "$threads" -f "$alt" output/indexes/"$output_name"

if [ "$input_type" == 'fasta' ]; then
    hisat2 -x output/indexes/"$output_name" \
           -f -p "$threads" --no-softclip --no-spliced-alignment \
           -k 1 --no-unal -U "$source" -S output/strict_aln.sam --norc \
           2>output/strict_aln.stderr
elif [ "$input_type" == 'fastq_u' ]; then
    hisat2 -x output/indexes/"$output_name" \
           -q -p "$threads" --no-softclip \
           -k 1 --no-unal -U "$readsu" -S output/strict_aln.sam \
           2>output/strict_aln.stderr
elif [ "$input_type" == 'fastq_p' ]; then
    hisat2 -x output/indexes/"$output_name" \
           -q -p "$threads" --no-softclip \
           -k 1 --no-unal -1 "$reads1" -2 "$reads2" -S output/strict_aln.sam \
           2>output/strict_aln.stderr
fi


error=$(samtools sort -l 9 -O BAM -@ "$threads" -o output/strict_aln.bam output/strict_aln.sam 2>&1)
error=$(echo "$error" | grep -v 'merging from')
if [ -n "$error" ]; then
    echo "$error"
fi

samtools index output/strict_aln.bam
rm output/strict_aln.sam

samtools depth -aa -@ "$threads" -o output/depth.tsv output/strict_aln.bam

#I am commenting the below because I want to try a full chromosome
#stretcher instead of my weird bin stuff which requires extra annoying
#logic. This is possible now because stretcher is cool

#I used to use the needle function, but I think the stretcher program
#is just better. It seems to take much less RAM, it can even align
#whole Leishmania chromosomes. I think I'll still run it in parallel
#using my bin technique, since this will be almost mandatory on larger
#genomes, but my bin size can be whatever I set it to
lim=${#ref_ends[@]}
lim=$((lim - 1)) #because the above extracts length, but I use 0 indexed
mkdir needle_files
>needle_files/errors.txt
seq 0 "$lim" | while read idx; do
    ((x=x%threads)); ((x++==0)) && wait #Makes it run in parallel
    (
        if echo "${breaks[@]}" | grep -qw "$idx"; then
            stretcher -asequence alt/"${contigs[$idx]}" \
                   -bsequence ref/"${contigs[$idx]}" \
                   -sbegin1 1 -send1 "${alt_ends[$idx]}" \
                   -sbegin2 1 -send2 "${ref_ends[$idx]}" \
                   -aformat3 fasta -sid1 alt -sid2 ref \
                   -outfile needle_files/"$idx" \
                   -gapopen 50 -gapextend 15 \
                   2>>needle_files/errors.txt
        else
            stretcher -asequence alt/"${contigs[$idx]}" \
                   -bsequence ref/"${contigs[$idx]}" \
                   -sbegin1 "$((${alt_ends[$idx - 1]} + 1))" \
                   -send1 "${alt_ends[$idx]}" \
                   -sbegin2 "$((${ref_ends[$idx - 1]} + 1))" \
                   -send2 "${ref_ends[$idx]}" \
                   -aformat3 fasta -sid1 alt -sid2 ref \
                   -outfile needle_files/"$idx" \
                   -gapopen 50 -gapextend 15 \
                   2>>needle_files/errors.txt
        fi
    ) &
done
