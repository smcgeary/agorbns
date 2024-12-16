# !/bin/sh
#
# Pipeline for analyzing RNA-Seq data from Bioo NEXTflex library prep
#

# check for the correct number of arguments
if [ "$#" -ne 5 ]; then
    echo -e "-h print this screen and exit"
    printf "Usage: `basename $0` <fastq file> <STAR index> <annotation file> <input directory> <output directory>";
    exit
fi


fastq_file="$1"
star_index="$2"
annot_file="$3"
temp_dir="$4"
out_dir="$5"

filename=$(basename $1 "$fullfile")
filename="${filename%.*}"

mkdir -p "$temp_dir" &&
mkdir -p "$out_dir" &&

# check existence of input files and directories
if [ ! -f "$fastq_file" ] ;
then
   echo -e "Cannot find fastq file $fastq_file";
   exit
fi
if [ ! -d "$star_index" ] ;
then
   echo -e "Cannot find directory $star_index";
   exit
fi
if [ ! -f "$annot_file" ] ;
then
   echo -e "Cannot find exon annotation file $annot_file";
   exit
fi
if [ ! -d "$temp_dir" ] ;
then
   echo -e "Cannot find directory $temp_dir";
   exit
fi
if [ ! -d "$out_dir" ] ;
then
   echo -e "Cannot find directory $out_dir";
   exit
fi
echo "$annot_file" &&
# reverse complement sequences - required for Bioo NEXTflex library prep
zcat "$fastq_file" | fastx_reverse_complement -z -o "$temp_dir/${filename}_revcomp.fastq.gz" &&

echo "Done with zcat/reverse complement." &&
# align to STAR
# flags: only take unique reads, only allow (readlength/0.04) mismatches, ...
STAR --genomeDir "$star_index" --runThreadN "10" --outFilterMultimapNmax "1" --outFilterMismatchNoverLmax "0.04" --outFilterIntronMotifs "RemoveNoncanonicalUnannotated" --outSJfilterReads "Unique" --outReadsUnmapped "Fastx" --readFilesCommand zcat --readFilesIn "$temp_dir/${filename}_revcomp.fastq.gz" --outFileNamePrefix "$temp_dir/${filename}_revcomp_star_" > "$temp_dir/${filename}_revcomp_star_logFile.txt" &&
echo "Done with STAR aglignment." &&


# convert SAM to BAM, then sort BAM file
samtools view -bS "$temp_dir/${filename}_revcomp_star_Aligned.out.sam" | samtools sort -T "$temp_dir/${filename}_TEMP2" -@ "24" -o "$temp_dir/${filename}_revcomp_star_sorted.bam" &&
echo "Done with samtools SAM to BAM, then sort." &&

## annotate
htseq-count "$temp_dir/${filename}_revcomp_star_sorted.bam" "$annot_file" -f "bam" -s "yes" -t "CDS" -i "transcript_id" -o "$temp_dir/${filename}_revcomp_star_sorted_counts.sam" > "$out_dir/${filename}_compiled.txt" &&
echo "Done with htseq-count." &&

# remove extraneous files
rm -rf "$temp_dir/${filename}_revcomp_star__STARtmp" &&
rm -f "$temp_dir/${filename}_revcomp_star_Aligned.out.sam" &&
rm "$temp_dir/${filename}_revcomp_star_sorted_counts.sam"

