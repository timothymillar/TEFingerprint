module load samtools/1.2
samtools view -F 0x800 -F 0x100 -F 4 -F 8 -b $1 > ${2}_mappedPairs.bam
samtools view -F 0x800 -F 0x100 -F 4 -f 8 $1 > ${2}_mappedRead.bam
samtools view -F 0x800 -F 0x100 -F 8 -f 4 $1 | awk '{print "@"$3"__"$1"\n"$10"\n+\n"$11}' > ${2}_danglers.fastq
samtools view -f 0x2 ${2}_mappedPairs.bam > ${2}_mappedPairsproper.bam
module unload samtools/1.2
