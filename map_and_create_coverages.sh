files=$@
for file in $files
do
    name=$(basename -s .fastq $file)
    echo $name
    echo Map
    bowtie2 -x /mnt/chr7/data/maciosz/ref/hg38_chr21_bowtie2_index/hg38_chr21 -U $file > ${name}_unsorted.bam 2> ${name}_mapping_statistics.txt
    echo Sort and index
    samtools sort ${name}_unsorted.bam > ${name}.bam 2> ${name}_sort.err
    samtools index ${name}.bam
    rm ${name}_unsorted.bam
    echo Calculate coverage
    bedtools genomecov -bg -ibam ${name}.bam > ${name}.bedgraph 2> ${name}_bedgraphs.err
    sort -k1,1 -k2,2n ${name}.bedgraph > tmp
    mv tmp ${name}.bedgraph
    echo Bedgraph to bigwig
    ~/software/bedGraphToBigWig ${name}.bedgraph /mnt/chr7/data/maciosz/ref/hg38_chr21_bowtie2_index/hg38_chr21.size ${name}.bigwig > ${name}2bigwig.out 2> ${name}2bigwig.err
done


