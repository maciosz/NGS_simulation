wzbogacenie=$1
zaszumienie=$2
# pozadane srednie pokrycie
# zakladajac ze mamy piki o ustalonej dlugosci na chr21

bam_z_pikami=$3
bam_z_szumami=../szum.bam
# dostosowac sciezki w razie potrzeby

width=$(basename -s .bam $bam_z_pikami | cut -f2 -d_)


ile_pikow=$(echo "$wzbogacenie*$width*5/2" | bc)
ile_szumow=$(echo "$zaszumienie*45000000 / 60" | bc)

echo $ile_pikow, $ile_szumow

name=$(basename -s .bam $bam_z_pikami)

nazwa=${name}_wzbogacenie_${wzbogacenie}_zaszumienie_${zaszumienie}

if [ -e ${nazwa}.bam ]
then
    echo Nothing to do, woo-hoo!
    exit 0
elif [ -e ../${name}_peaks_${ile_pikow}_szum_${ile_szumow}.bam ]
then
    cp ../${name}_peaks_${ile_pikow}_szum_${ile_szumow}.bam ${nazwa}.bam .
    echo Nothing to do, woo-hoo!
    exit 0
fi


if [ -e ../${name}_peaks_${ile_pikow}.bam ]
then
    echo peaks already sampled, cool
else
    echo Sampling peaks...
    bedtools sample -n $ile_pikow  -i $bam_z_pikami > tmp_${nazwa}.bam #2> err
    samtools sort tmp_${nazwa}.bam > ${name}_peaks_${ile_pikow}.bam
    rm tmp_${nazwa}.bam
    samtools index ${name}_peaks_${ile_pikow}.bam
    mv ${name}_peaks_${ile_pikow}.bam ../
fi


if [ -e ../szum_${ile_szumow}.bam ]
then
    echo background already sampled, cool
else
    echo Sampling background...
    bedtools sample -n $ile_szumow  -i $bam_z_szumami  > tmp_${nazwa}.bam #2> err
    samtools sort tmp_${nazwa}.bam > szum_${ile_szumow}.bam
    rm tmp_${nazwa}.bam
    samtools index szum_${ile_szumow}.bam
    mv szum_${ile_szumow}.bam* ../
#    mv szum_${ile_szumow}.bam* ../../mixed_bams/
fi

echo Merging
#samtools merge tmp_${nazwa}.bam ${name}_peaks_${ile_pikow}.bam ../../mixed_bams/szum_${ile_szumow}.bam 
samtools merge tmp_${nazwa}.bam ../${name}_peaks_${ile_pikow}.bam ../szum_${ile_szumow}.bam 
samtools sort tmp_${nazwa}.bam > ${nazwa}.bam
rm tmp_${nazwa}.bam
samtools index ${nazwa}.bam

echo Creating coverages
bedtools genomecov -bg -ibam ${nazwa}.bam > ${nazwa}.bedgraph
~/software/bedGraphToBigWig ${nazwa}.bedgraph /mnt/chr7/data/maciosz/ref/hg38_chr21_bowtie2_index/hg38_chr21.size ${nazwa}.bigwig


