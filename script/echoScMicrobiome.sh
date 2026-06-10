#!/bin/bash
# bash echoScMicrobiome_step1-4.sh sample.info project.id > jobs.sh

BASE=/path/to/your/folder
#BASE=/public/home/name/public#sample
Bacteria=$BASE/UHGG
software=$BASE/software
PPL=$BASE/01_Pipeline
python=$software/bin/python
htslib_path=$BASE/panhuaran/htslib
sendru_msr=$BASE/Sendru/prev6v2.msr
krakenDb=$BASE/kraken2_fungi

echo "#!/bin/bash"
echo "#SBATCH -J $2"
echo "#SBATCH -p cu"
echo "#SBATCH --nodes=1"
echo "#SBATCH --cpus-per-task=10"
echo "#SBATCH --mem=50G"
echo ""

p=`pwd`

echo "set -e"
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$htslib_path:$software/lib2"

# Step 1: Create directories and config.ini for each sample
cut -f 1-3 $1 | while read sample r1F1Path r2FqPath
do
    if [ ! -d "$sample" ]; then
        echo "mkdir $p/$sample"
    fi
    echo "sample=$sample" > $sample/config.ini
    echo "outdir=$p/$sample" >> $sample/config.ini
    echo "process=10" >> $sample/config.ini
    echo "nubeam_dedup=$software/nubeamdedup-master/Linux/nubeam-dedup" >> $sample/config.ini
    echo "kraken=$software/kraken2/kraken2" >> $sample/config.ini
    echo "braken=$PPL/microbiomePipe/est_abundance_240613.o" >> $sample/config.ini
    echo "filter_threshold=0" >> $sample/config.ini
    echo "krakenDb=$krakenDb" >> $sample/config.ini
done
echo ""

# Step 2: Remove adapters and run scMeta (taxonomy classification)
cut -f 1-3 $1 | while read sample r1FqPath r2FqPath
do
    echo "mkdir $p/$sample/clean_data"
    spaceR1FqPath=`echo $r1FqPath | tr ',' ' '`
    spaceR2FqPath=`echo $r2FqPath | tr ',' ' '`
    echo "$software/anchorbcm.o -1 <(zcat $spaceR1FqPath ) -2 <(zcat $spaceR2FqPath) -b $sendru_msr -o $p/${sample}/clean_data/${sample}_ "
    echo "cd $p/$sample"
    echo "echo '# Start scMeta for ${sample}: '\`date\`"
    echo "$python $PPL/microbiomePipe/scMeta.py --cfg config.ini"
    echo "echo '# Done:                  '\`date\`"
    echo "echo ''"
done
echo "cd $p"
echo ""

# Step 3: Select target genomes and build STAR index
cut -f 1-3 $1 | while read sample r1F1Path r2FqPath
do
    echo "tail -n +2 $sample/Result/${sample}_sc_taxonomy.report | awk -F '\t' '\$7 >= 0 && \$8 >= 0 {print \$3}' | sort -n | uniq > $sample/species.txt"
    if [ ! -d "$sample/ref" ]; then
        echo "mkdir $sample/ref"
    fi
    echo "echo '# Start extracting genome info for ${sample}: '\`date\`"
    echo "$python $PPL/microbiomePipe/selectGenomes.py -i $sample/species.txt -r $Bacteria/refseq_mapping.tsv -d $Bacteria/raw -o $sample/ref/genome"
    echo "echo '# Done:                                  '\`date\`"
    echo "$software/bin/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir $sample/ref/STAR_index --genomeFastaFiles $sample/ref/genome.fna --sjdbGTFfile $sample/ref/genome.gtf --sjdbOverhang 122 --sjdbGTFfeatureExon gene --sjdbGTFtagExonParentTranscript gene_id --sjdbGTFtagExonParentGene gene_id --genomeSAindexNbases 10 --limitGenomeGenerateRAM 120000000000"
done
echo -e "wait\n"

# Step 4: Align with STARsolo and generate gene expression matrix
cut -f 1-5 $1 | while read sample r1FqPath r2FqPath species expt
do
    # pre-index (when multiple fastq files are provided)
    echo "$software/bin/STAR --runThreadN 10 --soloType CB_UMI_Simple --soloCBwhitelist $PPL/starSoloPipe-preindex/20bp.whitelist.txt --genomeDir $sample/ref/STAR_index --soloCBstart 1 --soloCBlen 20 --soloUMIstart 21 --soloUMIlen 8 --readFilesIn $sample/clean_data/${sample}_2.fq $sample/clean_data/${sample}_1.fq --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS CB UB CR UR GX GN --alignSJoverhangMin 1000 --soloFeatures GeneFull --outFileNamePrefix $sample/${sample}_ --soloStrand Reverse --limitBAMsortRAM 120000000000 --soloCellFilter TopCells $expt --outFilterScoreMinOverLread 0.5"
    num=$(echo $r1FqPath | grep -o 'gz' | wc -l)
    if [ $num -eq 1 ]; then
        # ignore pre-index (single gzipped input)
        echo "$software/bin/STAR --runThreadN 10 --soloType CB_UMI_Simple --soloCBwhitelist $PPL/starSoloPipe-preindex/15bp.whitelist.txt --genomeDir $sample/ref/STAR_index --soloCBstart 1 --soloCBlen 15 --soloUMIstart 21 --soloUMIlen 8 --soloBarcodeReadLength 0 --readFilesIn $sample/clean_data/${sample}_2.fq $sample/clean_data/${sample}_1.fq --outSAMtype BAM SortedByCoordinate --outMultimapperOrder Random --runRNGseed 1 --outSAMattributes NH HI AS CB UB CR UR GX GN --alignSJoverhangMin 1000 --soloFeatures GeneFull --outFileNamePrefix $sample\"-nopi\"/${sample}-nopi_ --soloStrand Reverse --limitBAMsortRAM 120000000000 --soloCellFilter TopCells $expt --outFilterScoreMinOverLread 0.5"
    fi
done
echo ""
# Step 5: Functional gene matrix construction
echo "# Step 5: Functional gene matrix construction"
cut -f 1,2 $1 | while read sample r1FqPath; do
    echo "mkdir -p $p/\$sample/functional"
    # 获取正确的 solo 输出目录名
    echo "if [ \$(echo $r1FqPath | grep -o 'gz' | wc -l) -eq 1 ]; then"
    echo "  solo_dir=$sample-nopi/${sample}-nopi"
    echo "else"
    echo "  solo_dir=$sample/${sample}"
    echo "fi"
    echo "Rscript $PPL/functional_gene_matrix.R \\"
    echo "  --ref_genome $BASE/genomes-all_metadata_2.0.tsv \\"
    echo "  --gene_annot $BASE/df5_new_adPv.txt \\"
    echo "  --input_dir $p/\$solo_dir\_Solo.out/GeneFull/filtered \\"
    echo "  --output_base $p/\$sample/functional \\"
    echo "  --samples \$sample"
done