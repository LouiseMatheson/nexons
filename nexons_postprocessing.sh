#!/bin/bash

processing_script=$1
shift 1
nexdir=${@: -1}/
gtf=${@: -2:1}

nexdir=${nexdir/"//"/"/"}

outfile="nexons_output_refined.txt"
f="5"
m="2"
i="50"
x="30"
g="500"
e="10"
d="30"
p="0.1"
r="5"
c="1"
s="100"
mem="50G"

while getopts o:f:m:i:x:g:e:d:p:r:c:s:b: flag
do
    case "${flag}" in
        o) outfile=${OPTARG};;
        f) f=${OPTARG};;
        m) m=${OPTARG};;
        i) i=${OPTARG};;
        x) x=${OPTARG};;
        g) g=${OPTARG};;
        e) e=${OPTARG};;
        d) d=${OPTARG};;
        p) p=${OPTARG};;
        r) r=${OPTARG};;
        c) c=${OPTARG};;
        s) s=${OPTARG};;
        b) mem=${OPTARG};;
    esac
done


file_prefix=`echo $outfile | sed -e "s/....$//"`
file_suffix=`echo $outfile | sed -e "s/.*\(....\)$/\1/"`

outfile_ends="${file_prefix}_completeEnds$file_suffix"

echo "nexons postprocessing script: $processing_script"
echo "gtf file: $gtf"
echo "directory containing nexons files to process: $nexdir"
echo "output file: $outfile"
echo "ends output: $outfile_ends"
echo "flexibility: $f"
echo "min exons: $m"
echo "min intron length: $i"
echo "min exon length: $x"
echo "overlap gap: $g"
echo "ends flexibility: $e"
echo "ends flexibility2: $d"
echo "end proportion reported: $p"
echo "end supporting reads: $r"
echo "ncores per job: $c"
echo "ngenes per job: $s"
echo "memory per job: $mem"


if compgen -G "FINISHED_subset*log" > /dev/null; then
    echo "Files of format FINISHED_subset*.log already exist in current directory; please (re)move before running"
    exit 1
fi


if [ -f $outfile ]; then 
    echo "Output filename already exists; exiting"
    exit 1
else 
    echo "Creating output file"
    printf "File\tGeneID\tGeneName\tChr\tStrand\tSplicePattern\tTranscriptID\tCount\tStarts\tEnds\tMergedIsoforms\tshortIntron_nIsoform\tshortIntron_count\tshortIntron_nUnique\tshortIntron_minWidth\tshortIntron_maxWidth\tMissingExonIntron_TranscriptID\tMissingExons\tMissingExonWidths\tMissingFirstLastExon\tMissingIntrons\tMissingIntronWidths\tOverlapping\tOverlappingGene\tSpliceMatch\tReadLocation\tn_junctions\tTruncationOrigin\tTruncationStartEndScore\tTruncationScoreMax\tTruncationScoreMin\tTruncationAnnoScore\tTruncationAnnoScoreMax\tTruncationAnnoScoreMin\tTxEndsWithinTrunc\tTxEndsWithinAnnoTrunc\tTruncationType\tTruncationTypeAnno\tStartEndProportions\tDefinedStart\tDefinedEnd\tTranscriptProportions\n" > ${file_prefix}_subset0$file_suffix
fi


if [ -f $outfile_ends ]; then
    echo "Ends output filename already exists; exiting"
    exit 1
else
    echo "Creating ends output file"
    printf "GeneName\tTranscriptID\tStart\tEnd\tFile\tSplicePattern\tCount\tTranscriptEndMatch\tProportion\n" > ${file_prefix}_completeEnds_subset0$file_suffix
fi


njobs=$(((`find $nexdir -maxdepth 1 -iregex ".*txt" | wc -l`-1)/$s+1))

echo "Files split between ${njobs} jobs"

if [ -d "logfiles_nexons_postprocess" ]; then
    echo "Warning: logfiles_nexons_postprocess directory already exists; existing files may be overwritten"
else
    mkdir logfiles_nexons_postprocess
fi


for job in `seq 1 $njobs`;
do
    mkdir -p "${nexdir}subset$job";
    find $nexdir -maxdepth 1 -iregex ".*txt" | head -n $s | xargs -i mv "{}" "${nexdir}subset$job"
    ssub -c $c --mem $mem -o logfiles_nexons_postprocess/${file_prefix}_subset${job}.log -j subset$job --email Rscript --vanilla $processing_script -o "${file_prefix}_subset${job}$file_suffix" -q "${file_prefix}_completeEnds_subset${job}$file_suffix" -f $f -m $m -i $i -x $x -g $g -e $e -d $d -p $p -r $r -c $c -z "FINISHED_subset${job}.log" $gtf ${nexdir}subset$job
done

for job in `seq 1 $njobs`;
do
    while ! [ -f "FINISHED_subset${job}.log" ]; 
    do 
        echo "waiting for subset ${job}..."
        sleep 1m
    done
done

rm FINISHED_subset*log

cat ${file_prefix}_subset*$file_suffix > $outfile
cat ${file_prefix}_completeEnds_subset*$file_suffix > $outfile_ends


rm ${file_prefix}_subset*$file_suffix
rm ${file_prefix}_completeEnds_subset*$file_suffix

mv ${nexdir}subset*/* $nexdir
rm -r ${nexdir}subset*

echo "Finished!"
