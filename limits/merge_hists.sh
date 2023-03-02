n_lines=$(wc -l < ../config/procs.txt)
#years=(2016 2017 2018)
years=(2016)
histPath="hists"
echo "number of lines" $n_lines
mkdir -p /vols/cms/hsfar/hists_merged_taus

for ((i=1;i<=n_lines;i++)); do
    proc=$(awk "NR == $i" ../config/procs.txt)
    echo $i, $proc
    for year in "${years[@]}"; do
        echo $year
        files=($histPath/${proc}*_${year}*.root)
        if [ -e "${files[0]}" ]; then
            #rm hists_merged/${proc}_${year}.root
            hadd -f hists_merged/${proc}_${year}.root $histPath/${proc}*_${year}*.root
        else
            echo "skip"
        fi
    done
done

for year in "${years[@]}"; do
    rm hists_merged_taus/${year}.root
    hadd -f /vols/cms/hsfar/hists_merged_taus/${year}.root /vols/cms/hsfar/hists_merged_taus/*_${year}.root

done

