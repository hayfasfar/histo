n_lines=$(wc -l < ../config/procs.txt)
years=(2016)
#years=($1)
histPath="hists_09Mar23"
echo "number of lines" $n_lines
mkdir -p hists_merged

for ((i=1;i<=n_lines;i++)); do
    proc=$(awk "NR == $i" ../config/procs.txt)
    for year in "${years[@]}"; do
        echo $i/$n_lines ${proc} $year
        files=($histPath/${proc}*_${year}*.root)
        if [ -e "${files[0]}" ]; then
            #rm hists_merged/${proc}_${year}.root
            if [ ! -f hists_merged/${proc}_${year}.root ]; then
                hadd -f hists_merged/${proc}_${year}.root $histPath/${proc}*_${year}*.root
            else
                echo "skip - output exists"
            fi
        else    
            echo "skip - input missing"
        fi
    done
done

for year in "${years[@]}"; do
    #rm hists_merged2/${year}.root
    hadd -f hists_merged/${year}.root hists_merged/*_${year}.root
done

