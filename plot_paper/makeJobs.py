outputPath = "/vols/cms/mkomm/HNL/histo/plot_paper/hists_met_jlFrac"

exportVars = {}
exportVars['OUTPATH'] = outputPath


jobArrayCfg = []
for year in [2016,2017,2018]:
    for plotType in [
        "tagger_CR_boosted","tagger_CR_resolved",
        "tagger_SR_boosted","tagger_SR_resolved",
        #"bdt_SR",
        "mllj_SR",
        "dR_SR",
        #"mllj_SRcuts"
        #"mll_central","mll_fwd"
    ]:
        opts = [
            "--plot "+plotType,
            "--year "+str(year),
            "--syst nominal",
            "-o $OUTPATH/hist"
        ]

        jobArrayCfg.append(opts)

        for syst in [
            'jesTotal','jer','unclEn','pu','muEff','track'
        ]:
            for systVar in ['Up','Down']:
                opts = [
                    "--plot "+plotType,
                    "--year "+str(year),
                    "--syst "+syst+systVar,
                    "-o $OUTPATH/hist"
                ]

                jobArrayCfg.append(opts)
        
                    



logPath = outputPath+"/log"
submitFile = open("submit_hists.sh", 'w')

submitFile.write('''#!/bin/bash
#$ -cwd
#$ -q hep.q
#$ -l h_rt=02:55:00
#$ -t 1-'''+str(len(jobArrayCfg))+'''
#$ -e '''+logPath+'''/log.$TASK_ID.err
#$ -o '''+logPath+'''/log.$TASK_ID.out

hostname
date
source ~/.bashrc
cd /vols/cms/mkomm/HNL/histo
source env.sh
cd /vols/cms/mkomm/HNL/histo/plot_paper
''')

for k,v in exportVars.items():
    submitFile.write("export "+k+"="+v+"\n")

submitFile.write("JOBS=(\n")
submitFile.write(" \"pseudo job\"\n")
for jobCfg in jobArrayCfg:
    submitFile.write(" \"")
    for opt in jobCfg:
        submitFile.write(" "+opt)
    submitFile.write("\"\n")
submitFile.write(")\n")

submitFile.write('''
echo ${JOBS[$SGE_TASK_ID]}
python makeHists.py ${JOBS[$SGE_TASK_ID]}
date
''')

submitFile.close()

print ("Created ",len(jobArrayCfg)," tasks created")
