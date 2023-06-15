import os

categories=[
    "mumu_OS", 
    "mumu_SS",
    "ee_OS",
    "ee_SS",
    "mue_OS",
    "mue_SS",
    "emu_OS",
    "emu_SS",
    ]

ncats = len(categories)

regions = ["A", "B", "C", "D"]
nregions = len(regions)

procsbkg = []
procsHNL = []
procsData = []
years = ["2016"]#, "2017", "2018"]

import yaml
with open("../config/samples.yml") as samples_file:
    samples_dict = yaml.load(samples_file, Loader=yaml.FullLoader)

for l in samples_dict:
    if ("HNL" in l) and ("ntau" in l):
        procsHNL.append(l)
    #elif l == "muon" or l == "electron": 
    #    procsData.append(l)
    else:
        continue
        procsbkg.append(l)

nprocsbkg = len(procsbkg)
nprocsHNL = len(procsHNL)
nprocsData = len(procsData)

header = """Executable   = exec_histo.sh
Output       = %s/log.$(Process).out
Error        = %s/log.$(Process).err
Log          = /dev/null

Request_CPUs   = 1
Request_Memory = 2.5GB
+RequestRuntime = 10440
Should_Transfer_Files = NO

Requirements = ( OpSysAndVer == "CentOS7" )

queue arguments from (
"""
dagJobs = []
for year in years:
    jobs = []

    '''
    for proc in procsbkg:
        for region in regions:
            for category in categories:
                jobs.append(f'--proc {proc} --category {category} --region {region} --year {year}')
    '''
    
    for couplingRange in [range(1,17),range(17,34),range(34,51),range(51,68)]:
        optCoupling = ""
        for coupling in couplingRange:
            optCoupling += "--coupling %i "%coupling
        for proc in procsHNL:
            for region in regions : 
                for category in categories:
                    jobs.append(f'--proc {proc} --category {category} --region {region} --year {year} {optCoupling}')
    
    for proc in procsData:
        for region in regions:
            for category in categories:
                if proc == "muon" and category.startswith("mu"):
                    jobs.append(f'--data --proc {proc} --category {category} --region {region} --year {year}')
                elif proc == "electron" and category.startswith("e"):
                    jobs.append(f'--data --proc {proc} --category {category} --region {region} --year {year}')
    
    njobs = len(jobs)

    print(f"The number of jobs is {njobs}.")

    chunk_size = 4500
    for i,ichunk in enumerate(range(0, len(jobs), chunk_size)): 
        logFolder = f"log/{year}/chunk{i}"
        try:
            os.makedirs(logFolder)
        except:
            pass
    
        jobsPerChunk = jobs[ichunk:ichunk+chunk_size]

        dagJobs.append(f"subHists_{year}_{i}.condor")
        with open(f"subHists_{year}_{i}.condor", "w") as f:
            print("   - make job with "+str(len(jobsPerChunk)))
            f.write(header%(logFolder,logFolder))
            for job in jobsPerChunk:
                f.write(f'"--suffix _chunk{i} '+job+'"\n')
            
            f.write(")\n")
            
with open("dag.condor","w") as f:
    for i,dagJob in enumerate(dagJobs):
        f.write(f"Job J{i+1} {dagJob}\n")
    for i in range(len(dagJobs)-1):
        f.write(f"PARENT J{i+1} CHILD J{i+2}\n")
    for i in range(len(dagJobs)):
        f.write(f"RETRY J{i+1} 10\n")    
            
