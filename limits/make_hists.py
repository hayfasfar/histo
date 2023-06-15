import ROOT
import os
import sys
import yaml
import json
import argparse
from histo import Process, Sample

ROOT.gErrorIgnoreLevel = ROOT.kError


EMPTY_BIN_FILL = 1e-7
EMPTY_BIN_ERROR = 1e-8

def write_hist(hist_nano: ROOT.TH1D, category_dict: dict, name: str, isMC: bool=True):
    """
    Make a new histogram for categorised events passing the cuts and write it to file

    Args:
        hist_nano     : nanoAOD histogram
        category_dict : dictionary of category cuts
        name : category name
        isMC : flag for MC/data
    """

    index_new = 0
    hist_limits = ROOT.TH1D("hist", "hist", len(category_dict), 0, len(category_dict))
    hist_limits.SetDirectory(0)
    print(f"Writing hist {name}")
    for index, category in category_dict.items():
        index_new += 1
        hist_content = hist_nano.GetBinContent(hist_nano.FindBin(index))
        hist_error = hist_nano.GetBinError(hist_nano.FindBin(index))
        print(index, index_new, category, hist_content)
        hist_limits.SetBinContent(index_new, hist_content)
        hist_limits.SetBinError(index_new, hist_error)
        hist_limits.GetXaxis().SetBinLabel(index_new, category)

    if isMC:
        fill_zero_bins(hist_limits)
    hist_limits.SetTitle(name)
    hist_limits.SetName(name)
    hist_limits.SetDirectory(root_file)
    hist_limits.Write()


def remove_neg_entries(hist: ROOT.TH1D):
    """ Removes negative and/or small entries from a histogram """

    alpha = 1. - 0.6827
    upErr = ROOT.Math.gamma_quantile_c(alpha/2,1,1)
    avgWeight = hist.Integral()/hist.GetEntries() if hist.GetEntries()>0 else -1
    #print "weight",avgWeight
    for ibin in range(hist.GetNbinsX()):
        c = hist.GetBinContent(ibin+1)
        if c<10**-4:
            hist.SetBinContent(ibin+1,10**-3)
            #note: in case of 0 entries the uncertainty is also small
            #(this is not the case with negative events)
            if hist.GetBinError(ibin+1)<10**-4 and avgWeight>0:
                #set uncertainties for empy bins
                #https://twiki.cern.ch/twiki/bin/viewauth/CMS/PoissonErrorBars
                hist.SetBinError(ibin+1,upErr*avgWeight)
            else:
                hist.SetBinError(ibin+1,10**-4)
        #print "bin%2i, %.1f+-%.1f (+-%.1f%%)"%(ibin,c,hist.GetBinError(ibin+1),100.*hist.GetBinError(ibin+1)/c if c>0 else -1)

def fill_zero_bins(hist):
    for ibin in range(hist.GetNbinsX()):
        c = hist.GetBinContent(ibin+1)
        if c<EMPTY_BIN_FILL:
            hist.SetBinContent(ibin+1,EMPTY_BIN_FILL)
            hist.SetBinError(ibin+1,EMPTY_BIN_ERROR)

def mass_cut(delta_m:float=10., region:str="D", syst:str="nominal", single_lepton=False) -> str:
    """ 
    Returns mass window cut repending on region

    Args:
    delta_m : mW+-delta_m used to define signal region band
    region: A, B, C, D
    syst: systematic variation
    single_lepton: deprecated

    Returns:
    mass cut str 
    """
    
    mW=80.
    m_upper = mW + (- delta_m)
    m_lower = mW - (- delta_m)

    #m_upper = mW + ( delta_m)
    #m_lower = mW - ( delta_m)

    if region not in ["A", "B", "C", "D"]:
        raise ValueError(f"Invalid region {region} selected")

    if region == "A" or region == "C":
        #return  f'({syst}_m_llj<{m_lower}) '
        return  f'(({syst}_m_llj>{m_upper} and {syst}_m_llj < 110. ) or( {syst}_m_llj<{m_lower} and {syst}_m_llj > 50.) ) '
    elif region == "B" or region == "D":
        return  f'({syst}_m_llj<{m_upper} and {syst}_m_llj>{m_lower}) '
    else:
        return ""

def tagger_cut(tagger_threshold: float, lower_threshold: float=0.1, region:str="D", syst:str="nominal") -> str:

    """ Returns tagger score cut repending on region
    Args:
    tagger_threshold : define signal region
    lower_threshold: lower tagger score threshold
    region: A, B, C, D
    syst: systematic variation

    Returns:
    tagger cut str 
    """
    
    if region not in ["A", "B", "C", "D"]:
        raise ValueError(f"Invalid region {region} selected")

    if tagger_threshold < lower_threshold:
        raise ValueError("Inconsistent tagger thresholds")

    if region == "A" or region == "B":
        return f'(tagger_score_{syst} > {lower_threshold} and tagger_score_{syst} < {tagger_threshold})'
    elif region == "C" or region == "D":
        return f'(tagger_score_{syst} > {tagger_threshold})'
        
def tagger_compound_variable(syst:str="nominal", single_lepton=False) -> str:
    """ Compound tagger variable to facilitate categorisation """
    # if single_lepton:
    #     return f"hnlJet_{syst}_llpdnnx_ratio_LLP_Q"
    # else:
    
    return f"({syst}_dR_l2j>0.4 and {syst}_dR_l2j<1.3)*hnlJet_{syst}_llpdnnx_ratio_LLP_Q+\
        ({syst}_dR_l2j<0.4)*hnlJet_{syst}_llpdnnx_ratio_LLP_QMU*subleadingLeptons_isMuon[0]+\
        +({syst}_dR_l2j<0.4)*hnlJet_{syst}_llpdnnx_ratio_LLP_QE*subleadingLeptons_isElectron[0]"
        
def make_hists(process, systematics_shapes, systematics_rates, cut_nominal, category_variable_nominal, thresholds, region,  coupling=None):

    hists = {}

    def make_hist(process, category_variable, thresholds, weight, cut, region,  syst="nominal"):
        
        thresholds_merged_prompt = thresholds["prompt"]["merged"]
        thresholds_merged_medium = thresholds["medium"]["merged"]
        thresholds_merged_displaced = thresholds["displaced"]["merged"]
        thresholds_resolved_prompt = thresholds["prompt"]["resolved"]
        thresholds_resolved_medium = thresholds["medium"]["resolved"]
        thresholds_resolved_displaced = thresholds["displaced"]["resolved"]
    
        
        lower_threshold_merged = 0.2
        lower_threshold_resolved = 0.1
      
        mass_cut_merged_prompt = mass_cut(delta_m=thresholds_merged_prompt[1], region=region, syst=syst)
        mass_cut_merged_medium = mass_cut(delta_m=thresholds_merged_medium[1], region=region, syst=syst)
        mass_cut_merged_displaced = mass_cut(delta_m=thresholds_merged_displaced[1], region=region, syst=syst)

        mass_cut_resolved_prompt = mass_cut(delta_m=thresholds_resolved_prompt[1], region=region, syst=syst)
        mass_cut_resolved_medium = mass_cut(delta_m=thresholds_resolved_medium[1], region=region, syst=syst)
        mass_cut_resolved_displaced = mass_cut(delta_m=thresholds_resolved_displaced[1], region=region, syst=syst)

        #print(f"Tagger thresholds, merged: {threshold_merged}, resolved: {threshold_resolved}")
        #print(f"DeltaM, merged: {deltam_merged}, resolved: {deltam_resolved}")
        '''
        mass_cut_merged_prompt = mass_cut(10., region=region, syst=syst)
        mass_cut_merged_medium = mass_cut(10. , region=region, syst=syst)
        mass_cut_merged_displaced = mass_cut(10. , region=region, syst=syst)

        mass_cut_resolved_prompt = mass_cut(10. , region=region, syst=syst)
        mass_cut_resolved_medium = mass_cut(10. , region=region, syst=syst)
        mass_cut_resolved_displaced = mass_cut(10. , region=region, syst=syst)
        ''' 
         
        tagger_cut_merged_prompt = tagger_cut(thresholds_merged_prompt[0], region=region, syst=syst , lower_threshold=lower_threshold_merged)
        tagger_cut_merged_medium = tagger_cut(thresholds_merged_medium[0], region=region, syst=syst , lower_threshold=lower_threshold_merged)
        tagger_cut_merged_displaced = tagger_cut(thresholds_merged_displaced[0], region=region, syst=syst , lower_threshold=lower_threshold_merged)

        tagger_cut_resolved_prompt = tagger_cut(thresholds_resolved_prompt[0], region=region, syst=syst, lower_threshold=lower_threshold_resolved)
        tagger_cut_resolved_medium = tagger_cut(thresholds_resolved_medium[0], region=region, syst=syst , lower_threshold=lower_threshold_resolved)
        tagger_cut_resolved_displaced = tagger_cut(thresholds_resolved_displaced[0], region=region, syst=syst , lower_threshold=lower_threshold_resolved)
       
       
        cut += f" and ("
        cut += f"({category_variable}==1 and {mass_cut_merged_prompt} and {tagger_cut_merged_prompt} ) or "
        cut += f"({category_variable}==2 and {mass_cut_merged_medium} and {tagger_cut_merged_medium} ) or "
        cut += f"({category_variable}==3 and {mass_cut_merged_displaced} and {tagger_cut_merged_displaced}) or "
        cut += f"({category_variable}==4 and {mass_cut_resolved_prompt} and {tagger_cut_resolved_prompt} ) or "
        cut += f"({category_variable}==5 and {mass_cut_resolved_medium} and {tagger_cut_resolved_medium}) or "
        cut += f"({category_variable}==6 and {mass_cut_resolved_displaced} and {tagger_cut_resolved_displaced} )"
        cut += f")"

        #print(syst, category_variable, cut, weight)

        #print (category_variable,cut,weight)
        hist_nano = process.Histo1D((category_variable, category_variable, 6, 0.5, 6.5), category_variable, cut=cut, weight=weight)
        hist_nano = hist_nano.Clone()
        hist_nano.SetDirectory(0)
        return hist_nano
        
        
    if systematics_rates is not None:
    # variations with constant shape but changing weight
        for syst, abrv in systematics_rates.items():
            for variation in ["Up", "Down"]:
                if "HNL" in process.name:
                    name = f"{process.name}_coupling_{coupling}_{abrv}{year}{variation}"
                    weight = f"weightHNL_{coupling}_{abrv}{variation}"
                else:
                    name = f"{process.name}_{abrv}{year}{variation}"
                    weight = f"weight_{abrv}{variation}"

                hists[name] = make_hist(process, category_variable_nominal, thresholds, weight, cut_nominal, region, syst="nominal")

    if systematics_shapes is not None:
        for syst in systematics_shapes:
            # add name for variations
            variations = ["Up", "Down"] if "nominal" not in syst else [""]
            syst_hack = syst.replace("nominal", "")
            for variation in variations:
                if "HNL" in process.name:
                    if "nominal" in syst:
                        name = f"{process.name}_coupling_{coupling}"
                    else:
                        name = f"{process.name}_coupling_{coupling}_{syst_hack}{year}{variation}"

                    weight = f"weightNominalHNL_{coupling}"
                else:
                    if "nominal" in syst:
                        name = f"{process.name}"
                    else:
                        name = f"{process.name}_{syst_hack}{year}{variation}"
                    weight = "weightNominal"

                syst_var_name = f"{syst}{variation}"

                # Systematic variation -- replace nominal by systematic in all cuts
                cut = cut_nominal.replace("nominal", syst_var_name)
                cut = cut.replace("nselectedJets_unclEnUp", "nselectedJets_nominal") #Hack!
                cut = cut.replace("nselectedJets_unclEnDown", "nselectedJets_nominal") #Hack!
                category_variable = category_variable_nominal.replace("nominal", syst_var_name)
                # read in hist from nanoAOD friends
                hists[name] = make_hist(process, category_variable, thresholds, weight, cut, region, syst=syst_var_name)

    return hists

# make histograms per year, process
parser = argparse.ArgumentParser()
parser.add_argument("--year",default="2016")
parser.add_argument("--proc", default="wjets")
parser.add_argument("--category", default="mumu_OS")
parser.add_argument("--region", default="D")

#parser.add_argument("--ntuple_path", default="/nfs/dust/cms/user/mkomm/HNL/ntuples/19Jan23")
#parser.add_argument("--output_path", default="/nfs/dust/cms/user/mkomm/HNL/histo/limits/hists_19Jan23")
#parser.add_argument("--ntuple_path", default="/vols/cms/hsfar/nanoAOD_friends/21Sep20")
#parser.add_argument("--output_path", default="/vols/cms/hsfar/hists")
parser.add_argument("--ntuple_path", default="/nfs/dust/cms/user/mkomm/HNL/ntuples/09Mar23")
parser.add_argument("--output_path", default="/nfs/dust/cms/user/mkomm/HNL/histo/limits/hists_09Mar23")

#parser.add_argument("--ntuple_path", default="/vols/cms/hsfar/nanoAOD_friends/19Jan23")
#parser.add_argument("--ntuple_path", default="/vols/cms/hsfar/nanoAOD_friends/09Mar23")
#parser.add_argument("--output_path", default="/vols/cms/mkomm/HNL/histo/limits/hists_09Mar23")

parser.add_argument("--data", action="store_true", default=False)
parser.add_argument("--test", action="store_true", dest="oneFile", default=False)
parser.add_argument("--couplings", default=[], action='append', type=int, dest="couplings")
parser.add_argument("--suffix", dest="suffix", default='')

args = parser.parse_args()
print(vars(args))



year = args.year
proc = args.proc

category_name = args.category
ntuple_path = os.path.join(f"{args.ntuple_path}", year)
region = args.region
oneFile = args.oneFile
isData = args.data
isMC = not isData
output_path = args.output_path
suffix = args.suffix


outputFile = os.path.join(output_path, f"{proc}_{args.category}_{region}_{year}{suffix}.root")
if os.path.exists(outputFile):
    outputExists = True
    rootFile = ROOT.TFile(outputFile)
    if (rootFile and rootFile.IsZombie()):
        print ("ROOT file is zombie: ",outputFile)
        outputExists = False
    else:
        print ("ROOT file not zombie: ",outputFile)

    if (outputExists and not rootFile.Get(category_name+"_"+region)):
        print ("category not found: ",category_name+"_"+region,outputFile)
        outputExists = False
    else:
        print ("category found: ",category_name+"_"+region)
    
    if (outputExists and isMC and proc.find("HNL")>=0):
        for coupling in args.couplings:
            if not rootFile.Get(category_name+"_"+region+"/HNL_coupling_%i"%coupling):
                outputExists = False
                print ("HNL not found: ", category_name+"_"+region+"/HNL_coupling_%i"%coupling,outputFile)
                break
            else:
                print ("HNL found: ", category_name+"_"+region+"/HNL_coupling_%i"%coupling)
    if (outputExists and isData):
        if not rootFile.Get(category_name+"_"+region+"/data"):
            outputExists = False
            print ("Data not found: ", category_name+"_"+region+"/data",outputFile)
        else:
            print ("Data found: ", category_name+"_"+region+"/data")
    rootFile.Close()
    if (outputExists):
        print("output file exists -> skip")
        sys.exit(0)

with open("../config/samples.yml") as samples_file:
    samples_dict = yaml.load(samples_file, Loader=yaml.FullLoader)
    subprocesses = samples_dict[proc]

#####################################
### Various configurations go here

# Systematic uncertainties
systematics_rates = {}

systematics_rates["IsoMuTrigger_weight_trigger"] = "trigger"
systematics_rates["tightMuons_weight_iso"] = "tight_muon_iso"
systematics_rates["tightMuons_weight_id"] = "tight_muon_id"
systematics_rates["tightMuons_weight_reco"] = "tight_muon_reco"
systematics_rates["tightElectrons_weight_id"] = "tight_electron_id"
systematics_rates["tightElectrons_weight_reco"] = "tight_electron_reco"
systematics_rates["puweight"] = "pu"

systematics_rates["(nominal_dR_l2j<0.4)*1.+(nominal_dR_l2j>0.4)*hnlJet_track_weight"] = "tagger_q"
systematics_rates["(nominal_dR_l2j>0.4)*1.+((nominal_dR_l2j<0.4)*subleadingLeptons_isElectron[0])+((nominal_dR_l2j<0.4)*subleadingLeptons_isMuon[0])*hnlJet_track_weight"] = "tagger_qmu"
systematics_rates["(nominal_dR_l2j>0.4)*1.+((nominal_dR_l2j<0.4)*subleadingLeptons_isMuon[0])+((nominal_dR_l2j<0.4)*subleadingLeptons_isElectron[0])*hnlJet_track_weight"] = "tagger_qe"

systematics_rates["pdf"] = "pdf"
systematics_rates["scale_shapeonly"] = "scale"
systematics_rates["looseMuons_weight_id"] = "loose_muon_id"
systematics_rates["looseMuons_weight_reco"] = "loose_muon_reco"
systematics_rates["looseElectrons_weight_id"] = "loose_electron_id"
systematics_rates["lepton2_track"] = "resolvedLepton_track_reco"
systematics_shapes = ["nominal", "jesTotal", "jer", "unclEn"]

#systematics_rates = {}
#systematics_shapes = ["nominal"]
####################################
if len(args.couplings)==0:
    # couplings to consider
    #couplings = range(1, 68)
    #couplings = [ 1]#, 2, 7, 12, 47, 52]
    couplings = [1,2,7,12,47,52,67]
    #couplings = [12]
    print ("Using default couplings: "+str(couplings))
else:
    couplings = list(set(args.couplings)) #remove duplications
    print ("Using couplings from options: "+str(couplings))


category_file = '../config/categories_2l_inclusive.json'
threshold_file = f'../config/coordsBestThresholds_{year}.json'

with open(category_file, 'r') as fp:
    categories_2l = json.load(fp)

with open(threshold_file, 'r') as fp:
    threshold_dict = json.load(fp)

category_cut = categories_2l[category_name]["varexp"]
thresholds = {}
thresholds["prompt"] = threshold_dict[category_name+"_prompt"]
thresholds["medium"]= threshold_dict[category_name+"_medium"]
thresholds["displaced"]= threshold_dict[category_name+"_displaced"]


dilepton_category_dict = {}
dilepton_category_dict[1] = "ql, dxysig<3"
dilepton_category_dict[2] = "ql, 3<dxysig<10"
dilepton_category_dict[3] = "ql, dxysig>10"
dilepton_category_dict[4] = "q, dxysig<3"
dilepton_category_dict[5] = "q, 3<dxysig<10"
dilepton_category_dict[6] = "q, dxysig>10"

# Process configuration
if "HNL" in proc:
    if "_ntau_" in proc:
        process = Process("HNL", proc)
        for subprocess in subprocesses: 
            process.Add(Sample(subprocess, ntuple_path, ["{}-{}".format(subprocess, year)], year=year, limits=True))
    else:
        process = Process("HNL", proc)
        process.Add(Sample(proc, ntuple_path, ["{}-{}".format(proc, year)], year=year, limits=True))
else:
    process = Process(proc, proc)
    subprocesses = subprocesses[int(year)]
    for sample_name, sample_list in subprocesses.items():
        print(sample_name)
        sample = Sample(sample_name, ntuple_path, sample_list, year=year, oneFile=oneFile, isMC=isMC)
        process.Add(sample)


# Event weights: MC only
if year=="2016":

    if category_name == "mumu_OS":
        fractionThresPromptBoosted = 3.
        fractionThresDisplacedBoosted = 3.
    else: 
        fractionThresPromptBoosted = 3.
        fractionThresDisplacedBoosted = 4.

elif year == "2017":
    if category_name == "mumu_OS":
        fractionThresPromptBoosted = 2.
        fractionThresDisplacedBoosted = 2.
    else:
        fractionThresPromptBoosted = 2.
        fractionThresDisplacedBoosted = 4.
elif year == "2018":
    if category_name == "mumu_OS":
        fractionThresPromptBoosted = 2.
        fractionThresDisplacedBoosted = 2.
    else: 
        fractionThresPromptBoosted = 2.
        fractionThresDisplacedBoosted = 4.

print ("boosted prompt frac threshold: ",fractionThresPromptBoosted)

if isMC:
    for syst, abrv in systematics_rates.items():
        for variation in ["Up", "Down"]:
            if "HNL" in process.name:
                for coupling in couplings:
                    process.Define("weightHNL_{}_{}{}".format(coupling, abrv, variation), "weightNominalHNL_{}/({})*({})".format(coupling, syst+"_nominal", syst+"_"+variation.lower()))
            else:
                process.Define("weight_{}{}".format(abrv, variation), "weightNominal/({})*({})".format(syst+"_nominal", syst+"_"+variation.lower()))
    for syst in systematics_shapes:
        for variation in ["Up", "Down"]:
            if "nominal" in syst:
                variation = ""
            # Define resolved and merged categories & tagger score variable & mass cuts
            #process.Define(f"category_{syst}_index", f"1.*({syst}_dR_l2j<0.4) \
                                                    #+ 2.*({syst}_dR_l2j>0.4 and {syst}_dR_l2j<1.3)")
            
            

            process.Define(f"category_{syst+variation}_index", f"1.*({syst+variation}_dR_l2j<0.4 and subleadingLeptons_dxysig[0] < 3. and hnlJet_{syst+variation}_ptorig/subleadingLeptons_pt[0] < {fractionThresPromptBoosted}  and subleadingLeptons_dxyErr[0]< 0.05 ) \
                                                    + 2.*({syst+variation}_dR_l2j<0.4 and subleadingLeptons_dxysig[0]>3. and subleadingLeptons_dxysig[0]<10. and hnlJet_{syst+variation}_ptorig/subleadingLeptons_pt[0] < 4.)\
                                                    + 3.*({syst+variation}_dR_l2j<0.4 and subleadingLeptons_dxysig[0]>10. and hnlJet_{syst+variation}_ptorig/subleadingLeptons_pt[0] < 4.) \
                                                    + 4.*({syst+variation}_dR_l2j>0.4 and {syst+variation}_dR_l2j<1.3 and subleadingLeptons_dxysig[0]< 3. and hnlJet_{syst+variation}_ptorig/subleadingLeptons_pt[0] < 4.) \
                                                    + 5.*({syst+variation}_dR_l2j>0.4 and {syst+variation}_dR_l2j<1.3 and subleadingLeptons_dxysig[0]>3. and subleadingLeptons_dxysig[0]<10. and hnlJet_{syst+variation}_ptorig/subleadingLeptons_pt[0] < 4.) \
                                                    + 6.*({syst+variation}_dR_l2j>0.4 and {syst+variation}_dR_l2j<1.3  and subleadingLeptons_dxysig[0]>10. and hnlJet_{syst+variation}_ptorig/subleadingLeptons_pt[0] < 4.) ")
           

            process.Define(f"tagger_score_{syst+variation}", tagger_compound_variable(syst+variation, single_lepton=False))
             
           
else:

   
    process.Define(f"category_nominal_index", f"1.*(nominal_dR_l2j<0.4 and subleadingLeptons_dxysig[0]<3.  and hnlJet_nominal_ptorig/subleadingLeptons_pt[0] < {fractionThresPromptBoosted}  and subleadingLeptons_dxyErr[0]< 0.05  ) \
                                                + 2.*(nominal_dR_l2j<0.4 and subleadingLeptons_dxysig[0]>3. and subleadingLeptons_dxysig[0]<10. and hnlJet_nominal_ptorig/subleadingLeptons_pt[0] < 4. )\
                                                + 3.*(nominal_dR_l2j<0.4 and subleadingLeptons_dxysig[0]>10.  and hnlJet_nominal_ptorig/subleadingLeptons_pt[0] < 4.) \
                                                + 4.*(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3 and subleadingLeptons_dxysig[0]<3. and hnlJet_nominal_ptorig/subleadingLeptons_pt[0] < 4. ) \
                                                + 5.*(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3 and subleadingLeptons_dxysig[0]>3. and subleadingLeptons_dxysig[0]<10. and hnlJet_nominal_ptorig/subleadingLeptons_pt[0] < 4.) \
                                                + 6.*(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3 and subleadingLeptons_dxysig[0]>10. and hnlJet_nominal_ptorig/subleadingLeptons_pt[0] < 4.)")
                                                
   
    
    process.Define(f"tagger_score_nominal", tagger_compound_variable(syst="nominal", single_lepton=False))


# create root file with nominal value histogram and various systematic variations
# to be used with Combine Harvester


category_dict = dilepton_category_dict
category_variable_nominal = "category_nominal_index"

histsList = []

coupling = 0
while coupling < 67:
    
    # different scenarios
    # Need to calculate yield per coupling
    if "HNL" in process.name:
        coupling += 1
        if coupling not in couplings:
            continue
        print (f"Processing coupling {coupling}")
    else:
        coupling = 68
    
    if isMC:
        hists = make_hists(process, systematics_shapes, systematics_rates, category_cut, category_variable_nominal, thresholds, region,  coupling=coupling)
        for name, hist in hists.items():
            hist.SetDirectory(0)
            histsList.append({"hist":hist, "name": name, "isMC": True})
            #write_hist(hist, category_dict, name, isMC=True)
    else:
        hists = make_hists(process, ["nominal"], None, category_cut, category_variable_nominal, thresholds, region,  coupling=coupling)
        for name, hist in hists.items():
            hist.SetDirectory(0)
            histsList.append({"hist":hist, "name": "data", "isMC": False})
            #write_hist(hist, category_dict, "data", isMC=False)

root_file = ROOT.TFile.Open(outputFile, "RECREATE")
print("The category name and cut are:", category_name, category_cut)
root_file.cd()
root_file.mkdir(category_name+"_"+region)
root_file.cd(category_name+"_"+region)

if "HNL" in proc:
    for coupling in couplings:
        histNominal = list(filter(lambda x: x['name']=="HNL_coupling_"+str(coupling), histsList))[0]['hist']
        print (coupling)
        for syst in list(systematics_rates.values())+systematics_shapes:
            if syst=="nominal":
                continue
            histUp = list(filter(lambda x: x['name']=="HNL_coupling_"+str(coupling)+"_"+syst+str(year)+"Up", histsList))[0]['hist']
            histDown = list(filter(lambda x: x['name']=="HNL_coupling_"+str(coupling)+"_"+syst+str(year)+"Down", histsList))[0]['hist']

            for ibin in range(histNominal.GetNbinsX()):
                c = histNominal.GetBinContent(ibin+1)
                u = histUp.GetBinContent(ibin+1)
                d = histDown.GetBinContent(ibin+1)
                if c<2.*EMPTY_BIN_FILL:
                    u = c
                    d = c
                else:
                    if u<2.*EMPTY_BIN_FILL:
                        u = c
                    if d<2.*EMPTY_BIN_FILL:
                        d = c
                histUp.SetBinContent(ibin+1,u)
                histDown.SetBinContent(ibin+1,d)

for histDict in histsList:
    #histDict['hist'].SetDirectory(root_file)
    write_hist(histDict['hist'], category_dict, histDict['name'], isMC=histDict['isMC'])
root_file.Close()




