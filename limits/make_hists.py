import ROOT
import os
import yaml
import json
import argparse
from histo import Process, Sample

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
        remove_neg_entries(hist_limits)
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


def mass_cut(delta_m:float=5., region:str="D", syst:str="nominal", single_lepton=False) -> str:
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
    m_upper = mW + delta_m
    m_lower = mW - delta_m

    if region not in ["A", "B", "C", "D"]:
        raise ValueError(f"Invalid region {region} selected")

    if region == "B" or region == "C":
        return  f'({syst}_m_llj>{m_upper} or {syst}_m_llj<{m_lower}) '
    elif region == "A" or region == "D":
        return  f'({syst}_m_llj<{m_upper} and {syst}_m_llj>{m_lower}) '
    else:
        return ""

def tagger_cut(tagger_threshold: float, lower_threshold: float=0.25, region:str="D", syst:str="nominal") -> str:
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
        return f'tagger_score_{syst} > {lower_threshold} and tagger_score_{syst} < {tagger_threshold}'
    elif region == "C" or region == "D":
        return f'tagger_score_{syst} > {tagger_threshold}'


def tagger_compound_variable(syst:str="nominal", single_lepton=False) -> str:
    """ Compound tagger variable to facilitate categorisation """
    # if single_lepton:
    #     return f"hnlJet_{syst}_llpdnnx_ratio_LLP_Q"
    # else:
    return f"({syst}_dR_l2j>0.4 and {syst}_dR_l2j<1.3)*hnlJet_{syst}_llpdnnx_ratio_LLP_Q+\
        ({syst}_dR_l2j<0.4)*hnlJet_{syst}_llpdnnx_ratio_LLP_QMU*subleadingLeptons_isMuon[0]+\
        +({syst}_dR_l2j<0.4)*hnlJet_{syst}_llpdnnx_ratio_LLP_QE*subleadingLeptons_isElectron[0]"

def make_hists(process, systematics_shapes, systematics_rates, cut_nominal, category_variable_nominal, thresholds, region, coupling=None):

    hists = {}

    def make_hist(process, category_variable, thresholds, weight, cut, region, syst="nominal"):
        threshold_merged, deltam_merged = thresholds["merged"]
        threshold_resolved, deltam_resolved = thresholds["resolved"]
        #print(f"Tagger thresholds, merged: {threshold_merged}, resolved: {threshold_resolved}")
        #print(f"DeltaM, merged: {deltam_merged}, resolved: {deltam_resolved}")
        mass_cut_merged = mass_cut(delta_m=deltam_merged, region=region, syst=syst)
        mass_cut_resolved = mass_cut(delta_m=deltam_resolved, region=region, syst=syst)
        tagger_cut_merged = tagger_cut(threshold_merged, region=region, syst=syst)
        tagger_cut_resolved = tagger_cut(threshold_resolved, region=region, syst=syst)

        cut += f"and ( ({category_variable}==1 and {mass_cut_merged} and {tagger_cut_merged}) or ({category_variable}==2 and {mass_cut_resolved} and {tagger_cut_resolved}) )"
        print(syst, category_variable, cut, weight)

        hist_nano = process.Histo1D((category_variable, category_variable, 2, 0.5, 2.5), category_variable, cut=cut, weight=weight)
        hist_nano = hist_nano.Clone()
        return hist_nano
    if systematics_rates is not None:
    # variations with constant shape but changing weight
        for syst, abrv in systematics_rates.items():
            for variation in ["Up", "Down"]:
                if "HNL" in process.name:
                    name = f"{process.name}_coupling_{coupling}_{abrv}{variation}"
                    weight = f"weightHNL_{coupling}_{abrv}{variation}"
                else:
                    name = f"{process.name}_{abrv}{variation}"
                    weight = f"weight_{abrv}{variation}"
                hists[name] = make_hist(process, category_variable_nominal, thresholds, weight, cut_nominal, region, syst="nominal")

    if systematics_shapes is not None:

        for syst in systematics_shapes:
            if "HNL" in process.name:
                name = f"{process.name}_coupling_{coupling}"
                weight = f"weightNominalHNL_{coupling}"
            else:
                name = process.name
                weight = "weightNominal"

            # add name for variations
            if syst != "nominal":
                name += f"_{syst}"

            # Systematic variation -- replace nominal by systematic in all cuts
            cut = cut_nominal.replace("nominal", syst)
            cut = cut.replace("nselectedJets_unclEnUp", "nselectedJets_nominal") #Hack!
            cut = cut.replace("nselectedJets_unclEnDown", "nselectedJets_nominal") #Hack!
            category_variable = category_variable_nominal.replace("nominal", syst)
            # read in hist from nanoAOD friends
            hists[name] = make_hist(process, category_variable, thresholds, weight, cut, region, syst=syst)

    return hists

# make histograms per year, process
parser = argparse.ArgumentParser()
parser.add_argument("--year",default="2016")
parser.add_argument("--proc", default="wjets")
parser.add_argument("--category", default="mumu_OS_displaced")
parser.add_argument("--region", default="D")
parser.add_argument("--ntuple_path", default="/vols/cms/vc1117/LLP/nanoAOD_friends/HNL/28May21")
parser.add_argument("--output_path", default="hists")
parser.add_argument("--data", action="store_true", default=False)
parser.add_argument("--test", action="store_true", dest="oneFile", default=False)

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
systematics_rates["tightElectrons_weight_id"] = "tight_electron_id"
systematics_rates["tightElectrons_weight_reco"] = "tight_electron_reco"
systematics_rates["looseElectrons_weight_reco"] = "loose_electron_reco"
systematics_rates["puweight"] = "pu"

systematics_shapes = ["nominal", "jesTotalUp", "jesTotalDown", "jerUp", "jerDown", "unclEnUp", "unclEnDown"]

####################################

# couplings to consider
couplings = range(2, 68)
couplings = [2, 7, 12, 47, 52]

category_file = '../config/categories_2l.json'
threshold_file = f'../config/coordsBestThresholds_{year}.json'

with open(category_file, 'r') as fp:
    categories_2l = json.load(fp)

with open(threshold_file, 'r') as fp:
    threshold_dict = json.load(fp)

category_cut = categories_2l[category_name]["varexp"]
thresholds = threshold_dict[category_name]

dilepton_category_dict = {}
dilepton_category_dict[1] = "ql" # Merged
dilepton_category_dict[2] = "q" # Resolved

# Process configuration
if "HNL" in proc:
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
if isMC:
    for syst, abrv in systematics_rates.items():
        for variation in ["Up", "Down"]:
            if "HNL" in process.name:
                for coupling in couplings:
                    process.Define("weightHNL_{}_{}{}".format(coupling, abrv, variation), "weightNominalHNL_{}/{}*{}".format(coupling, syst+"_nominal", syst+"_"+variation.lower()))
            else:
                process.Define("weight_{}{}".format(abrv, variation), "weightNominal/{}*{}".format(syst+"_nominal", syst+"_"+variation.lower()))
    for syst in systematics_shapes:
        # Define resolved and merged categories & tagger score variable & mass cuts
        process.Define(f"category_{syst}_index", f"1.*({syst}_dR_l2j<0.4) \
                                                 + 2.*({syst}_dR_l2j>0.4 and {syst}_dR_l2j<1.3)")
        process.Define(f"tagger_score_{syst}", tagger_compound_variable(syst, single_lepton=False))

else:
    process.Define("category_nominal_index", "1.*(nominal_dR_l2j<0.4) \
                                             + 2.*(nominal_dR_l2j>0.4 and nominal_dR_l2j<1.3)")
    process.Define(f"tagger_score_nominal", tagger_compound_variable(syst="nominal", single_lepton=False))


# create root file with nominal value histogram and various systematic variations
# to be used with Combine Harvester
root_file = ROOT.TFile.Open(os.path.join(output_path, f"{proc}_{args.category}_{region}_{year}.root"), "RECREATE")
print("The category name and cut are:", category_name, category_cut)
root_file.cd()
root_file.mkdir(category_name+"_"+region)
root_file.cd(category_name+"_"+region)

category_dict = dilepton_category_dict
category_variable_nominal = "category_nominal_index"

coupling = 1
while coupling < 67:
    # different scenarios
    # Need to calculate yield per coupling
    if "HNL" in process.name:
        coupling += 1
        if coupling not in couplings:
            continue
    else:
        coupling = 68

    if isMC:
        hists = make_hists(process, systematics_shapes, systematics_rates, category_cut, category_variable_nominal, thresholds, region, coupling=coupling)
        for name, hist in hists.items():
            write_hist(hist, category_dict, name, isMC=True)
    else:
        hists = make_hists(process, ["nominal"], None, category_cut, category_variable_nominal, thresholds, region, coupling=coupling)
        for name, hist in hists.items():
            write_hist(hist, category_dict, "data", isMC=False)

root_file.Close()