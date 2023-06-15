import ROOT
import os
import json
from histo.lumi import lumi

def find_xsec(path, xsecs):
    for key, val in xsecs.items():
        if key in path:
            return val

# This class prepares a given sample by scaling to int. luminosity
class Sample:
    """
    Makes a "sample" object which is stored as a ROOT RDataFrame
    Calculate event weights based on cross section and yields 
    """
    def __init__(self, name, ntuple_path, paths, isMC=True, year="2016", cut=None, limits=False, oneFile=False):
        if os.path.exists("/vols/cms/LLP"):
            with open(os.path.join("/vols/cms/LLP/yields_201117", year, "eventyields.json")) as json_file:
                yields = json.load(json_file)
            with open(os.path.join("/vols/cms/LLP/yields_201117", year, "eventyieldsHNL.json")) as json_file:
                self.yieldsHNL = json.load(json_file)        
            with open(os.path.join("/vols/cms/LLP/yields_230309", year, "eventyieldsHNL.json")) as json_file:
                self.yieldsHNL.update(json.load(json_file))

            with open("/vols/cms/LLP/filterTable.json") as json_file:
                gen_filter = json.load(json_file)
            with open("/vols/cms/LLP/filterLPairTable.json") as json_file:
                gen_filter.update(json.load(json_file))
                
            with open("/vols/cms/LLP/gridpackLookupTable.json") as lookup_table_file:
                lookup_table = json.load(lookup_table_file)
            with open("/vols/cms/LLP/xsec.json") as xsec_file:
                xsecs = json.load(xsec_file)
                
        if os.path.exists("/nfs/dust/cms/user/mkomm/HNL"):
            with open(os.path.join("/nfs/dust/cms/user/mkomm/HNL/LLP/yields_201117", year, "eventyields.json")) as json_file:
                yields = json.load(json_file)
            with open(os.path.join("/nfs/dust/cms/user/mkomm/HNL/LLP/yields_201117", year, "eventyieldsHNL.json")) as json_file:
                self.yieldsHNL = json.load(json_file)      
            with open(os.path.join("/nfs/dust/cms/user/mkomm/HNL/LLP/yields_230309", year, "eventyieldsHNL.json")) as json_file:
                self.yieldsHNL.update(json.load(json_file))  
            with open("/nfs/dust/cms/user/mkomm/HNL/LLP/filterTable.json") as json_file:
                gen_filter = json.load(json_file)
            with open("/nfs/dust/cms/user/mkomm/HNL/LLP/filterLPairTable.json") as json_file:
                gen_filter.update(json.load(json_file))
                
            with open("/nfs/dust/cms/user/mkomm/HNL/LLP/gridpackLookupTable.json") as lookup_table_file:
                lookup_table = json.load(lookup_table_file)
            with open("/nfs/dust/cms/user/mkomm/HNL/LLP/xsec.json") as xsec_file:
                xsecs = json.load(xsec_file)  
                            
        self.name = name
        self.file_list = ROOT.std.vector('string')()
        self.sum_weight = 0
        self.isMC = isMC
        if oneFile:
            counter = 0
        for path in paths:
            for f in os.listdir(os.path.join(ntuple_path, path)):
                self.file_list.push_back(os.path.join(ntuple_path, path, f))
                if oneFile:
                    counter +=1
                    if counter > 20:
                        break
            if self.isMC:
                if "HNL" not in name:
                    self.sum_weight += yields[path]["weighted"]
                else:
                    self.sum_weight=1. #dummy value 
        self.rdf = ROOT.RDataFrame("Friends", self.file_list)
        count = self.rdf.Count().GetValue()
        print (name,count)
        #if count > 0:
        if cut is not None:
            self.rdf = self.rdf.Define("sample_cut", cut)
            self.rdf = self.rdf.Filter("sample_cut")
        #selected = self.rdf.Count().GetValue()

        #print("RDF {} has entries {}/{}".format(name, selected, count))

        if self.isMC:
            if "HNL" in name:
                if not limits:
                    if "ntau" in name:
                        lu_infos = lookup_table[name.replace('ntau0', 'all').replace('ntau1', 'all').replace('ntau2', 'all')]['weights'][str(int(coupling))]
                    elif "pt20" in name:
                        lu_infos = lookup_table[name.replace('pt20', 'all')]['weights'][str(int(coupling))]
                    else:
                        lu_infos = lookup_table[name]['weights'][str(int(coupling))]
                    xsec = lu_infos['xsec']['nominal']
                else:
                    xsec = 1.
            else:
                xsec = find_xsec(path, xsecs)

            #*\
            self.rdf = self.rdf.Define("weightNominal", f"IsoMuTrigger_weight_trigger_nominal*\
                tightMuons_weight_iso_nominal*tightMuons_weight_id_nominal*tightMuons_weight_reco_nominal*\
                tightElectrons_weight_id_nominal*tightElectrons_weight_reco_nominal*\
                looseMuons_weight_reco_nominal*looseMuons_weight_id_nominal*\
                looseElectrons_weight_id_nominal*\
                puweight_nominal*genweight*hnlJet_track_weight_nominal*lepton2_track_nominal*\
                {lumi[year]}*1000.0*{xsec}/{self.sum_weight}") #NB: weight will be 1. for HNL
            

            #self.rdf = self.rdf.Define("weightNominalCorrectedUp", "weightNominal*hnlJet_track_weight_adapted_nominal")

            if "HNL" in name:
                for coupling in range(1, 68):
                    if ("pt20" in name) or ("ntau" in name):
                        filtereff = gen_filter[name]['weights'][str(coupling)]['eff']
                        weightNormSum = self.yieldsHNL[name+"-"+str(year)]['LHEWeights_coupling_'+str(coupling)]
                        if filtereff<1e-8 and weightNormSum<1e-8:
                            filtereff = 0.0
                            weightNormSum = 1.0
                        self.rdf = self.rdf.Define("weightNominalHNL_{}".format(coupling), f"{filtereff}*weightNominal*LHEWeights_coupling_{coupling}/{weightNormSum}")
                    else:
                        weightNormSum = self.yieldsHNL[name+"-"+str(year)]['LHEWeights_coupling_'+str(coupling)]
                        self.rdf = self.rdf.Define("weightNominalHNL_{}".format(coupling), f"weightNominal*LHEWeights_coupling_{coupling}/{weightNormSum}")

        else:
            self.rdf = self.rdf.Define("weightNominal", "1")
