import ROOT
import math
import random
import os
import sys
import numpy as np
from histo import style

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetErrorX(0)
ROOT.gStyle.SetPaperSize(7.0*1.35,6.7*1.35)

colors = []
    
def newColorRGB(red,green,blue):
    color=ROOT.TColor(red,green,blue,1.0)
    colors.append(color)
    return color

newColorRGB.colorindex = 301

mcStyles = {
    'topbkg': {'legend': "t#bar{t} / t", 'fill': newColorRGB(1.,0.8,0.0)},
    'wjets': {'legend': "W+jets", 'fill': newColorRGB(0.33,0.75,0.35)},
    'dyjets':{'legend': "Z/#gamma*+jets", 'fill': newColorRGB(0.3,0.75,0.95)},
    'vgamma': {'legend': "V#gamma*+jets", 'fill': newColorRGB(0.73,0.25,0.96)},
    'qcd': {'legend': "Multijet", 'fill': newColorRGB(0.85,0.85,0.85)},
    'nonisoqcd': {'legend': "Multijet", 'fill': newColorRGB(0.85,0.85,0.85)}
}

for name,mcStyle in mcStyles.items():
    fillColor = mcStyle['fill']
    mcStyle['line'] = newColorRGB(
        fillColor.GetRed()*0.6,
        fillColor.GetGreen()*0.6,
        fillColor.GetBlue()*0.6,
    )

#bins: 0=underflow, 1=empty, 
# 2=mumu OS, 3=ee OS, 4=mue OS, 5=emu OS
# 6=mumu SS, 7=ee SS, 8=mue SS, 9=emu SS

def combineAll(h):
    return h.ProjectionX(h.GetName()+"tot")

def combineOS(h):
    return h.ProjectionX(h.GetName()+"tot",2,5)

def combineSS(h):
    return h.ProjectionX(h.GetName()+"tot",6,9)
    
    
def smooth(h,factor=0.2,n=1):
    hsmooth = h.Clone(h.GetName()+"smooth")
    for it in range(n):
        for ibin in range(2,h.GetNbinsX()-2):
            c1 = h.GetBinContent(ibin)
            e1 = h.GetBinError(ibin)
            
            c2 = h.GetBinContent(ibin+1)
            e2 = h.GetBinError(ibin+1)
            
            c3 = h.GetBinContent(ibin+2)
            e3 = h.GetBinError(ibin+2)
            
            c = factor*c1+factor*c3+(1.-2*factor)*c2
            
            e = 0.
            if c1>1e-3:
                e += factor*e1/c1*c
            if c3>1e-3:
                e += factor*e3/c3*c
            if c2>1e-3:
                e += (1.-2*factor)*e2/c2*c
            
            hsmooth.SetBinContent(ibin+1,c)
            hsmooth.SetBinError(ibin+1,e)
            h.SetBinContent(ibin+1,c)
            h.SetBinError(ibin+1,e)
    return h

class Plot():
    def __init__(self,
        plot,
        title,
        combine = combineAll,
        extraTitles=[],
        logy=False,
        yRange= None,
        yspace=1.2,
        unit="",
        binRange = [0,1],
        rebin = 1,
        oneLegend=True,
        outputSuffix = "",
        header="(ee,#kern[-0.5]{ }e#mu,#kern[-0.5]{ }#mue,#kern[-0.5]{ }#mu#mu)#kern[-0.2]{ }+#kern[-0.2]{ }jets",
        procs = ['topbkg','wjets','dyjets','vgamma','nonisoqcd'],
        showData=True,
        path='/vols/cms/mkomm/HNL/histo/plot_paper/hists_met/hist'
    ):
        self.plot = plot
        self.title = title
        self.combine = combine
        self.extraTitles = extraTitles
        self.logy = logy
        self.yspace = yspace
        self.yRange = yRange
        self.unit = unit
        self.binRange = binRange
        self.rebin = rebin
        self.oneLegend = oneLegend
        self.outputSuffix = outputSuffix
        self.header = header
        self.procs = procs
        self.showData = showData
        self.path = path
        
        self.signals = []

    def addSignal(self, name, scale, legend, style=[2,2,ROOT.kRed+1]):
        self.signals.append({
            'name':name,
            'scale': scale,
            'legend': legend,
            'style': style
        })
        
    def getMCSignal(self,signal,syst='nominal'):
        mcHistSum = None
        for year in [2016,2017,2018]:
            filepath = self.path+"_"+str(year)+"_"+self.plot+"_"+syst+".root"
            f = ROOT.TFile(filepath)
            if f==None:
                raise Exception("Cannot open file '%s'"%filepath)
                
            hist = f.Get(signal['name'])
            hist = self.combine(hist)
            hist.SetDirectory(0)
            if self.rebin>1:
                hist.Rebin(self.rebin)
            hist.SetFillStyle(0)
            hist.SetLineWidth(signal['style'][0])
            hist.SetLineStyle(signal['style'][1])
            hist.SetLineColor(signal['style'][2])
    
            if mcHistSum==None:
                mcHistSum = hist.Clone(self.plot+str(random.random())+"sum")
                mcHistSum.SetDirectory(0)
            else:
                mcHistSum.Add(hist)
            f.Close()
        mcHistSum.Scale(signal['scale'])
        return mcHistSum
        
    def getMCStackSumDict(self,syst='nominal'):
        mcHistDict = {}
        mcHistSum = None
        for year in [2016,2017,2018]:
            filepath = self.path+"_"+str(year)+"_"+self.plot+"_"+syst+".root"
            f = ROOT.TFile(filepath)
            if f==None:
                raise Exception("Cannot open file '%s'"%filepath)
                
            for i,sample in enumerate(['topbkg','wjets','dyjets','vgamma','qcd','nonisoqcd']):
                hist = f.Get(sample)
                hist = self.combine(hist)
                hist.SetDirectory(0)
                
                hist.SetFillColor(mcStyles[sample]['fill'].GetNumber())
                hist.SetLineColor(mcStyles[sample]['line'].GetNumber())
                hist.SetLineWidth(2)
                
                if sample not in mcHistDict.keys():
                    mcHistDict[sample] = hist.Clone(self.plot+sample+str(random.random()))
                    mcHistDict[sample].SetDirectory(0)
                else:
                    mcHistDict[sample].Add(hist)
                    
                
            f.Close()
            
        mcHistDict['nonisoqcd'].Scale(mcHistDict['qcd'].Integral()/mcHistDict['nonisoqcd'].Integral())
        
        smooth(mcHistDict['qcd'],0.2,5)
        #smooth(mcHistDict['nonisoqcd'],0.2,5)
        
        if self.rebin>1:
            for hist in mcHistDict.values():
                hist.Rebin(self.rebin)
        
        mcStack = ROOT.THStack()
        for sample in self.procs:
            mcStack.Add(mcHistDict[sample])
            if mcHistSum==None:
                mcHistSum = mcHistDict[sample].Clone(self.plot+str(random.random())+"sum")
                mcHistSum.SetDirectory(0)
            else:
                mcHistSum.Add(mcHistDict[sample])
                
        return mcStack,mcHistDict,mcHistSum
        
    def getMC(self):
        mcStackNominal, mcHistDictNominal, mcHistSumNominal = self.getMCStackSumDict()
        
        totalUnc2 = np.zeros(mcHistSumNominal.GetNbinsX())
        for syst in ['jesTotal','jer','unclEn','pu','muEff','track']:
            _, _, mcHistSumSystUp = self.getMCStackSumDict(syst+"Up")
            _, _, mcHistSumSystDown = self.getMCStackSumDict(syst+"Down")
            for ibin in range(mcHistSumNominal.GetNbinsX()):
                nom = mcHistSumNominal.GetBinContent(ibin+1)
                upDiff = math.fabs(mcHistSumSystUp.GetBinContent(ibin+1)-nom)
                downDiff = math.fabs(mcHistSumSystDown.GetBinContent(ibin+1)-nom)
                diff = 0.5*(upDiff+downDiff) #take average
                totalUnc2[ibin]+=diff**2
                
        for ibin in range(mcHistSumNominal.GetNbinsX()):
            mcHistSumNominal.SetBinError(ibin+1,math.sqrt(mcHistSumNominal.GetBinError(ibin+1)**2+totalUnc2[ibin])) #add mcstat err
                
        return mcStackNominal, mcHistDictNominal, mcHistSumNominal
        
    def getData(self):
        dataHistSum = None
        for year in [2016,2017,2018]:
            filepath = self.path+"_"+str(year)+"_"+self.plot+"_nominal.root"
            f = ROOT.TFile(filepath)
            if f==None:
                raise Exception("Cannot open file '%s'"%filepath)
                
            for sample in ['electron','muon']:
                hist = f.Get(sample)
                hist = self.combine(hist)
                hist.SetDirectory(0)
                hist.SetMarkerStyle(20)
                hist.SetMarkerSize(1.5)
                hist.SetLineColor(ROOT.kBlack)
                hist.SetMarkerColor(ROOT.kBlack)
                if dataHistSum==None:
                    dataHistSum = hist.Clone(self.plot+str(random.random())+"sum")
                    dataHistSum.SetDirectory(0)
                else:
                    dataHistSum.Add(hist)
            f.Close()
            
        if self.rebin>1:
            dataHistSum.Rebin(self.rebin)
            
        return dataHistSum

    def __call__(self):
        
        cv = style.makeCanvas("cv"+str(random.random()),700,670)
        cv.Divide(1,2,0,0)
        cv.GetPad(1).SetPad(0.0, 0.0, 1.0, 1.0)
        cv.GetPad(1).SetFillStyle(4000)
        cv.GetPad(2).SetPad(0.0, 0.00, 1.0,1.0)
        cv.GetPad(2).SetFillStyle(4000)

        cvxmin=0.13
        cvxmax=0.96
        cvymin=0.13
        cvymax=0.92
        resHeight=0.35

        for i in range(1,3):
            #for the canvas:
            cv.GetPad(i).SetBorderMode(0)
            cv.GetPad(i).SetGridx(False)
            cv.GetPad(i).SetGridy(False)


            #For the frame:
            cv.GetPad(i).SetFrameBorderMode(0)
            cv.GetPad(i).SetFrameBorderSize(1)
            cv.GetPad(i).SetFrameFillColor(0)
            cv.GetPad(i).SetFrameFillStyle(0)
            cv.GetPad(i).SetFrameLineColor(1)
            cv.GetPad(i).SetFrameLineStyle(1)
            cv.GetPad(i).SetFrameLineWidth(1)

            # Margins:
            cv.GetPad(i).SetLeftMargin(cvxmin)
            cv.GetPad(i).SetRightMargin(1-cvxmax)

            # For the Global title:
            cv.GetPad(i).SetTitle("")

            # For the axis:
            cv.GetPad(i).SetTickx(1)  # To get tick marks on the opposite side of the frame
            cv.GetPad(i).SetTicky(1)

            # Change for log plots:
            cv.GetPad(i).SetLogx(0)
            cv.GetPad(i).SetLogy(0)
            cv.GetPad(i).SetLogz(0)

        cv.GetPad(2).SetTopMargin(1-cvymax)
        cv.GetPad(2).SetBottomMargin(resHeight)
        cv.GetPad(1).SetTopMargin(1-resHeight)
        cv.GetPad(1).SetBottomMargin(cvymin)

        cv.cd(2)
        
        
        mcStackNominal, mcHistDictNominal, mcHistSumNominal = self.getMC()
        
        for sample in mcHistDictNominal.keys():
            print ("   %s: %.1f"%(sample,mcHistDictNominal[sample].Integral()))
        
        if self.unit!="":
            xaxisTitle = self.title+" ("+self.unit+")"
            
            binwidth = mcHistSumNominal.GetBinWidth(1)

            yaxisTitle = "Events / %.0f %s"%(binwidth,self.unit)

        else:
            xaxisTitle = self.title
            yaxisTitle = "Events / bins"
        
        ymax = mcHistSumNominal.GetMaximum()

        if self.logy:
            ymin = 0.07
            ymax = 10**(self.yspace*math.log10(max([ymax,1.0])))
        
            if self.yRange is not None:
                ymin = self.yRange[0]
                ymax = self.yRange[1]
        
            axis=ROOT.TH2F(
                "axis"+str(random.randint(0,99999)),";;"+yaxisTitle,
                50,self.binRange[0],self.binRange[1],
                50,ymin,ymax
            )
        else:
            ymin = 0.0
            ymax = self.yspace*ymax
            
            if self.yRange is not None:
                ymin = self.yRange[0]
                ymax = self.yRange[1]
        
            axis=ROOT.TH2F(
                "axis"+str(random.randint(0,99999)),";;"+yaxisTitle,
                50,self.binRange[0],self.binRange[1],
                50,ymin,ymax
            )

        axis.GetXaxis().SetLabelSize(0)
        axis.GetXaxis().SetTitle("")
        axis.GetXaxis().SetTickLength(0.015/(1-cv.GetPad(2).GetLeftMargin()-cv.GetPad(2).GetRightMargin()))
        axis.GetYaxis().SetTickLength(0.015/(1-cv.GetPad(2).GetTopMargin()-cv.GetPad(2).GetBottomMargin()))
        #axis.GetYaxis().SetNoExponent(True)
        axis.Draw("AXIS")
        
        if self.logy:
            cv.GetPad(2).SetLogy(1)
            

        dataHist = self.getData()
        
        print ("   %s: %.1f"%('data',dataHist.Integral()))
        
        '''
        if 'nonisoqcd' in mcHistDictNominal.keys():
            totalMCSum = mcHistSumNominal.Integral()
            totalDataSum = dataHist.Integral()
            sumQCD = mcHistDictNominal['nonisoqcd'].Integral()
            
            scaleQCD = (totalDataSum-(totalMCSum-sumQCD))/sumQCD
            print ("   %s: %.3f"%('qcd scale',scaleQCD))
            mcHistDictNominal['nonisoqcd'].Scale(scaleQCD)
        '''
            

        mcStackNominal.Draw("HISTSame")

        rootObj=[]

        for ibin in range(mcHistSumNominal.GetNbinsX()):
            c = mcHistSumNominal.GetBinCenter(ibin+1)
            w = mcHistSumNominal.GetBinWidth(ibin+1)
            m = mcHistSumNominal.GetBinContent(ibin+1)
            err = mcHistSumNominal.GetBinError(ibin+1)
            if c>self.binRange[1] or c<self.binRange[0]:
                continue
            if m>0.0:
                box = ROOT.TBox(c-0.5*w,m-err,c+0.5*w,m+err)
                box.SetFillStyle(3345)
                box.SetLineColor(ROOT.kGray+1)
                box.SetFillColor(ROOT.kGray)
                box.SetLineWidth(2)
                rootObj.append(box)
                box.Draw("SameF")
                '''
                box2 = ROOT.TBox(c-0.5*w,m-err,c+0.5*w,m+err)
                box2.SetFillStyle(0)
                box2.SetLineColor(ROOT.kGray+1)
                box2.SetFillColor(ROOT.kGray)
                box2.SetLineWidth(2)
                rootObj.append(box2)
                box2.Draw("SameL")
                '''
        
        if self.showData:
            dataHist.Draw("HISTPESame")

        isignalEntry=0
        if len(self.signals)>0:
            
            legendSignal = ROOT.TLegend(cvxmin+0.025,cvymax-0.07-0.048*(len(self.extraTitles)),cvxmin+0.3,cvymax-0.07-0.048*(len(self.extraTitles)+sum(map(lambda x: len(x['legend']),self.signals))),"","NDC")
            legendSignal.SetBorderSize(0)
            legendSignal.SetFillStyle(0)
            legendSignal.SetTextFont(43)
            legendSignal.SetTextSize(25)
                
            for signal in self.signals:
                histSignal = self.getMCSignal(signal)
                rootObj.append(histSignal)
                for ientry,legendEntry in enumerate(signal['legend']):
                    if ientry==0:
                        legendSignal.AddEntry(histSignal,legendEntry,"L")
                    else:
                        legendSignal.AddEntry("","","")
                        pText = ROOT.TPaveText(cvxmin+0.035,cvymax-0.075-0.048*(len(self.extraTitles)+isignalEntry),cvxmin+0.035,cvymax-0.075-0.048*(len(self.extraTitles)+isignalEntry),"NDC")
                        rootObj.append(pText)
                        pText.SetBorderSize(0)
                        pText.SetFillStyle(0)
                        pText.SetTextFont(43)
                        pText.SetTextSize(23)
                        pText.SetTextAlign(13)
                        pText.AddText(legendEntry)
                        pText.Draw()
                    isignalEntry+=1
                histSignal.Draw("HISTSame")

                
        ROOT.gPad.RedrawAxis()

        if len(self.signals)>0:
            legendSignal.Draw()
            
        
        pCMS=ROOT.TPaveText(cvxmin+0.025,cvymax-0.025,cvxmin+0.025,cvymax-0.025,"NDC")
        pCMS.SetFillColor(ROOT.kWhite)
        pCMS.SetBorderSize(0)
        pCMS.SetTextFont(63)
        pCMS.SetTextSize(27)
        pCMS.SetTextAlign(13)
        pCMS.AddText("CMS")
        pCMS.Draw("Same")

        pLumi=ROOT.TPaveText(cvxmax,0.94,cvxmax,0.94,"NDC")
        pLumi.SetFillColor(ROOT.kWhite)
        pLumi.SetBorderSize(0)
        pLumi.SetTextFont(43)
        pLumi.SetTextSize(30)
        pLumi.SetTextAlign(31)
        if self.header!="":
            pLumi.AddText(self.header+", 138#kern[-0.5]{ }fb#lower[-0.7]{#scale[0.7]{-1}} (13 TeV)")
        else:
            pLumi.AddText("137#kern[-0.5]{ }fb#lower[-0.7]{#scale[0.7]{-1}} (13 TeV)")
        pLumi.Draw("Same")
        
        
        if self.oneLegend:
            legend = ROOT.TLegend(cvxmax-0.21,cvymax-0.02,cvxmax-0.02,cvymax-0.02-0.048*7,"","NDC")
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextFont(43)
            legend.SetTextSize(25)
            
            legend.AddEntry(dataHist,"Data","P")
            for sample in self.procs[:2]:
                legend.AddEntry(mcHistDictNominal[sample],mcStyles[sample]['legend'],"F")
            for sample in self.procs[2:]:
                legend.AddEntry(mcHistDictNominal[sample],mcStyles[sample]['legend'],"F")
                   
            legend.AddEntry(box,"Unc.","F")
            legend.Draw()
            
        else:
            legend1 = ROOT.TLegend(cvxmax-0.39,cvymax-0.02,cvxmax-0.21,cvymax-0.02-0.048*3,"","NDC")
            legend1.SetBorderSize(0)
            legend1.SetFillStyle(0)
            legend1.SetTextFont(43)
            legend1.SetTextSize(25)
            
            legend2 = ROOT.TLegend(cvxmax-0.21,cvymax-0.02,cvxmax-0.02,cvymax-0.02-0.048*4,"","NDC")
            legend2.SetBorderSize(0)
            legend2.SetFillStyle(0)
            legend2.SetTextFont(43)
            legend2.SetTextSize(25)
        
            legend1.AddEntry(dataHist,"Data","P")
            for sample in self.procs[:2]:
                legend1.AddEntry(mcHistDictNominal[sample],mcStyles[sample]['legend'],"F")
            for sample in self.procs[2:]:
                legend2.AddEntry(mcHistDictNominal[sample],mcStyles[sample]['legend'],"F")
                
                
            legend2.AddEntry(box,"Unc.","F")
            
            
            legend1.Draw()
            legend2.Draw()



        
        for i,text in enumerate(self.extraTitles):
            pText = ROOT.TPaveText(cvxmin+0.025,cvymax-0.08-0.048*i,cvxmin+0.025,cvymax-0.08-0.048*i,"NDC")
            rootObj.append(pText)
            pText.SetBorderSize(0)
            pText.SetFillStyle(0)
            pText.SetTextFont(43)
            pText.SetTextSize(25)
            pText.SetTextAlign(13)
            pText.AddText(text)
            pText.Draw()
            
        cv.cd(1)
        
        axisRes=ROOT.TH2F(
            "axisRes"+str(random.randint(0,99999)),";"+xaxisTitle+";Data/Pred.",
            50,self.binRange[0],self.binRange[-1],50,0.1,1.9)#0.55,1.45)#0.1,1.9)
        axisRes.GetYaxis().SetNdivisions(406)
        axisRes.GetXaxis().SetTickLength(0.020/(1-cv.GetPad(1).GetLeftMargin()-cv.GetPad(1).GetRightMargin()))
        axisRes.GetYaxis().SetTickLength(0.015/(1-cv.GetPad(1).GetTopMargin()-cv.GetPad(1).GetBottomMargin()))

        axisRes.Draw("AXIS")
        '''
        xaxisTitlePave = ROOT.TMathText()
        xaxisTitlePave.SetNDC()
        xaxisTitlePave.SetTextFont(43)
        xaxisTitlePave.SetTextAlign(31);
        xaxisTitlePave.SetTextSize(32);
        xaxisTitlePave.DrawMathText(cvxmax,0.03,)
        '''
        
        if self.showData:
            dataRes = dataHist.Clone(dataHist.GetName()+str(random.random())+"res")

        for ibin in range(mcHistSumNominal.GetNbinsX()):
            c = mcHistSumNominal.GetBinCenter(ibin+1)
            w = mcHistSumNominal.GetBinWidth(ibin+1)
            m = mcHistSumNominal.GetBinContent(ibin+1)
            if c>self.binRange[1] or c<self.binRange[0]:
                continue
                
            if self.showData:
                if m>1e-3:
                    dataRes.SetBinContent(ibin+1,dataHist.GetBinContent(ibin+1)/m)
                    dataRes.SetBinError(ibin+1,dataHist.GetBinError(ibin+1)/m)

            if m>0.0:
                h = min(mcHistSumNominal.GetBinError(ibin+1)/m,0.899)
                box = ROOT.TBox(c-0.5*w,1-h,c+0.5*w,1+h)
                box.SetFillStyle(3345)
                box.SetLineColor(ROOT.kGray+1)
                box.SetFillColor(ROOT.kGray)
                box.SetLineWidth(2)
                rootObj.append(box)
                box.Draw("SameF")
                box2 = ROOT.TBox(c-0.5*w,1-h,c+0.5*w,1+h)
                box2.SetFillStyle(0)
                box2.SetLineColor(ROOT.kGray+1)
                box2.SetFillColor(ROOT.kGray)
                box2.SetLineWidth(2)
                rootObj.append(box2)
                box2.Draw("SameL")
                
        if self.showData:
            dataRes.Draw("HISTPESame")
                
        
        axisLine = ROOT.TF1("axisLine"+str(random.randint(0,99999)),"1",self.binRange[0],self.binRange[1])
        axisLine.SetLineColor(ROOT.kBlack)
        axisLine.SetLineWidth(1)
        axisLine.Draw("SameL")
        
        ROOT.gStyle.SetLineWidth(1)
        ROOT.gPad.RedrawAxis()


        hidePave=ROOT.TPaveText(cvxmin-0.06,resHeight-0.029,cvxmin-0.005,resHeight+0.028,"NDC")
        hidePave.SetFillColor(ROOT.kWhite)
        hidePave.SetFillStyle(1001)
        hidePave.Draw("Same")

                
        
        cv.Update()      
        cv.Print(self.plot+self.outputSuffix+".png")
        cv.Print(self.plot+self.outputSuffix+".pdf")
        
        pPreliminary=ROOT.TPaveText(cvxmin+0.025+0.09,cvymax-0.025,cvxmin+0.025+0.09,cvymax-0.025,"NDC")
        pPreliminary.SetFillColor(ROOT.kWhite)
        pPreliminary.SetBorderSize(0)
        pPreliminary.SetTextFont(53)
        pPreliminary.SetTextSize(27)
        pPreliminary.SetTextAlign(13)
        pPreliminary.AddText("Preliminary")
        pPreliminary.Draw("Same")
        
        cv.Update()      
        cv.Print(self.plot+self.outputSuffix+"_PAS.pdf")
        
        '''
        os.system("gs -sDEVICE=pdfwrite -dDEVICEWIDTHPOINTS=%i -dDEVICEHEIGHTPOINTS=%i -dPDFFitPage -o %s %s"%(
            7.0*1.35*28.3465,6.5*1.35*28.3465, #cm to points
            self.plot+self.outputSuffix+".pdf",
            self.plot+self.outputSuffix+".eps"
        ))
        '''
jSym = "j#lower[-0.2]{#scale[0.8]{*}}"
pqj = "P#lower[0.3]{#scale[0.7]{q}}#kern[-0.5]{ }("+jSym+")"
plj = "P#lower[0.3]{#scale[0.7]{#font[12]{l}}}#kern[-0.5]{ }("+jSym+")"
pj = "P#lower[0.3]{#scale[0.7]{q,#kern[-0.7]{ }#font[12]{l}}}#kern[-0.5]{ }("+jSym+")"

mll = "m#lower[0.3]{#scale[0.7]{#kern[-0.6]{ }#font[12]{ll}}}"
mllj = "m#lower[0.3]{#scale[0.7]{#kern[-0.6]{ }#font[12]{ll}"+jSym+"}}"

dR = "#Delta#kern[-0.25]{ }R(#font[12]{l}#lower[0.2]{#scale[0.8]{2}},#kern[-0.2]{ }"+jSym+")"

'''
bdtOS = Plot("bdt_SR","BDT score",combine=combineOS,binRange=[0,1],extraTitles=["SR OS"],yspace=1.7,outputSuffix="_OS")
bdtOS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e5, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{5}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
bdtOS()

bdtSS = Plot("bdt_SR","BDT score",combine=combineSS,binRange=[0,1],extraTitles=["SR SS"],yspace=1.8,outputSuffix="_SS")
bdtSS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e4, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{4}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
bdtSS()
'''
mlljOS = Plot("mllj_SR",mllj,combine=combineOS,binRange=[30,200], unit="GeV", rebin=5, extraTitles=["#font[62]{SR}, OS leptons"],yspace=1.8,outputSuffix="_OS")
mlljOS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e5, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{5}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
mlljOS()

mlljSS = Plot("mllj_SR",mllj,combine=combineSS,binRange=[30,200], unit="GeV", rebin=5, extraTitles=["#font[62]{SR}, SS leptons"],yspace=2.0,outputSuffix="_SS")
mlljSS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e4, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{4}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
mlljSS()

dROS = Plot("dR_SR",dR,combine=combineOS,binRange=[0,1.3], rebin=5, unit="", extraTitles=["#font[62]{SR}, OS leptons"],yRange=[1e3,1e9],outputSuffix="_OS",logy=True)
dROS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e5, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{5}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
dROS()

dRSS = Plot("dR_SR",dR,combine=combineSS,binRange=[0,1.3], rebin=5, unit="", extraTitles=["#font[62]{SR}, SS leptons"],yRange=[1e3,1e9],outputSuffix="_SS",logy=True)
dRSS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e4, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{4}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
dRSS()




tagger_SR_boosted_OS = Plot("tagger_SR_boosted",plj,combine=combineOS,binRange=[0,1],logy=True, rebin=10, extraTitles=["#font[62]{SR}, OS leptons, boosted "+jSym],yspace=1.8,outputSuffix="_OS")
tagger_SR_boosted_OS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e2, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{2}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
tagger_SR_boosted_OS()

tagger_SR_resolved_OS = Plot("tagger_SR_resolved",pqj,combine=combineOS,binRange=[0,1],logy=True, rebin=10, extraTitles=["#font[62]{SR}, OS leptons, resolved "+jSym],yspace=1.8,outputSuffix="_OS")
tagger_SR_resolved_OS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e2, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{2}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
tagger_SR_resolved_OS()

tagger_SR_boosted_SS = Plot("tagger_SR_boosted",plj,combine=combineSS,binRange=[0,1],logy=True, rebin=10, extraTitles=["#font[62]{SR}, SS leptons, boosted "+jSym],yspace=1.8,outputSuffix="_SS")
tagger_SR_boosted_SS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e2, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{2}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
tagger_SR_boosted_SS()

tagger_SR_resolved_SS = Plot("tagger_SR_resolved",pqj,combine=combineSS,binRange=[0,1],logy=True, rebin=10, extraTitles=["#font[62]{SR}, SS leptons, resolved "+jSym],yspace=1.8,outputSuffix="_SS")
tagger_SR_resolved_SS.addSignal("HNL_majorana_pt20_ctau1p0e00_massHNL10p0_Vall1p177e-03_all", 1e2, ["Majorana HNL (#times10#lower[-0.7]{#scale[0.7]{2}})","m#lower[0.3]{#scale[0.7]{N}}#kern[-0.2]{ }=#kern[-0.25]{ }10#kern[-0.1]{ }GeV, c#tau#lower[0.3]{#scale[0.7]{0}}#kern[-0.2]{ }=#kern[-0.25]{ }1#kern[-0.1]{ }mm"], style=[3,2,ROOT.kRed+1])
tagger_SR_resolved_SS()

tagger_CR_boosted = Plot("tagger_CR_boosted",pj,binRange=[0,1],logy=True, rebin=10, extraTitles=["#font[62]{CR}, OS#kern[-0.2]{ }+#kern[-0.2]{ }SS leptons, boosted "+jSym],yspace=1.7,showData=True)
tagger_CR_boosted()

tagger_CR_resolved = Plot("tagger_CR_resolved",pj,binRange=[0,1],logy=True, rebin=10, extraTitles=["#font[62]{CR}, OS#kern[-0.2]{ }+#kern[-0.2]{ }SS leptons, resolved "+jSym],yspace=1.7,showData=True)
tagger_CR_resolved()

'''
mll_fwd = Plot("mll_fwd",pj,combine=combineOS,binRange=[10,40],logy=False,yspace=0.2)
mll_fwd()
'''



