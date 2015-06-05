from os import sys
from ROOT import TFile, TTree, TH1, TCut, kRed, kYellow, kBlue, kOrange, kGreen, TF1, kMagenta, TMath, TH1F, TGraphAsymmErrors, Double
from math import sqrt



def gethessuncertainty(central, pdfsets):
    sumsqup   = 0
    sumsqdown = 0
    for i in range(int(len(pdfsets)/2)):
        sumsqup   = sumsqup   + pow(max(pdfsets[i*2] - central, max(pdfsets[i*2+1] - central,0)),2)
        sumsqdown = sumsqdown + pow(max(central - pdfsets[i*2], max(central - pdfsets[i*2+1],0)),2)
    return [sqrt(sumsqup ), sqrt(sumsqdown)]

def get_rel_graph(name, g, c):
    x1 = Double(0.0)
    y1 = Double(0.0)
    x2 = Double(0.0)
    y2 = Double(0.0)

    g2 = g.Clone(name)
    if "TGraph" in g.ClassName():
        for i in range(g.GetN()):
            g.GetPoint(i, x1, y1)
            c.GetPoint(i, x2, y2)
            g2.SetPoint(i, x1, y1 - y2)
    else:
        for i in range(g.GetNbinsX()):
            c.GetPoint(i, x2, y2)
            g2.SetBinContent(i+1, g2.GetBinContent(i+1) - y2)
    return g2


def xsectoEvts(hist, scale, L = 0.001, evts = False, syst = None):
    #convert a cross-section distribution to number of events with stat uncertainty
    #note that if cross-section is in pb, 0.001 pb^{-1} corresponds to a 1 fb^{-1}
    #if evts is false just return same histogram with stat err, else return expected events
    outhist = hist.Clone(hist.GetTitle()+"_L_"+str(L)+"_S_"+str(scale))
    for i in range(hist.GetNbinsX()):
        nevts = hist.GetBinContent(i) * L * scale
        err = sqrt(nevts)
        if syst_hist:
            if (isinstance(syst,float)):
               systerr = syst
            else:
               systerr = (syst_hist.GetBinError(i)/syst_hist.GetBinContent(i)) * nevts
            err = sqrt(err*err + systerr*systerr)
        if evts:
            outhist.SetBinContent(i, nevts)
            outhist.SetBinError(i, err)
        else:
            outhist.SetBinError(outhist.GetBinContent(i)*(stat_err/nevts))
    return outhist


def GetSpreadHist(name, hists, addErr = False):
    if len(hists) == 0: return
    if len(hists) == 1: return hists[0]
    hist = hists[0].Clone(name)
    wdiffs = [0.0 for i in range(hist.GetNbinsX())]
    for i in range(hist.GetNbinsX()):
        for whist in hists:
            diff = abs(hist.GetBinContent(i+1) - whist.GetBinContent(i+1))
            if diff > wdiffs[i]:
                wdiffs[i] = diff
        if addErr: hist.SetBinError(i+1, sqrt(pow(hist.GetBinError(i+1),2) + pow(wdiffs[i],2)))  
        else: hist.SetBinError(i+1, wdiffs[i])
    return hist

def GetSpreadGraph(name, hists, addErr = False):
    if len(hists) == 0: return
    if len(hists) == 1: return hists[0]
    graph = TGraphAsymmErrors()
    hist = hists[0]
    wdiffsup   = [0.0 for i in range(hist.GetNbinsX())]
    wdiffsdown = [0.0 for i in range(hist.GetNbinsX())]
    for i in range(hist.GetNbinsX()):
        for whist in hists:
            diff = hist.GetBinContent(i+1) - whist.GetBinContent(i+1)
            if diff > 0 and abs(diff) > wdiffsdown[i]:
                wdiffsdown[i] = abs(diff)
            elif diff < 0 and abs(diff) > wdiffsup[i]:
                wdiffsup[i] = abs(diff)
        graph.SetPoint(i, hist.GetBinCenter(i+1), hist.GetBinContent(i+1))
        if addErr:
            graph.SetPointError(i,
                                hist.GetBinCenter(i+1) - hist.GetBinLowEdge(i+1), 
                                hist.GetXaxis().GetBinUpEdge(i+1) - hist.GetBinCenter(i+1),
                                sqrt(pow(hist.GetBinError(i+1),2) + pow(wdiffsdown[i],2)),
                                sqrt(pow(hist.GetBinError(i+1),2) + pow(wdiffsup[i],2))
                                )
        else: 
            graph.SetPointError(i,
                                hist.GetBinCenter(i+1) - hist.GetBinLowEdge(i+1), 
                                hist.GetXaxis().GetBinUpEdge(i+1) - hist.GetBinCenter(i+1),
                                wdiffsdown[i],
                                wdiffsup[i]
                                )
            
    return graph

def GetPDFHist(name, hists, addErr = False):
    if len(hists) == 0: return
    if len(hists) == 1: return hists[0]
    h = hists[0].Clone(hists[0].GetName()+"_pdferr")
    hist = hists[0]
    #wdiffsup   = [0.0 for i in range(hist.GetNbinsX())]
    #wdiffsdown = [0.0 for i in range(hist.GetNbinsX())]
    for i in range(hist.GetNbinsX()):
        sqsum = 0
        for whist in hists:
            sqsum += pow(hist.GetBinContent(i+1) - whist.GetBinContent(i+1),2)
        wdiff   = sqrt(sqsum / (len(hists) - 2))
        #h.SetBinContent(i+1, hist.GetBinContent(i+1))
        if addErr:
            h.SetBinError(i+1, sqrt(pow(hist.GetBinError(i+1),2) + pow(wdiff,2)))
        else: 
            h.SetBinError(i+1, wdiff)
    return h

def GetPDFGraph(name, hists, addErr = False):
    if len(hists) == 0: return
    if len(hists) == 1: return hists[0]
    graph = TGraphAsymmErrors()
    hist = hists[0]
    wdiffsup   = [0.0 for i in range(hist.GetNbinsX())]
    wdiffsdown = [0.0 for i in range(hist.GetNbinsX())]
    for i in range(hist.GetNbinsX()):
        sqsum = 0
        for whist in hists:
            sqsum += pow(hist.GetBinContent(i+1) - whist.GetBinContent(i+1),2)
        wdiffsup[i]   = sqrt(sqsum / (len(hists) - 2))
        wdiffsdown[i] = sqrt(sqsum / (len(hists) - 2))
        graph.SetPoint(i, hist.GetBinCenter(i+1), hist.GetBinContent(i+1))
        if addErr:
            graph.SetPointError(i,
                                hist.GetBinCenter(i+1) - hist.GetBinLowEdge(i+1), 
                                hist.GetXaxis().GetBinUpEdge(i+1) - hist.GetBinCenter(i+1),
                                sqrt(pow(hist.GetBinError(i+1),2) + pow(wdiffsdown[i],2)),
                                sqrt(pow(hist.GetBinError(i+1),2) + pow(wdiffsup[i],2))
                                )
        else: 
            graph.SetPointError(i,
                                hist.GetBinCenter(i+1) - hist.GetBinLowEdge(i+1),
                                hist.GetXaxis().GetBinUpEdge(i+1) - hist.GetBinCenter(i+1),
                                wdiffsdown[i],
                                wdiffsup[i]
                                )
    return graph



def gethessuncertainty(central, pdfsets):
    sumsqup   = 0
    sumsqdown = 0
    for i in range(int(len(pdfsets)/2)):
        sumsqup   = sumsqup   + pow(max(pdfsets[i*2] - central, max(pdfsets[i*2+1] - central,0)),2)
        sumsqdown = sumsqdown + pow(max(central - pdfsets[i*2], max(central - pdfsets[i*2+1],0)),2)
    return [sqrt(sumsqup ), sqrt(sumsqdown)]

def GetPDFHessGraph(name, hists, addErr = False, factor = 1.0):
    if len(hists) == 0: return
    if len(hists) == 1: return hists[0]
    graph = TGraphAsymmErrors()
    hist = hists[0]
    #wdiffsup   = [0.0 for i in range(hist.GetNbinsX())]
    #wdiffsdown = [0.0 for i in range(hist.GetNbinsX())]

    
    for i in range(hist.GetNbinsX()):
        #sqsumup   = 0
        #sqsumdown = 0
        #for j in range(len(hists)/2):
        #    sumsqup    = sqsumup   + pow(max(hists[j*2].GetBinContent(i+1) - hist.GetBinContent(i+1),
        #                                    max(hists[j*2+1].GetBinContent(i+1) - hist.GetBinContent(i+1),0)),2)
        #    sumsqdown  = sqsumdown + pow(max(hist.GetBinContent(i+1) - hists[j*2].GetBinContent(i+1),
        #                                    max(hist.GetBinContent(i+1) - hists[j*2+1].GetBinContent(i+1),0)),2)
        #wdiffsup[i]   = sqrt(sumsqup)
        #wdiffsdown[i] = sqrt(sumsqdown)
        wdiffs = gethessuncertainty(hist.GetBinContent(i+1), [h.GetBinContent(i+1) for h in hists])
        wdiffup   = wdiffs[0]
        wdiffdown = wdiffs[1]

        graph.SetPoint(i, hist.GetBinCenter(i+1), hist.GetBinContent(i+1))
        if addErr:
            graph.SetPointError(i,
                                hist.GetBinCenter(i+1) - hist.GetBinLowEdge(i+1), 
                                hist.GetXaxis().GetBinUpEdge(i+1) - hist.GetBinCenter(i+1),
                                sqrt(pow(hist.GetBinError(i+1),2) + pow(wdiffdown*factor,2)),
                                sqrt(pow(hist.GetBinError(i+1),2) + pow(wdiffup*factor,2))
                                )
        else: 
            graph.SetPointError(i,
                                hist.GetBinCenter(i+1) - hist.GetBinLowEdge(i+1),
                                hist.GetXaxis().GetBinUpEdge(i+1) - hist.GetBinCenter(i+1),
                                wdiffdown*factor,
                                wdiffup*factor
                                )
    return graph

def HistConsistency(hist):
    for i in range(hist.GetNbinsX() + 2):
        if hist.GetBinError(i) > hist.GetBinContent(i):
            hist.SetBinError(i, hist.GetBinContent(i))

def GetWeightHist(name, histA, histB):
    hist = TH1F(name, name, histA.GetNbinsX(), histA.GetBinLowEdge(1), histA.GetBinLowEdge(histA.GetXaxis().GetLast() + 1))
    if (histA.GetNbinsX() == histB.GetNbinsX() and 
        histA.GetBinLowEdge(0) == histB.GetBinLowEdge(0) and
        histA.GetBinLowEdge(histA.GetXaxis().GetLast() +1) == histB.GetBinLowEdge(histB.GetXaxis().GetLast() + 1)):
        for i in range(histA.GetNbinsX()+1):
            if (histB.GetBinContent(i) == 0):
                hist.SetBinContent(i, 1)
            else:
                hist.SetBinContent(i, (histA.GetBinContent(i)/histB.GetBinContent(i))*(histB.Integral()/histA.Integral()))
    return hist

def GetWeightHist2D(name, histA, histB):
    hist_norm = histA.Clone(name)
    for i in range(hist_norm.GetXaxis().GetNbins()):
        for j in range(hist_norm.GetYaxis().GetNbins()):
            a = histA.GetBinContent(i+1, j+1)
            b = histB.GetBinContent(i+1, j+1)
            if a == 0 or b == 0 :
                hist_norm.SetBinContent(i+1, j+1, 1)
            else:
                hist_norm.SetBinContent(i+1, j+1, a/b)
    return hist_norm

def GetMCLumi(name, mcdict):
    mclist = mcdict[name]
    lumi   = mclist[2]/(mclist[0] * mclist[1])
    return lumi

def GetVPT(hist, step = 5):
    tot = hist.GetNbinsY() / step
    normhists = []
    for i in range(tot):
        binned = hist.ProjectionX(str(i), (step + 1) * (i), step * (i+1))
        binned.Scale(1 / binned.Integral())
        normhists += [binned]
    return normhists
            
