from os import sys
from ROOT import TFile, TTree, TH1, TCut, kRed, kYellow, kBlue, kOrange, kGreen, TF1, kMagenta, TMath, TH1F, TGraphAsymmErrors, Double, TCanvas
from math import sqrt


def getpdfuncertainty(central, pdfsets):
    sumsq = 0
    for pdf in pdfsets:
        sumsq = sumsq + pow(central - pdf,2)
    return sqrt(sumsq / (len(pdfsets) - 1))


def normalise_to_bins(name, histA, histB, bins = [0, -1]):
    histAA = histA.Clone(name+"A")
    histBB = histB.Clone(name+"B")
    histAA.Scale(1 / histAA.Integral(bins[0], bins[1]))
    histBB.Scale(1 / histBB.Integral(bins[0], bins[1]))
    c = TCanvas(name)
    histAA.SetLineColor(kBlue)
    histAA.Draw()
    histBB.SetLineColor(kRed)
    histBB.Draw("same")
    return [histAA, histBB, c]



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
            if b == 0 :
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
            

#classes for dealing with scale and pdf uncertainties

class asymm_details:
    def __init__(self, name, c1, scales1, c2, scales2, pdfs1 = None, pdfs2 = None):
        self.name = name
        self.tot_asymm        = (c1.Integral() - c2.Integral()) / (c1.Integral() + c2.Integral())
        self.tot_asymm_scales = [(p.Integral() - m.Integral()) / (p.Integral() + m.Integral()) for p, m in zip(scales1, scales2)]
        self.tot_asymm_scale_errhi = max(self.tot_asymm_scales) - self.tot_asymm
        self.tot_asymm_scale_errlo = self.tot_asymm - min(self.tot_asymm_scales)
        #differential plot with scale errors
        self.asymm = c1.GetAsymmetry(c2)
        self.asymm_scales = [p.GetAsymmetry(m) for p, m in zip(scales1, scales2)]
        self.asymm_scaleerrs        = GetSpreadGraph(name+"_scaleerrs", [self.asymm] + self.asymm_scales, addErr = True)  
        self.asymm_scaleerrs_nostat = GetSpreadGraph(name+"_scaleerrs_nostat", [self.asymm] + self.asymm_scales, addErr = False) 
        if pdfs1 != None and pdfs2 != None:
            self.tot_asymm_pdfs = [(p.Integral() - m.Integral())/(p.Integral() + m.Integral()) for p,m in zip(pdfs1, pdfs2)]
            self.tot_asymm_pdf_err = getpdfuncertainty(self.tot_asymm, self.tot_asymm_pdfs)
            self.asymm_pdfs = [p.GetAsymmetry(m) for p, m in zip(pdfs1, pdfs2)]
            self.asymm_pdferrs        = GetPDFGraph(name+"_pdferrs", [self.asymm] + self.asymm_pdfs, addErr = True)  
            self.asymm_pdferrs_nostat = GetPDFGraph(name+"_pdferrs_nostat", [self.asymm] + self.asymm_pdfs, addErr = False) 
        else:
            self.asymm_pdferrs = None
            self.asymm_pdferrs_nostat = None
            self.tot_asymm_pdf_err = 0
            

def get_ratio(name, c1, c2):
    r = c1.Clone(name+"_ratio")
    r.Sumw2()
    r.Divide(c2)
    return r

class ratio_details:
    def __init__(self, name, c1, scales1, c2, scales2, pdfs1 = None, pdfs2 = None):
        self.name = name
        self.tot_ratio        = c1.Integral()/ c2.Integral()
        self.tot_ratio_scales = [p.Integral() / m.Integral() for p, m in zip(scales1, scales2)]
        self.tot_ratio_scale_errhi = max(self.tot_ratio_scales) - self.tot_ratio
        self.tot_ratio_scale_errlo = self.tot_ratio - min(self.tot_ratio_scales)
        #differential plot with scale errors
        self.ratio = get_ratio(name, c1, c2)
        self.ratio_scales = [get_ratio(name+'_scale'+str(i), p, m) for p, m, i in zip(scales1, scales2, range(len(scales1)))]
        self.ratio_scaleerrs        = GetSpreadGraph(name+"_scaleerrs", [self.ratio] + self.ratio_scales, addErr = True)  
        self.ratio_scaleerrs_nostat = GetSpreadGraph(name+"_scaleerrs_nostat", [self.ratio] + self.ratio_scales, addErr = False) 
        if pdfs1 != None and pdfs2 != None:
            self.tot_ratio_pdfs = [p.Integral()/m.Integral() for p,m in zip(pdfs1, pdfs2)]
            self.tot_ratio_pdf_err = getpdfuncertainty(self.tot_ratio, self.tot_ratio_pdfs)
            self.ratio_pdfs = [get_ratio(name+'_pdfs'+str(i), p, m) for p, m, i in zip(pdfs1, pdfs2, range(len(pdfs1)))]
            self.ratio_pdferrs        = GetPDFGraph(name+"_pdferrs", [self.ratio] + self.ratio_pdfs, addErr = True)  
            self.ratio_pdferrs_nostat = GetPDFGraph(name+"_pdferrs_nostat", [self.ratio] + self.ratio_pdfs, addErr = False) 
        else:
            self.ratio_pdferrs = None
            self.ratio_pdferrs_nostat = None
            self.tot_ratio_pdf_err = 0

class details:
    def __init__(self, name, c, scales, pdfs = None):
        self.name = name
        self.tot_xsec  = c.Integral()
        self.tot_xsec_scales = [m.Integral() for m in scales]
        self.tot_xsec_scale_errhi = max(self.tot_xsec_scales) - self.tot_xsec
        self.tot_xsec_scale_errlo = self.tot_xsec - min(self.tot_xsec_scales)
        #declare normalised cross-sections
        self.norm_xsec = c.Clone(name+"_diff_norm")
        #declare norm cross-section for scale sets
        self.norm_xsec_scales = [h.Clone(h.GetName()+"_norm") for h in scales]
        #normalise
        for h in ([self.norm_xsec] + self.norm_xsec_scales):
            h.Scale(1/h.Integral())
        #differential plot with scale errors
        self.xsec_scaleerrs        = GetSpreadGraph(name+"_scaleerrs", [c] + scales, addErr = True)  
        self.xsec_scaleerrs_nostat = GetSpreadGraph(name+"_scaleerrs_nostat", [c] + scales, addErr = False) 
        #normalised differential plot with scale errors
        self.norm_xsec_scaleerrs        = GetSpreadGraph(name+"_norm_scaleerrs", [self.norm_xsec] + self.norm_xsec_scales, addErr = True)  
        self.norm_xsec_scaleerrs_nostat = GetSpreadGraph(name+"_norm_scaleerrs_nostat", [self.norm_xsec] + self.norm_xsec_scales, addErr = False) 
        if pdfs != None:
            self.tot_xsec_pdfs = [m.Integral() for m in pdfs]
            self.tot_xsec_pdf_err   = getpdfuncertainty(self.tot_xsec, self.tot_xsec_pdfs)
            self.xsec_pdferrs        = GetPDFGraph(name+"_pdferrs", [c] +  pdfs, addErr = True)  
            self.xsec_pdferrs_nostat = GetPDFGraph(name+"_pdferrs_nostat", [c] + pdfs, addErr = False) 
            #declare norm cross-section for pdf sets
            self.norm_xsec_pdfs = [h.Clone(h.GetName()+"_norm") for h in pdfs]
            #normalise
            for h in (self.norm_xsec_pdfs):
                h.Scale(1/h.Integral())
            self.xsec_pdferrs        = GetPDFGraph(name+"_pdferrs", [c] +  pdfs, addErr = True)  
            self.xsec_pdferrs_nostat = GetPDFGraph(name+"_pdferrs_nostat", [c] + pdfs, addErr = False) 
        else:
            self.tot_xsec_pdf_err = 0
            self.xsec_pdferrs = None
            self.xsec_pdferrs_nostat = None


