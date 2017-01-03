from os import sys
from ROOT import *
from math import sqrt
from Jawa import GetLumi

class Bunch:
        def __init__(self, **kwds):
                self.__dict__.update(kwds)

def normalise_to_min(hist):
    min_width = -1.0
    for i in range(hist.GetXaxis().GetNbins()):
        if hist.GetBinWidth(i+1) < min_width or min_width == -1.0:
            min_width = hist.GetBinWidth(i+1)
    for i in range(hist.GetXaxis().GetNbins()):
        if hist.GetBinWidth(i+1) != min_width:
            s = min_width/hist.GetBinWidth(i+1)
            hist.SetBinContent(i+1, hist.GetBinContent(i+1)*s)


def combine_graph_errs(name, pdf, scale, alphas = None):
    graph = pdf.Clone(name)
    for i in range(graph.GetN()):
        if alphas is not  None:
            errhi = sqrt(pow(pdf.GetErrorYhigh(i),2) + pow(alphas.GetErrorYhigh(i),2)) + scale.GetErrorYhigh(i)
            errlo = sqrt(pow(pdf.GetErrorYlow(i),2) + pow(alphas.GetErrorYlow(i),2))   + scale.GetErrorYlow(i)
        else:
            errhi = pdf.GetErrorYhigh(i) + scale.GetErrorYhigh(i)
            errlo = pdf.GetErrorYlow(i)  + scale.GetErrorYlow(i)
            
        graph.SetPointEYhigh(i, errhi)
        graph.SetPointEYlow(i, errlo)
    return graph

def print_diff(f, g):
    for k in f.GetListOfKeys():
        ha = f.Get(k.GetName())
        hb = g.Get(k.GetName())
        if isinstance(ha, TH1F) and isinstance(hb, TH1F):
            print "---------------------"
            print k.GetName()
            print "a", "b", "(a-b)/a", "ae", "be", "precision"
            for i in range(ha.GetNbinsX()):
                a = ha.GetBinContent(i+1)
                b = hb.GetBinContent(i+1)
                ae = ha.GetBinError(i+1)
                be = hb.GetBinError(i+1)
                print a, b, (a-b)/a, ae, be, sqrt(pow((1/a)*be,2) + pow(b/(a*a)*ae,2))
            for i in range(ha.GetNbinsX()):
                a = ha.GetBinContent(i+1)
                b = hb.GetBinContent(i+1)
                ae = ha.GetBinError(i+1)
                be = hb.GetBinError(i+1)
                print (a-b)/a
        

def GetAcceptance(name, numerator, denominator, common):
    asymm = numerator.Clone(name)
    asymm.Divide(denominator)
    #calculate error correctly to take into account common events
    # i.e. acceptance = (a+b)/(a+c)
    #df/da = (c-b)/(a+c)^2
    #df/db = (a+c)/(a+c)^2
    #df/dc = -(a+b)/(a+c)^2
    for i in range(asymm.GetNbinsX()):
        a = common.GetBinContent(i+1)
        b = numerator.GetBinContent(i+1) - a
        c = denominator.GetBinContent(i+1) - a
        print a, b, c
        err = sqrt(pow(sqrt(a)*(c - b)/pow(a+c,2),2) + pow(sqrt(b)*(a+c)/pow(a+c,2),2)
                   + pow(sqrt(c) * (a+b)/pow(a+c,2),2))
        asymm.SetBinError(i+1, err)
    return asymm


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
    def __init__(self, name, c1, scales1, c2, scales2, pdfs1 = None, pdfs2 = None, alphas1 = None, alphas2 = None):
        self.name = name
        self.tot_asymm        = (c1.Integral() - c2.Integral()) / (c1.Integral() + c2.Integral())
        self.tot_asymm_scales = [(p.Integral() - m.Integral()) / (p.Integral() + m.Integral()) for p, m in zip(scales1, scales2)]
        self.tot_asymm_scale_errhi = max(self.tot_asymm_scales) - self.tot_asymm
        self.tot_asymm_scale_errlo = self.tot_asymm - min(self.tot_asymm_scales)
        self.tot_asymm_pdf_err = 0.0
        self.tot_asymm_alphas_errhi = 0.0
        self.tot_asymm_alphas_errlo = 0.0
        #differential plot with scale errors
        self.asymm = c1.GetAsymmetry(c2)
        self.asymm_stat = c1.GetAsymmetry(c2)
        self.asymm_scales = [p.GetAsymmetry(m) for p, m in zip(scales1, scales2)]
        self.asymm_scaleerrs        = GetSpreadGraph(name+"_scaleerrs", [self.asymm] + self.asymm_scales, addErr = True)  
        self.asymm_scaleerrs_nostat = GetSpreadGraph(name+"_scaleerrs_nostat", [self.asymm] + self.asymm_scales, addErr = False) 
        if pdfs1 != None and pdfs2 != None:
            self.tot_asymm_pdfs = [(p.Integral() - m.Integral())/(p.Integral() + m.Integral()) for p,m in zip(pdfs1, pdfs2)]
            self.tot_asymm_pdf_err = getpdfuncertainty(self.tot_asymm, self.tot_asymm_pdfs)
            self.asymm_pdfs = [p.GetAsymmetry(m) for p, m in zip(pdfs1, pdfs2)]
            self.asymm_pdferrs        = GetPDFGraph(name+"_pdferrs", [self.asymm] + self.asymm_pdfs, addErr = True)  
            self.asymm_pdferrs_nostat = GetPDFGraph(name+"_pdferrs_nostat", [self.asymm] + self.asymm_pdfs, addErr = False) 
            self.asymm_toterrs        = combine_graph_errs(name+"_toterrs",
                                                           scale = self.asymm_scaleerrs_nostat,
                                                           pdf = self.asymm_pdferrs_nostat)
        if alphas1 != None and alphas2 != None:
            self.tot_asymm_alphas = [(p.Integral() - m.Integral())/(p.Integral() + m.Integral()) for p,m in zip(alphas1, alphas2)]
            self.tot_asymm_alphas_errhi = max(self.tot_asymm_alphas) - self.tot_asymm
            self.tot_asymm_alphas_errlo = self.tot_asymm - min(self.tot_asymm_alphas)
            self.asymm_alphas = [p.GetAsymmetry(m) for p, m in zip(alphas1, alphas2)]
            self.asymm_alphaserrs        = GetSpreadGraph(name+"_alphaserrs", [self.asymm] + self.asymm_alphas, addErr = True)  
            self.asymm_alphaserrs_nostat = GetSpreadGraph(name+"_alphaserrs_nostat", [self.asymm] + self.asymm_alphas, addErr = False) 
            self.asymm_toterrs        = combine_graph_errs(name+"_toterrs",
                                                           scale = self.asymm_scaleerrs_nostat,
                                                           pdf = self.asymm_pdferrs_nostat,
                                                           alphas = self.asymm_alphaserrs_nostat)
        else:
            self.asymm_alphaserrs = None
            self.asymm_alphaserrs_nostat = None            

        self.tot_asymm_errhi = sqrt(pow(self.tot_asymm_pdf_err,2) + pow(self.tot_asymm_alphas_errhi,2)) + self.tot_asymm_scale_errhi
        self.tot_asymm_errlo = sqrt(pow(self.tot_asymm_pdf_err,2) + pow(self.tot_asymm_alphas_errlo,2)) + self.tot_asymm_scale_errlo

    def write(self):
        TParameter(float)(self.name+"_tot", self.tot_asymm).Write()
        TParameter(float)(self.name+"_tot_scaleerr_hi", self.tot_asymm_scale_errhi).Write()
        TParameter(float)(self.name+"_tot_scaleerr_lo", self.tot_asymm_scale_errlo).Write()
        TParameter(float)(self.name+"_tot_alphaserr_hi", self.tot_asymm_alphas_errhi).Write()
        TParameter(float)(self.name+"_tot_alphaserr_lo", self.tot_asymm_alphas_errlo).Write()
        TParameter(float)(self.name+"_tot_pdferr", self.tot_asymm_pdf_err).Write()
        TParameter(float)(self.name+"_tot_errhi", self.tot_asymm_errhi).Write()
        TParameter(float)(self.name+"_tot_errlo", self.tot_asymm_errlo).Write()
        self.asymm_scaleerrs_nostat.Write(self.name+"_scaleerr")
        self.asymm_stat.Write(self.name+"_staterrs")
        if (self.asymm_pdferrs_nostat):
            self.asymm_pdferrs_nostat.Write(self.name+"_pdferrs")   
        if (self.asymm_alphaserrs_nostat):
            self.asymm_alphaserrs_nostat.Write(self.name+"_alphaserrs")    
        if (self.asymm_toterrs):
            self.asymm_toterrs.Write(self.name+"_toterrs")             
                

def get_ratio(name, c1, c2):
    r = c1.Clone(name+"_ratio")
    r.Sumw2()
    r.Divide(c2)
    return r

class ratio_details:
    def __init__(self, name, c1, scales1, c2, scales2, pdfs1 = None, pdfs2 = None, alphas1 = None, alphas2 = None):
        self.name = name
        self.tot_ratio        = c1.Integral()/ c2.Integral()
        self.tot_ratio_scales = [p.Integral() / m.Integral() for p, m in zip(scales1, scales2)]
        self.tot_ratio_scale_errhi = max(self.tot_ratio_scales) - self.tot_ratio
        self.tot_ratio_scale_errlo = self.tot_ratio - min(self.tot_ratio_scales)
        self.tot_ratio_alphas_errhi = 0.0
        self.tot_ratio_alphas_errlo = 0.0
        self.tot_ratio_pdf_err = 0.0
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
            self.ratio_toterrs        = combine_graph_errs(name+"_toterrs",
                                                           scale = self.ratio_scaleerrs_nostat,
                                                           pdf = self.ratio_pdferrs_nostat)
        else:
            self.ratio_pdferrs = None
            self.ratio_pdferrs_nostat = None
            self.tot_ratio_pdf_err = 0
            self.ratio_toterrs = None
        if alphas1 != None and alphas2 != None:
            self.tot_ratio_alphass = [p.Integral() / m.Integral() for p, m in zip(alphas1, alphas2)]
            self.tot_ratio_alphas_errhi = max(self.tot_ratio_alphass) - self.tot_ratio
            self.tot_ratio_alphas_errlo = self.tot_ratio - min(self.tot_ratio_alphass)
            #differential plot with alphas errors
            self.ratio_alphas = [get_ratio(name+'_alphas'+str(i), p, m) for p, m, i in zip(alphas1, alphas2, range(len(alphas1)))]
            self.ratio_alphaserrs        = GetSpreadGraph(name+"_alphaserrs", [self.ratio] + self.ratio_alphas, addErr = True)  
            self.ratio_alphaserrs_nostat = GetSpreadGraph(name+"_alphaserrs_nostat", [self.ratio] + self.ratio_alphas, addErr = False) 
        else:
            self.ratio_alphaserrs = None
            self.ratio_alphaserrs_nostat = None     
  
        self.ratio_toterrs        = combine_graph_errs(name+"_toterrs",
                                                       scale = self.ratio_scaleerrs_nostat,
                                                       pdf = self.ratio_pdferrs_nostat,
                                                       alphas = self.ratio_alphaserrs_nostat)
        self.tot_ratio_errhi = sqrt(pow(self.tot_ratio_pdf_err,2) + pow(self.tot_ratio_alphas_errhi,2)) + self.tot_ratio_scale_errhi
        self.tot_ratio_errlo = sqrt(pow(self.tot_ratio_pdf_err,2) + pow(self.tot_ratio_alphas_errlo,2)) + self.tot_ratio_scale_errlo

    def write(self):
        TParameter(float)(self.name+"_tot", self.tot_ratio).Write()
        TParameter(float)(self.name+"_tot_scaleerr_hi", self.tot_ratio_scale_errhi).Write()
        TParameter(float)(self.name+"_tot_scaleerr_lo", self.tot_ratio_scale_errlo).Write()
        TParameter(float)(self.name+"_tot_alphaserr_hi", self.tot_ratio_alphas_errhi).Write()
        TParameter(float)(self.name+"_tot_alphaserr_lo", self.tot_ratio_alphas_errlo).Write()
        TParameter(float)(self.name+"_tot_pdferr", self.tot_ratio_pdf_err).Write()
        TParameter(float)(self.name+"_tot_errhi", self.tot_ratio_errhi).Write()
        TParameter(float)(self.name+"_tot_errlo", self.tot_ratio_errlo).Write()
        self.ratio_scaleerrs_nostat.Write(self.name+"_scaleerr")
        if (self.ratio_pdferrs_nostat):
            self.ratio_pdferrs_nostat.Write(self.name+"_pdferr")
        if (self.ratio_alphaserrs_nostat):
            self.ratio_alphaserrs_nostat.Write(self.name+"_alphaserr")
        if (self.ratio_toterrs):
            self.ratio_toterrs.Write(self.name+"_toterr")

class details:
    def __init__(self, name, c, scales, pdfs = None, alphas = None):
        self.name = name
        self.tot_xsec  = c.Integral()
        self.tot_xsec_scales = [m.Integral() for m in scales]
        self.tot_xsec_scale_errhi = max(0, max(self.tot_xsec_scales) - self.tot_xsec)
        self.tot_xsec_scale_errlo = max(0, self.tot_xsec - min(self.tot_xsec_scales))
        #declare normalised cross-sections
        self.norm_xsec = c.Clone(name+"_diff_norm")
        #declare norm cross-section for scale sets
        self.norm_xsec_scales = [h.Clone(h.GetName()+"_norm") for h in scales]
        #normalise
        for h in ([self.norm_xsec] + self.norm_xsec_scales):
		if h.Integral() > 0:
			h.Scale(1/h.Integral())
        #differential plot with scale errorss
	self.xsec_staterrs         = c.Clone(name+"_staterrs")
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
		    if h.Integral() > 0:
			    h.Scale(1/h.Integral())
            self.xsec_pdferrs        = GetPDFGraph(name+"_pdferrs", [c] +  pdfs, addErr = True)  
            self.xsec_pdferrs_nostat = GetPDFGraph(name+"_pdferrs_nostat", [c] + pdfs, addErr = False) 
            self.xsec_toterrs        = combine_graph_errs(name+"_toterrs",
                                                          scale = self.xsec_scaleerrs_nostat,
                                                          pdf = self.xsec_pdferrs_nostat)
        else:
            self.tot_xsec_pdf_err = 0
            self.xsec_pdferrs = None
            self.xsec_pdferrs_nostat = None
            self.xsec_toterrs = None
        if alphas != None:
            self.tot_xsec_alphas = [m.Integral() for m in alphas]
            self.tot_xsec_alphas_errhi   = max(self.tot_xsec_alphas) - self.tot_xsec
            self.tot_xsec_alphas_errlo   = self.tot_xsec - min(self.tot_xsec_alphas)
            self.xsec_alphaserrs        = GetSpreadGraph(name+"_alphaserrs", [c] +  alphas, addErr = True)  
            self.xsec_alphaserrs_nostat = GetSpreadGraph(name+"_alphaserrs_nostat", [c] + alphas, addErr = False) 
            #declare norm cross-section for different alphas values
            self.norm_xsec_alphas = [h.Clone(h.GetName()+"_norm") for h in alphas]
            #normalise
            for h in (self.norm_xsec_alphas):
		    if h.Integral() > 0:
			    h.Scale(1/h.Integral())
            self.xsec_toterrs        = combine_graph_errs(name+"_toterrs",
                                                          scale = self.xsec_scaleerrs_nostat,
                                                          pdf = self.xsec_pdferrs_nostat,
                                                          alphas = self.xsec_alphaserrs_nostat)
        else:
            self.tot_xsec_alphas_errhi = 0
            self.tot_xsec_alphas_errlo = 0
            self.xsec_alphaserrs = None
            self.xsec_alphaserrs_nostat = None
        self.tot_xsec_errhi = sqrt(pow(self.tot_xsec_pdf_err,2) + pow(self.tot_xsec_alphas_errhi,2)) + self.tot_xsec_scale_errhi
        self.tot_xsec_errlo = sqrt(pow(self.tot_xsec_pdf_err,2) + pow(self.tot_xsec_alphas_errlo,2)) + self.tot_xsec_scale_errlo

    def write(self):
        TParameter(float)(self.name+"_tot", self.tot_xsec).Write()
        TParameter(float)(self.name+"_scaleerr_hi", self.tot_xsec_scale_errhi).Write()
        TParameter(float)(self.name+"_scaleerr_lo", self.tot_xsec_scale_errlo).Write()
        TParameter(float)(self.name+"_pdferr"     , self.tot_xsec_pdf_err).Write()
        TParameter(float)(self.name+"_alphaserr_hi", self.tot_xsec_alphas_errhi).Write()
        TParameter(float)(self.name+"_alphaserr_lo", self.tot_xsec_alphas_errlo).Write()
        TParameter(float)(self.name+"_toterr_hi", self.tot_xsec_errhi).Write()
        TParameter(float)(self.name+"_toterr_lo", self.tot_xsec_errlo).Write()
        self.xsec_staterrs.Write(self.name+"_staterrs")
        self.xsec_scaleerrs_nostat.Write(self.name+"_scaleerrs")
        #self.norm_xsec_scaleerrs_nostat.Write(self.name+"_norm_scaleerrs")
        if self.xsec_pdferrs:
            self.xsec_pdferrs_nostat.Write(self.name+"_pdferrs")
        if self.xsec_alphaserrs:
            self.xsec_alphaserrs.Write(self.name+"_alphaserrs")
        if self.xsec_toterrs:
            self.xsec_toterrs.Write(self.name+"_toterrs")


class wzdetails:
    def __init__(self, name, wp, wm, z, wp_scales, wm_scales, z_scales, wp_pdfs = None, wm_pdfs = None, z_pdfs = None, wp_alphas = None, wm_alphas = None, z_alphas=None):
        self.name = name
        self.wp = wp.Integral()
        self.wm = wm.Integral()
        self.z  = z.Integral()
        self.wp_scales = [p.Integral() for p in wp_scales]
        self.wm_scales = [m.Integral() for m in wm_scales ]
        self.z_scales  = [zz.Integral() for zz in z_scales]
        print (self.wp + self.wm)/self.z
        if wp_pdfs and wm_pdfs and z_pdfs:
            self.wp_pdfs   = [p.Integral() for p in wp_pdfs  if isinstance(p, TH1) ]
            self.wm_pdfs   = [m.Integral() for m in wm_pdfs  if isinstance(m, TH1) ]
            self.z_pdfs    = [zz.Integral() for zz in z_pdfs   if isinstance(zz, TH1)]
        
        if wp_alphas and wm_alphas and z_alphas:
            self.wp_alphas = [p.Integral() for p in wp_alphas]
            self.wm_alphas = [m.Integral() for m in wm_alphas]
            self.z_alphas  = [zz.Integral() for zz in z_alphas ]

        self.rwpm  = wp.Integral()/wm.Integral()
        self.rwz   = (wp.Integral() + wm.Integral())/z.Integral()
        print name, self.rwz, (wp.Integral() + wm.Integral())/z.Integral(), (self.wp + self.wm)/self.z
        self.rwpz  = wp.Integral()/z.Integral()
        self.rwmz  = wm.Integral()/z.Integral()
        self.rwpm_scales = [p.Integral()/m.Integral() for p,m in zip(wp_scales, wm_scales)]
        self.rwz_scales  = [(p.Integral() + m.Integral())/t.Integral() for p,m,t in zip(wp_scales, wm_scales, z_scales)]
        self.rwpz_scales = [p.Integral()/t.Integral() for p,m,t in zip(wp_scales, wm_scales, z_scales)]
        self.rwmz_scales = [m.Integral()/t.Integral() for p,m,t in zip(wp_scales, wm_scales, z_scales)]

        if wp_pdfs and wm_pdfs and z_pdfs:
            self.rwpm_pdfs   = [p.Integral()/m.Integral() for p,m in zip(wp_pdfs, wm_pdfs) if isinstance(p, TH1)]
            self.rwz_pdfs    = [(p.Integral() + m.Integral())/t.Integral() for p,m,t in zip(wp_pdfs, wm_pdfs, z_pdfs) if isinstance(p, TH1)]
            self.rwpz_pdfs   = [p.Integral()/t.Integral() for p,m,t in zip(wp_pdfs, wm_pdfs, z_pdfs) if isinstance(p, TH1)]
            self.rwmz_pdfs   = [m.Integral()/t.Integral() for p,m,t in zip(wp_pdfs, wm_pdfs, z_pdfs) if isinstance(p, TH1)]

        if wp_alphas and wm_alphas and z_alphas:
            self.rwpm_alphas   = [p.Integral()/m.Integral() for p,m in zip(wp_alphas, wm_alphas) if isinstance(p, TH1)]
            self.rwz_alphas    = [(p.Integral() + m.Integral())/t.Integral() for p,m,t in zip(wp_alphas, wm_alphas, z_alphas) if isinstance(p, TH1)]
            self.rwpz_alphas   = [p.Integral()/t.Integral() for p,m,t in zip(wp_alphas, wm_alphas, z_alphas) if isinstance(p, TH1)]
            self.rwmz_alphas   = [m.Integral()/t.Integral() for p,m,t in zip(wp_alphas, wm_alphas, z_alphas) if isinstance(p, TH1)]
        
        self.rwpm_scale_errhi = max(0, max(self.rwpm_scales) - self.rwpm)
        self.rwpm_scale_errlo = max(0, self.rwpm - min(self.rwpm_scales))

        if wp_pdfs and wm_pdfs and z_pdfs:
            self.rwpm_pdf_err = getpdfuncertainty(self.rwpm, self.rwpm_pdfs)
        else:
            self.rwpm_pdf_err = 0.0
        
        if wp_alphas and wm_alphas and z_alphas:
            self.rwpm_alphas_errhi = max(0, max(self.rwpm_alphas) - self.rwpm)
            self.rwpm_alphas_errlo = max(0, self.rwpm - min(self.rwpm_alphas))
        else:
            self.rwpm_alphas_errhi = 0.0
            self.rwpm_alphas_errlo = 0.0
        
        self.rwpm_tot_errhi = self.rwpm_scale_errhi + (sqrt(pow(self.rwpm_pdf_err,2) + pow(self.rwpm_alphas_errhi,2)))
        self.rwpm_tot_errlo = self.rwpm_scale_errlo + (sqrt(pow(self.rwpm_pdf_err,2) + pow(self.rwpm_alphas_errlo,2)))
        
        self.rwz_scale_errhi = max(0, max(self.rwz_scales) - self.rwz)
        self.rwz_scale_errlo = max(0, self.rwz - min(self.rwz_scales))
        if wp_pdfs and wm_pdfs and z_pdfs:
            self.rwz_pdf_err = getpdfuncertainty(self.rwz, self.rwz_pdfs)
        else:
            self.rwz_pdf_err = 0.0

        if wp_alphas and wm_alphas and z_alphas:
            self.rwz_alphas_errhi = max(0, max(self.rwz_alphas) - self.rwz)
            self.rwz_alphas_errlo = max(0, self.rwz - min(self.rwz_alphas))
        else:
            self.rwz_alphas_errhi = 0.0
            self.rwz_alphas_errlo = 0.0
        
        self.rwz_tot_errhi = self.rwz_scale_errhi + (sqrt(pow(self.rwz_pdf_err,2) + pow(self.rwz_alphas_errhi,2)))
        self.rwz_tot_errlo = self.rwz_scale_errlo + (sqrt(pow(self.rwz_pdf_err,2) + pow(self.rwz_alphas_errlo,2)))
        
        self.rwpz_scale_errhi = max(0, max(self.rwpz_scales) - self.rwpz)
        self.rwpz_scale_errlo = max(0, self.rwpz - min(self.rwpz_scales))
        if wp_pdfs and wm_pdfs and z_pdfs:
            self.rwpz_pdf_err = getpdfuncertainty(self.rwpz, self.rwpz_pdfs)
        else:
            self.rwpz_pdf_err = 0.0
        if wp_alphas and wm_alphas and z_alphas:
            self.rwpz_alphas_errhi = max(0, max(self.rwpz_alphas) - self.rwpz)
            self.rwpz_alphas_errlo = max(0, self.rwpz - min(self.rwpz_alphas))
        else:
            self.rwpz_alphas_errhi = 0.0
            self.rwpz_alphas_errlo = 0.0


        self.rwpz_tot_errhi = self.rwpz_scale_errhi + (sqrt(pow(self.rwpz_pdf_err,2) + pow(self.rwpz_alphas_errhi,2)))
        self.rwpz_tot_errlo = self.rwpz_scale_errlo + (sqrt(pow(self.rwpz_pdf_err,2) + pow(self.rwpz_alphas_errlo,2)))
        
        self.rwmz_scale_errhi = max(0, max(self.rwmz_scales) - self.rwmz)
        self.rwmz_scale_errlo = max(0, self.rwmz - min(self.rwmz_scales))
        if wp_pdfs and wm_pdfs and z_pdfs:
            self.rwmz_pdf_err = getpdfuncertainty(self.rwmz, self.rwmz_pdfs)
        else:
            self.rwmz_pdf_err = 0.0
        if wp_alphas and wm_alphas and z_alphas:
            self.rwmz_alphas_errhi = max(0, max(self.rwmz_alphas) - self.rwmz)
            self.rwmz_alphas_errlo = max(0, self.rwmz - min(self.rwmz_alphas))
        else:
            self.rwmz_alphas_errhi = 0.0
            self.rwmz_alphas_errlo = 0.0
            
        self.rwmz_tot_errhi = self.rwmz_scale_errhi + (sqrt(pow(self.rwmz_pdf_err,2) + pow(self.rwmz_alphas_errhi,2)))
        self.rwmz_tot_errlo = self.rwmz_scale_errlo + (sqrt(pow(self.rwmz_pdf_err,2) + pow(self.rwmz_alphas_errlo,2)))
        print name, self.rwz, (wp.Integral() + wm.Integral())/z.Integral()

    def write(self):
        TParameter(float)(self.name+"_rwpm", self.rwpm).Write()
        TParameter(float)(self.name+"_rwpm_scale_errhi" , self.rwpm_scale_errhi).Write()
        TParameter(float)(self.name+"_rwpm_scale_errlo" , self.rwpm_scale_errlo).Write()
        TParameter(float)(self.name+"_rwpm_alphas_errhi", self.rwpm_alphas_errhi).Write()
        TParameter(float)(self.name+"_rwpm_alphas_errlo", self.rwpm_alphas_errlo).Write()
        TParameter(float)(self.name+"_rwpm_pdferr"      , self.rwpm_pdf_err).Write()
        TParameter(float)(self.name+"_rwpm_tot_errhi"   , self.rwpm_tot_errhi).Write()
        TParameter(float)(self.name+"_rwpm_tot_errlo"   , self.rwpm_tot_errlo).Write()
        TParameter(float)(self.name+"_rwz"              , self.rwz).Write()
        TParameter(float)(self.name+"_rwz_scale_errhi"  , self.rwz_scale_errhi).Write()
        TParameter(float)(self.name+"_rwz_scale_errlo"  , self.rwz_scale_errlo).Write()
        TParameter(float)(self.name+"_rwz_alphas_errhi" , self.rwz_alphas_errhi).Write()
        TParameter(float)(self.name+"_rwz_alphas_errlo" , self.rwz_alphas_errlo).Write()
        TParameter(float)(self.name+"_rwz_tot_errhi"    , self.rwz_tot_errhi).Write()
        TParameter(float)(self.name+"_rwz_tot_errlo"    , self.rwz_tot_errlo).Write()
        TParameter(float)(self.name+"_rwz_pdferr"       , self.rwz_pdf_err).Write()
        TParameter(float)(self.name+"_rwpz"             , self.rwpz).Write()
        TParameter(float)(self.name+"_rwpz_scale_errhi" , self.rwpz_scale_errhi).Write()
        TParameter(float)(self.name+"_rwpz_scale_errlo" , self.rwpz_scale_errlo).Write()
        TParameter(float)(self.name+"_rwpz_alphas_errhi", self.rwpz_alphas_errhi).Write()
        TParameter(float)(self.name+"_rwpz_alphas_errlo", self.rwpz_alphas_errlo).Write()
        TParameter(float)(self.name+"_rwpz_pdferr"      , self.rwpz_pdf_err).Write()
        TParameter(float)(self.name+"_rwpz_tot_errhi"   , self.rwpz_tot_errhi).Write()
        TParameter(float)(self.name+"_rwpz_tot_errlo"   , self.rwpz_tot_errlo).Write()
        TParameter(float)(self.name+"_rwmz"             , self.rwmz).Write()
        TParameter(float)(self.name+"_rwmz_scale_errhi" , self.rwmz_scale_errhi).Write()
        TParameter(float)(self.name+"_rwmz_scale_errlo" , self.rwmz_scale_errlo).Write()
        TParameter(float)(self.name+"_rwmz_alphas_errhi", self.rwmz_alphas_errhi).Write()
        TParameter(float)(self.name+"_rwmz_alphas_errlo", self.rwmz_alphas_errlo).Write()
        TParameter(float)(self.name+"_rwmz_pdferr"      , self.rwmz_pdf_err).Write()
        TParameter(float)(self.name+"_rwmz_tot_errhi"   , self.rwmz_tot_errhi).Write()
        TParameter(float)(self.name+"_rwmz_tot_errlo"   , self.rwmz_tot_errlo).Write()
