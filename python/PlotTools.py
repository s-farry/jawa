import ROOT
from random import gauss
from types import IntType
from collections import OrderedDict
import time
from math import sqrt

#ROOT.gStyle.SetTextFont(80)
#ROOT.gStyle.SetLabelFont(80,"XYZ")
#ROOT.gStyle.SetTitleFont(80,"XYZ")
#ROOT.gStyle.SetTitleFont(80)
#ROOT.gStyle.SetOptStat(0)
#from lhcbStyle import SetStyle()

def GetTGraphSum(graph):
    x = ROOT.Double(0.0)
    y = ROOT.Double(0.0)
    ysum = 0
    for i in graph.GetN():
        graph.GetPoint(i, x, y)
        ysum += y
    return ysum

def split_list(list_items, split_idx):
    split_list = []
    #round up
    splits = round(len(list_items)/split_idx + 0.5)
    for i in range(splits):
        min_idx = i * split_idx
        max_idx = (i + 1) * split_idx
        max_idx = min(max_idx, len(list_items) -1)
        split_list += list_items[min_idx:max_idx]
        ymin = ROOT.Double(0.0)
        ymax = ROOT.Double(0.0)
        
        graph.GetPoint(0, xmin, ymin)
        graph.GetPoint(graph.GetN() - 1, xmax, ymax)
        ymin = ROOT.Double(0.0)
        ymax = ROOT.Double(0.0)
        
        graph.GetPoint(0, xmin, ymin)
        graph.GetPoint(graph.GetN() - 1, xmax, ymax)

def TGraph2Hists(name, graph):
    if isinstance(graph, ROOT.TGraph):
        xmin = ROOT.Double(0.0)
        xmax = ROOT.Double(0.0)
        ymin = ROOT.Double(0.0)
        ymax = ROOT.Double(0.0)
        
        graph.GetPoint(0, xmin, ymin)
        graph.GetPoint(graph.GetN() - 1, xmax, ymax)
        ymin = ROOT.Double(0.0)
        ymax = ROOT.Double(0.0)
        
        graph.GetPoint(0, xmin, ymin)
        graph.GetPoint(graph.GetN() - 1, xmax, ymax)
      
        xmin = xmin - (graph.GetErrorXlow(2))
        xmax = xmax + (graph.GetErrorXhigh(2))
      
        hist_ctr = ROOT.TH1F(name+"_ctr",name+"_ctr",graph.GetN(), xmin, xmax)
        hist_lo  = ROOT.TH1F(name+"_lo" ,name+"_lo" ,graph.GetN(), xmin, xmax)
        hist_hi  = ROOT.TH1F(name+"_hi" ,name+"_hi" ,graph.GetN(), xmin, xmax)
        for bin in range(graph.GetN()):
            valx = ROOT.Double(0.0)
            valy = ROOT.Double(0.0)
            graph.GetPoint(bin, valx, valy)
            valyhi = valy + graph.GetErrorYhigh(bin)
            valylo = valy - graph.GetErrorYlow(bin)
            hist_ctr.SetBinContent(bin+1,valy)
            hist_lo.SetBinContent(bin+1,valylo)
            hist_hi.SetBinContent(bin+1,valyhi)
            hist_ctr.SetBinError(bin+1,0)
            hist_lo.SetBinError(bin+1,0)
            hist_hi.SetBinError(bin+1,0)
        hist_ctr.SetBinContent(0, hist_ctr.GetBinContent(1))
        hist_lo.SetBinContent(0, hist_lo.GetBinContent(1))
        hist_hi.SetBinContent(0, hist_hi.GetBinContent(1))

        return [hist_ctr, hist_lo, hist_hi]
            

def SetROOTMarkerStyle(hist, marker):
    markerdict = {'o' : 20,
                  '^' : 22,
                  'v' : 23,
                  'x' : 5}
    if marker in markerdict:
        hist.SetMarkerStyle(markerdict[marker])
    else:
        hist.SetMarkerStyle(marker)
    hist.SetMarkerSize(0.7)

def SetROOTFillColor(hist, mcolor):
    colordict = {'red'     : ROOT.kRed,
                 'blue'    : ROOT.kBlue,
                 'yellow'  : ROOT.kYellow,
                 'magenta' : ROOT.kMagenta,
                 'green'   : ROOT.kGreen,
                 'black'   : ROOT.kBlack,
                 'orange'  : ROOT.kOrange}
    
    if mcolor in colordict:
       hist.SetFillColor(colordict[mcolor])
    elif mcolor is 'hatch':
       hist.SetFillStyle(3004)
    elif mcolor is 'None' or mcolor is 'none' or mcolor is None:
        hist.SetFillStyle(0)
    else:
       hist.SetFillColor(mcolor)

def SetROOTMarkerColor(hist, mcolor):
   colordict = {'red'     : ROOT.kRed,
                 'blue'    : ROOT.kBlue,
                 'yellow'  : ROOT.kYellow,
                 'magenta' : ROOT.kMagenta,
                 'green'   : ROOT.kGreen,
                 'black'   : ROOT.kBlack,
                 'orange'  : ROOT.kOrange}
   if mcolor in colordict:
        hist.SetMarkerColor(colordict[mcolor])
   else:
        hist.SetMarkerColor(mcolor)

def SetROOTLineColor(hist, mcolor):
   colordict = {'red'     : ROOT.kRed,
                'blue'    : ROOT.kBlue,
                'yellow'  : ROOT.kYellow,
                'magenta' : ROOT.kMagenta,
                'green'   : ROOT.kGreen,
                'black'   : ROOT.kBlack,
                'orange'  : ROOT.kOrange}
   if mcolor in colordict:
      hist.SetLineColor(colordict[mcolor])
      hist.SetMarkerColor(colordict[mcolor])
   else:
      hist.SetLineColor(mcolor)
      hist.SetMarkerColor(mcolor)

def MeVtoGeV(graph):
    i = graph.GetN()
    x = ROOT.Double(0.0)
    y = ROOT.Double(0.0)
    for j in range(0,i):
        graph.GetPoint(j,x,y)
        xerrhi = graph.GetErrorXhigh(j)
        xerrlo = graph.GetErrorXlow(j)
        yerrhi = graph.GetErrorYhigh(j)
        yerrlo = graph.GetErrorYlow(j)
        graph.SetPoint(j,x/1000.,y)
        graph.SetPointError(j,xerrlo/1000,xerrhi/1000,yerrlo,yerrhi)

def makeComp(grapha, graphb):
   if isinstance(grapha, ROOT.TGraph):
      xmin = ROOT.Double(0.0)
      xmax = ROOT.Double(0.0)
      ymin = ROOT.Double(0.0)
      ymax = ROOT.Double(0.0)
      
      grapha.GetPoint(0, xmin, ymin)
      grapha.GetPoint(grapha.GetN() - 1, xmax, ymax)
      
      xmin = xmin - (grapha.GetErrorXlow(2))
      xmax = xmax + (grapha.GetErrorXhigh(2))
      
      hist = ROOT.TH1F("hist","hist",grapha.GetN(), xmin, xmax)
      for bin in range(0, grapha.GetN()):
         val1x = ROOT.Double(0.0)
         val1y = ROOT.Double(0.0)
         val2x = ROOT.Double(0.0)
         val2y = ROOT.Double(0.0)
         grapha.GetPoint(bin, val1x, val1y)
         graphb.GetPoint(bin, val2x, val2y)
         hist.SetBinContent(bin+1,(val2y/val1y))
         val1err = max(grapha.GetErrorYhigh(bin), grapha.GetErrorYlow(bin))
         val2err = max(graphb.GetErrorYhigh(bin), graphb.GetErrorYlow(bin))
         hist.SetBinError(bin+1, val2y/val1y * numpy.sqrt((val1err/val1y)**2 + (val2err/val2y)**2))
      return hist
   elif ("TH1" in grapha.ClassName()):
      hist = grapha.Clone("CompGraph")
      for bin in range(0, hist.GetNbinsX() + 2):
         if graphb.GetBinContent(bin) > 0:
            hist.SetBinContent(bin, grapha.GetBinContent(bin)/graphb.GetBinContent(bin))
            hist.SetBinError(bin, grapha.GetBinError(bin)/graphb.GetBinContent(bin))
         else:
            hist.SetBinContent(bin, 0)
            hist.SetBinError(bin, 0)
      return hist

def makeNormErr(graph, name=''):
    if (isinstance(graph, ROOT.TGraph)):
        if name == '': name = graph.GetName()
        xmin = ROOT.Double(0.0)
        xmax = ROOT.Double(0.0)
        ymin = ROOT.Double(0.0)
        ymax = ROOT.Double(0.0)
        
        graph.GetPoint(0, xmin, ymin)
        graph.GetPoint(graph.GetN() - 1, xmax, ymax)
        
        xmin = xmin - (graph.GetErrorXlow(2))
        xmax = xmax + (graph.GetErrorXhigh(2))
        
        hist = ROOT.TH1F(name + "_normerr","hist",graph.GetN(), xmin, xmax)
        for bin in range(graph.GetN()):
            val1x = ROOT.Double(0.0)
            val1y = ROOT.Double(0.0)
            graph.GetPoint(bin, val1x, val1y)
            hist.SetBinContent(bin+1,1.0)
            valerr = max(graph.GetErrorYhigh(bin), graph.GetErrorYlow(bin))/val1y
            hist.SetBinError(bin+1, valerr)
        return hist
    elif (isinstance(graph, ROOT.TH1)):
        hist = graph.Clone(graph.GetName()+"_normerr")
        for bin in range(0, hist.GetNbinsX() + 2):
            if graph.GetBinContent(bin) != 0:
                hist.SetBinContent(bin, 1.0)
                hist.SetBinError(bin, abs(graph.GetBinError(bin)/graph.GetBinContent(bin)))
            else:
                hist.SetBinContent(bin, 1.0)
                hist.SetBinError(bin, 0)
        return hist
    elif isinstance(graph, list) and len(list) == 3:
        hist = graph[0].Clone(graph[0].GetName() + "_normerr")
        for bin in range(hist.GetNbinsX() + 2):
            if graph[0].GetBinContent(bin) != 0:
                hist.SetBinContent(bin, 1.0)
                bin_ctr = graph[0].GetBinContent(bin)
                bin_hi = graph[1].GetBinContent(bin)
                bin_lo = graph[2].GetBinContent(bin)
                hist.SetBinError(bin, max(bin_hi - bin_ctr, bin_ctr - bin_lo))
            else:
                hist.SetBinContent(bin, 1.0)
                hist.SetBinError(bin, 0)
        return hist

def makeNormErr2Graph(name, graph):
    if (isinstance(graph, ROOT.TGraph)):
        if name == '': name = graph.GetName()
        normgraph = graph.Clone(name+"_norm")

        for i in range(graph.GetN()):
            val1x = ROOT.Double(0.0)
            val1y = ROOT.Double(0.0)
            graph.GetPoint(i, val1x, val1y)
            normgraph.SetPoint(i, val1x, 1.0)
            if isinstance(graph, ROOT.TGraphAsymmErrors):
                if val1y != 0 :
                    normgraph.SetPointEYhigh(i, graph.GetErrorYhigh(i)/val1y)
                    normgraph.SetPointEYlow (i, graph.GetErrorYlow(i)/val1y)
                else:
                    normgraph.SetPointEYhigh(i,0)
                    normgraph.SetPointEYlow(i,0)
        return normgraph


def makeResidual(grapha, graphb):
   if ("TH1" in grapha.ClassName()):
      hist = grapha.Clone("CompGraph")
      for bin in range(0, hist.GetNbinsX() + 2):
         if graphb.GetBinContent(bin) > 0:
            hist.SetBinContent(bin, (grapha.GetBinContent(bin) - graphb.GetBinContent(bin))/grapha.GetBinError(bin))
            hist.SetBinError(bin, 0)
         else:
            hist.SetBinContent(bin, 0)
            hist.SetBinError(bin, 0)
      return hist

class StyleObj:
    def __init__(self, color = 1, fillcolor=None, linecolor=None, fill = False, fillstyle = 0, linestyle = 1):
        #make sure things are set correctly
        if fillstyle != 0  and fill == False:
            fill = True
        if fillcolor is not None:
            fill = True
        if fillstyle == 0 and fill:
            fillstyle = 1001

        self.linecolor = 1 # set default of black if it is also filled
        if fill == False : self.linecolor = color # if not filled, line is main color
        self.fillcolor = color

        if linecolor is not None: self.linecolor = color
        if fillcolor is not None: self.fillcolor = color

def SetStyle(h, s):
    h.SetFillStyle(s.fillstyle)
    h.SetLineStyle(s.linestyle)
    h.SetFillColor(s.fillcolor)
    h.SetLineColor(s.linecolor)

def getStyleObjs(colors, linecolors, fillcolors, fillstyles, linestyles):
    styles = []
    #make copies
    colors_cpy = copy.copy(colors)
    linecolors_cpy = copy.copy(linecolors)
    fillcolors_cpy = copy.copy(fillcolors)
    linestyles_cpy = copy.copy(linestyles)
    fillstyles_cpy = copy.copy(fillstyles)

class compObj:
    def __init__(self, name, plot, c ):
        self.name = name
        self.plot = plot
        self.comp = c
        self.drawOpt = self.comp.drawOpt
        if not isinstance(self.plot, ROOT.TH1) and not isinstance(self.plot, ROOT.TGraph):
            print "Not a plot, no comparison added"
            return
        if not isinstance(self.comp.plot, ROOT.TH1) and not isinstance(c.plot, ROOT.TGraph):
            print "Trying to compare to not a plot, no comparison added"
            return
        self.ratio     = self.comp.plot.Clone(self.name+'ratiocomp')
        self.diff      = self.comp.plot.Clone(self.name+'diffcomp')
        self.rationorm = self.comp.plot.Clone(self.name+'rationormcomp')
        self.diffnorm  = self.comp.plot.Clone(self.name+'diffnormcomp')
        #unpleasant hack, will do for now
        self.ratio.GetYaxis().SetTickLength(0.029 * 0.85/0.2)
        self.diff.GetYaxis().SetTickLength(0.029 * 0.85/0.2)
        self.ratio.GetYaxis().SetNdivisions(505)
        self.diff.GetYaxis().SetNdivisions(505)
        self.rationorm.GetYaxis().SetNdivisions(505)
        self.diffnorm.GetYaxis().SetNdivisions(505)
        if 'TH1' in self.ratio.ClassName() and self.ratio.GetLabelSize() == 0:
            self.ratio.SetLabelSize(0.05)
            self.rationorm.SetLabelSize(0.05)
        if 'TH1' in self.ratio.ClassName() and self.ratio.GetLabelOffset() == 999:
            self.ratio.SetLabelOffset(0.02)
            self.rationorm.SetLabelOffset(0.02)
        if 'TH1' in self.diff.ClassName() and self.diff.GetLabelSize() == 0:
            self.diff.SetLabelSize(0.05)
            self.diffnorm.SetLabelSize(0.05)
        if 'TH1' in self.diff.ClassName() and self.diff.GetLabelOffset() == 999:
            self.diff.SetLabelOffset(0.02)
            self.diffnorm.SetLabelOffset(0.02)
        if ('TH1' in self.ratio.ClassName() and 'TH1' in self.plot.ClassName()):
            totA = self.ratio.Integral()
            totB = self.plot.Integral()
            for i in range(self.ratio.GetNbinsX()):
                if self.plot.GetBinContent(i+1) != 0:
                    a    = self.ratio.GetBinContent(i+1)
                    b    = self.plot.GetBinContent(i+1)
                    aerr = self.ratio.GetBinError(i+1)
                    berr = self.plot.GetBinError(i+1)
                    self.ratio.SetBinContent(i+1, a/b)
                    self.ratio.SetBinError(  i+1 , sqrt(pow(aerr/b,2) + pow(a*berr/pow(b,2),2)))
                    self.rationorm.SetBinContent(i+1 , (totB/totA)*a/b )
                    self.rationorm.SetBinError(  i+1 , (totB/totA)*sqrt(pow(aerr/b,2) + pow(a*berr/pow(b,2),2)))
                else:
                    self.ratio.SetBinContent(i+1, 1.0 )
                    self.ratio.SetBinError(i+1, 0.0 )
                    self.rationorm.SetBinContent(i+1, 1.0 )
                    self.rationorm.SetBinError(i+1, 0.0 )
                self.diff.SetBinContent(i+1, self.diff.GetBinContent(i+1) - self.plot.GetBinContent(i+1))
                self.diffnorm.SetBinContent(i+1, self.diff.GetBinContent(i+1)/totA - self.plot.GetBinContent(i+1)/totB)
        elif ('TH1' in self.ratio.ClassName() and 'TGraph' in self.plot.ClassName()):
            totA = self.ratio.Integral()
            totB = GetTGraphSum(self.plot)
            if self.ratio.GetNbinsX() == self.plot.GetN():
                x = ROOT.Double(0.0)
                y = ROOT.Double(0.0)
                for i in range(self.plot.GetN()):
                    self.plot.GetPoint(i, x, y)
                    a    = self.ratio.GetBinContent(i+1)
                    aerr = self.ratio.GetBinError(i+1)
                    if y != 0:
                        self.ratio.SetBinContent(i+1, a/ y)
                        self.ratio.SetBinError(i+1, aerr/ y)
                        self.rationorm.SetBinContent(i+1, (totB/totA)*a/y)
                        self.rationorm.SetBinError(i+1  , (totB/totA)*aerr/y)
                    else:
                        self.ratio.SetBinContent(i+1, 1)
                        self.ratio.SetBinError(i+1, 0)
                        self.rationorm.SetBinContent(i+1, 1.0 )
                        self.rationorm.SetBinError(i+1, 0.0 )
                    self.diff.SetBinContent(i+1, self.diff.GetBinContent(i+1) - y)
        elif ('TGraph' in self.ratio.ClassName() and 'TH1F' in self.plot.ClassName()):
            totA = GetTGraphSum(self.ratio)
            totB = self.plot.Integral()
            if self.ratio.GetN() == self.plot.GetNbinsX():
                x = ROOT.Double(0.0)
                y = ROOT.Double(0.0)
                for i in range(self.ratio.GetN()):
                    b    = self.plot.GetBinContent(i+1)
                    berr = self.plot.GetBinError(i+1)
                    self.ratio.GetPoint(i, x, y)
                    yerrhi = self.ratio.GetErrorYhigh(i)
                    yerrlo = self.ratio.GetErrorYlow(i)
                    self.ratio.SetPoint(i, x, y / b)
                    self.ratio.SetPointEYhigh(i, yerrhi/b )
                    self.ratio.SetPointEYlow(i, yerrlo/b)
                    self.ratio.SetPoint(i, x, y*(totB/totA) / b )
                    self.ratio.SetPointEYhigh(i, yerrhi*(totB/totA)/b )
                    self.ratio.SetPointEYlow(i , yerrlo*(totB/totA)/b )
                    self.diff.SetPoint(i, x, y - b)
                    self.diffnorm.SetPoint(i, x, y/totA - b/totB)
                    self.diffnorm.SetPointEYhigh(i, yerrhi*(totB/totA)/b)
                    self.diffnorm.SetPointEYlow(i,  yerrlo*(totB/totA)/b)

        elif ('TGraph' in clone.ClassName() and 'TGraph' in self.plot.ClassName()):
            totA = GetTGraphSum(self.ratio)
            totB = GetTGraphSum(self.plot)
            if self.ratio.GetN() == self.plot.GetN():
                x1 = ROOT.Double(0.0)
                y1 = ROOT.Double(0.0)
                x2 = ROOT.Double(0.0)
                y2 = ROOT.Double(0.0)
                self.ratio.GetPoint(i, x1, y1)
                self.plot.GetPoint(i, x2, y2)
                self.ratio.SetPoint(i, x1, y1/ y2)
                self.ratio.SetPointEYhigh(i, self.ratio.GetErrorYhigh(i)/y2)
                self.ratio.SetPointEYlow(i, self.ratio.GetErrorYlow(i)/y2)
                self.diff.SetPoint(i, x1, y1 - y2)

                self.rationorm.SetPoint(i, x1, (totB/totA)*y1/ y2)
                self.rationorm.SetPointEYhigh(i, (totB/totA)*self.ratio.GetErrorYhigh(i)/y2)
                self.rationorm.SetPointEYlow(i, (totB/totA)*self.ratio.GetErrorYlow(i)/y2)
                self.diffnorm.SetPoint(i, x1, y1/totA - y2/totB)
        else:
            print "Not valid plots for comparison: ",clone.ClassName(), self.plot.ClassName()

    def Style(self, xlabel, xlims, ycomplims):
        for p in [self.ratio, self.diff, self.rationorm, self.diffnorm]:
            SetROOTFillColor(p, self.comp.color)
            p.SetFillStyle(self.comp.style)
            SetROOTLineColor(p, self.comp.linecolor)
            p.SetLineStyle(self.comp.linestyle)
            if self.comp.markerstyle != 0 : p.SetMarkerStyle(self.comp.markerstyle)
            if xlims is not None and len(xlims) == 2: p.GetXaxis().SetRangeUser(xlims[0], xlims[1])
            if ycomplims is not None and len(ycomplims) == 2: p.GetYaxis().SetRangeUser(ycomplims[0], ycomplims[1])
            p.GetXaxis().SetTitle(xlabel)
            p.GetYaxis().SetTitle("")
            p.GetYaxis().SetNdivisions(505)
            p.GetXaxis().SetLabelSize(0.05)

class PlotObj:
    def __init__(self,name, plot, color, style, linecolor, linestyle, markerstyle, drawOpt, xticks, compareAsRatio = False):
        self.name      = name
        self.plot      = plot
        self.color     = color
        self.linecolor = linecolor
        self.style     = style
        self.linestyle = linestyle
        self.markerstyle = markerstyle
        self.drawOpt   = drawOpt
        self.compareAsRatio = False
        self.friends     = []
        self.comparisons  = []
        self.xticks = xticks
        if isinstance(self.plot, ROOT.TGraph):
            self.split     = TGraph2Hists(self.name, self.plot)
            self.norm      = makeNormErr2Graph(self.name, self.plot)
            self.normsplit = TGraph2Hists(self.name+"_norm", self.norm)
        elif isinstance(self.plot, ROOT.TH1):
            self.split = []
            self.norm = self.plot.Clone(self.name+"_norm")
            self.norm.Scale(1/self.norm.Integral())
            self.normsplit = []
               
    def ApplyStyle(self, xlabel, ylabel, xlims, ylims, ynormlims, ycomplims, offset = 0.0):
        if 'TLine' in self.plot.ClassName():
            print "This is a Tline"
        elif 'TBox' in self.plot.ClassName():
            print "This is a TBox"
        else:
            for p in [self.plot] + self.split:
                SetROOTFillColor(p, self.color)
                p.SetFillStyle(self.style)
                SetROOTLineColor(p, self.linecolor)
                if self.markerstyle != 0 : p.SetMarkerStyle(self.markerstyle)
                p.SetLineStyle(self.linestyle)
                if xlims is not None and len(xlims) == 2: p.GetXaxis().SetRangeUser(xlims[0], xlims[1])
                if ylims is not None and len(ylims) == 2:
                    p.GetYaxis().SetRangeUser(ylims[0], ylims[1])
                p.GetXaxis().SetTitle(xlabel)
                p.GetYaxis().SetTitle(ylabel)
                if offset != 0.0: p.GetYaxis().SetTitleOffset(offset)
            for p in [self.norm] + self.normsplit:
                if p != None:
                    SetROOTFillColor(p, self.color)
                    p.SetFillStyle(self.style)
                    SetROOTLineColor(p, self.linecolor)
                    p.SetLineStyle(self.linestyle)
                    if self.markerstyle != 0 : p.SetMarkerStyle(self.markerstyle)
                    if xlims is not None and len(xlims) == 2: p.GetXaxis().SetRangeUser(xlims[0], xlims[1])
                    if ynormlims is not None and len(ynormlims) == 2: 
                        p.GetYaxis().SetRangeUser(ynormlims[0], ynormlims[1])
                    p.GetXaxis().SetTitle(xlabel)
                    #p.GetYaxis().SetTitle("Rel. Err")
                    p.GetYaxis().SetNdivisions(505)
                    p.GetXaxis().SetLabelSize(0.05)
            for p in self.comparisons:
                p.Style(xlabel, xlims, ycomplims)

    def AddComparison(self, c):
        self.friends     += [c]
        self.comparisons += [ compObj(self.name+'_'+str(len(self.comparisons)), self.plot, c) ]

    def Draw(self, extra=""):
        self.plot.Draw(self.drawOpt+extra)
    def DrawCtr(self, extra=""):
        self.split[0].Draw(self.drawOpt+extra)
    def DrawLo(self, extra=""):
        self.split[1].Draw(self.drawOpt+extra)
    def DrawHi(self, extra=""):
        self.split[2].Draw(self.drawOpt+extra)
    def DrawNorm(self, extra=""):
        self.norm.Draw(self.drawOpt+extra)
    def DrawNormLo(self, extra=""):
        self.normsplit[1].Draw(self.drawOpt+extra)
    def DrawNormHi(self, extra=""):
        self.normsplit[2].Draw(self.drawOpt+extra)


class Plot:
   def __init__(self, plots, name='default'):
      ROOT.gROOT.SetBatch(True)
      self.plots = plots
      self.name = name
      self.properties = {
          #outputb
          "filename" : "PlotTool", "location" : ".",
          #labels
          "xlabel" : "", "ylabel" : "", 'ylowlabel' : "",
          #legend
          "legend": True, "labels" : [], "legwidth" : 0, "legheight" : 0, "legloc" : 2,
          #x and y axes limits
          "xlims" : None, "ylims" : None, "yloglims" : None, "ycomplims" : None, "ynormlims" : None,
          #plot scale
          "logscale" : False, 'logx' : False, "normalised" : False,
          "leglims" : [], "split" : False, "splitstyles" : [1001,3001,3001],
          "includeFit" : False,
          "includeNormErr" : False, "title" : "",
          #colors
          "colors" : None, "linecolors" : None, "markercolors": None,
          #styles
          "fillstyles" : 1001, "linestyles" : 1, "markerstyles" : 0, # colors and styles
          "extraObjs" : [],
          "drawOpts" : '', "drawOpts2" : None, "yoffset" : 0.0,
          "style" : True, "forcestyle" : False,
          #
          'canvas_size' : [700, 500],
          #
          'xticks' : 505,
          #comparison plots
          'toCompare' : {}, 'compareAsRatio' : True
          }

   def getPlotObjs(self, hide_label = False):
      #get colors and styles
      colors       = self.properties['colors']
      linecolors   = self.properties['linecolors']
      fillstyles   = self.properties['fillstyles']
      linestyles   = self.properties['linestyles']
      markerstyles = self.properties['markerstyles']
      compareAsRatio = self.properties['compareAsRatio']
      plots        = self.plots
      objs = []
      
      size  = len(plots)
      
      #if no colors are specified use existing hist colors
      if not colors:
          colors = [hist.GetLineColor() for hist in plots]
      # if just one instance make it an array
      elif isinstance(colors,str):
          colors = [self.properties['colors'] for i in range(size)]
      if isinstance(linecolors,str):
          linecolors = [self.properties['linecolors'] for i in range(size)]
      elif not linecolors:
          linecolors = colors
      #if just one option is specified for styles use it for all
      if isinstance(fillstyles, int):
          fillstyles     = [fillstyles for i in range(size)]
      if isinstance(linestyles, int):
          linestyles = [linestyles for i in range(size)]
      if isinstance(markerstyles, int):
          markerstyles = [markerstyles for i in range(size)]

      drawOpts = self.properties['drawOpts']
      drawOpts2 = self.properties['drawOpts2']
      if isinstance(drawOpts,str):
          drawOpts = [drawOpts for j in range(size)]
      if not drawOpts2:
          drawOpts2 = drawOpts
      elif isinstance(drawOpts2,str):
          drawOpts2 = [drawOpts2 for k in range(size)]

      xlabel      = self.properties['xlabel']
      ylabel      = self.properties['ylabel']
      xlims       = self.properties['xlims']
      ylims       = self.properties['ylims']
      normlims    = self.properties['ynormlims']
      complims    = self.properties['ycomplims']
      if complims is None : complims = [0.8, 1.2]
      offset      = self.properties['yoffset']
      xticks = self.properties['xticks']
      
      for plot,color,linecolor,fillstyle,linestyle,markerstyle,drawOpt,i in zip(plots, colors, linecolors, fillstyles, linestyles, markerstyles, drawOpts, range(size)):
          p = PlotObj(self.properties['filename']+"_"+str(color)+"_"+str(fillstyle)+"_"+str(i), plot,color,fillstyle, linecolor,linestyle, markerstyle, drawOpt, xticks)
          p.ApplyStyle(xlabel, ylabel, xlims, ylims, normlims, complims, offset)    
          objs += [p]
      for c in self.properties['toCompare'].keys():
          for c2 in self.properties['toCompare'][c]:
              if isinstance(c,int) and isinstance(c2,int) and c >= 0 and c < len(objs) and c2 >= 0 and c2 < len(objs):
                  objs[c].AddComparison(objs[c2])

      for o in objs:
          o.ApplyStyle(xlabel, ylabel, xlims, ylims, normlims, complims, offset)    
          if hide_label and (isinstance(o.plot, ROOT.TH1) or isinstance(o.plot, ROOT.TGraph)):
              for h in [o.plot ] + o.split:
                  h.GetXaxis().SetLabelOffset(999)
                  h.GetXaxis().SetLabelSize(0)
                  h.GetXaxis().SetTitle("")
      return objs



   def getPads(self, n, large = False, x1 = 0.005, x2 = 0.995, y0 = 0.05, y1 = 0.32, y2 = 0.995, bm = 0.14, y3 = 0.46):
       pads = []
       #x1 = 0.005
       #x2 = 0.995
       #y0 = 0.05
       #y1 = 0.32
       #y2 = 0.995
       #bm = 0.12
       #y3 = 0.46

       if n > 0 :
           pads += [ROOT.TPad("upperPad", "upperPad", x1, y0, x2, y2)]
       if n > 1 :
           pads += [ROOT.TPad("lowerPad", "lowerPad", x1, y0, x2, y2)]
       if n > 2 :
           pads += [ROOT.TPad("middlPad", "middlPad", x1, y0, x2, y2)]
       #another pad covering the whole bottom for axes
       if n > 1 :
           pads += [ROOT.TPad("bottomPad", "bottomPad", x1, y0, x2, y2)]

       if n == 1:
           pads[0].SetBottomMargin(bm)
       if n == 2:
           pads[0].SetBottomMargin(y1-y0)
           pads[1].SetBottomMargin(bm)
           pads[1].SetTopMargin(y2-(y1-y0))
           if large:
               pads[0].SetBottomMargin(y3 - y0)
               pads[1].SetTopMargin(y2 - (y3 - y0))
       if n == 3:
           pads[0].SetBottomMargin(y3-y0)
           pads[2].SetBottomMargin(y1-y0)
           pads[1].SetBottomMargin(bm)
           pads[2].SetTopMargin(y2 - (y3 - y0))
           pads[1].SetTopMargin(y2 - (y1 - y0))
           pads[-1].SetTopMargin(y2-(y3-y0))
           pads[-1].SetBottomMargin(bm)
                                                        
       for pad in pads: pad.SetFillStyle(0)
                                                        
       return pads
         
   def GetLegend(self):
       leglims   = self.properties['leglims']
       legloc    = self.properties['legloc']
       legheight = self.properties['legheight']
       legwidth  = self.properties['legwidth']
       labels    = self.properties['labels']
       title     = self.properties['title']
       drawOpts      = self.properties['drawOpts']
       markerstyles  = self.properties['markerstyles']

       styles = []
       if isinstance(drawOpts,str):
           drawOpts = [drawOpts for j in range(len(labels))]
       if isinstance(markerstyles,int):
           markerstyles = [drawOpts for j in range(len(labels))]
       for d,m in zip(drawOpts,markerstyles):
           if 'p' in d or m != 0:
               styles += ['lep']
           elif 'f' in d:
               styles += ['f']
           else:
               styles += ['l']
                      
          

       if (len(leglims) != 4):
           legheight = 0.15
           legwidth  = 0.5
           if (legloc == 1):#top left
               leglims = [0.1, 0.9 - legheight, legwidth, 0.9]
           elif (legloc == 2):#top right
               leglims = [0.9 - legwidth, 0.9 - legheight, 0.9, 0.9]
           elif (legloc == 3):#bottom right
               leglims = [0.9 - legwidth, 0.1, 0.9, 0.1 + legheight]
           elif (legloc == 4):#bottom left
               leglims = [0.1, 0.1, 0.1 + legwidth, 0.1 + legheight]
       
       leg = ROOT.TLegend(leglims[0],leglims[1],leglims[2],leglims[3])
       leg.SetFillColor(10)
       leg.SetFillStyle(0)
       for hist, label, style in zip(self.plots, labels, styles):
           leg.AddEntry(hist, label, style)
       return leg
                                                       


   def drawROOT(self, filename = None):
      plots  = self.plots
      if filename == None: filename = self.properties['filename']

      includeNormErr = self.properties['includeNormErr']

      canvas_size = self.properties['canvas_size']

      c = ROOT.TCanvas(self.name+"_Canvas", self.name+"Canvas", canvas_size[0], canvas_size[1])
      if self.properties['forcestyle'] : 
          ROOT.gROOT.ForceStyle(True)
          c.UseCurrentStyle()
      self.pads = []
      x1 = 0.005
      x2 = 0.995
      y0 = 0.05
      y1 = 0.32
      y2 = 0.995
      bm = 0.12
    
      if includeNormErr or len(self.properties['toCompare'].keys()) == 1:
          self.pads = Plot.getPads(self,2, y1 = 0.38, x1 = 0.005)
      elif len(self.properties['toCompare'].keys()) == 2:
          self.pads = Plot.getPads(self, 3)
      else:
          self.pads = Plot.getPads(self, 1)
        

      self.pads[-1].Draw("epl")
      for i in range(len(self.pads)-1): self.pads[i].Draw("epl")

      if (self.properties["logscale"]):
          self.pads[0].SetLogy()
      if (self.properties["logx"])    : 
          for p in self.pads: 
              p.SetLogx()

      hide_label = False
      if len(self.pads) > 1: hide_label = True
      
      self.plotObjs = self.getPlotObjs(hide_label = hide_label)

      self.pads[0].cd()
      for i,p in enumerate(self.plotObjs):
          drawOpt = ''
          if isinstance(p.plot, ROOT.TGraph) and not isinstance(p.plot, ROOT.RooCurve): drawOpt += 'p'
          if i == 0 and isinstance(p.plot, ROOT.TGraph) and not isinstance(p.plot, ROOT.RooCurve): drawOpt += 'A'
          if i > 0: drawOpt += 'same'
          if self.properties['logx'] and (isinstance(p.plot, ROOT.TH1) or isinstance(p.plot, ROOT.TGraph)): 
              p.plot.GetXaxis().SetNdivisions(self.properties['xticks'])

          if self.properties['normalised']: 
              #p.plot.Scale(1/p.plot.Integral())
              p.DrawNorm(drawOpt)
          else : 
              p.Draw(drawOpt)

      leg = self.GetLegend()
      if self.properties['legend']: leg.Draw("same")
      for obj in self.properties["extraObjs"]: obj.Draw("same")

      if len(self.properties['toCompare']) > 1:
          self.pads[-1].cd()
          h = ROOT.TH1F()
          h.GetXaxis().SetLabelOffset(999)
          h.GetXaxis().SetLabelSize(0)
          h.GetXaxis().SetTitle("")
          h.GetYaxis().SetTitle('Ratio')
          h.GetYaxis().SetLabelOffset(999)
          h.SetFillStyle(0)
          h.GetYaxis().CenterTitle()
          h.Draw()

      if len(self.properties['toCompare']) > 0:
          for i,k in enumerate(self.properties['toCompare'].keys()):
              self.pads[i+1].cd()
              p = self.plotObjs[k]
              for j,pp in enumerate(p.comparisons):
                  drawOpt = ''
                  if i == 0 and j == 0 and (isinstance(pp.diff, ROOT.TGraph)): drawOpt += 'A'
                  if i > 0 or j > 0 : drawOpt += 'same'
                  if len(self.pads) == 3 and i == 0:
                      pp.diff.GetYaxis().SetTitle('Diff.')
                      pp.ratio.GetYaxis().SetTitle('Ratio')
                      pp.ratio.GetYaxis().CenterTitle()
                      pp.diff.GetYaxis().SetTitleOffset(1.2)
                      pp.diffnorm.GetYaxis().SetTitle('Diff.')
                      pp.rationorm.GetYaxis().SetTitle('Ratio')
                      pp.rationorm.GetYaxis().CenterTitle()
                      pp.diffnorm.GetYaxis().SetTitleOffset(1.2)
                  if len(self.pads) == 4 and i != 0 :
                      pp.ratio.GetXaxis().SetLabelOffset(999)
                      pp.ratio.GetXaxis().SetLabelSize(0)
                      pp.ratio.GetXaxis().SetTitle("")
                      pp.diff.GetXaxis().SetLabelOffset(999)
                      pp.diff.GetXaxis().SetLabelSize(0)
                      pp.diff.GetXaxis().SetTitle("")

                      pp.rationorm.GetXaxis().SetLabelOffset(999)
                      pp.rationorm.GetXaxis().SetLabelSize(0)
                      pp.rationorm.GetXaxis().SetTitle("")
                      pp.diffnorm.GetXaxis().SetLabelOffset(999)
                      pp.diffnorm.GetXaxis().SetLabelSize(0)
                      pp.diffnorm.GetXaxis().SetTitle("")
                  if self.properties['compareAsRatio']:
                      if self.properties['normalised']:
                          pp.rationorm.Draw(pp.drawOpt+drawOpt)
                      else:
                          pp.ratio.Draw(pp.drawOpt+drawOpt)
                  else:
                      if self.properties['normalised']:
                          pp.diffnorm.Draw(pp.drawOpt+drawOpt)
                      else:
                          pp.diff.Draw(pp.drawOpt+drawOpt)
                  
      if '.pdf' not in filename and '.eps' not in filename and '.C' not in filename:
          filenames = [filename + '.pdf', filename+'.eps', filename+'.C', filename+'.png']
      else: filenames = [filename]
      location = self.properties["location"]
      
      #for p in self.plotObjs: 
      #    if isinstance(p, ROOT.TH1F): p.Draw("SAMEAXIS")
      for pad in self.pads:
          pad.RedrawAxis()
      for fName in filenames:
          c.Print(location + '/' + fName) ## Include this line to save image
      return c

      '''
      #get colors and styles
      colors = self.properties['colors']
      linecolors = self.properties['linecolors']
      styles = self.properties['styles']
      linestyles = self.properties['linestyles']
      

      #if no colors are specified use existing hist colors
      if not colors:
         colors = [hist.GetLineColor() for hist in plots]
      elif isinstance(colors,str):
          colors = [self.properties['colors'] for i in range(len(plots))]
      if isinstance(linecolors,str):
          linecolors = [self.properties['linecolors'] for i in range(len(plots))]
      elif not linecolors:
          linecolors = colors
      #if just one option is specified for styles use it for all
      if isinstance(styles, int):
          styles = [styles for i in range(len(plots))]
      if isinstance(linestyles, int):
          linestyles = [linestyles for i in range(len(plots))]

      for hist, color,style,linecolor,linestyle in zip(plots, colors,styles, linecolors, linestyles):
         if self.properties['normalised']:
            hist.Scale(1/hist.Integral())
         if color == 'None':
             hist.SetFillStyle(0)
         else:
             SetROOTFillColor(hist,color)
         SetROOTLineColor(hist,linecolor)
         hist.SetFillStyle(style)
         hist.SetLineStyle(linestyle)

      dummyhist = ROOT.TH1F("dummy","dummy",10,0,1)
      dummyhist.SetLineColor(10)
      dummyhist.SetFillColor(10)
      dummyhist.SetMarkerColor(10)

      leglims   = self.properties['leglims']
      legloc    = self.properties['legloc']
      legheight = self.properties['legheight']
      legwidth  = self.properties['legwidth']
      labels    = self.properties['labels']
      title     = self.properties['title']

      if (len(leglims) != 4):
         legheight = 0.15
         legwidth  = 0.5
         if (legloc == 1):#top left
            leglims = [0.1, 0.9 - legheight, legwidth, 0.9]
         elif (legloc == 2):#top right
            leglims = [0.9 - legwidth, 0.9 - legheight, 0.9, 0.9]
         elif (legloc == 3):#bottom right
            leglims = [0.9 - legwidth, 0.1, 0.9, 0.1 + legheight]
         elif (legloc == 4):#bottom left
            leglims = [0.1, 0.1, 0.1 + legwidth, 0.1 + legheight]


      leg = Plot.GetLegend(self)
      #ROOT.TLegend(leglims[0],leglims[1],leglims[2],leglims[3])
      #leg.SetFillColor(10)
      #leg.SetFillStyle(0)
      #if (title != ""): leg.AddEntry(dummyhist, title)
      #for i, hist, label in zip(range(len(plots)),plots, labels):
      #        leg.AddEntry(hist, label, "l")

      #ylimit = plot.GetMaximum()*2
      #plt.ylim(100,ylimit)

      xlabel     = self.properties['xlabel']
      ylabel     = self.properties['ylabel']
      xlims      = self.properties['xlims']
      ylims      = self.properties['ylims']

      if includeNormErr:
         upperPad.cd()
         if (self.properties["logscale"]): upperPad.SetLogy()
      
      i = 0
      drawOpts = self.properties['drawOpts']
      drawOpts2 = self.properties['drawOpts2']
      if isinstance(drawOpts,str):
          drawOpts = [drawOpts for j in range(len(plots))]
      if not drawOpts2:
          drawOpts2 = drawOpts
      elif isinstance(drawOpts2,str):
          drawOpts2 = [drawOpts2 for k in range(len(plots))]
      for hist, drawOpt in zip(plots, drawOpts):
         if i == 0:
             if xlims is not None: hist.GetXaxis().SetLimits(xlims[0],xlims[1])
             if ylims is not None: hist.GetYaxis().SetRangeUser(ylims[0],ylims[1])
             hist.GetXaxis().SetTitle(xlabel)
             hist.GetYaxis().SetTitle(ylabel)
             if ('TGraph' in hist.ClassName()):
                 hist.Draw("A"+drawOpt)
             else:
                 hist.Draw(drawOpt)
         else:
             hist.Draw(drawOpt+"same")
         i = i + 1
      if self.properties['legend']: leg.Draw("same")
      for obj in self.properties["extraObjs"]: obj.Draw("same")

      if includeNormErr:
         lowerPad.cd()
         plots[0].GetXaxis().SetLabelSize(0);
         i = 0
         normerrplots = [makeNormErr(plot,str(ii)) for plot,ii in zip(plots,range(len(plots)))]
         normerrplots[0].GetXaxis().SetLimits(2,4.5)
         for errplot,color,style,drawOpt in zip(normerrplots,colors,styles,drawOpts2):
             SetROOTLineColor(errplot,color)
             SetROOTFillColor(errplot,color)
             errplot.GetXaxis().SetTitle(xlabel)
             errplot.GetYaxis().SetTitle("Re. Err.")
             errplot.GetYaxis().SetNdivisions(505)
             errplot.GetXaxis().SetLabelSize(0.05)
             errplot.SetFillStyle(style)
             if xlims is not None :
                 errplot.GetXaxis().SetRangeUser(xlims[0], xlims[1])
             if (self.properties['normerrlims']):
                 normerrlims = self.properties['normerrlims']
                 errplot.GetYaxis().SetRangeUser(normerrlims[0],normerrlims[1])
             else: errplot.GetYaxis().SetRangeUser(0.8,1.2)
             if i == 0:
                 errplot.Draw(drawOpt)
             else:
                 errplot.Draw(drawOpt+"same")
             i = i + 1

      if '.pdf' not in filename and '.eps' not in filename:
              filenames = [filename + '.pdf', filename + '.eps']
      else:
          filenames = [filename]
      location = self.properties["location"]

      if includeNormErr:
         lowerPad.Update()
         upperPad.Update()
         c.Update()
      c.RedrawAxis()
      if self.properties['forcestyle'] : 
          ROOT.gROOT.ForceStyle(True)
          c.UseCurrentStyle()
   
      for fname in filenames:
          c.Print(location + '/' + fname) ## Include this line to save image
      return c
      '''
   
   def ShiftLegX(self, x):
       oldleg = self.properties['leglims']
       self.properties['leglims'] = [oldleg[0] + x, oldleg[1], oldleg[2] + x, oldleg[3]]
   def ShiftLegY(self, y):
       oldleg = self.properties['leglims']
       self.properties['leglims'] = [oldleg[0], oldleg[1] + y, oldleg[2], oldleg[3] + y]
   def setProps(self, props):
      for key in props:
         self.properties[key] = props[key]
   def setProp(self, keyname, prop):
       key_exists = False
       for key in self.properties:
           if key == keyname:
               self.properties[key] = prop
               key_exists = True
       if key_exists == False: print "Did not find property: ",keyname
   def getProp(self, keyname):
       if keyname in self.properties.keys():
           return self.properties[keyname]
       else:
           'Property '+keyname+' does not exist'
           return
       #print 'set ',keyname,' to ',prop
               

