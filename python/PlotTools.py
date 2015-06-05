import ROOT
from random import gauss
from types import IntType
from collections import OrderedDict
import time

#ROOT.gStyle.SetTextFont(80)
#ROOT.gStyle.SetLabelFont(80,"XYZ")
#ROOT.gStyle.SetTitleFont(80,"XYZ")
#ROOT.gStyle.SetTitleFont(80)
#ROOT.gStyle.SetOptStat(0)
#from lhcbStyle import SetStyle()

def split_list(list_items, split_idx):
    split_list = []
    #round up
    splits = round(len(list_items)/split_idx + 0.5)
    for i in range(splits):
        min_idx = i * split_idx
        max_idx = (i + 1) * split_idx
        max_idx = min(max_idx, len(list_items) -1)
        split_list += list_items[min_idx:max_idx]

def TGraph2Hists(name, graph):
    if ("TGraph" in graph.ClassName()):
        xmin = ROOT.Double(0.0)
        xmax = ROOT.Double(0.0)
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
      print "Filing manually for color: ",mcolor
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
   if ("TGraph" in grapha.ClassName()):
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
            normgraph.SetPointEYhigh(i, graph.GetErrorYhigh(i)/val1y)
            normgraph.SetPointEYlow (i, graph.GetErrorYlow(i)/val1y)
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

class PlotObj:
    def __init__(self,name, plot, color, style, linecolor, linestyle, drawOpt):
        self.name      = name
        self.plot      = plot
        self.color     = color
        self.linecolor = linecolor
        self.style     = style
        self.linestyle = linestyle
        self.drawOpt   = drawOpt
        if ('TGraph' in self.plot.ClassName()):
            self.split     = TGraph2Hists(self.name, self.plot)
            self.norm      = makeNormErr2Graph(self.name, self.plot)
            self.normsplit = TGraph2Hists(self.name+"_norm", self.norm)
        else:
            self.split = []
            self.norm = None
            self.normsplit = []
               
    def ApplyStyle(self, xlabel, ylabel, xlims, ylims, ynormlims, offset = 0.0):
        if 'TLine' in self.plot.ClassName():
            print "This is a Tline"
        else:
            for p in [self.plot] + self.split:
                SetROOTFillColor(p, self.color)
                p.SetFillStyle(self.style)
                SetROOTLineColor(p, self.linecolor)
                p.SetLineStyle(self.linestyle)
                p.GetXaxis().SetRangeUser(xlims[0], xlims[1])
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
                    p.GetXaxis().SetRangeUser(xlims[0], xlims[1])
                    p.GetYaxis().SetRangeUser(ynormlims[0], ynormlims[1])
                    p.GetXaxis().SetTitle(xlabel)
                    #p.GetYaxis().SetTitle("Rel. Err")
                    p.GetYaxis().SetNdivisions(505)
                    p.GetXaxis().SetLabelSize(0.05)

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
   def __init__(self, plots):
      self.plots = plots
      self.properties = {"xlabel" : "", "ylabel" : "", "filename" : "", # labels
                         "logscale": False, "legwidth" : 0, "legheight" : 0, "legloc" : 2, "legend" : True,
                         "leglims" : [], "split" : False, "splitstyles" : [1001,3001,3001],
                         "labels" : [], "includeFit" : False,
                         "xlims" : None, "ylims" : None,
                         "includeNormErr" : False, "title" : "",
                         "colors" : None, "linecolors" : None, "styles" : 1001, "linestyles" : 1, # colors and styles
                         "normerrlims" : None, 
                         "mplcolors" : [], "mpllabels" : [], "mplxlabel" : "", "mplylabel" : "",
                         'mpltitle' : "", 'normalised' : False,
                         "location" : "/Users/sfarry/WTau/Graphics/PDF", "extraObjs" : [],
                         "drawOpts" : 'e1', "drawOpts2" : None, "yoffset" : 0.0
                         }

   def getPlotObjs(self, hide_label = False):
      #get colors and styles
      colors      = self.properties['colors']
      linecolors  = self.properties['linecolors']
      styles      = self.properties['styles']
      linestyles  = self.properties['linestyles']
      plots       = self.plots
      objs = []
      
      size  = len(plots)
      
      #if no colors are specified use existing hist colors
      if not colors:
          colors = [hist.GetLineColor() for hist in plots]
      elif isinstance(colors,str):
          colors = [self.properties['colors'] for i in range(size)]
      if isinstance(linecolors,str):
          linecolors = [self.properties['linecolors'] for i in range(size)]
      elif not linecolors:
          linecolors = colors
      #if just one option is specified for styles use it for all
      if isinstance(styles, int):
          styles     = [styles for i in range(size)]
      if isinstance(linestyles, int):
          linestyles = [linestyles for i in range(size)]

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
      normerrlims = self.properties['normerrlims']
      if normerrlims is None : normerrlims = [0.8, 1.2]
      offset      = self.properties['yoffset']


      for plot,color,linecolor,style,linestyle,drawOpt,i in zip(plots, colors, linecolors, styles, linestyles,drawOpts, range(size)):
          print "initialising object with linestyle "+str(linestyle)
          p = PlotObj(self.properties['filename']+"_"+str(color)+"_"+str(style)+"_"+str(i), plot,color,style, linecolor,linestyle, drawOpt)
          p.ApplyStyle(xlabel, ylabel, xlims, ylims, normerrlims, offset)    
          if hide_label:
              for h in [p.plot ] + p.split:
                  h.GetXaxis().SetLabelOffset(999)
                  h.GetXaxis().SetLabelSize(0)
                  h.GetXaxis().SetTitle("")
          objs += [p]
      return objs



   def getPads(self, n, large = False, x1 = 0.005, x2 = 0.995, y0 = 0.05, y1 = 0.32, y2 = 0.995, bm = 0.12, y3 = 0.46):
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
                                                        
       for pad in pads: pad.SetFillStyle(0)
                                                        
       return pads
         
   def GetLegend(self):
       leglims   = self.properties['leglims']
       legloc    = self.properties['legloc']
       legheight = self.properties['legheight']
       legwidth  = self.properties['legwidth']
       labels    = self.properties['labels']
       title     = self.properties['title']
       drawOpts  = self.properties['drawOpts']


       styles = []
       if isinstance(drawOpts,str):
           drawOpts = [drawOpts for j in range(len(labels))]
       for d in drawOpts:
           if 'p' in d:
               styles += ['lep']
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

      c = ROOT.TCanvas()
      if includeNormErr:
         x1 = 0.005
         x2 = 0.995
         y0 = 0.05
         y1 = 0.32
         y2 = 0.995
         bm = 0.12
         upperPad = ROOT.TPad("upperPad", "upperPad", x1, y0, x2, y2)
         lowerPad = ROOT.TPad("lowerPad", "lowerPad", x1, y0, x2, y2)
         upperPad.Draw("epl")
         lowerPad.Draw("epl")
         upperPad.SetBottomMargin(y1-y0)
         lowerPad.SetBottomMargin(bm)
         lowerPad.SetFillStyle(0)
         lowerPad.SetTopMargin(y2 - (y1 - y0))

      if (self.properties["logscale"]): c.SetLogy()

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
      for fname in filenames:
          c.Print(location + '/' + fname) ## Include this line to save image
      return c
   
   def ShiftLegX(self, x):
       oldleg = self.properties['leglims']
       self.properties['leglims'] = [oldleg[0] + x, oldleg[1], oldleg[2] + x, oldleg[3]]
   def ShiftLegY(self, y):
       oldleg = self.properties['leglims']
       self.properties['leglims'] = [oldleg[0], oldleg[1] + y, oldleg[2], oldleg[3] + y]
   def setProps(self, props):
      for key in props:
         self.properties[key] = props[key]


class WJetPlot(Plot):
    def drawROOT(self, filename = None):
        plots  = self.plots
        includeNormErr = self.properties['includeNormErr']
        if filename == None: filename = self.properties['filename']

        if len(plots) != 4:
            print "We need four plots for the WJet Plot"
            return

        c = ROOT.TCanvas("c1", "c1", 600, 600)

        pads = Plot.getPads(self,3)
        
        for pad in pads: pad.Draw("epl")

        if (self.properties["logscale"]): c.SetLogy()

        #get plot objs
        plotObjs = Plot.getPlotObjs(self, hide_label = True)

        pads[0].cd()

        leg = Plot.GetLegend(self)
      
        k = 0
        for plot in plotObjs:
            if k == 0: 
                plot.DrawLo("][")
            else:
                plot.DrawLo("same][")
            plot.DrawHi("same][")
            k = k+1

        if self.properties['legend']: leg.Draw("same")
        for obj in self.properties["extraObjs"]: obj.Draw("same")

        pads[1].cd()
        i = 0
        for plot in plotObjs[2:4]:
            if i == 0:
                plot.DrawNormLo("][")
            else:
                plot.DrawNormLo("same][")
            plot.DrawNormHi("same][")
            i = i + 1
        pads[2].cd()
        i = 0
        for plot in plotObjs[0:2]:
            plot.normsplit[0].GetXaxis().SetLabelOffset(999)
            plot.normsplit[0].GetXaxis().SetLabelSize(0)
            plot.normsplit[0].GetXaxis().SetTitleSize(0)
            plot.normsplit[0].SetTitle("")
            plot.normsplit[1].GetXaxis().SetLabelOffset(999)
            plot.normsplit[1].GetXaxis().SetLabelSize(0)
            plot.normsplit[1].GetXaxis().SetTitleSize(0)
            plot.normsplit[1].SetTitle("")
            if i == 0:
                plot.DrawNormLo("][")
            else:
                plot.DrawNormLo("same][")
            plot.DrawNormHi("same][")
            i = i + 1        

        if '.pdf' not in filename and '.eps' not in filename:
            filenames = [filename + '.pdf', filename+'.eps']
        else: filenames = [filename]
        location = self.properties["location"]

        for pad in pads: pad.Update()
        c.Update()

        for fName in filenames:
            c.Print(location + '/' + fName) ## Include this line to save image
        return c

class WAsymmPlot(Plot):
    def __init__(self, plothi, plotlo):
        self.PlotHi = plothi
        self.PlotLo = plotlo
        self.properties = {
            'legend' : True,
            'filename' : 'file',
            'logscale' : False,
            'extraObjs' : [],
            'location' : None
            }
        for prop in ['legheight', 'legwidth','labels', 'title', 'drawOpts', 'leglims', 'legloc']:
            self.properties[prop] = plothi.properties[prop]

    def drawROOT(self, filename = None):
        plotshi    = self.PlotHi.getPlotObjs(hide_label = True)
        plotslo    = self.PlotLo.getPlotObjs()

        if filename == None: filename = self.properties['filename']

        c = ROOT.TCanvas("c1"+str(time.time()), "c1", 800, 700)

        pads = Plot.getPads(self,2, y1 = 0.38)
        
        for pad in pads: pad.Draw("epl")

        if (self.properties["logscale"]): c.SetLogy()

        #get plot objs

        pads[0].cd()

        k = 0
        for plot in plotshi:
            if k == 0: 
                if ('TGraph' in plot.plot.ClassName()):
                    plot.Draw("A")
                else:
                    plot.Draw()
            else:
                plot.Draw("same")
            k = k+1
        
        for obj in self.PlotHi.properties['extraObjs']: obj.Draw("same")

        for obj in self.properties["extraObjs"]: obj.Draw("same")

        pads[1].cd()
        i = 0
        plotslo[0].plot.GetYaxis().SetNdivisions(505)
        plotslo[0].plot.GetXaxis().SetLabelSize(0.05)

        for plot in plotslo:
            if i == 0:
                plot.Draw("][")
            else:
                plot.Draw("same][")
            i = i + 1

        for obj in self.PlotLo.properties['extraObjs']: obj.Draw("same")

        if '.pdf' not in filename and '.eps' not in filename:
            filenames = [filename + '.pdf', filename+'.eps']
        else: filenames = [filename]
        location = self.properties["location"]

        for pad in pads: pad.Update()
        c.Update()

        for fName in filenames:
            c.Print(location + '/' + fName) ## Include this line to save image
        return c

class StackPlot:
   def __init__(self, h, stack, fit = None):
      self.hist       = h
      self.stack      = stack
      self.fit        = fit
      self.properties = {"xlabel" : "", "ylabel" : "", "filename" : "",
                         "logscale": False, "legwidth" : 0, "legheight" : 0, "legloc" : 2,
                         "leglims" : [], "labels" : [], "WithLegend": True, "title" : "", "includeFit" : False,
                         "includeResidual" : False, "title" : "", "colors" : [],
                         "mplcolors" : [], "mpllabels" : [], "mplxlabel" : "", "mplylabel" : "",
                         'mpltitle' : "",
                         "location" : "/Users/sfarry/WTau/Graphics/PDF", "extraObjs" : []
                         ,"stacked" : True, "stackstyles" : [3004, 3005, 3006, 3002, 3003, 3944]
                         }

   def drawROOT(self, filename = None):
      plot  = self.hist
      stack = self.stack
      fit   = self.fit
      
      includeResidual = self.properties['includeResidual']

      c = ROOT.TCanvas()
      if includeResidual:
         x1 = 0.005
         x2 = 0.995
         y0 = 0.05
         y1 = 0.32
         y2 = 0.995
         bm = 0.12
         upperPad = ROOT.TPad("upperPad", "upperPad", x1, y0, x2, y2)
         lowerPad = ROOT.TPad("lowerPad", "lowerPad", x1, y0, x2, y2)
         upperPad.Draw("epl")
         lowerPad.Draw("epl")
         upperPad.SetBottomMargin(y1-y0)
         lowerPad.SetBottomMargin(bm)
         lowerPad.SetFillStyle(0)
         lowerPad.SetTopMargin(y2 - (y1 - y0))

      if (self.properties["logscale"]): c.SetLogy()
      histlist      = stack.GetHists()
      colors = []

      if len(self.properties["colors"]) > 0:
         colors = self.properties["colors"]
      else:
         colors = [hist.GetFillColor() for hist in histlist]

      stacked     = self.properties["stacked"]
      stackstyles = self.properties["stackstyles"]

      for hist, color, style in zip(histlist, colors, stackstyles):
         if stacked:
            SetROOTFillColor(hist, color)
            hist.SetLineColor(ROOT.kBlack)
         else:
            SetROOTLineColor(hist,color)
            SetROOTFillColor(hist, color)
            hist.SetFillStyle(style)

      if (fit): fit.SetLineColor(ROOT.kRed)

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

      if (self.properties['WithLegend']):
          leg = ROOT.TLegend(leglims[0],leglims[1],leglims[2],leglims[3])
          leg.SetFillColor(10)
          leg.SetFillStyle(0)
          leg.AddEntry(plot, "2012 Data")
          for hist, label in zip(histlist, labels):
              leg.AddEntry(hist, label, 'f')

      #ylimit = plot.GetMaximum()*2
      #plt.ylim(100,ylimit)

      xlabel     = self.properties['xlabel']
      ylabel     = self.properties['ylabel']
      includeFit = self.properties['includeFit']

      if includeResidual: upperPad.cd()
      if stacked:
         stack.Draw("hist")
         stack.GetXaxis().SetTitle(xlabel)
         stack.GetYaxis().SetTitle(ylabel)
      else:
         plot.GetXaxis().SetTitle(xlabel)
         plot.GetYaxis().SetTitle(ylabel)
         plot.Draw("e1")
         for hist in histlist:
            hist.Draw("histsame")
            
      if includeFit :
         fit.Draw("same")
      plot.Draw("E1same")

      if self.properties['WithLegend']:
          leg.Draw("same")
      for obj in self.properties["extraObjs"]: obj.Draw("same")

      if includeResidual: 
         lowerPad.cd()
         #stack.GetXaxis().SetLabelOffset(999);
         stack.GetXaxis().SetLabelSize(0);
         com = makeResidual(plot, fit)
         com.SetFillColor(ROOT.kBlack)
         com.GetYaxis().SetRangeUser(-3.9,3.9)
         com.GetXaxis().SetTitle(xlabel)
         com.GetYaxis().SetTitle("Residual")
         com.GetYaxis().SetNdivisions(505)
         com.GetXaxis().SetLabelSize(0.05)
         com.Draw("e1")
         com.GetXaxis().SetLabelSize(0);
      #if includeResidual:
      #else:
      if filename is None:
         filename = self.properties["filename"]
      location = self.properties["location"]

      if includeResidual:
         lowerPad.Update()
         upperPad.Update()
         c.Update()

      #plot.Delete()
      #fit.Delete()
      #stack.Delete()
      if '.pdf' not in filename and '.eps' not in filename and '.png' not in filename:
          filename = filename+'.pdf'
      c.Print(location + '/' + filename) ## Include this line to save image
      return c

   def setProps(self, props):
      for key in props:
         self.properties[key] = props[key]
