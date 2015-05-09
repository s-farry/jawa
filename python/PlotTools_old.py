import ROOT
from random import gauss
from types import IntType
from collections import OrderedDict

#ROOT.gStyle.SetTextFont(80)
#ROOT.gStyle.SetLabelFont(80,"XYZ")
#ROOT.gStyle.SetTitleFont(80,"XYZ")
#ROOT.gStyle.SetTitleFont(80)
#ROOT.gStyle.SetOptStat(0)
#from lhcbStyle import SetStyle()

 
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

def makeNormErr(graph):
   if ("TGraph" in graph.ClassName()):
      xmin = ROOT.Double(0.0)
      xmax = ROOT.Double(0.0)
      ymin = ROOT.Double(0.0)
      ymax = ROOT.Double(0.0)
      
      graph.GetPoint(0, xmin, ymin)
      graph.GetPoint(graph.GetN() - 1, xmax, ymax)
      
      xmin = xmin - (graph.GetErrorXlow(2))
      xmax = xmax + (graph.GetErrorXhigh(2))
      
      hist = ROOT.TH1F("hist","hist",graph.GetN(), xmin, xmax)
      for bin in range(0, graph.GetN()):
         val1x = ROOT.Double(0.0)
         val1y = ROOT.Double(0.0)
         graph.GetPoint(bin, val1x, val1y)
         hist.SetBinContent(bin+1,1.0)
         valerr = max(graph.GetErrorYhigh(bin), graph.GetErrorYlow(bin))/val1y
         hist.SetBinError(bin+1, valerr)
      return hist
   elif ("TH1" in graph.ClassName()):
      hist = graph.Clone("CompGraph")
      for bin in range(0, hist.GetNbinsX() + 2):
         if graph.GetBinContent(bin) != 0:
            hist.SetBinContent(bin, 1.0)
            hist.SetBinError(bin, abs(graph.GetBinError(bin)/graph.GetBinContent(bin)))
         else:
            hist.SetBinContent(bin, 1.0)
            hist.SetBinError(bin, 0)
      return hist

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


class Plot:
   def __init__(self, plots):
      self.plots = plots
      self.properties = {"xlabel" : "", "ylabel" : "", "filename" : "",
                         "logscale": False, "legwidth" : 0, "legheight" : 0, "legloc" : 2,
                         "leglims" : [], "labels" : [], "title" : "", "includeFit" : False,
                         "xlims" : None, "ylims" : None,
                         "includeNormErr" : False, "title" : "", "colors" : [],
                         "normerrlims" : None, "styles" : [3004, 3005, 3006, 3002, 3003, 3944],
                         "mplcolors" : [], "mpllabels" : [], "mplxlabel" : "", "mplylabel" : "",
                         'mpltitle' : "", 'normalised' : False,
                         "location" : "/Users/sfarry/WTau/Graphics/PDF", "extraObjs" : [],
                         "drawOpts" : 'e5'
                         }

   def drawROOT(self, filename = None):
      plots  = self.plots
      
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
      colors = []

      if len(self.properties["colors"]) > 0:
         colors = self.properties["colors"]
      else:
         colors = [hist.GetLineColor() for hist in plots]

      styles = self.properties['styles']

      for hist, color,style in zip(plots, colors,styles):
         if self.properties['normalised']:
            hist.Scale(1/hist.Integral())
         SetROOTLineColor(hist,color)
         SetROOTFillColor(hist,color)
         hist.SetFillStyle(style)

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

      
      drawOpts = self.properties['drawOpts']
      if isinstance(drawOpts,str):
          drawOpts = [drawOpts for p in range(len(plots))]

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
      if len(title) > 0:
          leg.AddEntry(dummyhist, title)
      for hist, label,drawOpt in zip(plots, labels,drawOpts):
          if ('h' in drawOpt):
              leg.AddEntry(hist, label, "l")
          elif ('e' in drawOpt):
              leg.AddEntry(hist, label, 'pl')
          else:
              leg.AddEntry(hist,label)

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

      for hist,drawOpt in zip(plots,drawOpts):
         if i == 0:
            if ('TGraph' in hist.ClassName()):
                if xlims is not None: hist.GetXaxis().SetRangeUser(xlims[0],xlims[1])
                if ylims is not None: hist.GetYaxis().SetRangeUser(ylims[0],ylims[1])
                hist.GetXaxis().SetTitle(xlabel)
                hist.GetYaxis().SetTitle(ylabel)
                #hist.SetFillStyle(3002)
                hist.Draw("A"+drawOpt)
            else:
                hist.SetXTitle(xlabel)
                hist.SetYTitle(ylabel)
                hist.Draw(drawOpt)
         else:
            if ('TGraph' in hist.ClassName()):
                #hist.SetFillStyle(3002)
                hist.Draw(drawOpt+"same")
            else:
                hist.Draw(drawOpt+"same")
         i = i + 1
      leg.Draw("same")
      for obj in self.properties["extraObjs"]: obj.Draw("same")

      if includeNormErr: 
         lowerPad.cd()
         plots[0].GetXaxis().SetLabelSize(0);
         i = 0
         normerrplots = [makeNormErr(plot) for plot in plots]
         normerrplots[0].GetXaxis().SetLimits(2,4.5)
         for errplot,color,style,drawOpt in zip(normerrplots,colors,styles,drawOpts):
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
               #errplot.Draw("e1")
               errplot.Draw(drawOpt)
            else:
               #errplot.Draw("e1same")
               errplot.Draw(drawOpt+"same")
            i = i + 1
            #plots.GetXaxis().SetLabelSize(0);
      #if includeNormErr:
      #else:
      if filename is None:
         filename = self.properties["filename"]
      location = self.properties["location"]

      if includeNormErr:
         lowerPad.Update()
         upperPad.Update()
         c.Update()

      #plot.Delete()
      #fit.Delete()
      #stack.Delete()
      c.Print(location + '/' + filename) ## Include this line to save image
      return c
   
    
   def setProps(self, props):
      for key in props:
         self.properties[key] = props[key]


class StackPlot:
   def __init__(self, h, stack, fit = None):
      self.hist       = h
      self.stack      = stack
      self.fit        = fit
      self.properties = {"xlabel" : "", "ylabel" : "", "filename" : "",
                         "logscale": False, "legwidth" : 0, "legheight" : 0, "legloc" : 2,
                         "leglims" : [], "labels" : [], "title" : "", "includeFit" : False,
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

      leg = ROOT.TLegend(leglims[0],leglims[1],leglims[2],leglims[3])
      leg.SetFillColor(10)
      leg.AddEntry(dummyhist, title)
      leg.AddEntry(plot, "2012 Data")
      for hist, label in zip(histlist, labels):
         leg.AddEntry(hist, label, 'fl')

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
      c.Print(location + '/' + filename) ## Include this line to save image
      return c

   def setProps(self, props):
      for key in props:
         self.properties[key] = props[key]
