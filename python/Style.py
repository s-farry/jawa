import ROOT


black=1
red=2
green=3
blue=4
yellow=5 
magenta=6
cyan=7
purple=9

def SetLHCbStyle():
    lhcbStyle= ROOT.TStyle("lhcbStyle","Standard LHCb plots style")

    lhcbFont = 132
    lhcbWidth = 2

    lhcbStyle.SetFrameBorderMode(0)
    lhcbStyle.SetCanvasBorderMode(0)
    lhcbStyle.SetPadBorderMode(0)
    lhcbStyle.SetPadColor(0)
    lhcbStyle.SetCanvasColor(0)
    lhcbStyle.SetStatColor(0)
    lhcbStyle.SetPalette(1)
    
    lhcbStyle.SetPaperSize(20,26)
    lhcbStyle.SetPadTopMargin(0.05)
    lhcbStyle.SetPadRightMargin(0.05) # increase for colz plots
    lhcbStyle.SetPadBottomMargin(0.16)
    lhcbStyle.SetPadLeftMargin(0.14)
    
    lhcbStyle.SetTextFont(lhcbFont)
    lhcbStyle.SetTextSize(0.08)
    lhcbStyle.SetLabelFont(lhcbFont,"x")
    lhcbStyle.SetLabelFont(lhcbFont,"y")
    lhcbStyle.SetLabelFont(lhcbFont,"z")
    lhcbStyle.SetTitleFont(lhcbFont,"x")
    lhcbStyle.SetTitleFont(lhcbFont,"y")
    lhcbStyle.SetTitleFont(lhcbFont,"z")
    lhcbStyle.SetLabelSize(0.05,"x")
    lhcbStyle.SetLabelSize(0.05,"y")
    lhcbStyle.SetLabelSize(0.05,"z")
    #lhcbStyle.SetTitleFont(lhcbFont)
    lhcbStyle.SetTitleSize(0.06,"x")
    lhcbStyle.SetTitleSize(0.06,"y")
    lhcbStyle.SetTitleSize(0.06,"z")
    
    lhcbStyle.SetLineWidth(lhcbWidth)
    lhcbStyle.SetFrameLineWidth(lhcbWidth)
    lhcbStyle.SetHistLineWidth((lhcbWidth))
    lhcbStyle.SetFuncWidth(lhcbWidth)
    lhcbStyle.SetGridWidth(lhcbWidth/2)
    lhcbStyle.SetLineStyleString(2,"[12 12]") # postscript dashes
    lhcbStyle.SetMarkerStyle(0)
    lhcbStyle.SetMarkerSize(1)
    
    lhcbStyle.SetLabelOffset(0.02)
    
    
    lhcbStyle.SetOptStat(0)  
    #lhcbStyle.SetOptStat("emr")  # show only nent -e , mean - m , rms -r
    
    lhcbStyle.SetStatFormat("6.3g") # specified as c printf options
    lhcbStyle.SetOptTitle(0)
    lhcbStyle.SetOptFit(0)
    
    lhcbStyle.SetStatBorderSize(0)
    lhcbStyle.SetStatFont(lhcbFont)
    lhcbStyle.SetStatFontSize(0.05)
    lhcbStyle.SetStatX(0.9)
    lhcbStyle.SetStatY(0.9)
    lhcbStyle.SetStatW(0.25)
    lhcbStyle.SetStatH(0.15)
    # put tick marks on top and RHS of plots
    lhcbStyle.SetPadTickX(1)
    lhcbStyle.SetPadTickY(1)

    # histogram divisions: only 5 in x to avoid label overlaps
    lhcbStyle.SetNdivisions(1005,"x")
    lhcbStyle.SetNdivisions(510,"y")
    #Set legends
    lhcbStyle.SetLegendBorderSize(0)
    lhcbStyle.SetLegendFillColor(10)
    lhcbStyle.SetLegendFont(lhcbFont)

    #lhcbStyle.SetPadGridX(True)
    #lhcbStyle.SetPadGridY(True)
    lhcbStyle.SetGridStyle(3)
    #define style for text
    lhcbLabel = ROOT.TText()
    lhcbLabel.SetTextFont(lhcbFont)
    lhcbLabel.SetTextColor(1)
    lhcbLabel.SetTextSize(0.04)
    lhcbLabel.SetTextAlign(12)    
    # define style of latex text
    lhcbLatex = ROOT.TLatex()
    lhcbLatex.SetTextFont(lhcbFont)
    lhcbLatex.SetTextColor(1)
    lhcbLatex.SetTextSize(0.04)
    lhcbLatex.SetTextAlign(12)

    # set this style
    ROOT.gROOT.SetStyle("lhcbStyle")
    ROOT.gROOT.ForceStyle()


#
# based on a style file from BaBar
#

#..BABAR style from RooLogon.C in workdir

def SetAtlasStyle():
    atlasStyle= ROOT.TStyle("ATLAS","Atlas style");

    # use plain black on white colors
    icol=0
    atlasStyle.SetFrameBorderMode(icol)
    atlasStyle.SetCanvasBorderMode(icol)
    atlasStyle.SetPadBorderMode(icol)
    atlasStyle.SetPadColor(icol)
    atlasStyle.SetCanvasColor(icol)
    atlasStyle.SetStatColor(icol)
    #atlasStyle.SetFillColor(icol)

    # set the paper & margin sizes
    atlasStyle.SetPaperSize(20,26)
    atlasStyle.SetPadTopMargin(0.05)
    atlasStyle.SetPadRightMargin(0.05)
    atlasStyle.SetPadBottomMargin(0.16)
    atlasStyle.SetPadLeftMargin(0.12)

    # use large fonts
    #Int_t font=72
    font=42
    tsize=0.05
    atlasStyle.SetTextFont(font)


    atlasStyle.SetTextSize(tsize)
    atlasStyle.SetLabelFont(font,"x")
    atlasStyle.SetTitleFont(font,"x")
    atlasStyle.SetLabelFont(font,"y")
    atlasStyle.SetTitleFont(font,"y")
    atlasStyle.SetLabelFont(font,"z")
    atlasStyle.SetTitleFont(font,"z")

    atlasStyle.SetLabelSize(tsize,"x")
    atlasStyle.SetTitleSize(tsize,"x")
    atlasStyle.SetLabelSize(tsize,"y")
    atlasStyle.SetTitleSize(tsize,"y")
    atlasStyle.SetLabelSize(tsize,"z")
    atlasStyle.SetTitleSize(tsize,"z")
    
    atlasStyle.SetLabelOffset(tsize/4,"x")
    atlasStyle.SetLabelOffset(tsize/4,"y")
    
    
    #use bold lines and markers
    atlasStyle.SetMarkerStyle(20)
    atlasStyle.SetMarkerSize(1.2)
    atlasStyle.SetHistLineWidth(2)
    atlasStyle.SetLineStyleString(2,"[12 12]") # postscript dashes

    #get rid of X error bars and y error bar caps
    #atlasStyle.SetErrorX(0.001)

    #do not display any of the standard histogram decorations
    atlasStyle.SetOptTitle(0)
    #atlasStyle.SetOptStat(1111)
    atlasStyle.SetOptStat(0)
    #atlasStyle.SetOptFit(1111)
    atlasStyle.SetOptFit(0)

    # put tick marks on top and RHS of plots
    atlasStyle.SetPadTickX(1)
    atlasStyle.SetPadTickY(1)

    ROOT.gROOT.SetStyle("Plain")
    
    #gStyle.SetPadTickX(1)
    #gStyle.SetPadTickY(1)

    ROOT.gROOT.SetStyle("ATLAS")
    ROOT.gROOT.ForceStyle()

