from analysis import *
from ROOT import *

dpm = "root://hepgrid11.ph.liv.ac.uk///dpm/ph.liv.ac.uk/home/lhcb/electroweak/"

tune = Expr("muminus_ipyubs")
eta   = Var("ETA", "muminus_ETA", 10, 2, 4.5)

dataf = TFile.Open(dpm + "ZMuMu.2012.root")
mcf   = TFile.Open(dpm + "ZMuMu.Zg_mumu.MC2012.root")

datat = Tree("data", dataf.Get("ZMuMu/DecayTree"))
mct   = Tree("mc"  , mcf.Get("TupleOS/DecayTree"))
zmumu = Expr("boson_M > 60000 && boson_M < 120000 && muplus_ETA > 2 && muplus_ETA < 4.5 && muminus_ETA > 2 && muminus_ETA < 4.5 && muminus_TRACK_PCHI2 > 0.01 && muplus_TRACK_PCHI2 > 0.01 && sqrt(muminus_PERR2)/muminus_P < 0.1 && sqrt(muplus_PERR2)/muplus_P < 0.1")

ipx = Expr("muminus_ipxubs")
ipy = Expr("muminus_ipyubs")
ipz = Expr("muminus_ipzubs")
iptx = Expr("muminus_iptxubs")
ipty = Expr("muminus_iptyubs")


tuner_ipx = Tune("ipx_tune", datat, mct, ipx, eta, zmumu.GetExpr())

a = [[25, 15, -5], [15, 18, 0], [-5, 0, 11] ]




cov = datat.getCorrelationMatrix([ipx, ipy, ipz, iptx, ipty], zmumu)

tuner_ipx = Tune("ipx_tune", datat, mct, ipx, eta, zmumu.GetExpr())
tuner_ipx.tune()
tuner_ipx.SaveToFile()

tuner_ipy = Tune("ipy_tune", datat, mct, ipy, eta, zmumu.GetExpr())
tuner_ipy.tune()
tuner_ipy.SaveToFile()

tuner_ipz = Tune("ipz_tune", datat, mct, ipz, eta, zmumu.GetExpr())
tuner_ipz.tune()
tuner_ipz.SaveToFile()

tuner_iptx = Tune("iptx_tune", datat, mct, iptx, eta, zmumu.GetExpr())
tuner_iptx.tune()
tuner_iptx.SaveToFile()

tuner_ipty = Tune("ipty_tune", datat, mct, ipty, eta, zmumu.GetExpr())
tuner_ipty.tune()
tuner_ipty.SaveToFile()

