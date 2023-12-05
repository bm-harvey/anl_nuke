from scipy import stats
import numpy as numpy
import matplotlib.pyplot as plt
import seaborn as sns

from ROOT import TH1D, TCanvas, gStyle

data: list[float] = stats.norm.rvs(0, 1, size=2000)


hist = TH1D("hist", "Binned", 100, -3, 3)
gStyle.SetOptStat(0)
hist.SetLineColor(1)
hist.SetLineWidth(2)

for x in data:
    hist.Fill(x)

hist_3maw = hist.Clone()
for bin in range(2, hist_3maw.GetNbinsX() - 1):
    new_value = (hist.GetBinContent(bin - 1) + hist.GetBinContent(bin) + hist.GetBinContent(bin + 1))/3
    hist_3maw.SetBinContent(bin, new_value)


canvas = TCanvas("c", "c", 800, 800)

hist.Draw("E")
canvas.Print("fig/gauss.png")

hist.Smooth()
hist.Draw("E")
canvas.Print("fig/gauss_smoothed.png")

hist_3maw.Draw("E")
canvas.Print("fig/gauss_3maw.png")
