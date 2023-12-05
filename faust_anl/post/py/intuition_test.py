import random
import numpy as np
from scipy import stats
from ROOT import TCanvas, TGraph, TFile, TH1D, TTree
from ROOT import gStyle, gROOT
gROOT.SetBatch(1)
print("done with imports")

hist_counter = 0
bins = 300


def main():
    file_7 = TFile(
        'D:\\tamu_data\\exp\\si28_c_35\\anl\\filt_7a_std\\real\\relative_energy_7a.root', 'READ')
    tree_7 = file_7.Get('T')

    file_3 = TFile(
        'D:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_std\\real\\relative_energy_3a.root', 'READ')
    tree_3 = file_3.Get('T')

    file_4 = TFile(
        'D:\\tamu_data\\exp\\si28_c_35\\anl\\filt_4a_std\\real\\relative_energy_4a.root', 'READ')
    tree_4 = file_4.Get('T')

    data_dict = extract_data(tree_7, tree_3, tree_4)

    keys = [key for key, _ in data_dict.items()]
    random.shuffle(keys)

    # ouput file to write the answer key to
    answer_key = open("fig/answer_key.txt", "w")

    for key in keys:
        draw_hists(key, data_dict[key], answer_key)


def extract_data(tree_7a: TTree, tree_3a: TTree, tree_4a: TTree) -> dict[str, list[float]]:
    """
    Extracts data from the input tree and creates a dictionary
    which maps a name to a data series. The extracted data is a combination
    of data which we all agree is smooth, which some think have clear
    structure and some are just sampled from a smooth distribution.
    """

    data = {}

    data["raw_3a"] = []

    data["raw_4a"] = []
    data["raw_4a_thru_2be"] = []

    data["raw_7a"] = []
    data["7a_thru_1be"] = []
    data["7a_thru_2be"] = []
    data["7a_thru_hoyle"] = []
    data["7a_thru_2be_r"] = []
    data["7a_thru_c_3m"] = []
    data["7a_thru_c_3m_ds_cut"] = []

    data["beta_4_9"] = stats.beta.rvs(4, 9, size=3142)
    data["beta_5_8"] = stats.beta.rvs(5, 8, size=1618)
    data["beta_5_11"] = stats.beta.rvs(5, 11, size=2718)
    data["beta_3_9"] = stats.beta.rvs(3, 9, size=2718)
    data["beta_7_8"] = stats.beta.rvs(7, 8, size=2718)

    for entry in tree_3a:
        e_rel = entry.e_rel
        data["raw_3a"].append(e_rel)

        if len(data["raw_3a"]) == 1_000_000:
            break

    for entry in tree_4a:
        e_rel = entry.e_rel
        data["raw_4a"].append(e_rel)

        if entry.be8gs_count == 2:
            data["raw_4a_thru_2be"].append(e_rel)

        if len(data["raw_4a"]) == 1000000:
            break

    for entry in tree_7a:
        e_rel = entry.e_rel
        sphericity = entry.sphericity
        coplanarity = entry.coplanarity

        # distance to ds line is a little involved
        x = .25 * (3 - np.sqrt(3) * coplanarity + sphericity)
        y = np.sqrt(3) * (1 - x)
        d_ds = np.sqrt((x - sphericity)**2 + (y - coplanarity) ** 2)

        data["raw_7a"].append(e_rel)

        if entry.be8gs_count == 2:
            data["7a_thru_2be"].append(e_rel)
            if sphericity < 0.41 and coplanarity < 0.18:
                data["7a_thru_2be_r"].append(e_rel)

        if entry.be8gs_count == 1:
            data["7a_thru_1be"].append(e_rel)

        if entry.hoyle_count >= 1:
            data["7a_thru_hoyle"].append(e_rel)

        if entry.c12_3m_count >= 1:
            data["7a_thru_c_3m"].append(e_rel)
            if 0.5 < d_ds and d_ds < 0.6:
                data["7a_thru_c_3m_ds_cut"].append(e_rel)

    print(len(data["7a_thru_c_3m"]))
    print(len(data["7a_thru_c_3m_ds_cut"]))

    data["7a_thru_1be_scaled"] = data["7a_thru_1be"][0:1000]

    return data


def draw_hists(name: str, data: list[float], answer_key) -> int:
    # create histogram of the actual data
    hist, min_val, max_val = make_hist(data)
    integral = len(data)
    print(hist.GetMean())
    bin_width = hist.GetBinWidth(1)

    # fit the histogram with a beta distribution
    #   - no physical reason but works really really well
    a, b, loc, scale = stats.beta.fit(data)

    # create a list of histograms. The zeroeth is the actual data.
    # the remaining `num_sampled` are sampled from the fit above
    num_sampled = 9
    hists = [hist]

    for _ in range(num_sampled):
        hists.append(make_beta_sampled_hist(a, b, loc, scale, len(data), min_val, max_val))

    # randomize the order in which histograms are accessed for plotting.
    hist_idx = np.arange(0, len(hists), 1)
    random.shuffle(hist_idx)

    # prepare a canvas to draw on
    can = TCanvas("c", "c", 400, 800)
    gStyle.SetOptStat(0)
    vertical_shift = hists[0].GetMaximum() * .75

    # plot all of the histograms
    for count, idx in enumerate(hist_idx):
        hist = hists[idx]
        hist.SetTitle("")
        if idx == 0:
            answer_key.write(f"\n{name}\n")
            answer_key.write(f"Answer: {count}\n")
            answer_key.write(f"Count: {hist.Integral()}\n")

        offset = (count + 1) * vertical_shift
        for bin in range(1, hist.GetNbinsX() + 1):
            err = np.sqrt(hist.GetBinContent(bin))
            hist.SetBinContent(bin, hist.GetBinContent(bin) + offset)
            hist.SetBinError(bin, err)

        hist.SetLineColor(int(count) % 4 + 1)
        hist.Draw("e, same")

        hist.SetNdivisions(000, "X")
        hist.GetXaxis().SetTitle("Relative Energy (Arb. Units)")

    hists[hist_idx[0]].GetYaxis().SetRangeUser(
        0., hists[hist_idx[-1]].GetMaximum() * 1.02)

    # Plot the fit so it can be seen is being sampled
    graph = TGraph()
    for x in np.arange(min_val, max_val, .05):
        graph.AddPoint(x, stats.beta.pdf(x, a, b, loc, scale) * bin_width * integral)
    graph.SetLineColor(435)
    graph.SetLineWidth(2)
    graph.Draw("l, same")

    # output
    can.Print(f"fig/{name}.png")
    can.Print(f"fig/{name}.pdf")


def make_kde_graph(hist: TH1D, kde) -> TGraph:
    bin_width = hist.GetBinWidth(1)
    integral = hist.Integral()
    kde_graph = TGraph()
    for x in np.linspace(hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax(), 1000):
        kde_graph.SetPoint(kde_graph.GetN(), x, kde.pdf(x)[
                           0] * integral * bin_width)
    return kde_graph


def get_bins(size: int) -> int:
    if size > 100_000:
        return 500
    elif size > 1000:
        return 200
    else:
        return 150


def make_hist(data: list[float]) -> (TH1D, float, float):
    min_val = min(data)
    max_val = max(data)
    buffer = 0  # .02 * (max_val - min_val)
    hist = TH1D(hist_name(), 'hist', get_bins(len(data)), min_val - buffer, max_val + buffer)
    max_val = max(data)
    min_val = min(data)

    for entry in data:
        hist.Fill(entry)

    return (hist, min_val, max_val)


def make_beta_sampled_hist(a: float, b: float, location: float, scale: float, count: int, min_val: float, max_val: float) -> TH1D:
    data = stats.beta.rvs(a, b, location, scale, size=count)
    hist = TH1D(hist_name(), 'hist', get_bins(len(data)), min_val, max_val)
    for x in data:
        hist.Fill(x)
    return hist


def make_sampled_hist(kde, num: int, min_val: int, max_val: int) -> TH1D:
    data = kde.resample()[0]
    hist = TH1D(hist_name(), 'hist', get_bins(len(data)), min_val, max_val)
    for val in data:
        hist.Fill(val)
    return hist


def make_kde(data: list[float]):
    kde = stats.gaussian_kde(data)
    return kde


def hist_name() -> str:
    global hist_counter
    hist_counter += 1
    return f"hist_{hist_counter}"


if __name__ == '__main__':
    main()
