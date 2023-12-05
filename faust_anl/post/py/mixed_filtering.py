import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import json
from flatten_json import flatten
import ROOT as rt
if __name__ == "__main__":
    real_file_path = "D:tamu_data\\exp\\si28_c_35\\anl\\2a_rndm_phi\\real\\erel.json"
    mixed_file_path = "D:tamu_data\\exp\\si28_c_35\\anl\\2a_rndm_phi\\mixed\\erel.json"

    with open(real_file_path) as real_file:
        with open(mixed_file_path) as mixed_file:
            real_data = json.load(real_file)
            mixed_data = json.load(mixed_file)

            real_data = pd.DataFrame((flatten(record, '_')
                                     for record in real_data))
            mixed_data = pd.DataFrame((flatten(record, '_')
                                       for record in mixed_data))

            print(real_data)
            real_data.to_csv("real.csv")
            print(mixed_data)
            mixed_data.to_csv("mixed.csv")

            """ E rel. plotting"""
            h_real_erel = rt.TH1D(
                "h_real_erel", ";E_{rel} (MeV);", 1_000, -5, 45)
            h_mixed_erel = rt.TH1D(
                "h_mixed_erel", ";E_{rel} (MeV);", 1_000, -5, 45)

            for idx in range(real_data.shape[0]):
                h_real_erel.Fill(
                    real_data.iloc[idx]["e_rel"]
                )

            for idx in range(mixed_data.shape[0]):
                h_mixed_erel.Fill(
                    mixed_data.iloc[idx]["e_rel"]
                )

            canvas = rt.TCanvas()

            h_real_erel.SetLineColor(rt.kBlack)
            h_mixed_erel.SetLineColor(rt.kRed)

            h_real_erel.Scale(1/h_real_erel.Integral())
            h_mixed_erel.Scale(.95/h_mixed_erel.Integral())

            h_real_erel.Draw("")
            h_mixed_erel.Draw("same")
            canvas.Print("e_rel.png")

            """ Pair plotting """
            # ax = sns.pairplot(real_data.iloc[:1_000, :], kind="hist")
            # plt.savefig("real.png")

            # ax = sns.pairplot(
            # mixed_data.iloc[:1_000, :], kind="hist")
            # plt.savefig("mixed.png")
