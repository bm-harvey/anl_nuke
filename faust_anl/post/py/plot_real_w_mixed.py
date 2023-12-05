from ROOT import TFile
# import ROOT as rt
import os


if __name__ == "__main__":
    mixer_name = "7a_erel"
    real_data_path = os.path("D:\\tamu_data\\exp\\si28_c_35\\anl").join(
        mixer_name).join("real\\erel.root")
    mixed_data_path = os.path("D:\\tamu_data\\exp\\si28_c_35\\anl").join(
        mixer_name).join("mixed\\erel.root")

    real_file = TFile(real_data_path)
    mixed_file = TFile(real_data_path)

    real_file.ls()
    mixed_file.ls()
