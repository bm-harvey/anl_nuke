#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>

#include <algorithm>
#include <filesystem>
#include <memory>

int main() {
    auto be_0 = "7a_not_through_be8gs";
    auto be_1 = "7a_through_be8gs";
    auto be_2 = "7a_through_2be8gs";
    auto be_all = "7a_std";

    auto path_0 = std::filesystem::path("D:\\tamu_data\\exp\\si28_c_35\\anl")
        / be_0 / "real\\relative_energy_7a.root";

    auto path_1 = std::filesystem::path("D:\\tamu_data\\exp\\si28_c_35\\anl")
        / be_1 / "real\\relative_energy_7a.root";

    auto path_2 = std::filesystem::path("D:\\tamu_data\\exp\\si28_c_35\\anl")
        / be_2 / "real\\relative_energy_7a.root";

    auto path_all = std::filesystem::path("D:\\tamu_data\\exp\\si28_c_35\\anl")
        / be_all / "real\\relative_energy_7a.root";

    auto file_0 = std::make_unique<TFile>(path_0.string().c_str());
    auto file_1 = std::make_unique<TFile>(path_1.string().c_str());
    auto file_2 = std::make_unique<TFile>(path_2.string().c_str());
    auto file_all = std::make_unique<TFile>(path_all.string().c_str());

    auto hist_0 = file_0->Get<TH1D>("source_particles_rel_ke_hist");
    auto hist_1 = file_1->Get<TH1D>("source_particles_rel_ke_hist");
    auto hist_2 = file_2->Get<TH1D>("source_particles_rel_ke_hist");
    auto hist_all = file_all->Get<TH1D>("source_particles_rel_ke_hist");

    hist_2->Rebin(4);
    hist_1->Rebin(2);

    hist_0->Scale(1. / hist_0->Integral());
    hist_all->Scale(1. / hist_all->Integral());
    hist_1->Scale(.5 / hist_1->Integral());
    hist_2->Scale(.25 / hist_2->Integral());

    auto canvas = std::make_shared<TCanvas>("canvas", "", 800, 800);
    gStyle->SetOptStat(0);

    hist_0->Draw("e");
    hist_1->Draw("e, same");
    hist_2->Draw("e, same");
    hist_all->Draw("e, same");

    auto legend = std::make_shared<TLegend>(.5, .7, .9, .85);
    legend->AddEntry(hist_0, "0 ^{8}Be g.s.");
    legend->AddEntry(hist_1, "1 ^{8}Be g.s.");
    legend->AddEntry(hist_2, "2 ^{8}Be g.s.");
    legend->AddEntry(hist_all, "all");
    legend->Draw("same");

    hist_0->GetYaxis()->SetRangeUser(0, .008);
    hist_0->GetXaxis()->SetRangeUser(0, 40);
    hist_0->SetLineColor(kBlack);

    hist_0->SetTitle("Source frame kinetic energy of intermediate particles");

    hist_1->SetLineColor(kGreen);

    hist_2->SetLineColor(kBlue);

    hist_all->SetLineColor(kRed);

    canvas->Print("fig/source_frame_cm_ke.png");
}
