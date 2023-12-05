#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <memory>

#include "formatting.hpp"

void preliminary_formatting(
    TH1D* hist,
    int style_idx,
    const char* title,
    double low,
    double high);

TH1D* generate_normalized(
    TH1D* standard,
    TH1D* to_be_scaled,
    double low,
    double high);

double integrate_hist(TH1D* hist, double low, double high);

struct TreeBranches {
    double e_rel {};
    double inner_angle {};
};

int main() {
    auto filter_name = "filt_7a_std";
    //auto filter_name = "filt_3a_min";
    //auto mixer_name = "mixed_rndm_phi";
    auto mixer_name = "mixed_mix_7a_unique";
    //auto mixer_name = "mixed_mix_3a";
    //auto mixer_name = "mixed_rndm_phi_faust";
    auto system = std::string("7a");
    double scaling_region_low = 20.;
    double scaling_region_up = 100.;
    double range_low = -2.;
    double range_up = 140.;

    auto file_name = std::string("relative_energy_" + system + ".root");

    auto r_path = std::filesystem::path("K:\\tamu_data\\exp\\si28_c_35\\anl")
    / filter_name / "real" / file_name;

    auto m_path = std::filesystem::path("K:\\tamu_data\\exp\\si28_c_35\\anl")
    / filter_name / mixer_name / file_name;

    //auto r_path = std::filesystem::path("K:\\tamu_data\\exp\\c12_si_35\\anl")
        /// filter_name / "real" / file_name;

    //auto m_path = std::filesystem::path("K:\\tamu_data\\exp\\c12_si_35\\anl")
        /// filter_name / mixer_name / file_name;

    auto r_file = std::make_shared<TFile>(r_path.string().c_str());
    auto m_file = std::make_shared<TFile>(m_path.string().c_str(), "read");

    auto e_rel_real = r_file->Get<TH1D>("e_rel_hist");
    auto e_rel_mixed = m_file->Get<TH1D>("e_rel_hist");

    preliminary_formatting(e_rel_real, 1, "Real", range_low, range_up);
    preliminary_formatting(e_rel_mixed, 2, "Real", range_low, range_up);

    TH1D* e_rel_mixed_scaled = generate_normalized(
        e_rel_real,
        e_rel_mixed,
        scaling_region_low,
        scaling_region_up);

    auto residuals = (TH1D*)e_rel_real->Clone();
    residuals->Add(e_rel_mixed_scaled, -1.);
    residuals->GetYaxis()->SetTitle("Residuals");

    auto std_residual = std::make_shared<TGraph>();
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    for (int bin = 0; bin < e_rel_real->GetNbinsX(); ++bin) {
        double bin_center = e_rel_real->GetBinCenter(bin);
        double bin_content = e_rel_real->GetBinContent(bin);
        double background = e_rel_mixed_scaled->GetBinContent(bin);

        double error = std::sqrt(bin_content + background);

        double diff = bin_content - background;
        double std_res = diff / error;

        if (std_res < min) {
            min = std_res;
        }
        if (std_res > max) {
            max = std_res;
        }

        std_residual->SetPoint(bin, bin_center, std_res);
    }

    double range = max - min;

    auto canvas_2panel =
        std::make_unique<TCanvas>("canvas", "canvas", 800, 600);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetOptTitle(0);
    gStyle->SetFrameLineWidth(2);

    canvas_2panel->Divide(1, 2, 0, 0);
    canvas_2panel->cd(1);
    gPad->SetLeftMargin(0.15f);
    gPad->SetRightMargin(0.02f);
    gPad->SetTicks(1, 1);
    residuals->SetLineColor(kBlack);
    auto baseline = std::make_shared<TLine>(range_low, 0., range_up, 0.);
    baseline->SetLineColor(kRed);
    baseline->SetLineWidth(4);
    residuals->Draw("hist");
    baseline->Draw("same");
    canvas_2panel->cd(2);
    gPad->SetBottomMargin(0.15f);
    gPad->SetRightMargin(0.02f);
    gPad->SetLeftMargin(0.15f);
    gPad->SetTicks(1, 1);
    auto h_blank = std::make_shared<TH1D>(
        "blank",
        "blank",
        e_rel_real->GetNbinsX(),
        range_low,
        range_up);
    preliminary_formatting(h_blank.get(), 0, "", range_low, range_up);
    h_blank->SetStats(false);
    h_blank->SetTitle("#sigma");
    h_blank->GetXaxis()->SetTitle("E_{rel} [MeV]");
    h_blank->GetYaxis()->SetTitle("#sigma");
    h_blank->GetYaxis()->SetRangeUser(min - range * .1, max + range * .1);
    h_blank->Draw();
    std_residual->SetMarkerStyle(20);
    std_residual->SetMarkerSize(.5);
    std_residual->Draw("P");
    baseline->Draw("same");

    canvas_2panel->SaveAs("fig/real_w_mixed_hist_residuals.png");

    auto canvas_raw = std::make_unique<TCanvas>("canvas", "canvas", 800, 600);
    canvas_raw->SetBottomMargin(0.15f);
    canvas_raw->SetLeftMargin(0.15f);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetOptTitle(0);
    gStyle->SetFrameLineWidth(2);
    e_rel_real->Draw("E");
    e_rel_real->GetXaxis()->SetTitle("E_{rel} [MeV]");
    canvas_raw->SaveAs("fig/real_w_mixed_hist_raw.png");
    e_rel_mixed_scaled->Draw("E SAME");
    canvas_raw->SaveAs("fig/real_w_mixed_hist_w_mixed.png");



    for (int bin = 1; bin <= e_rel_real->GetNbinsX(); ++bin){
        e_rel_real->SetBinContent(bin, sqrt(e_rel_real->GetBinContent(bin)));
    }
    e_rel_real->GetYaxis()->SetTitle("#sqrt{Yield}"); 
    e_rel_real->SetMarkerSize(.5);
    e_rel_real->SetMarkerStyle(20);

    e_rel_real->Draw("P");
    canvas_raw->SaveAs("fig/real_sqrt.png");



}

void preliminary_formatting(
    TH1D* hist,
    int style_idx,
    const char* title,
    double low,
    double high) {
    hist->SetLineColor(style_idx);
    hist->SetMarkerColor(style_idx);
    hist->SetLineWidth(2);
    hist->SetMarkerSize(2);

    for (int i = 0; i < 4; ++i) {
        hist->RebinX(2);
    }

    hist->GetXaxis()->SetNdivisions(108);
    hist->GetYaxis()->SetNdivisions(105);

    hist->SetTitle(title);

    hist->GetXaxis()->SetLabelSize(0.04f);
    hist->GetYaxis()->SetLabelSize(0.04f);

    hist->GetYaxis()->SetTitleOffset(0.95f);
    hist->GetXaxis()->SetTitleOffset(0.95f);
    hist->GetXaxis()->SetTitleSize(0.06f);
    hist->GetYaxis()->SetTitleSize(0.06f);
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    //hist->GetYaxis()->ChangeLabel(1, -1, 0, -1, -1, -1, "");
    //hist->GetXaxis()->ChangeLabel(-1, -1, 0, -1, -1, -1, "");
    hist->GetXaxis()->SetLabelFont(112);
    hist->GetYaxis()->SetLabelFont(112);
    hist->GetXaxis()->SetTitleFont(112);
    hist->GetYaxis()->SetTitleFont(112);
    //hist->GetYaxis()->ChangeLabel(-1, -1, 0, -1, -1, -1, "");

    hist->GetXaxis()->SetRangeUser(low, high);
}

TH1D* generate_normalized(
    TH1D* standard,
    TH1D* to_be_scaled,
    double low,
    double high) {
    auto result = (TH1D*)to_be_scaled->Clone();

    double norm_standard = integrate_hist(standard, low, high);
    double norm = integrate_hist(result, low, high);

    result->Scale(norm_standard / norm);

    return result;
}

double integrate_hist(TH1D* hist, double low, double high) {
    int bin_low = hist->GetXaxis()->FindBin(low);
    int bin_high = hist->GetXaxis()->FindBin(high);

    return hist->Integral(bin_low, bin_high);
}
