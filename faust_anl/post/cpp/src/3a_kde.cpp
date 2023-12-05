

#include <Math/PdfFuncMathCore.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TTree.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#include "formatting.hpp"

TGraphErrors*
generate_correlation_func(TH1D* real, TH1D* mix_scaled, double offset = 0) {
    auto* corr = new TGraphErrors();
    for (int bin = 1; bin <= real->GetNbinsX(); ++bin) {
        auto e_rel = real->GetBinCenter(bin);
        auto real_yield = real->GetBinContent(bin);
        auto mixed_yield = mix_scaled->GetBinContent(bin);
        auto corr_val = real_yield / mixed_yield;
        auto corr_err = corr_val
            * sqrt(pow(real->GetBinError(bin) / real_yield, 2)
                   + pow(mix_scaled->GetBinError(bin) / mixed_yield, 2));

        if (real_yield > 2 && mixed_yield > 2) {
            corr->SetPoint(bin - 1, e_rel, corr_val + offset);
            corr->SetPointError(bin - 1, 0, corr_err);
        }
    }
    return corr;
}

double kernel(double x, double h) {
    return exp(-x * x / 2 / h / h) / sqrt(2 * M_PI) / h;
}

TGraph* generate_kde(
    TTree* tree,
    double bandwidth,
    int max_points = std::numeric_limits<int>::max()) {
    double e_rel = 0;
    auto* result = new TGraph();
    for (double x = 0; x < 100; x += 1) {
        result->AddPoint(x, 0);
    }

    tree->SetBranchAddress("e_rel", &e_rel);
    int num_entries = std::min((int)tree->GetEntries(), max_points);
    double norm = std::max(1., double(num_entries) / max_points);

    for (int i = 0; i < (int)num_entries; ++i) {
        tree->GetEntry(i);
        for (int j = 0; j < result->GetN(); ++j) {
            auto x = result->GetX()[j];
            auto y = result->GetY()[j];
            auto k = kernel(x - e_rel, bandwidth);
            result->SetPoint(j, x, (y + k * norm));
        }

        if (i % 10'000 == 0) {
            std::cout << i << std::endl;
        }
    }
    return result;
}

int main() {
    formatting::set_style();

    auto file_hhh_real = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\real\\relative_energy_3h.root");
    auto file_hha_real = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\real\\relative_energy_1a_2h.root");
    auto file_haa_real = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\real\\relative_energy_2a_1h.root");
    auto file_aaa_real = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\real\\relative_energy_3a.root");

    //auto file_hhh_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_rndm_phi_faust\\relative_energy_3h.root");
    //auto file_hha_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_rndm_phi_faust\\relative_energy_1a_2h.root");
    //auto file_haa_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_rndm_phi_faust\\relative_energy_2a_1h.root");
    //auto file_aaa_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_rndm_phi_faust\\relative_energy_3a.root");

    //auto file_hhh_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_rndm_phi\\relative_energy_3h.root");
    //auto file_hha_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_rndm_phi\\relative_energy_1a_2h.root");
    //auto file_haa_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_rndm_phi\\relative_energy_2a_1h.root");
    //auto file_aaa_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_rndm_phi\\relative_energy_3a.root");

    auto file_hhh_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_mix_3h_unique\\relative_energy_3h.root");
    auto file_hha_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_mix_1a_2h_unique\\relative_energy_1a_2h.root");
    auto file_haa_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_mix_2a_1h_unique\\relative_energy_2a_1h.root");
    auto file_aaa_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_mix_3a_unique\\relative_energy_3a.root");

    //auto file_hhh_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_mix_3h\\relative_energy_3h.root");
    //auto file_hha_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_mix_1a_2h\\relative_energy_1a_2h.root");
    //auto file_haa_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_mix_2a_1h\\relative_energy_2a_1h.root");
    //auto file_aaa_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_mix_3a\\relative_energy_3a.root");

    double scale_low = 10;
    double scale_high = 100;

    auto hhh_r_tree = file_hhh_real->Get<TTree>("T");
    auto hha_r_tree = file_hha_real->Get<TTree>("T");
    auto haa_r_tree = file_haa_real->Get<TTree>("T");
    auto aaa_r_tree = file_aaa_real->Get<TTree>("T");

    auto hhh_m_tree = file_hhh_mixed->Get<TTree>("T");
    auto hha_m_tree = file_hha_mixed->Get<TTree>("T");
    auto haa_m_tree = file_haa_mixed->Get<TTree>("T");
    auto aaa_m_tree = file_aaa_mixed->Get<TTree>("T");

    auto hist_hhh_r = file_hhh_real->Get<TH1D>("e_rel_hist");
    auto hist_hha_r = file_hha_real->Get<TH1D>("e_rel_hist");
    auto hist_haa_r = file_haa_real->Get<TH1D>("e_rel_hist");
    auto hist_aaa_r = file_aaa_real->Get<TH1D>("e_rel_hist");

    auto hist_hhh_m = file_hhh_mixed->Get<TH1D>("e_rel_hist");
    auto hist_hha_m = file_hha_mixed->Get<TH1D>("e_rel_hist");
    auto hist_haa_m = file_haa_mixed->Get<TH1D>("e_rel_hist");
    auto hist_aaa_m = file_aaa_mixed->Get<TH1D>("e_rel_hist");

    auto hist_aaa_m_norm = formatting::generate_normalized(
        hist_aaa_r,
        hist_aaa_m,
        scale_low,
        scale_high);
    auto hist_haa_m_norm = formatting::generate_normalized(
        hist_haa_r,
        hist_haa_m,
        scale_low,
        scale_high);
    auto hist_hha_m_norm = formatting::generate_normalized(
        hist_hha_r,
        hist_hha_m,
        scale_low,
        scale_high);
    auto hist_hhh_m_norm = formatting::generate_normalized(
        hist_hhh_r,
        hist_hhh_m,
        scale_low,
        scale_high);

    double bandwidth = 2.;
    //auto hhh_kde_r = generate_kde(hhh_r_tree, bandwidth);
    auto hha_kde_r = generate_kde(hha_r_tree, bandwidth);
    auto haa_kde_r = generate_kde(haa_r_tree, bandwidth);
    //auto aaa_kde_r = generate_kde(aaa_r_tree, bandwidth);

    //auto hhh_kde_m = generate_kde(hhh_m_tree, bandwidth);
    auto hha_kde_m = generate_kde(hha_m_tree, bandwidth);
    auto haa_kde_m = generate_kde(haa_m_tree, bandwidth);
    //auto aaa_kde_m = generate_kde(aaa_m_tree, bandwidth);

    for (int pnt = 0; pnt < hha_kde_r->GetN(); ++pnt) {
        auto x = hha_kde_r->GetX()[pnt];
        auto y = hha_kde_r->GetY()[pnt];
        hha_kde_r->SetPoint(
            pnt,
            x,
            y * hist_hha_m_norm->Integral() / hist_hha_m_norm->Integral());
    }

    auto corrected_mixed = new TGraph();
    auto corrected_residuals = new TGraph();
    //for (double e_rel = .1; e_rel < 100; e_rel += 1) {
    for (int bin = 1; bin <= hist_aaa_r->GetNbinsX(); ++bin) {
        auto e_rel = hist_aaa_r->GetBinCenter(bin);
        auto hha_corr = hha_kde_r->Eval(e_rel) / hha_kde_m->Eval(e_rel);

        auto haa_corr = haa_kde_r->Eval(e_rel) / haa_kde_m->Eval(e_rel);

        auto aaa_corr =
            hist_aaa_r->GetBinContent(bin) / hist_aaa_m->GetBinContent(bin);

        auto aaa_corr_coulomb = haa_corr + (haa_corr - hha_corr);

        auto real_aaa = hist_aaa_r->GetBinContent(bin);
        auto mixed_aaa = real_aaa / (aaa_corr + 1. - aaa_corr_coulomb);
        corrected_mixed->AddPoint(e_rel, mixed_aaa);
        if (real_aaa > 0) {
            corrected_residuals->AddPoint(
                e_rel,
                (real_aaa - mixed_aaa) / sqrt(real_aaa));
        }
    }

    auto canvas_residuals = new TCanvas("canvas", "", 800, 600);
    corrected_residuals->SetMarkerStyle(20);
    corrected_residuals->SetMarkerSize(.5);
    auto slate = TH2D("slate", "", 100, -5, 60, 100, -50, 400);
    slate.Draw();
    auto line = TLine(0, 0, 60, 0);
    line.Draw("same");
    corrected_residuals->Draw("pl, same");
    canvas_residuals->SaveAs("corrected_residuals.png");
    canvas_residuals->SaveAs("corrected_residuals.pdf");

    auto canvas_raw = new TCanvas("canvas_raw", "", 800, 1600);
    canvas_raw->SetLogy();
    canvas_raw->Divide(1, 4, 0, 0);
    canvas_raw->cd(1);
    hist_hhh_r->Draw("ep");
    hist_hhh_m_norm->Draw("ep same");
    canvas_raw->cd(2);
    hist_hha_r->Draw("ep");
    hha_kde_r->Draw("pl same");
    hha_kde_r->SetMarkerStyle(20);
    hha_kde_r->SetMarkerSize(.5);

    //hha_kde_m->Draw("l same");
    hist_hha_m_norm->Draw("ep same");
    canvas_raw->cd(3);
    hist_haa_r->Draw("ep");
    hist_haa_m_norm->Draw("ep same");
    canvas_raw->cd(4);
    hist_aaa_r->Draw("ep");
    hist_aaa_m_norm->Draw("ep same");
    corrected_mixed->Draw("l same");
    corrected_mixed->SetLineColor(kBlack);
    canvas_raw->SaveAs("canvas_raw.png");
}
