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

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "formatting.hpp"

#define PRINTLN std::cout << __LINE__ << std::endl;

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

    //auto file_hhh_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_mix_3h_unique\\relative_energy_3h.root");
    //auto file_hha_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_mix_1a_2h_unique\\relative_energy_1a_2h.root");
    //auto file_haa_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_mix_2a_1h_unique\\relative_energy_2a_1h.root");
    //auto file_aaa_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_mix_3a_unique\\relative_energy_3a.root");

    auto file_hhh_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_mix_3h\\relative_energy_3h.root");
    auto file_hha_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_mix_1a_2h\\relative_energy_1a_2h.root");
    auto file_haa_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_mix_2a_1h\\relative_energy_2a_1h.root");
    auto file_aaa_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_mix_3a\\relative_energy_3a.root");

    double scale_low = 40;
    double scale_high = 100;

    auto hist_hhh_r = file_hhh_real->Get<TH1D>("e_rel_hist");
    auto hist_hha_r = file_hha_real->Get<TH1D>("e_rel_hist");
    auto hist_haa_r = file_haa_real->Get<TH1D>("e_rel_hist");
    auto hist_aaa_r = file_aaa_real->Get<TH1D>("e_rel_hist");

    auto hist_hhh_m = file_hhh_mixed->Get<TH1D>("e_rel_hist");
    auto hist_hha_m = file_hha_mixed->Get<TH1D>("e_rel_hist");
    auto hist_haa_m = file_haa_mixed->Get<TH1D>("e_rel_hist");
    auto hist_aaa_m = file_aaa_mixed->Get<TH1D>("e_rel_hist");

    for (int i = 0; i < 3; ++i) {
        hist_hhh_r->Rebin(2);
        hist_hha_r->Rebin(2);
        hist_haa_r->Rebin(2);
        hist_aaa_r->Rebin(2);
        hist_hhh_m->Rebin(2);
        hist_hha_m->Rebin(2);
        hist_haa_m->Rebin(2);
        hist_aaa_m->Rebin(2);
    }

    formatting::preliminary_formatting(hist_hhh_r, kOrange + 1);
    formatting::preliminary_formatting(hist_hha_r, kViolet + 1);
    formatting::preliminary_formatting(hist_haa_r, kGreen + 1);
    formatting::preliminary_formatting(hist_aaa_r, kCyan + 1);

    formatting::preliminary_formatting(hist_hhh_m, kOrange + 3);
    formatting::preliminary_formatting(hist_hha_m, kViolet + 3);
    formatting::preliminary_formatting(hist_haa_m, kGreen + 3);
    formatting::preliminary_formatting(hist_aaa_m, kCyan + 3);

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
    //
    // FIT ALL WITH fit_func_2 DISTRIBUTIONS
    //
    auto* fit_func_2 = new TF1(
        "fit_func_2",
        "[0]*ROOT::Math::beta_pdf(x * [4] , [1], [2])",
        //"([0] + [1]*x + [2] * x * x)/([3] + [4] * x + [5] * x * x)",
        //"pol10",
        //"[0] * x*x*x *exp(-x*x/[1])",
        //"[0]*pow(x,[2]) * exp(-x * x/ 2 / [1] / [1])",
        0,
        80);
    fit_func_2->SetParameters(1, 4, 2.5, 0, 1. / 11, 0, 0);
    fit_func_2->SetParLimits(0, -100, 1000000);
    fit_func_2->SetParLimits(1, 0, 100);
    fit_func_2->SetParLimits(2, 0, 100);
    fit_func_2->SetParLimits(3, -100, 100);
    fit_func_2->SetParLimits(4, -100, 100);
    fit_func_2->SetParLimits(5, -100, 100);

    //hist_hhh_m_norm->Fit(fit_func_2, "R");
    //hist_hhh_r->Fit(fit_func_2, "R");
    //hist_hha_m_norm->Fit(fit_func_2, "R");
    //hist_hha_r->Fit(fit_func_2, "R");
    //hist_haa_m_norm->Fit(fit_func_2, "R");
    //hist_haa_r->Fit(fit_func_2, "R");
    //hist_aaa_m_norm->Fit(fit_func_2, "R");
    //hist_aaa_r->Fit(fit_func_2, "R");
    //hist_hhh_m_norm->Fit(fit_func_2, "R");

    hist_hhh_r->Fit(fit_func_2, "RL");
    hist_hha_m_norm->Fit(fit_func_2, "RL");
    hist_hha_r->Fit(fit_func_2, "RL");
    hist_haa_m_norm->Fit(fit_func_2, "RL");
    hist_haa_r->Fit(fit_func_2, "RL");
    hist_aaa_m_norm->Fit(fit_func_2, "RL");
    hist_aaa_r->Fit(fit_func_2, "RL");

    //hist_hhh_m_norm->Fit("pol9", "R");
    //hist_hhh_r->Fit("pol9", "R");
    //hist_hha_m_norm->Fit("pol9", "R");
    //hist_hha_r->Fit("pol9", "R");
    //hist_haa_m_norm->Fit("pol9", "R");
    //hist_haa_r->Fit("pol9", "R");
    //hist_aaa_m_norm->Fit("pol9", "R");
    //hist_aaa_r->Fit("pol9", "R");

    auto* fit_func_2_hhh_m_norm = hist_hhh_m_norm->GetFunction("fit_func_2");
    auto* fit_func_2_hhh_r = hist_hhh_r->GetFunction("fit_func_2");
    auto* fit_func_2_hha_m_norm = hist_hha_m_norm->GetFunction("fit_func_2");
    auto* fit_func_2_hha_r = hist_hha_r->GetFunction("fit_func_2");
    auto* fit_func_2_haa_m_norm = hist_haa_m_norm->GetFunction("fit_func_2");
    auto* fit_func_2_haa_r = hist_haa_r->GetFunction("fit_func_2");
    auto* fit_func_2_aaa_m_norm = hist_aaa_m_norm->GetFunction("fit_func_2");
    auto* fit_func_2_aaa_r = hist_aaa_r->GetFunction("fit_func_2");

    hist_hhh_m_norm->GetXaxis()->SetRangeUser(0, 80);
    hist_hhh_r->GetXaxis()->SetRangeUser(0, 80);
    hist_hha_m_norm->GetXaxis()->SetRangeUser(0, 80);
    hist_hha_r->GetXaxis()->SetRangeUser(0, 80);
    hist_haa_m_norm->GetXaxis()->SetRangeUser(0, 80);
    hist_haa_r->GetXaxis()->SetRangeUser(0, 80);
    hist_aaa_m_norm->GetXaxis()->SetRangeUser(0, 80);
    hist_aaa_r->GetXaxis()->SetRangeUser(0, 80);

    auto corrected_mixed = new TGraph();
    auto corrected_residuals = new TGraph();
    //for (double e_rel = .1; e_rel < 100; e_rel += 1) {
    for (int bin = 1; bin < hist_aaa_r->GetNbinsX(); ++bin) {
        auto e_rel = hist_aaa_r->GetBinCenter(bin);
        auto hha_corr =
            fit_func_2_hha_r->Eval(e_rel) / fit_func_2_hha_m_norm->Eval(e_rel);
        auto haa_corr =
            fit_func_2_haa_r->Eval(e_rel) / fit_func_2_haa_m_norm->Eval(e_rel);
        auto aaa_corr = hist_aaa_r->GetBinContent(bin)
            / hist_aaa_m_norm->GetBinContent(bin);

        auto aaa_corr_coulomb = haa_corr + haa_corr - hha_corr;

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
    auto slate = TH2D("slate", "", 100, 0, 60, 100, -50, 400);
    slate.Draw();
    auto line = TLine(0, 0, 100, 0);
    line.Draw("same");
    corrected_residuals->Draw("pl, same");

    canvas_residuals->SaveAs("fig/canvas_residuals.png");

    auto canvas_raw = new TCanvas("canvas_raw", "", 800, 1600);
    canvas_raw->SetLogy();
    canvas_raw->Divide(1, 4, 0, 0);
    canvas_raw->cd(1);
    gPad->SetLogy();
    hist_hhh_r->Draw("ep");
    fit_func_2_hhh_r->Draw("same");
    hist_hhh_m_norm->Draw("ep same");
    fit_func_2_hhh_m_norm->Draw("same");
    canvas_raw->cd(2);
    gPad->SetLogy();
    hist_hha_r->Draw("ep");
    fit_func_2_hha_r->Draw("same");
    hist_hha_m_norm->Draw("ep same");
    fit_func_2_hha_m_norm->Draw("same");
    canvas_raw->cd(3);
    gPad->SetLogy();
    hist_haa_r->Draw("ep");
    fit_func_2_haa_r->Draw("same");
    hist_haa_m_norm->Draw("ep same");
    fit_func_2_haa_m_norm->Draw("same");
    canvas_raw->cd(4);
    gPad->SetLogy();
    hist_aaa_r->Draw("ep");
    fit_func_2_aaa_r->Draw("same");
    hist_aaa_m_norm->Draw("ep same");
    fit_func_2_aaa_m_norm->Draw("same");
    corrected_mixed->Draw("l same");
    corrected_mixed->SetLineColor(kBlack);
    canvas_raw->SaveAs("fig/canvas_raw.png");

    //
    //
    //
    //
    //
    auto* corr_hhh = generate_correlation_func(hist_hhh_r, hist_hhh_m_norm);
    auto* corr_hha = generate_correlation_func(hist_hha_r, hist_hha_m_norm);
    auto* corr_haa = generate_correlation_func(hist_haa_r, hist_haa_m_norm);
    auto* corr_aaa = generate_correlation_func(hist_aaa_r, hist_aaa_m_norm);

    auto* fit_func = new TF1(
        "fit_func",
        //"exp(-[0] * x) * (1 - [1]/x) + (1-exp(-[0] * x)) * ([2] + [3]*x + [4]*x*x)",
        //"exp(-[0] * x) * (1 - [1]/x) + (1-exp(-[0] * x)) * ([2] + [3]*x + [4]*x*x)",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1]))- 1 ) * (exp(-x*[4]) + [2]*x*x + [3])",
        "exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1]))- 1 ) * ([2] + [3]*x + [4]*x*x )",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1]))- 1 ) * (1 + [3]*x + [4]*x*x )",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 ) * ([3]/x + x*[2] )",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 )",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 ) + (1 + [3]*x )",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 ) * (1 + [3] * cosh((x-[4]) *[2]))",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 ) * (1 + [2] * x * log(x) )",

        2,
        100);
    fit_func->SetParameter(0, 1.2);
    fit_func->SetParameter(1, 1.2);
    fit_func->SetParameter(2, -3);
    fit_func->SetParameter(3, -3);
    fit_func->SetParameter(4, 1);
    //fit_func->SetParameter(5, 0);

    //fit_func->SetParLimits(0, 0, 0);
    fit_func->SetParLimits(0, 0, 5);
    fit_func->SetParLimits(1, 0, 10);
    fit_func->SetParLimits(2, -5, 10);
    fit_func->SetParLimits(3, -5, 10);
    fit_func->SetParLimits(4, -5, 10);
    //fit_func->SetParLimits(5, -100, 100);

    corr_hhh->Fit(fit_func, "R");
    corr_hha->Fit(fit_func, "R");
    corr_haa->Fit(fit_func, "R");

    auto* fit_hhh = (TF1*)corr_hhh->GetFunction("fit_func");
    auto* fit_hha = (TF1*)corr_hha->GetFunction("fit_func");
    auto* fit_haa = (TF1*)corr_haa->GetFunction("fit_func");

    auto* aaa_correction = (TF1*)fit_func->Clone("fit_func");

    for (int param_idx = 0; param_idx < aaa_correction->GetNpar();
         ++param_idx) {
        auto* params = new TGraphErrors();
        //params->SetPoint(0, 0, fit_hhh->GetParameter(param_idx));
        //params->SetPointError(0, 0, fit_hhh->GetParError(param_idx));

        params->SetPoint(1, 1, fit_hha->GetParameter(param_idx));
        params->SetPointError(1, 0, fit_hha->GetParError(param_idx));

        params->SetPoint(2, 2, fit_haa->GetParameter(param_idx));
        params->SetPointError(2, 0, fit_haa->GetParError(param_idx));

        params->Fit("pol1", "Q");

        auto correction = params->GetFunction("pol1")->Eval(3);

        aaa_correction->SetParameter(param_idx, correction);
    }

    aaa_correction->SetLineColor(kBlack);

    auto* canvas = new TCanvas("canvas", "", 1000, 1000);

    corr_aaa->SetLineColor(kCyan + 3);
    corr_aaa->SetMarkerColor(kCyan + 3);
    corr_aaa->SetMarkerStyle(20);
    corr_aaa->SetMarkerSize(.5);
    corr_aaa->SetLineWidth(1);

    aaa_correction->SetLineColor(kBlack);
    aaa_correction->SetLineWidth(3);

    fit_haa->SetLineColor(kGreen + 1);
    fit_hha->SetLineColor(kViolet + 1);
    fit_hhh->SetLineColor(kOrange + 1);

    auto* slate_corrected = new TH2D(
        "slate_corrected",
        ";E_{rel} [MeV];Real / Mixed",
        100,
        0,
        100,
        100,
        0,
        5);
    formatting::preliminary_formatting(slate_corrected);

    slate_corrected->Draw();
    corr_aaa->Draw("EP SAME");
    fit_hhh->Draw("SAME");
    fit_hha->Draw("SAME");
    fit_haa->Draw("SAME");
    aaa_correction->Draw("SAME");
    corr_aaa->Fit("fit_func");
    canvas->SaveAs("fig/aaa_correction.png");

    auto* residuals = (TH1D*)hist_aaa_r->Clone("residuals");
    auto* residuals_2 = (TH1D*)hist_aaa_r->Clone("residuals");
    formatting::preliminary_formatting(residuals);
    residuals->SetTitle(";E_{rel} [MeV];(Real - Mixed) / #sqrt{Real}");

    for (int bin = 1; bin <= residuals->GetNbinsX(); ++bin) {
        double e_rel = residuals->GetBinCenter(bin);
        double real = residuals->GetBinContent(bin);

        auto corrected_correlation =
            corr_aaa->Eval(e_rel) + 1. - aaa_correction->Eval(e_rel);

        double mixed = real / corrected_correlation;

        if (real == 0) {
            residuals->SetBinContent(bin, 0);
        } else {
            residuals->SetBinContent(bin, (real - mixed) / sqrt(real));
        }

        mixed = hist_aaa_m->GetBinContent(bin) * aaa_correction->Eval(e_rel);

        if (real == 0) {
            residuals_2->SetBinContent(bin, 0);
        } else {
            residuals_2->SetBinContent(bin, (real - mixed) / sqrt(real));
        }
    }

    corr_hhh->Draw("EAP");
    corr_hhh->SetLineColor(kOrange + 1);
    corr_hhh->SetMarkerColor(kOrange + 1);
    corr_hhh->GetFunction("fit_func")->SetLineColor(kOrange + 4);
    corr_hhh->GetFunction("fit_func")->Draw("same");
    canvas->SaveAs("fig/hhh.png");

    corr_hha->Draw("EAP");
    corr_hha->SetLineColor(kViolet + 1);
    corr_hha->SetMarkerColor(kViolet + 1);
    corr_hha->GetFunction("fit_func")->SetLineColor(kViolet + 4);
    corr_hha->GetFunction("fit_func")->Draw("same");
    canvas->SaveAs("fig/hha.png");

    corr_haa->Draw("EAP");
    corr_haa->SetLineColor(kGreen + 1);
    corr_haa->SetMarkerColor(kGreen + 1);
    corr_haa->GetFunction("fit_func")->SetLineColor(kGreen + 4);
    corr_haa->GetFunction("fit_func")->Draw("same");
    canvas->SaveAs("fig/haa.png");

    canvas->Clear();
    slate_corrected->Draw();
    slate_corrected->GetYaxis()->SetRangeUser(-5, 200);

    //residuals->Draw("p");
    corrected_residuals->Draw("p same");
    formatting::preliminary_formatting(residuals);
    residuals->GetXaxis()->SetRangeUser(0, 50);
    canvas->SaveAs("fig/aaa_residuals.png");
}
