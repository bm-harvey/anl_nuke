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

    auto file_3h_real = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\real\\relative_energy_3h.root");
    auto file_1a_2h_real = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\real\\relative_energy_1a_2h.root");
    auto file_2a_1h_real = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\real\\relative_energy_2a_1h.root");
    auto file_3a_real = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\real\\relative_energy_3a.root");

    //auto file_3h_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_rndm_phi\\relative_energy_3h.root");
    //auto file_1a_2h_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_rndm_phi\\relative_energy_1a_2h.root");
    //auto file_2a_1h_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_rndm_phi\\relative_energy_2a_1h.root");
    //auto file_3a_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_rndm_phi\\relative_energy_3a.root");

    //auto file_3h_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_mix_3h_unique\\relative_energy_3h.root");
    //auto file_1a_2h_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_mix_1a_2h_unique\\relative_energy_1a_2h.root");
    //auto file_2a_1h_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_mix_2a_1h_unique\\relative_energy_2a_1h.root");
    //auto file_3a_mixed = TFile::Open(
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_mix_3a_unique\\relative_energy_3a.root");

    auto file_3h_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3h_min\\mixed_mix_3h\\relative_energy_3h.root");
    auto file_1a_2h_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_1a_2h_min\\mixed_mix_1a_2h\\relative_energy_1a_2h.root");
    auto file_2a_1h_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\mixed_mix_2a_1h\\relative_energy_2a_1h.root");
    auto file_3a_mixed = TFile::Open(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\mixed_mix_3a\\relative_energy_3a.root");

    //
    // GET THE HISTOGRAMS
    //
    auto hist_3h_real = file_3h_real->Get<TH1D>("e_rel_hist");
    auto hist_1a_2h_real = file_1a_2h_real->Get<TH1D>("e_rel_hist");
    auto hist_2a_1h_real = file_2a_1h_real->Get<TH1D>("e_rel_hist");
    auto hist_3a_real = file_3a_real->Get<TH1D>("e_rel_hist");
    auto hist_3h_mixed = file_3h_mixed->Get<TH1D>("e_rel_hist");
    auto hist_1a_2h_mixed = file_1a_2h_mixed->Get<TH1D>("e_rel_hist");
    auto hist_2a_1h_mixed = file_2a_1h_mixed->Get<TH1D>("e_rel_hist");
    auto hist_3a_mixed = file_3a_mixed->Get<TH1D>("e_rel_hist");

    //
    // REBIN
    //
    hist_3h_real->RebinX(16);
    hist_3h_mixed->RebinX(16);

    hist_1a_2h_real->RebinX(16);
    hist_1a_2h_mixed->RebinX(16);

    hist_2a_1h_real->RebinX(16);
    hist_2a_1h_mixed->RebinX(16);

    hist_3a_real->RebinX(16);
    hist_3a_mixed->RebinX(16);
    //
    // PRELIMINARY FORMATTING
    //
    formatting::preliminary_formatting(hist_3a_real, kCyan + 3);
    formatting::preliminary_formatting(hist_2a_1h_real, kGreen + 3);
    formatting::preliminary_formatting(hist_1a_2h_real, kMagenta + 3);
    formatting::preliminary_formatting(hist_3h_real, kRed + 3);
    formatting::preliminary_formatting(hist_3a_mixed, kCyan + 1);
    formatting::preliminary_formatting(hist_2a_1h_mixed, kGreen + 1);
    formatting::preliminary_formatting(hist_1a_2h_mixed, kMagenta + 1);
    formatting::preliminary_formatting(hist_3h_mixed, kRed + 1);

    //
    // SCALING PARAMETERS
    //
    double scaling_low = 40;
    double scaling_high = 100;

    //
    // NORMALIZE THE MIXED HISTOGRAMS
    //
    auto hist_3h_mix_scaled = formatting::generate_normalized(
        hist_3h_real,
        hist_3h_mixed,
        scaling_low,
        scaling_high);
    auto hist_1a_2h_mix_scaled = formatting::generate_normalized(
        hist_1a_2h_real,
        hist_1a_2h_mixed,
        scaling_low,
        scaling_high);
    auto hist_2a_1h_mix_scaled = formatting::generate_normalized(
        hist_2a_1h_real,
        hist_2a_1h_mixed,
        scaling_low,
        scaling_high);
    auto hist_3a_mix_scaled = formatting::generate_normalized(
        hist_3a_real,
        hist_3a_mixed,
        scaling_low,
        scaling_high);

    auto canvas = new TCanvas("canvas", "", 800, 600);
    canvas->SetLogy();

    hist_3h_real->GetXaxis()->SetRangeUser(-5, 100);
    hist_1a_2h_real->GetXaxis()->SetRangeUser(-5, 100);
    hist_2a_1h_real->GetXaxis()->SetRangeUser(-5, 100);
    hist_3a_real->GetXaxis()->SetRangeUser(-5, 100);

    hist_3a_real->Draw("ep same");
    hist_2a_1h_real->Draw("ep same");
    hist_1a_2h_real->Draw("ep same");
    hist_3h_real->Draw("ep same");

    hist_3h_mix_scaled->Draw("l same");
    hist_1a_2h_mix_scaled->Draw("l same");
    hist_2a_1h_mix_scaled->Draw("l same");
    hist_3a_mix_scaled->Draw("l same");

    hist_3a_real->SetTitle(";E_{rel} [MeV];Yield");

    auto legend = new TLegend(0.1, 0.9, 0.9, 0.95);
    legend->SetNColumns(4);
    legend->AddEntry(hist_3h_mixed, "hhh", "lp");
    legend->AddEntry(hist_1a_2h_mixed, "ahh", "lp");
    legend->AddEntry(hist_2a_1h_mixed, "aah", "lp");
    legend->AddEntry(hist_3a_mixed, "aaa", "lp");
    legend->SetTextFont(112);
    legend->Draw("same");
    canvas->SaveAs("3a.png");

    auto* corr_3h =
        generate_correlation_func(hist_3h_real, hist_3h_mix_scaled, 3);
    auto* corr_1a_2h =
        generate_correlation_func(hist_1a_2h_real, hist_1a_2h_mix_scaled, 2);
    auto* corr_2a_1h =
        generate_correlation_func(hist_2a_1h_real, hist_2a_1h_mix_scaled, 1);
    auto* corr_3a =
        generate_correlation_func(hist_3a_real, hist_3a_mix_scaled, 0);

    formatting::preliminary_formatting(corr_3h, kRed + 3);
    formatting::preliminary_formatting(corr_1a_2h, kMagenta + 3);
    formatting::preliminary_formatting(corr_2a_1h, kGreen + 3);
    formatting::preliminary_formatting(corr_3a, kCyan + 3);

    auto corr_canvas = new TCanvas("corr_canvas", "", 800, 600);

    auto corr_slate = new TH2D(
        "corr_slate",
        ";E_{rel} [MeV];Real / Mixed + Offset",
        100,
        0,
        100,
        100,
        0,
        5);

    formatting::preliminary_formatting(corr_slate);
    corr_slate->Draw();
    corr_2a_1h->Draw("ep same");
    corr_1a_2h->Draw("ep same");
    corr_3h->Draw("ep same");
    corr_3a->Draw("ep same");
    auto corr_legend = new TLegend(0.1, 0.9, 0.9, 0.95);

    corr_legend->SetNColumns(4);
    corr_legend->AddEntry(corr_3h, "hhh + 3", "lp");
    corr_legend->AddEntry(corr_1a_2h, "ahh + 2", "lp");
    corr_legend->AddEntry(corr_2a_1h, "aah + 1", "lp");
    corr_legend->AddEntry(corr_3a, "aaa + 0", "lp");
    corr_legend->SetTextFont(112);
    corr_legend->Draw("same");
    corr_canvas->SaveAs("corr.png");
    gPad->SetLogx();
    corr_canvas->SaveAs("corr_logx.png");

    //
    //
    //
    //

    auto* fit_func = new TF1(
        "fit_func",
        "exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 ) * ( [4]  + [5] * x + [3]*x*x)",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 )",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 ) + (1 + [3]*x )",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 ) * (1 + [3] * cosh(x *[2]))",
        //"exp(-[0]/[1]) * ((1 + exp([0]/[1])) / (1+exp(([0] - x)/[1])) - 1 ) * ([2] * x * log(x) )",

        0,
        100);
    fit_func->SetParameter(0, .5);
    fit_func->SetParameter(1, .5);
    fit_func->SetParameter(2, 0);
    fit_func->SetParameter(3, 0);
    fit_func->SetParameter(4, 0);
    fit_func->SetParameter(5, 0);
    fit_func->SetParLimits(0, 0, 10);
    fit_func->SetParLimits(1, -5, 10);
    fit_func->SetParLimits(2, -5, 10);
    fit_func->SetParLimits(3, -5, 10);
    fit_func->SetParLimits(4, -5, 10);
    fit_func->SetParLimits(5, -100, 100);

    //auto* fit_func = new TF1("fit_func", "pol5", 0, 100);
    //fit_func->SetParameter(0, 0);
    //fit_func->SetParLimits(0, 0, 0);

    auto* corr_1a_2h_unshifted = (TGraphErrors*)corr_1a_2h->Clone();
    for (int pnt = 0; pnt < corr_1a_2h->GetN(); ++pnt) {
        double x {}, y {};
        corr_1a_2h->GetPoint(pnt, x, y);
        corr_1a_2h_unshifted->SetPoint(pnt, x, y - 2);
    }
    corr_1a_2h_unshifted->Fit(fit_func, "R");
    auto* fit_1a_2h_corr = corr_1a_2h_unshifted->GetFunction("fit_func");
    fit_1a_2h_corr->SetLineColor(kMagenta + 4);
    fit_1a_2h_corr->SetLineWidth(4);

    auto* corr_2a_1h_unshifted = (TGraphErrors*)corr_2a_1h->Clone();
    for (int pnt = 0; pnt < corr_2a_1h->GetN(); ++pnt) {
        double x {}, y {};
        corr_2a_1h->GetPoint(pnt, x, y);
        corr_2a_1h_unshifted->SetPoint(pnt, x, y - 1);
    }

    corr_2a_1h_unshifted->Fit(fit_func, "R");
    auto* fit_2a_1h_corr = corr_2a_1h_unshifted->GetFunction("fit_func");
    fit_2a_1h_corr->SetLineColor(kGreen + 1);
    fit_2a_1h_corr->SetLineWidth(2);

    auto* corr_3h_unshifted = (TGraphErrors*)corr_3h->Clone();
    for (int pnt = 0; pnt < corr_3h->GetN(); ++pnt) {
        double x {}, y {};
        corr_3h->GetPoint(pnt, x, y);
        corr_3h_unshifted->SetPoint(pnt, x, y - 3);
    }

    corr_3h_unshifted->Fit(fit_func, "R");
    auto* fit_3h_corr = corr_3h_unshifted->GetFunction("fit_func");
    fit_3h_corr->SetLineColor(kRed + 1);
    fit_3h_corr->SetLineWidth(2);

    auto* aaa_correction = (TF1*)fit_func->Clone();
    aaa_correction->SetLineColor(kCyan + 1);
    aaa_correction->SetLineWidth(5);

    for (int param_idx = 0; param_idx < aaa_correction->GetNpar();
         ++param_idx) {
        auto* params = new TGraphErrors();
        params->SetPoint(0, 0, fit_3h_corr->GetParameter(param_idx));
        params->SetPointError(0, 0, fit_3h_corr->GetParError(param_idx));

        params->SetPoint(1, 1, fit_1a_2h_corr->GetParameter(param_idx));
        params->SetPointError(1, 0, fit_1a_2h_corr->GetParError(param_idx));

        params->SetPoint(2, 2, fit_2a_1h_corr->GetParameter(param_idx));
        params->SetPointError(2, 0, fit_2a_1h_corr->GetParError(param_idx));

        params->Fit("pol1", "Q");

        auto correction = params->GetFunction("pol1")->Eval(3);

        aaa_correction->SetParameter(param_idx, correction);
    }

    auto* corr_2a_1h_corrected = (TGraphErrors*)corr_2a_1h->Clone();
    corr_2a_1h_corrected->SetMarkerColor(kGreen + 4);
    corr_2a_1h_corrected->SetLineColor(kGreen + 4);

    auto* corr_3a_corrected = (TGraphErrors*)corr_3a->Clone();
    corr_3a_corrected->SetMarkerColor(kCyan + 4);
    corr_3a_corrected->SetLineColor(kCyan + 4);

    for (int pnt = 0; pnt < corr_3a->GetN(); ++pnt) {
        auto e_rel = corr_3a->GetX()[pnt];
        auto corr_3a_val = corr_3a->GetY()[pnt];

        auto corrected = corr_3a_val + 1 - aaa_correction->Eval(e_rel);
        double err = 0;

        corr_3a_corrected->SetPoint(pnt, e_rel, corrected);
        corr_3a_corrected->SetPointError(pnt, 0, err);
    }

    double x_max = 100;
    auto* canvas_corrected = new TCanvas("canvas_corrected", "", 800, 600);
    auto* slate_corrected = new TH2D(
        "slate_corrected",
        ";E_{rel} [MeV];Real / Mixed",
        100,
        0,
        x_max,
        100,
        0,
        5);

    formatting::preliminary_formatting(slate_corrected);
    slate_corrected->Draw();

    //corr_1a_2h_unshifted->Draw("ep ");
    //corr_1a_2h_unshifted->SetLineColor(kMagenta - 10);
    //corr_1a_2h_unshifted->SetMarkerColor(kMagenta - 10);

    //corr_2a_1h_unshifted->Draw("ep same");
    //fit_2a_1h_corr->Draw("l same");
    //corr_2a_1h_unshifted->SetLineColor(kGreen - 10);
    //corr_2a_1h_unshifted->SetMarkerColor(kGreen - 10);

    //corr_3h_unshifted->Draw("ep same");
    //corr_3h_unshifted->SetLineColor(kRed -10);
    //corr_3h_unshifted->SetMarkerColor(kRed -10);

    corr_3a->Draw("ep same");

    //corr_3a_corrected->Draw("ep same");
    //corr_3a_corrected->SetLineColor(kBlack);
    //corr_3a_corrected->SetMarkerColor(kBlack);

    fit_1a_2h_corr->SetLineColor(kMagenta + 1);
    fit_1a_2h_corr->SetLineWidth(2);

    fit_3h_corr->Draw("l same");
    fit_1a_2h_corr->Draw("l same");
    fit_2a_1h_corr->Draw("l same");

    //corr_2a_1h_corrected->Draw("ep same");

    auto* line = new TLine(0, 1, x_max, 1);
    //line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->SetLineColor(kGray + 2);
    line->Draw("same");

    aaa_correction->Draw("l same");
    aaa_correction->SetLineColor(kCyan + 4);
    aaa_correction->SetMarkerColor(kCyan + 4);

    canvas_corrected->SaveAs("corr_corrected.png");

    gPad->SetLogx();
    slate_corrected->GetXaxis()->SetRangeUser(.001, 40);
    canvas_corrected->SaveAs("corr_corrected_logx.png");

    gPad->SetLogx(0);
    x_max = 20;
    slate_corrected->GetXaxis()->SetRangeUser(0, x_max);
    line = new TLine(0, 1, x_max, 1);

    canvas_corrected->SaveAs("corr_corrected_zoom.png");

    auto* corrected_mixed = (TH1D*)hist_3a_mix_scaled->Clone();
    std::cout << "\n\n"
              << __LINE__ << " | " << corrected_mixed->GetNbinsX() << "\n"
              << std::endl;
    for (int bin = 1; bin < corrected_mixed->GetNbinsX(); ++bin) {
        corrected_mixed->SetBinContent(
            bin,
            corrected_mixed->GetBinContent(bin)
                * aaa_correction->Eval(corrected_mixed->GetBinCenter(bin)));
    }

    auto* residuals = (TH1D*)hist_3a_real->Clone();
    for (int bin = 1; bin <= hist_3a_real->GetNbinsX(); ++bin) {
        residuals->SetBinContent(
            bin,
            //real
            hist_3a_real->GetBinContent(bin)
                //corrected mixed
                - hist_3a_real->GetBinContent(bin)
                    / corr_3a_corrected->Eval(hist_3a_real->GetBinCenter(bin)));
    }

    auto* residuals_uncorrected = (TH1D*)hist_3a_real->Clone();
    residuals_uncorrected->Add(hist_3a_mixed, -1);

    auto* canvas_residuals = new TCanvas("canvas_residuals", "", 800, 600);
    //residuals->RebinX(2);
    residuals->SetLineColor(kBlack);
    residuals->Draw("hist");
    //residuals_uncorrected->Draw("same hist");
    residuals->GetXaxis()->SetRangeUser(0, 80);
    canvas_residuals->SaveAs("residuals.png");
    residuals->GetXaxis()->SetRangeUser(0, 20);
    canvas_residuals->SaveAs("residuals_zoom.png");
}
