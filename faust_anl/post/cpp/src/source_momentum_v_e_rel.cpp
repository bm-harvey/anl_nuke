
#include <Math/PdfFuncMathCore.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TH2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TMath.h>
#include <cmath>
#include <iostream>
#include <numeric>

#include "formatting.hpp"

int main() {
    formatting::set_style();

    auto real_file_name = std::string(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_7a_std\\real\\relative_energy_7a.root");

    auto* real_file = new TFile(real_file_name.c_str(), "READ");
    auto* real_hist_2d =
        (TH2D*)real_file->Get("source_particles_rel_p_w_only_2_8begs_v_e_rel");

    auto* canvas = new TCanvas("canvas", "canvas", 1200, 800);

    //
    real_hist_2d->SetTitle(
        "2 {}^{8}Be(gs) && 0 {}^{12}C(0+) && 0 {}^{12}C(3-);E_{rel} [MeV];C.M. Momenta [MeV / c]");
    formatting::preliminary_formatting(real_hist_2d);
    real_hist_2d->RebinX(4);
    real_hist_2d->RebinY(8);
    real_hist_2d->GetXaxis()->SetRangeUser(0, 125);
    auto* func = new TF1("func", "sqrt(2 * 3728 * x / 7)", 0, 125);
    real_hist_2d->Draw("col");
    func->Draw("same");
    func->SetLineColor(kBlack);

    canvas->SaveAs("source_momentum_v_e_rel.pdf");

    auto* hist = real_file->Get<TH1D>("source_particles_rel_p_w_only_2_8begs");
    hist->SetTitle(
        "2 {}^{8}Be(gs) && 0 {}^{12}C(0+) && 0 {}^{12}C(3-);C.M. Momenta [MeV / c];Yield");
    formatting::preliminary_formatting(hist);
    hist->RebinX(4);
    //auto* std_residuals = (TH1D*)hist->Clone();
    hist->Fit("pol9", "");

    canvas->SaveAs("source_momenta_2Begs.pdf");

    auto* std_residuals = new TGraph();
    auto* blank = new TH2D("blank", "blank", 100, 0, 550, 100, -5, 5);
    blank->SetTitle(
        "2 {}^{8}Be(gs) && 0 {}^{12}C(0+) && 0 {}^{12}C(3-);C.M. Momenta [MeV / c];std. residuals");
    
    auto* gauss_bkgnd = new TF2("gaus_bkgnd", "exp(- y * y / 2)", 0, 550, -5, 5);
    auto* gauss_bkgnd_hist = new TH2D("gauss_bkngd_hist", "gauss_bkngd_hist", 500, 0, 550, 500, -5, 5);
    for (int x = 1; x <= gauss_bkgnd_hist->GetNbinsX(); ++x) {
        for (int y = 1; y <= gauss_bkgnd_hist->GetNbinsY(); ++y) {
            auto x_val = gauss_bkgnd_hist->GetXaxis()->GetBinCenter(x); 
            auto y_val = gauss_bkgnd_hist->GetYaxis()->GetBinCenter(y);

            gauss_bkgnd_hist->Fill(x_val, y_val, gauss_bkgnd->Eval(x, y));
        }
    }

    //gauss_bkgnd_hist->SetContour(100);


    formatting::preliminary_formatting(blank);
    std_residuals->GetYaxis()->SetTitle("#sigma");
    auto* fit = hist->GetFunction("pol9");
    for (int bin = 0; bin < hist->GetNbinsX(); ++bin) {
        auto p = hist->GetBinCenter(bin);
        auto y = hist->GetBinContent(bin);
        auto y_fit = fit->Eval(p);
        auto sigma = hist->GetBinError(bin);

        if (y > 2) {
            std_residuals->AddPoint(p, (y - y_fit) / sigma);
        }
    }


    std_residuals->SetMarkerSize(2);
    std_residuals->SetMarkerStyle(20);
    blank->Draw();
    gauss_bkgnd_hist->Draw("same, col");
    //gauss_bkgnd_hist->Draw("col");
    std_residuals->Draw("P, same");

    canvas->SaveAs("source_momenta_2Begs_sigma.pdf");
    canvas->SaveAs("source_momenta_2Begs_sigma.png");
}
