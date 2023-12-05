#include <Rtypes.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVirtualPad.h>

#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <memory>

double integrate_hist(TH1D* hist, double low, double high) {
    int bin_low = hist->GetXaxis()->FindBin(low);
    int bin_high = hist->GetXaxis()->FindBin(high);

    return hist->Integral(bin_low, bin_high);
}

void preliminary_formatting(TH2D* hist) {
    gStyle->SetTitleFont(112, "please_work");
    hist->SetLineWidth(2);
    hist->SetMarkerSize(2);

    hist->GetXaxis()->SetNdivisions(108);
    hist->GetYaxis()->SetNdivisions(108);

    hist->SetTitle("");

    hist->GetXaxis()->SetLabelSize(0.03);
    hist->GetYaxis()->SetLabelSize(0.03);

    hist->GetYaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetTitleOffset(1.1);
    hist->GetXaxis()->SetTitleSize(0.03);
    hist->GetYaxis()->SetTitleSize(0.03);
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();
    hist->SetTitleOffset(.5, "main");

    hist->GetXaxis()->SetLabelFont(112);
    hist->GetYaxis()->SetLabelFont(112);
    hist->GetXaxis()->SetTitleFont(112);
    hist->GetYaxis()->SetTitleFont(112);
    hist->SetTitleFont(112);
    hist->SetLabelFont(112);
}

TH1D* generate_normalized(
    TH1D* standard,
    TH1D* to_be_scaled,
    double low,
    double high) {
    auto result = (TH1D*)to_be_scaled->Clone();

    double norm_standard = integrate_hist(standard, low, high);
    double norm = integrate_hist(result, low, high);

    result->Scale(
        norm_standard * standard->GetBinWidth(1) / norm
        / to_be_scaled->GetBinWidth(1));

    return result;
}

int main() {
    auto file_path = std::string(
        "K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_3a_min\\real\\relative_energy_after_scrambling_3a.root");
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2a_1h_min\\real\\relative_energy_after_scrambling_2a_1h.root");
    //"K:\\tamu_data\\exp\\si28_c_35\\anl\\filt_2d_strict\\real\\relative_energy_after_scrambling_2d.root");

    auto file = std::make_unique<TFile>(file_path.c_str(), "READ");

    //TH1D* mixed = file->Get<TH1D>("e_rel_mixed");
    //TH1D* real = file->Get<TH1D>("e_rel_real");
    TH2D* e_rel_real_v_rndm_phi = file->Get<TH2D>("e_rel_real_v_rndm_phi");
    auto e_rel_real_v_rndm_phi_clone = (TH2D*)e_rel_real_v_rndm_phi->Clone();

    auto real = e_rel_real_v_rndm_phi->ProjectionY("real");
    //e_rel_real_v_rndm_phi->Smooth();

    //e_rel_real_v_rndm_phi->Smooth();
    std::vector<double> ratios;
    auto norm_start = 20.;
    int bin_high = e_rel_real_v_rndm_phi->GetXaxis()->FindBin(20.);

    double corrective_factor = 0;
    for (int itr = 0; itr < 1; ++itr) {
        auto n_bins = e_rel_real_v_rndm_phi->GetNbinsX();

        //for (int bin = 1; bin <= e_rel_real_v_rndm_phi->GetNbinsX(); ++bin) {
        for (int bin = bin_high; bin >= 1; --bin) {
            auto slice_y =
                e_rel_real_v_rndm_phi->ProjectionY("slice_y", bin, bin);
            auto slice_x =
                e_rel_real_v_rndm_phi->ProjectionX("slice_x", bin, bin);

            if (slice_x->Integral(bin, n_bins) == 0
                || slice_y->Integral(bin, n_bins) == 0) {
                continue;
            }

            //auto ratio = slice_y->Integral(0, bin) / slice_x->Integral(0, bin);
            auto ratio =
                slice_y->Integral(bin, n_bins) / slice_x->Integral(bin, n_bins);

            ratios.push_back(ratio);

            if (corrective_factor == 0) {
                corrective_factor = 1. / ratio;
                std::cout << corrective_factor << std::endl;
                //corrective_factor = 1.;
            } else {
                slice_x->Scale(ratio * corrective_factor);

                //for (int bin_x = 1; bin_x <= slice_x->GetNbinsX(); ++bin_x) {
                //e_rel_real_v_rndm_phi->SetBinContent(
                //bin_x,
                //bin,
                //slice_x->GetBinContent(bin_x));
                //}
                for (int bin_x = 1; bin_x <= slice_x->GetNbinsX(); ++bin_x) {
                    e_rel_real_v_rndm_phi->SetBinContent(
                        bin_x,
                        bin,
                        slice_x->GetBinContent(bin_x));
                }
            }

            //mixed_rescaled->Add(slice_x);
        }
    }

    for (int bin_x = 1; bin_x <= e_rel_real_v_rndm_phi->GetNbinsX(); ++bin_x) {
        for (int bin_y = 1; bin_y < bin_x; ++bin_y) {
            //auto new_content =
            //e_rel_real_v_rndm_phi->GetBinContent(bin_y, bin_x);
            //e_rel_real_v_rndm_phi->SetBinContent(bin_x, bin_y, new_content);

            //auto new_content =
            //e_rel_real_v_rndm_phi->GetBinContent(bin_x, bin_y)
            //+ e_rel_real_v_rndm_phi->GetBinContent(bin_y, bin_x);
            //new_content /= 2;

            //e_rel_real_v_rndm_phi->SetBinContent(bin_y, bin_x, new_content);
            //e_rel_real_v_rndm_phi->SetBinContent(bin_x, bin_y, new_content);
        }
    }
    /*
    auto bkgnd = (TH1D*)real->Clone();
    for (int idx = 0; idx < ratios.size(); ++idx) {
        bkgnd->SetBinContent(idx + 1, real->GetBinContent(idx + 1) * ratios[idx]);
    }
    */
    auto bkgnd = e_rel_real_v_rndm_phi->ProjectionX("bkgnd");
    //auto bkgnd = e_rel_real_v_rndm_phi->ProjectionY("bkgnd");
    bkgnd = generate_normalized(real, bkgnd, norm_start, 100);

    auto canvas = std::make_unique<TCanvas>("canvas", "canvas", 800, 800);
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(255);

    canvas->cd();
    real->SetLineColor(kBlack);
    bkgnd->SetLineColor(kRed);
    real->GetXaxis()->SetRangeUser(0., 60.);
    bkgnd->GetXaxis()->SetRangeUser(0., 60.);
    real->Draw();
    bkgnd->Draw("SAME");
    canvas->SaveAs("bkgnd.png");

    auto line = std::make_unique<TLine>(0., 0., 60., 60.);
    line.get()->SetLineColor(kBlack);
    line.get()->SetLineWidth(2);

    preliminary_formatting(e_rel_real_v_rndm_phi);
    e_rel_real_v_rndm_phi->SetTitle("Post-Filter");
    e_rel_real_v_rndm_phi->Draw("COL");
    e_rel_real_v_rndm_phi->GetZaxis()->SetRangeUser(0., 1'000.);
    e_rel_real_v_rndm_phi->GetXaxis()->SetRangeUser(0., 60.);
    e_rel_real_v_rndm_phi->GetYaxis()->SetRangeUser(0., 60.);
    line->Draw("same");
    //gPad->SetLogz();
    canvas->SaveAs("e_rel_real_v_rndm_phi_filt.png");

    preliminary_formatting(e_rel_real_v_rndm_phi_clone);
    e_rel_real_v_rndm_phi_clone->Draw("COL");
    e_rel_real_v_rndm_phi_clone->SetTitle(
        "Pre-Filter;E_{rel} [MeV] (Scrambled); E_{rel} [MeV] (Real)");
    e_rel_real_v_rndm_phi_clone->GetXaxis()->SetRangeUser(0., 60.);
    e_rel_real_v_rndm_phi_clone->GetYaxis()->SetRangeUser(0., 60.);
    e_rel_real_v_rndm_phi_clone->GetZaxis()->SetRangeUser(0., 1'000.);
    line->Draw("same");
    //gPad->SetLogz();
    canvas->SaveAs("e_rel_real_v_rndm_phi.png");
}
