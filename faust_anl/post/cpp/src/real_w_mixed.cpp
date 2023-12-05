#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>

#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <memory>

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
    auto filter_name = "filt_3a_strict";
    auto mixer_name = "mixed_mix_3a_unique";
    auto system = std::string("3a"); 
    double scaling_region_low = 30.;
    double scaling_region_up = 50.;

    auto file_name = std::string("relative_energy_" + system + ".root");

    auto r_path = std::filesystem::path("K:\\tamu_data\\exp\\si28_c_35\\anl")
        / filter_name / "real" /file_name;

    auto m_path = std::filesystem::path("K:\\tamu_data\\exp\\si28_c_35\\anl")
        / filter_name / mixer_name / file_name;

    auto r_file = std::make_shared<TFile>(r_path.string().c_str());
    auto m_file = std::make_shared<TFile>(m_path.string().c_str(), "read");

    r_file->ls();
    m_file->ls();

    auto* r_tree = r_file->Get<TTree>("T");
    auto* m_tree = m_file->Get<TTree>("T");

    TreeBranches r_data;
    TreeBranches m_data;

    r_tree->SetBranchAddress("e_rel", &r_data.e_rel);
    m_tree->SetBranchAddress("e_rel", &m_data.e_rel);

    r_tree->SetBranchAddress("inner_angle_avg", &r_data.inner_angle);
    m_tree->SetBranchAddress("inner_angle_avg", &m_data.inner_angle);

    auto e_rel_r =
        std::make_shared<TH1D>("r_e_rel", ";E_{rel} (MeV);", 4'096, -5, 200);
    auto e_rel_m =
        std::make_shared<TH1D>("m_e_rel", ";E_{rel} (MeV);", 4'096, -5, 200);
    auto e_rel_m_resampled = std::make_shared<TH1D>(
        "m_resamled_e_rel",
        ";E_{rel} (MeV);",
        4'096,
        -5,
        80);

    auto inner_angle_r = std::make_shared<TH1D>(
        "r_inner_angle",
        ";Average of inner angles (Deg);",
        1'024,
        0.,
        360.);
    auto inner_angle_m = std::make_shared<TH1D>(
        "m_inner_angle",
        ";Average of inner angles (Deg);",
        1'024,
        0.,
        360.);
    auto inner_angle_m_resampled = std::make_shared<TH1D>(
        "m_resamled_inner_angle",
        ";Average of inner angles (Deg);",
        1'024,
        0.,
        360);

    for (unsigned int idx = 0; idx < (unsigned int)r_tree->GetEntries();
         ++idx) {
        r_tree->GetEntry(idx);
        e_rel_r->Fill(r_data.e_rel);
        inner_angle_r->Fill(r_data.inner_angle);
    }

    for (unsigned int idx = 0; idx < (unsigned int)m_tree->GetEntries();
         ++idx) {
        m_tree->GetEntry(idx);
        e_rel_m->Fill(m_data.e_rel);
        inner_angle_m->Fill(m_data.inner_angle);
    }

    auto* angle_ratio = (TH1D*)inner_angle_r->Clone();
    angle_ratio->Divide(inner_angle_m.get());

    for (unsigned int idx = 0; idx < (unsigned int)m_tree->GetEntries();
         ++idx) {
        m_tree->GetEntry(idx);
        if (gRandom->Uniform(0.1)
            < angle_ratio->Interpolate(m_data.inner_angle)) {
            e_rel_m_resampled->Fill(m_data.e_rel);
            inner_angle_m_resampled->Fill(m_data.inner_angle);
        }
    }

    //

    preliminary_formatting(e_rel_r.get(), 1, mixer_name, 0., 120);
    preliminary_formatting(e_rel_m.get(), 2, mixer_name, 0., 120);
    preliminary_formatting(e_rel_m_resampled.get(), 4, mixer_name, 0., 120.);

    preliminary_formatting(inner_angle_r.get(), 1, mixer_name, 0., 50.);
    preliminary_formatting(inner_angle_m.get(), 2, mixer_name, 0, 50.);
    preliminary_formatting(
        inner_angle_m_resampled.get(),
        4,
        mixer_name,
        0.,
        50.);

    auto* scaled_erel_m =
        generate_normalized(e_rel_r.get(), e_rel_m.get(), 0, 200);
    auto* scaled_erel_m_resampled =
        generate_normalized(e_rel_r.get(), e_rel_m_resampled.get(), 0, 200);

    auto canvas = std::make_shared<TCanvas>("canvas", "", 800, 800);

    e_rel_r->Draw("e");
    scaled_erel_m->Draw("l, same");
    //scaled_erel_m_resampled->Draw("e, same");
    canvas->Print("fig/relative_energy.png");

    inner_angle_m->Scale(1. / inner_angle_m->Integral());
    inner_angle_m_resampled->Scale(1. / inner_angle_m_resampled->Integral());
    inner_angle_r->Scale(1. / inner_angle_r->Integral());
    inner_angle_r->Draw("e");
    inner_angle_m->Draw("e, same");
    //inner_angle_m_resampled->Draw("e, same");
    canvas->Print("fig/inner_angle.png");

    auto* e_rel_r_minus_m = (TH1D*)e_rel_r->Clone();
    auto* neg_scaled_e_rel_m = (TH1D*)scaled_erel_m->Clone();
    neg_scaled_e_rel_m->Scale(-1.);
    e_rel_r_minus_m->Add(neg_scaled_e_rel_m);

    e_rel_r_minus_m->Draw();
    canvas->Print("fig/relative_energy_subtracted.png");

    auto norm_residuals = std::make_unique<TGraph>();
    auto x_axis = e_rel_r_minus_m->GetXaxis();
    for (int idx = 0; idx < x_axis->GetNbins(); ++idx) {
        auto x = x_axis->GetBinCenter(idx);
        auto y = e_rel_r_minus_m->GetBinContent(idx);

        if (y != 0) {
            y /= sqrt(fabs(y));
            norm_residuals->AddPoint(x, y);
        }
    }

    norm_residuals->SetLineWidth(2);
    norm_residuals->SetMarkerSize(2);
    norm_residuals->SetMarkerStyle(8);
    norm_residuals->Draw("P");
    canvas->Print("fig/std_residuals.png");
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

    while (hist->GetNbinsX() > 512 / 4) {
        //while (hist->GetNbinsX() > 1024) {
        hist->RebinX(2);
    }

    hist->SetTitle(title);

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
