#include <Math/PdfFuncMathCore.h>
#include <Rtypes.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TVirtualPad.h>

#include "formatting.hpp"

#define PRINTLN std::cout << __LINE__ << std::endl;

TH1D* generate_normalized(
    TH1D* standard,
    TH1D* to_be_scaled,
    double low,
    double high);

double integrate_hist(TH1D* hist, double low, double high);

int main() {
    formatting::set_style();

    auto file = TFile::Open(
        "/data/sjygroup/sjy20/bmharvey/acs/c12_si_35/anl/filt_3a_min/real/aaa_thru_be8gs.root");

    auto file_m = TFile::Open(
        "/data/sjygroup/sjy20/bmharvey/acs/c12_si_35/anl/filt_3a_min/mixed_mix_3a_unique/aaa_thru_be8gs.root");

    auto* e_ex = file->Get<TH1D>("e_ex");
    auto* e_ex_m = file_m->Get<TH1D>("e_ex");
    auto* e_ex_thru_be8_gs = file->Get<TH1D>("e_ex_thru_be8gs");
    auto* e_ex_thru_be8_2p = file->Get<TH1D>("e_ex_thru_be_2p");
    auto* e_ex_thru_be8_both = file->Get<TH1D>("e_ex_thru_both");
    auto* e_ex_2 = file->Get<TH1D>("e_ex_2a");
    auto* e_ex_thru_0p_theta = file->Get<TH2D>("ex_3_vs_theta_through_0p");
    auto* e_ex_thru_2p_theta = file->Get<TH2D>("ex_3_vs_theta_through_2p");
    auto* dalitz_0p = file->Get<TH2D>("dalitz_thru_c12_0p");
    auto* dalitz_all = file->Get<TH2D>("dalitz");
    auto* dalitz_low = file->Get<TH2D>("dalitz_thru_c12_low");
    auto* dalitz_3m = file->Get<TH2D>("dalitz_thru_c12_3m");

    e_ex_thru_be8_2p->Add(e_ex_thru_be8_both, -1.);

    formatting::preliminary_formatting(e_ex);
    formatting::preliminary_formatting(e_ex_m, kViolet + 2);

    formatting::preliminary_formatting(e_ex_thru_be8_gs, kRed + 2);
    formatting::preliminary_formatting(e_ex_thru_be8_2p, kBlue + 2);
    formatting::preliminary_formatting(e_ex_2, kBlack);

    e_ex->RebinX(4);
    e_ex_m->RebinX(4);
    e_ex_thru_be8_gs->RebinX(4);
    e_ex_thru_be8_2p->RebinX(4);
    e_ex_2->RebinX(4);

    e_ex->GetXaxis()->SetRangeUser(5, 40);
    //e_ex_m->GetXaxis()->SetRangeUser(5, 40);
    e_ex_thru_be8_gs->GetXaxis()->SetRangeUser(0, 60);
    e_ex_thru_be8_2p->GetXaxis()->SetRangeUser(0, 60);

    auto* legend = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend->SetTextFont(112);
    legend->AddEntry(e_ex, "All", "p");
    legend->AddEntry(e_ex_thru_be8_gs, "w/ ^{8}Be(0^{+})", "p");
    legend->AddEntry(e_ex_thru_be8_2p, "w/ ^{8}Be(2^{+})", "p");

    auto canvas = new TCanvas("canvas", "canvas", 800, 800);
    canvas->SetTicks(1, 1);

    e_ex->Draw("e1");
    e_ex->GetXaxis()->SetTitle("3#alpha Reconstructed E^{*}[MeV]");
    e_ex_thru_be8_gs->Draw("e1 same");
    e_ex_thru_be8_2p->Draw("e1 same");
    legend->Draw("same");
    canvas->SaveAs("decay_paths.pdf");
    canvas->SaveAs("decay_paths.png");

    e_ex_2->Draw("e1");
    e_ex_2->GetXaxis()->SetRangeUser(-2, 20);

    auto* gate_0p = new TBox(-.2, 0, 0.3, 130000);
    gate_0p->SetFillColorAlpha(kRed + 2, 0.3f);
    gate_0p->Draw("f same");

    auto* gate_2p = new TBox(2, 0, 4, 40000);
    gate_2p->SetFillColorAlpha(kBlue + 2, 0.3f);
    //gate_2p->SetFillStyle(1);
    gate_2p->Draw("f same");

    auto* legend_2 = new TLegend(0.6, 0.7, 0.8, 0.8);
    legend_2->SetTextFont(112);
    legend_2->AddEntry(gate_0p, "^{8}Be(0^{+}) gate", "f");
    legend_2->AddEntry(gate_2p, "^{8}Be(2^{+}) gate", "f");
    legend_2->Draw("same");

    e_ex_2->GetXaxis()->SetTitle("2#alpha Reconstructed E^{*}[MeV]");
    canvas->SaveAs("decay_paths_2.pdf");
    canvas->SaveAs("decay_paths_2.png");

    e_ex->GetXaxis()->SetRangeUser(0, 80);
    e_ex->Draw("e1");
    auto* scaled_m = generate_normalized(e_ex, e_ex_m, 40, 100);
    scaled_m->Draw("e1 same");
    canvas->SaveAs("c12_w_bkgnd_m.pdf");
    canvas->SaveAs("c12_w_bkgnd_m.png");

    formatting::preliminary_formatting(e_ex_thru_0p_theta);
    formatting::preliminary_formatting(e_ex_thru_2p_theta);
    e_ex_thru_0p_theta->RebinX(4);

    e_ex_thru_0p_theta->Draw("col");
    e_ex_thru_0p_theta->GetYaxis()->SetRangeUser(5, 30);
    e_ex_thru_0p_theta->GetYaxis()->SetTitle("#alpha#alpha#alpha E^{*} [MeV]");
    e_ex_thru_0p_theta->SetTitle("Through ^{8}Be(0^{+})");
    auto* energy_exis = e_ex_thru_0p_theta->GetYaxis();
    auto bin_low = energy_exis->FindBin(10.3);
    auto bin_high = energy_exis->FindBin(11.);
    auto c12_1m_theta =
        e_ex_thru_0p_theta->ProjectionX("c12_1m", bin_low, bin_high);

    bin_low = energy_exis->FindBin(6);
    bin_high = energy_exis->FindBin(8.);
    auto c12_0p_theta =
        e_ex_thru_0p_theta->ProjectionX("hoyle", bin_low, bin_high);
    formatting::preliminary_formatting(c12_1m_theta);
    formatting::preliminary_formatting(c12_0p_theta, kRed + 2);

    canvas->SaveAs("c12_w_bkgnd_m_theta_0p.pdf");
    canvas->SaveAs("c12_w_bkgnd_m_theta_0p.png");
    c12_0p_theta->Scale(1. / c12_0p_theta->Integral());
    c12_1m_theta->Scale(1. / c12_1m_theta->Integral());
    c12_0p_theta->RebinX(2);
    c12_1m_theta->RebinX(2);
    c12_1m_theta->SetMarkerSize(2);
    c12_0p_theta->SetMarkerSize(2);
    c12_1m_theta->SetLineWidth(2);
    c12_0p_theta->SetLineWidth(2);
    c12_0p_theta->Draw("e1l");
    c12_0p_theta->GetYaxis()->SetTitle("Normalized Yield");
    c12_1m_theta->Draw("e1l, same");

    auto* legend_3 = new TLegend(0.3, 0.2, 0.45, 0.3);
    legend_3->SetTextFont(112);
    legend_3->AddEntry(c12_0p_theta, "^{12}C(0^{+})", "lp");
    legend_3->AddEntry(c12_1m_theta, "^{12}C(1^{-})", "lp");

    legend_3->Draw("same");
    canvas->SaveAs("theta.png");
    canvas->SaveAs("theta.pdf");

    e_ex_thru_2p_theta->Draw("col");
    e_ex_thru_2p_theta->GetYaxis()->SetRangeUser(5, 40);
    e_ex_thru_2p_theta->GetYaxis()->SetTitle("#alpha#alpha#alpha E^{*} [MeV]");
    e_ex_thru_2p_theta->SetTitle("Through ^{8}Be(2^{+})");
    //e_ex_thru_2p_theta->RebinY(2);
    canvas->SaveAs("c12_w_bkgnd_m_theta_2p.pdf");
    canvas->SaveAs("c12_w_bkgnd_m_theta_2p.png");

    formatting::preliminary_formatting(dalitz_all);
    formatting::preliminary_formatting(dalitz_0p);
    formatting::preliminary_formatting(dalitz_3m);

    gPad->SetLeftMargin(.15);
    //gPad->SetLogz();
    dalitz_0p->SetTitle(
        "^{12}C(0^{+}_{2})#rightarrow#alpha#alpha#alpha;#sqrt{3}(#varepsilon_{2} - #varepsilon_{1});2#varepsilon_{3} - #varepsilon_{2} - #varepsilon_{1}");
    //dalitz_0p->Rebin2D(8);
    //dalitz_0p->RebinY(4);
    //dalitz_0p->Rebin2D(4);
    //dalitz_0p->Smooth();
    //dalitz_0p->GetYaxis()->SetRangeUser(-.6, .6);
    //dalitz_0p->GetXaxis()->SetRangeUser(-.6, .6);
    //
    dalitz_0p->SetMarkerSize(.25);
    dalitz_0p->SetMarkerColor(kBlack);
    dalitz_0p->SetMarkerStyle(20);

    auto* cirle = new TEllipse(0, 0, 1.);
    auto* cirle_in = new TEllipse(0, 0, 0.1);
    //auto * cirle  = new TEllipse(0, 0, 1./ 3.);
    cirle->SetFillStyle(0);
    cirle->SetLineColor(kBlack);
    cirle->SetLineWidth(5);

    cirle_in->SetFillStyle(0);
    cirle_in->SetLineColor(kBlack);
    cirle_in->SetLineWidth(5);

    dalitz_0p->Draw();
    //dalitz_0p->Draw("col");
    //dalitz_0p->SetMinimum(-1);
    cirle->Draw("same");
    cirle_in->Draw("same");
    canvas->SaveAs("dalitz_0p.pdf");
    canvas->SaveAs("dalitz_0p.png");

    dalitz_3m->SetTitle(
        "^{12}C(3^{-})#rightarrow#alpha#alpha#alpha;#sqrt{3}(#varepsilon_{2} - #varepsilon_{1});2#varepsilon_{3} - #varepsilon_{2} - #varepsilon_{1}");
    //dalitz_3m->Rebin2D(8);
    //dalitz_3m->RebinY(4);
    //dalitz_3m->GetYaxis()->SetRangeUser(-, .6);
    //dalitz_3m->GetXaxis()->SetRangeUser(-.6, .6);

    dalitz_3m->Draw();
    dalitz_3m->SetMarkerSize(.05);
    dalitz_3m->SetMarkerColor(kBlack);
    dalitz_3m->SetMarkerStyle(20);
    cirle->Draw("same");
    canvas->SaveAs("dalitz_3m.pdf");
    canvas->SaveAs("dalitz_3m.png");

    dalitz_all->SetTitle(
        "^{12}C#rightarrow#alpha#alpha#alpha;#frac{#sqrt{3}}{2}(#varepsilon_{2} - #varepsilon_{1});#frac{1}{2}(2#varepsilon_{3} - #varepsilon_{2} - #varepsilon_{1})");
    dalitz_all->Rebin2D(8);
    dalitz_all->RebinY(4);
    //dalitz_all->GetYaxis()->SetRangeUser(-.6, .6);
    //dalitz_all->GetXaxis()->SetRangeUser(-.6, .6);

    dalitz_all->Draw("col");
    cirle->Draw("same");
    canvas->SaveAs("dalitz_all.pdf");
    canvas->SaveAs("dalitz_all.png");

    formatting::preliminary_formatting(dalitz_low);
    dalitz_low->SetTitle(
        "^{12}C(low)#rightarrow#alpha#alpha#alpha;#frac{#sqrt{3}}{2}(#varepsilon_{2} - #varepsilon_{1});#frac{1}{2}(2#varepsilon_{3} - #varepsilon_{2} - #varepsilon_{1})");
    dalitz_low->Rebin2D(8);
    dalitz_low->RebinY(4);
    //dalitz_low->GetYaxis()->SetRangeUser(-.6, .6);
    //dalitz_low->GetXaxis()->SetRangeUser(-.6, .6);

    dalitz_low->Draw("col");
    cirle->Draw("same");
    canvas->SaveAs("dalitz_low.pdf");
    canvas->SaveAs("dalitz_low.png");
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
