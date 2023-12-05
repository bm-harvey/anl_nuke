#include <iostream>
#include "formatting.hpp"

namespace formatting {
void preliminary_formatting(TH1D* hist, int style_idx) {
    hist->SetTitleFont(112, "hist");
    hist->SetLineColor(style_idx);
    hist->SetMarkerColor(style_idx);
    hist->SetLineWidth(1);
    hist->SetMarkerSize(.5);
    hist->SetMarkerStyle(20);

    hist->GetXaxis()->SetNdivisions(108);
    hist->GetYaxis()->SetNdivisions(105);

    hist->GetXaxis()->SetLabelSize(0.03f);
    hist->GetYaxis()->SetLabelSize(0.03f);

    //hist->GetXaxis()->SetTitleOffset(1.1f);
    hist->GetYaxis()->SetTitleOffset(1.7f);
    hist->GetXaxis()->SetTitleSize(0.03f);
    hist->GetYaxis()->SetTitleSize(0.03f);
    hist->GetXaxis()->CenterTitle();
    hist->GetYaxis()->CenterTitle();

    hist->GetXaxis()->SetLabelFont(112);
    hist->GetYaxis()->SetLabelFont(112);
    hist->GetXaxis()->SetTitleFont(112);
    hist->GetYaxis()->SetTitleFont(112);
    hist->GetYaxis()->SetMaxDigits(3);
    hist->GetXaxis()->SetMaxDigits(3);
}
void preliminary_formatting(TGraph* graph, int style_idx) {
    graph->SetLineColor(style_idx);
    graph->SetMarkerColor(style_idx);
    graph->SetLineWidth(1);
    graph->SetMarkerSize(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(.5);

    graph->GetXaxis()->SetNdivisions(108);
    graph->GetYaxis()->SetNdivisions(105);

    graph->GetXaxis()->SetLabelSize(0.03f);
    graph->GetYaxis()->SetLabelSize(0.03f);

    //hist->GetXaxis()->SetTitleOffset(1.1f);
    graph->GetYaxis()->SetTitleOffset(1.5f);
    graph->GetXaxis()->SetTitleSize(0.03f);
    graph->GetYaxis()->SetTitleSize(0.03f);
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->CenterTitle();

    graph->GetXaxis()->SetLabelFont(112);
    graph->GetYaxis()->SetLabelFont(112);
    graph->GetXaxis()->SetTitleFont(112);
    graph->GetYaxis()->SetTitleFont(112);

}
void preliminary_formatting(TH2D* hist) {
    gStyle->SetTitleFont(112, "hist");
    hist->SetLineWidth(2);
    hist->SetMarkerSize(2);

    hist->GetXaxis()->SetNdivisions(108);
    hist->GetYaxis()->SetNdivisions(108);

    hist->SetTitle("");

    hist->GetXaxis()->SetLabelSize(0.03f);
    hist->GetYaxis()->SetLabelSize(0.03f);

    hist->GetYaxis()->SetTitleOffset(2.5f);
    hist->GetXaxis()->SetTitleOffset(1.1f);
    hist->GetXaxis()->SetTitleSize(0.03f);
    hist->GetYaxis()->SetTitleSize(0.03f);
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

double integrate_hist(TH1D* hist, double low, double high) {
    int bin_low = hist->GetXaxis()->FindBin(low);
    int bin_high = hist->GetXaxis()->FindBin(high);

    return hist->Integral(bin_low, bin_high);
}

TH1D* generate_normalized(
    TH1D* standard,
    TH1D* to_be_scaled,
    double low,
    double high) {
    auto result = (TH1D*)to_be_scaled->Clone();

    double norm_standard = integrate_hist(standard, low, high);
    double norm = integrate_hist(result, low, high);

    std::cout << norm << std::endl;
    std::cout << norm_standard << std::endl;

    result->Scale(
    norm_standard * standard->GetBinWidth(1) / norm
    / to_be_scaled->GetBinWidth(1));

    return result;
}

void set_style(){
    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(255);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
}

}  // namespace formatting
