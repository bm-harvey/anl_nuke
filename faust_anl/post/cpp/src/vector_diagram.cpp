#include <Rtypes.h>
#include <TArrow.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TStyle.h>

#include "formatting.hpp"

struct Point {
    double x;
    double y;
};

void fmt_arrow(TArrow* arrow, int color = kBlack) {
    arrow->SetLineColor(color);
    arrow->SetLineWidth(5);
    arrow->SetArrowSize(0.015f);
    arrow->SetAngle(40);
    arrow->SetFillColor(color);
}

void fmt_txt(TPaveText* txt, int color = kBlack) {
    txt->SetTextFont(112);
    txt->SetTextSize(0.02f);
    txt->SetTextColor(color);
    txt->SetFillColor(kWhite);
}

int main() {
    gStyle->SetOptStat(0);

    auto* canvas = new TCanvas("canvas", "canvas", 1200, 1200);
    gPad->SetFrameLineWidth(2);
    auto* slate = new TH2D("slate", ";;", 10, -1, .8, 10, -.6, .8);
    slate->Draw();
    formatting::preliminary_formatting(slate);
    slate->GetXaxis()->SetNdivisions(0);
    slate->GetYaxis()->SetNdivisions(0);

    auto* line = new TLine(0, 0, 0.5, 0.5);
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    //line->Draw();

    auto origin = Point {-.4, 0};
    auto a1 = Point {0.22f, 0.3f};
    auto a2 = Point {-0.05f, -0.5f};
    auto a3 = Point {0.5f, -0.01f};

    auto* arrow_a1 = new TArrow(origin.x, origin.y, a1.x, a1.y, 0.015f, "|>");
    auto* txt_a1 = new TPaveText(0.08f, 0.3f, 0.08f, 0.3f, "NB");
    auto* arrow_a2 = new TArrow(origin.x, origin.y, a2.x, a2.y, 0.015f, "|>");
    auto* txt_a2 = new TPaveText(-0.2f, -0.45f, -0.2f, -0.45f, "NB");
    auto* arrow_a3 = new TArrow(origin.x, origin.y, a3.x, a3.y, 0.015f, "|>");
    auto* txt_a3 = new TPaveText(0.5f, 0.025f, 0.5f, 0.025f, "NB");

    auto* txt_tgt =
        new TPaveText(origin.x, origin.y + .15, origin.x, origin.y + .15, "NB");
    txt_tgt->AddText("Si Target");

    auto* arrow_beam = new TArrow(-.75, 0, origin.x, origin.y, 0.015f, "|>");
    auto* txt_beam = new TPaveText(-.77, 0.05, -.77, 0.05, "NB");
    txt_beam->AddText("^{12}C Beam");

    auto c12 = Point {(a1.x + a2.x + a3.x) / 3, (a1.y + a2.y + a3.y) / 3.};
    auto* arrow_12c =
        new TArrow(origin.x, origin.y, c12.x, c12.y, 0.025f, "-|>-");

    auto* arrow_a1_c = new TArrow(c12.x, c12.y, a1.x, a1.y, 0.015f, "|>");
    auto* arrow_a2_c = new TArrow(c12.x, c12.y, a2.x, a2.y, 0.015f, "|>");
    auto* arrow_a3_c = new TArrow(c12.x, c12.y, a3.x, a3.y, 0.015f, "|>");

    auto* txt_c12 = new TPaveText(-0.15, -0.08, -0.15, -0.08, "NB");
    txt_c12->AddText("#vec{v}_{^{12}C^{*}}");

    auto* txt_a1_c = new TPaveText(0.33f, 0.2f, 0.33f, 0.2f, "NB");
    txt_a1_c->AddText("#vec{v}_{#alpha_{1},C^{*}}");
    auto* txt_a2_c = new TPaveText(0.1f, -0.45f, 0.1f, -0.45f, "NB");
    txt_a2_c->AddText("#vec{v}_{#alpha_{2},C^{*}}");
    auto* txt_a3_c = new TPaveText(0.5f, -0.06f, 0.5f, -0.06f, "NB");
    txt_a3_c->AddText("#vec{v}_{#alpha_{3},C^{*}}");

    auto* txt_eqn = new TPaveText(-0.9f, 0.35f, -0.2f, 0.75f);
    txt_eqn->AddText(
        "#vec{v}_{#alpha_{i},C^{*}} = #vec{v}_{#alpha_{i},lab} - #vec{v}_{^{12}C^{*}}");
    txt_eqn->AddText(
        "E_{rel.} = #sum_{i} #frac{1}{2}m_{#alpha}v_{#alpha_{i}}^{2}");
    txt_eqn->AddText("E^{*} = E_{rel.} - Q");

    auto* target = new TEllipse(origin.x, origin.y, 0.04f, 0.1f);
    target->SetFillStyle(0);
    target->SetLineColor(kBlue + 1);
    target->SetLineWidth(2);
    target->Draw();

    auto* legend = new TLegend(0.1, 0.9, 0.9, 0.95);
    legend->SetNColumns(3);
    legend->SetTextFont(112);
    legend->SetLineWidth(2);
    legend->AddEntry(arrow_beam, "Experimental Setup", "l");
    legend->AddEntry(arrow_a1, "Measured", "l");
    legend->AddEntry(arrow_a1_c, "Inferred", "l");
    legend->Draw("same");


    txt_a1->AddText("#vec{v}_{#alpha_{1},lab}");
    txt_a2->AddText("#vec{v}_{#alpha_{2},lab}");
    txt_a3->AddText("#vec{v}_{#alpha_{3},lab}");

    fmt_arrow(arrow_a1);
    fmt_arrow(arrow_a2);
    fmt_arrow(arrow_a3);
    fmt_arrow(arrow_a1_c, kRed + 1);
    fmt_arrow(arrow_a2_c, kRed + 1);
    fmt_arrow(arrow_a3_c, kRed + 1);
    fmt_arrow(arrow_beam, kBlue + 1);
    fmt_arrow(arrow_12c, kRed + 1);
    arrow_12c->SetLineStyle(2);
    arrow_12c->SetLineWidth(3);

    fmt_txt(txt_a1);
    fmt_txt(txt_a2);
    fmt_txt(txt_a3);
    fmt_txt(txt_eqn);
    fmt_txt(txt_c12, kRed + 1);
    fmt_txt(txt_tgt, kBlue + 1);
    fmt_txt(txt_beam, kBlue + 1);
    fmt_txt(txt_a1_c, kRed + 1);
    fmt_txt(txt_a2_c, kRed + 1);
    fmt_txt(txt_a3_c, kRed + 1);

    arrow_a1->Draw();
    arrow_a2->Draw();
    arrow_a3->Draw();
    arrow_a1_c->Draw();
    arrow_a2_c->Draw();
    arrow_a3_c->Draw();
    arrow_12c->Draw();
    arrow_beam->Draw();

    txt_a1->Draw();
    txt_a2->Draw();
    txt_a3->Draw();
    txt_tgt->Draw();
    txt_beam->Draw();
    txt_c12->Draw();
    txt_a1_c->Draw();
    txt_a2_c->Draw();
    txt_a3_c->Draw();
    txt_eqn->Draw();

    canvas->SaveAs("vector_diagram.png");
}
