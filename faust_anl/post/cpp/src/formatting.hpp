#pragma once

#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>

namespace formatting {
void preliminary_formatting(TH2D* hist);
void preliminary_formatting(TH1D* hist, int idx = 1);
void preliminary_formatting(TGraph* graph, int idx = 1);

double integrate_hist(TH1D* hist, double low, double high);

void set_style();

TH1D* generate_normalized(
    TH1D* standard,
    TH1D* to_be_scaled,
    double low,
    double high);

}  // namespace formatting
