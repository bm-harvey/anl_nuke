#pragma once
#include <memory>
#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "rust/cxx.h"

namespace file {
class RtFile {
  public:
    RtFile(const char* name, const char* opt) {
        file_ = std::make_shared<TFile>(name, opt);
    }
    void close() const;
    void write() const;
    void cd() const;

  private:
    std::shared_ptr<TFile> file_ {};
};

auto open(rust::string name) -> std::shared_ptr<RtFile>;
auto create(rust::string name) -> std::shared_ptr<RtFile>;
};  // namespace file

namespace h1d {
class RtH1D {
  public:
    RtH1D(
        rust::string name,
        rust::string title,
        int bins,
        double low,
        double high) {
        hist_ = std::make_shared<TH1D>(
            name.c_str(),
            title.c_str(),
            bins,
            low,
            high);
    }

    void fill(double value) const {
        hist_->Fill(value);
    }
    void weighted_fill(double value, double weight) const {
        hist_->Fill(value, weight);
    }
    void write() const {
        hist_->Write();
    }

  private:
    std::shared_ptr<TH1D> hist_ {};
};

auto new_h1d(
    rust::string name,
    rust::string title,
    int64_t bins,
    double low,
    double high) -> std::shared_ptr<RtH1D>;
}  // namespace h1d

namespace h2d {
class RtH2D {
  public:
    RtH2D(
        rust::string name,
        rust::string title,
        int bins_x,
        double low_x,
        double high_x,
        int bins_y,
        double low_y,
        double high_y) {
        hist_ = std::make_shared<TH2D>(
            name.c_str(),
            title.c_str(),
            bins_x,
            low_x,
            high_x,
            bins_y,
            low_y,
            high_y);
    }

    void fill(double value_x, double value_y) const {
        hist_->Fill(value_x, value_y);
    }
    void weighted_fill(double value_x, double value_y, double weight) const {
        hist_->Fill(value_x, value_y, weight);
    }
    void write() const {
        hist_->Write();
    }

  private:
    std::shared_ptr<TH2D> hist_ {};
};

auto new_h2d(
    rust::string name,
    rust::string title,
    int64_t bins_x,
    double low_x,
    double high_x,
    int64_t bins_y,
    double low_y,
    double high_y) -> std::shared_ptr<RtH2D>;
}  // namespace h2d

namespace tree {
class RtIntBranch {
  public:
    void set_value(int64_t value) {
        data_ = value;
    }
    int64_t* address() {
        return &data_;
    }

  private:
    int64_t data_ {};
};
class RtBranch {
  public:
    void set_value(double value) {
        data_ = value;
    }
    double* address() {
        return &data_;
    }

  private:
    double data_ {};
};

class RtTree {
  public:
    RtTree(rust::string name, rust::string title) {
        tree_ = std::make_shared<TTree>(name.c_str(), title.c_str());
    }
    
    std::unique_ptr<RtBranch> make_branch(rust::string name) const {
        auto branch = std::make_unique<RtBranch>();

        tree_->Branch(name.c_str(), branch->address());

        return branch;
    }
    std::unique_ptr<RtIntBranch> make_int_branch(rust::string name) const {
        auto branch = std::make_unique<RtIntBranch>();

        tree_->Branch(name.c_str(), branch->address());

        return branch;
    }

    void fill() const {
        tree_->Fill();
    }
    void write() const {
        tree_->Write();
    }

  private:
    std::shared_ptr<TTree> tree_ {};
};

auto new_tree(rust::string name, rust::string title) -> std::shared_ptr<RtTree>;
}  // namespace tree
