#include "my_bridge.hpp"

#include <iostream>
#include <memory>
#include <string>

namespace file {
std::shared_ptr<RtFile> open(rust::string name) {
    return std::make_shared<RtFile>(name.c_str(), "read");
}

std::shared_ptr<RtFile> create(rust::string name) {
    return std::make_shared<RtFile>(name.c_str(), "recreate");
}

void RtFile::write() const {
    file_->Write();
}
void RtFile::cd() const {
    file_->cd();
}
void RtFile::close() const {
    file_->Close();
}
};  // namespace file

namespace h1d {

std::shared_ptr<RtH1D> new_h1d(
    rust::string name,
    rust::string title,
    int64_t bins,
    double low,
    double high) {
    return std::make_shared<RtH1D>(name, title, bins, low, high);
}
};  // namespace h1d

namespace h2d {

std::shared_ptr<RtH2D> new_h2d(
    rust::string name,
    rust::string title,
    int64_t bins_x,
    double low_x,
    double high_x,
    int64_t bins_y,
    double low_y,
    double high_y) {
    return std::make_shared<
        RtH2D>(name, title, bins_x, low_x, high_x, bins_y, low_y, high_y);
}

};  // namespace h2d
    //
namespace tree {

auto new_tree(rust::string name, rust::string title)
    -> std::shared_ptr<RtTree> {
    return std::make_shared<RtTree>(name.c_str(), title.c_str());
}
}  // namespace tree
