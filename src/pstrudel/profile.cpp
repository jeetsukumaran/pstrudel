#include <colugo/utility.hpp>
#include "profile.hpp"

namespace pstrudel {

Profile::Profile(Profile::InterpolationMethod interpolation_method)
        : interpolation_method_(interpolation_method) {
}

void Profile::build_subprofiles() {
    for (auto & raw_data_set : this->raw_data_) {
        const std::string& raw_data_id = raw_data_set.first;
        unsigned long num_interpolated_points = this->num_interpolated_points_[raw_data_id];
        SubprofileDataType & raw_data_vector = raw_data_set.second;

        this->subprofiles_[raw_data_id] = SubprofileDataType();
        auto & subprofile = this->subprofiles_[raw_data_id];

        unsigned long raw_data_size = raw_data_vector.size();
        COLUGO_ASSERT(raw_data_size > 1);
        if (num_interpolated_points < raw_data_size) {
            colugo::colugo_abort("Error building subprofile '",
                    raw_data_id,
                    "': Number of requested interpolated points (",
                    num_interpolated_points,
                    ") less than raw data size (",
                    raw_data_size,
                    ")");
        }

        unsigned long default_bin_size = static_cast<unsigned long>(num_interpolated_points / (raw_data_size));
        // if (default_bin_size <= 0) {
        //     std::cerr << "*************\n";
        //     std::cerr << num_interpolated_points << "/" << raw_data_size << " = " << default_bin_size << std::endl;
        //     std::cerr << "*************\n";
        // }
        COLUGO_ASSERT(default_bin_size > 0);
        std::vector<unsigned long> bin_sizes(raw_data_size, default_bin_size);

        // due to rounding error, default bin size may not be enough
        // this hacky algorithm adds additional points to bins to make up the
        // difference, trying to distribute the additional points along the
        // length of the line
        unsigned long diff = num_interpolated_points - (default_bin_size * (raw_data_size));
        if (diff > 0) {
            double dv = static_cast<double>(diff) / (raw_data_size);
            double cv = 0;
            for (auto & bin : bin_sizes) {
                if (diff <= 1) {
                    break;
                }
                cv += dv;
                if (cv >= 1) {
                    bin += 1;
                    diff -= 1;
                    cv = cv - 1.0;
                }
            }
        }

        subprofile.reserve(num_interpolated_points + 1);
        if (this->interpolation_method_ == Profile::InterpolationMethod::STAIRCASE) {
            for (auto & rd : raw_data_vector) {
                this->interpolate_flat(subprofile, rd, default_bin_size);
            }
        } else {
            this->interpolate_linear(subprofile, 0, 0, raw_data_vector[0], bin_sizes[0]);
            for (unsigned long bin_idx=0; bin_idx < raw_data_size - 1; ++bin_idx) {
                this->interpolate_linear(subprofile, bin_idx, raw_data_vector[bin_idx], raw_data_vector[bin_idx+1], bin_sizes[bin_idx], num_interpolated_points);
            }
            // this->interpolate_linear(subprofile, raw_data_size-1, raw_data_vector[raw_data_size-2], raw_data_vector.back(), bin_size);
        }
        subprofile.push_back(raw_data_vector.back());
    }
}

const Profile::SubprofileDataType& Profile::get_subprofile(const std::string& subprofile_id) {
    COLUGO_ASSERT(this->subprofiles_.find(subprofile_id) != this->subprofiles_.end());
    return this->subprofiles_[subprofile_id];
}

double Profile::calc_subprofile_distance(const Profile& other, const std::string& subprofile_id) const {
    return std::sqrt(this->sum_of_squares(other, subprofile_id));
}

double Profile::calc_distance(const Profile& other) const {
    COLUGO_ASSERT(this->subprofiles_.size() == other.subprofiles_.size());
    double ss = 0.0;
    for (auto & subprofile_pair : this->subprofiles_) {
        const std::string& subprofile_id = subprofile_pair.first;
        ss += this->sum_of_squares(other, subprofile_id);
    }
    return std::sqrt(ss);
}

void Profile::clear() {
    this->raw_data_.clear();
    this->subprofiles_.clear();
}

void Profile:: poll_max_num_interpolated_points(std::map<const std::string, unsigned long>& num_interpolated_points) {
    for (auto & npp_iter : this->num_interpolated_points_) {
        if (this->num_interpolated_points_[npp_iter.first] > num_interpolated_points[npp_iter.first]) {
            num_interpolated_points[npp_iter.first] = this->num_interpolated_points_[npp_iter.first];
        }
    }
}

void Profile:: set_num_interpolated_points(std::map<const std::string, unsigned long>& num_interpolated_points) {
    for (auto & npp_iter : this->num_interpolated_points_) {
        this->num_interpolated_points_[npp_iter.first] = num_interpolated_points[npp_iter.first];
    }
}

void Profile::interpolate_flat(SubprofileDataType& subprofile, double value, unsigned long num_points) {
    while (num_points > 0) {
        subprofile.push_back(value);
        --num_points;
    }
}
void Profile::interpolate_linear(SubprofileDataType& subprofile, double x1, double y1, double y2, unsigned long num_points, unsigned long max_points) {
    // pstrudel_print(std::cout, "\n---\n", "x1 = ", x1, ", y1 = ", y1, ", y2 = ", y2, "; n=", num_points, "\n");
    COLUGO_ASSERT(num_points > 0);
    double slope = (y2 - y1)/num_points;
    COLUGO_ASSERT(slope >= 0);
    for (unsigned long xi = 0; xi < num_points; ++xi) {
        if (max_points > 0 && subprofile.size() >= max_points) {
            return;
        }
        subprofile.push_back(slope * xi + y1);
        // pstrudel_print(std::cout, xi, ": ", subprofile.back(), "\n");
    }
}

double Profile::sum_of_squares(const Profile& other, const std::string& subprofile_id) const {
    auto sp_pair1 = this->subprofiles_.find(subprofile_id);
    auto sp_pair2 = other.subprofiles_.find(subprofile_id);
    if (sp_pair1 == this->subprofiles_.end()) {
        colugo::colugo_abort("Subprofile '", subprofile_id, "' not found");
        exit(1);
    }
    if (sp_pair2 == other.subprofiles_.end()) {
        colugo::colugo_abort("Subprofile '", subprofile_id, "' not found");
        exit(1);
    }
    auto & subprofile1 = sp_pair1->second;
    auto & subprofile2 = sp_pair2->second;
    unsigned long sp1_size = subprofile1.size();
    unsigned long sp2_size = subprofile2.size();
    unsigned long sp1_idx = 0;
    unsigned long sp2_idx = 0;
    if (sp1_size > sp2_size) {
        sp1_idx = sp1_size - sp2_size;
    }
    if (sp2_size > sp1_size) {
        sp2_idx = sp2_size - sp1_size;
    }
    double ss = 0.0;
    unsigned long num_comparisons = 0;
    while (sp1_idx < sp1_size && sp2_idx < sp2_size) {
        ss += std::pow(subprofile1[sp1_idx] - subprofile2[sp2_idx], 2);
        ++sp1_idx;
        ++sp2_idx;
        ++num_comparisons;
    }
    // std::fstream x;
    // x.open("x.txt", std::fstream::out | std::fstream::app);
    // x << num_comparisons << std::endl;
    return ss;
}

} // namespace pstrudel
