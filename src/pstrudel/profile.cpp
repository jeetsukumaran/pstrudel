#include <algorithm>
#include <colugo/utility.hpp>
#include "profile.hpp"

namespace pstrudel {

//////////////////////////////////////////////////////////////////////////////
// Profile (public)

Profile::Profile(unsigned long fixed_size,
        Profile::InterpolationMethod interpolation_method)
    : fixed_size_(fixed_size)
      , interpolation_method_(interpolation_method)
      , last_profile_comparison_size_(0) {
}

void Profile::clear() {
    this->raw_data_.clear();
    this->interpolated_profiles_.clear();
}

std::vector<double> & Profile::get_profile(unsigned long profile_size) {
    if (profile_size == 0) {
        profile_size = this->raw_data_.size();
    }
    auto v1_iter = this->interpolated_profiles_.find(profile_size);
    if (v1_iter == this->interpolated_profiles_.end()) {
        this->build_interpolated_profile(profile_size);
    }
    return this->interpolated_profiles_[profile_size];
}

double Profile::get_distance(Profile & other, bool normalize_by_profile_size) {
    unsigned long profile_size = this->get_profile_comparison_size(other);
    double dist = this->calc_distance(other, profile_size);
    if (normalize_by_profile_size) {
        return dist/profile_size;
    } else {
        return dist;
    }
}

//////////////////////////////////////////////////////////////////////////////
// Profile (protected/private)

double Profile::calc_distance(Profile & other, unsigned long profile_size) {
    auto & v1 = this->get_profile(profile_size);
    auto & v2 = other.get_profile(profile_size);
    this->last_profile_comparison_size_ = profile_size;
    other.last_profile_comparison_size_ = profile_size;
    return Profile::calc_euclidean_distance(v1, v2);
}


unsigned long Profile::get_profile_comparison_size(Profile & other) {
    unsigned long profile_size = 0;
    if (this->fixed_size_ > 0 && other.fixed_size_ > 0) {
        // fixed size profiles
        if (this->fixed_size_ != other.fixed_size_) {
            colugo::colugo_abort("Comparing two profiles locked to different sizes: ",
                    this->fixed_size_, " and ", other.fixed_size_);
            exit(1);
        }
        profile_size = this->fixed_size_;
    } else if (this->fixed_size_ == 0 && other.fixed_size_ > 0) {
        if (other.fixed_size_ < this->raw_data_.size()) {
            colugo::colugo_abort("Cannot interpolate points in current profile:",
                    " raw data size is ", this->raw_data_.size(),
                    " but comparison requires ", other.fixed_size_,
                    " because other profile is locked to fixed size");
            exit(1);
        }
        profile_size = other.fixed_size_;
    } else if (this->fixed_size_ > 0 && other.fixed_size_ == 0) {
        if (this->fixed_size_ < other.raw_data_.size()) {
            colugo::colugo_abort("Cannot interpolate points in other profile:",
                    " other raw data size is ", other.raw_data_.size(),
                    " but comparison requires ", this->fixed_size_,
                    " because this profile is locked to fixed size");
            exit(1);
        }
        profile_size = this->fixed_size_;
    } else {
        profile_size = std::max(this->raw_data_.size(), other.raw_data_.size());
    }
    return profile_size;
}

void Profile::build_interpolated_profile(unsigned long profile_size) {
    std::vector<double> & interpolated_profile = this->interpolated_profiles_[profile_size];
    interpolated_profile.clear();
    unsigned long raw_data_size = this->raw_data_.size();
    COLUGO_ASSERT(raw_data_size > 1);
    if (profile_size < raw_data_size) {
        colugo::colugo_abort("Error interpolating points in profile ",
                ": Number of requested interpolated points (",
                profile_size,
                ") less than raw data size (",
                raw_data_size,
                ")");
    }
    unsigned long default_bin_size = static_cast<unsigned long>(profile_size / (raw_data_size));
    COLUGO_ASSERT(default_bin_size > 0);
    std::vector<unsigned long> bin_sizes(raw_data_size, default_bin_size);

    // due to rounding error, default bin size may not be enough
    // this hacky algorithm adds additional points to bins to make up the
    // difference, trying to distribute the additional points along the
    // length of the line
    unsigned long diff = profile_size - (default_bin_size * (raw_data_size));
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

    interpolated_profile.reserve(profile_size + 1);
    if (this->interpolation_method_ == Profile::InterpolationMethod::STAIRCASE) {
        for (auto & rd : this->raw_data_) {
            this->interpolate_flat(interpolated_profile, rd, default_bin_size);
        }
    } else {
        this->interpolate_linear(interpolated_profile, 0, 0, this->raw_data_[0], bin_sizes[0]);
        for (unsigned long bin_idx=0; bin_idx < raw_data_size - 1; ++bin_idx) {
            this->interpolate_linear(interpolated_profile, bin_idx, this->raw_data_[bin_idx], this->raw_data_[bin_idx+1], bin_sizes[bin_idx], profile_size);
        }
    }
    interpolated_profile.push_back(this->raw_data_.back());
}

void Profile::interpolate_flat(std::vector<double> & interpolated_profile, double value, unsigned long num_points) {
    while (num_points > 0) {
        interpolated_profile.push_back(value);
        --num_points;
    }
}
void Profile::interpolate_linear(std::vector<double> & interpolated_profile, double x1, double y1, double y2, unsigned long num_points, unsigned long max_points) {
    // pstrudel_print(std::cout, "\n---\n", "x1 = ", x1, ", y1 = ", y1, ", y2 = ", y2, "; n=", num_points, "\n");
    COLUGO_ASSERT(num_points > 0);
    double slope = (y2 - y1)/num_points;
    COLUGO_ASSERT(slope >= 0);
    for (unsigned long xi = 0; xi < num_points; ++xi) {
        if (max_points > 0 && interpolated_profile.size() >= max_points) {
            return;
        }
        interpolated_profile.push_back(slope * xi + y1);
        // pstrudel_print(std::cout, xi, ": ", interpolated_profile.back(), "\n");
    }
}

double Profile::calc_euclidean_distance(std::vector<double> & v1, std::vector<double> & v2) {
    unsigned long v1_size = v1.size();
    unsigned long v2_size = v2.size();
    unsigned long v1_idx = 0;
    unsigned long v2_idx = 0;
    if (v1_size > v2_size) {
        v1_idx = v1_size - v2_size;
    }
    if (v2_size > v1_size) {
        v2_idx = v2_size - v1_size;
    }
    double ss = 0.0;
    while (v1_idx < v1_size && v2_idx < v2_size) {
        ss += std::pow(v1[v1_idx] - v2[v2_idx], 2);
        ++v1_idx;
        ++v2_idx;
    }
    return std::sqrt(ss);
}










// void Profile::build_subprofiles() {
//     for (auto & raw_data_set : this->raw_data_) {
//         const std::string& raw_data_id = raw_data_set.first;
//         unsigned long num_interpolated_points = this->num_interpolated_points_[raw_data_id];
//         std::vector<double> & raw_data_vector = raw_data_set.second;

//         this->subprofiles_[raw_data_id] = std::vector<double>();
//         auto & subprofile = this->subprofiles_[raw_data_id];

//         unsigned long raw_data_size = raw_data_vector.size();
//         COLUGO_ASSERT(raw_data_size > 1);
//         if (num_interpolated_points < raw_data_size) {
//             colugo::colugo_abort("Error building subprofile '",
//                     raw_data_id,
//                     "': Number of requested interpolated points (",
//                     num_interpolated_points,
//                     ") less than raw data size (",
//                     raw_data_size,
//                     ")");
//         }

//         unsigned long default_bin_size = static_cast<unsigned long>(num_interpolated_points / (raw_data_size));
//         // if (default_bin_size <= 0) {
//         //     std::cerr << "*************\n";
//         //     std::cerr << num_interpolated_points << "/" << raw_data_size << " = " << default_bin_size << std::endl;
//         //     std::cerr << "*************\n";
//         // }
//         COLUGO_ASSERT(default_bin_size > 0);
//         std::vector<unsigned long> bin_sizes(raw_data_size, default_bin_size);

//         // due to rounding error, default bin size may not be enough
//         // this hacky algorithm adds additional points to bins to make up the
//         // difference, trying to distribute the additional points along the
//         // length of the line
//         unsigned long diff = num_interpolated_points - (default_bin_size * (raw_data_size));
//         if (diff > 0) {
//             double dv = static_cast<double>(diff) / (raw_data_size);
//             double cv = 0;
//             for (auto & bin : bin_sizes) {
//                 if (diff <= 1) {
//                     break;
//                 }
//                 cv += dv;
//                 if (cv >= 1) {
//                     bin += 1;
//                     diff -= 1;
//                     cv = cv - 1.0;
//                 }
//             }
//         }

//         subprofile.reserve(num_interpolated_points + 1);
//         if (this->interpolation_method_ == Profile::InterpolationMethod::STAIRCASE) {
//             for (auto & rd : raw_data_vector) {
//                 this->interpolate_flat(subprofile, rd, default_bin_size);
//             }
//         } else {
//             this->interpolate_linear(subprofile, 0, 0, raw_data_vector[0], bin_sizes[0]);
//             for (unsigned long bin_idx=0; bin_idx < raw_data_size - 1; ++bin_idx) {
//                 this->interpolate_linear(subprofile, bin_idx, raw_data_vector[bin_idx], raw_data_vector[bin_idx+1], bin_sizes[bin_idx], num_interpolated_points);
//             }
//             // this->interpolate_linear(subprofile, raw_data_size-1, raw_data_vector[raw_data_size-2], raw_data_vector.back(), bin_size);
//         }
//         subprofile.push_back(raw_data_vector.back());
//     }
// }

// const Profile::std::vector<double>& Profile::get_subprofile(const std::string& subprofile_id) {
//     COLUGO_ASSERT(this->subprofiles_.find(subprofile_id) != this->subprofiles_.end());
//     return this->subprofiles_[subprofile_id];
// }

// double Profile::calc_subprofile_distance(const Profile& other, const std::string& subprofile_id) const {
//     return std::sqrt(this->sum_of_squares(other, subprofile_id));
// }

// double Profile::calc_distance(const Profile& other) const {
//     COLUGO_ASSERT(this->subprofiles_.size() == other.subprofiles_.size());
//     double ss = 0.0;
//     for (auto & subprofile_pair : this->subprofiles_) {
//         const std::string& subprofile_id = subprofile_pair.first;
//         ss += this->sum_of_squares(other, subprofile_id);
//     }
//     return std::sqrt(ss);
// }

// void Profile::clear() {
//     this->raw_data_.clear();
//     this->subprofiles_.clear();
// }

// void Profile:: poll_max_num_interpolated_points(std::map<const std::string, unsigned long>& num_interpolated_points) {
//     for (auto & npp_iter : this->num_interpolated_points_) {
//         if (this->num_interpolated_points_[npp_iter.first] > num_interpolated_points[npp_iter.first]) {
//             num_interpolated_points[npp_iter.first] = this->num_interpolated_points_[npp_iter.first];
//         }
//     }
// }

// void Profile:: set_num_interpolated_points(std::map<const std::string, unsigned long>& num_interpolated_points) {
//     for (auto & npp_iter : this->num_interpolated_points_) {
//         this->num_interpolated_points_[npp_iter.first] = num_interpolated_points[npp_iter.first];
//     }
// }


} // namespace pstrudel
