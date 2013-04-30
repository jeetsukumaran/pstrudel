#ifndef PSTRUDEL_PROFILE_PROFILE_HPP
#define PSTRUDEL_PROFILE_PROFILE_HPP

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
#include <map>
#include <iterator>
#include <colugo/utility.hpp>

namespace pstrudel {

typedef double                                  ProfileMetricValueType;
typedef std::vector<ProfileMetricValueType>     ProfileMetricVectorType;

class Profile {
    public:

        enum class InterpolationMethod {
            STAIRCASE,
            PIECEWISE_LINEAR,
        };

    public:

        // if fixed_size == 0, then profiles are dynamically interpolated when comparing
        // profiles to the maximum size of the profiles being compared
        Profile(unsigned long fixed_size=0,
                Profile::InterpolationMethod interpolation_method=Profile::InterpolationMethod::PIECEWISE_LINEAR);
        Profile(const Profile & other);
        Profile & operator=(const Profile & other);

        template <class iterator>
        void set_data(
                iterator src_begin,
                iterator src_end,
                bool sort=true) {
            this->raw_data_.clear();
            this->interpolated_profiles_.clear();
            this->raw_data_.insert(this->raw_data_.end(), src_begin, src_end);
            if (this->fixed_size_ > 0 && (this->raw_data_.size() > this->fixed_size_)) {
                colugo::colugo_abort("Profile is fixed to ",
                        this->fixed_size_,
                        "but trying to add ",
                        this->raw_data_.size(),
                        " points as raw data");
                exit(1);
            }
            if (sort) {
                std::sort(this->raw_data_.begin(), this->raw_data_.end());
            }
            // // store raw data as profile for `n` points
            // auto & default_profile = this->interpolated_profiles_[this->raw_data_.size()];
            // default_profile.insert(this->raw_data_.end(), src_begin, src_end);
            // if (this->fixed_size_ > 0 && this->fixed_size_ != this->raw_data_.size()) {
            //     this->build_interpolated_profile(this->fixed_size_);
            // }
        }

        void clear();
        unsigned long data_size() const {
            return this->raw_data_.size();
        }
        ProfileMetricVectorType & get_profile(unsigned long profile_size=0);
        ProfileMetricValueType get_distance(Profile & other, bool normalize_by_profile_size=false);
        void fix_size(unsigned long size) {
            if (size > 0 && size < this->raw_data_.size()) {
                throw std::runtime_error("Cannot fix size to less than data size");
            }
            this->fixed_size_ = size;
        }
        void unfix_size() {
            this->fixed_size_ = 0;
        }

    private:
        ProfileMetricValueType calc_distance(Profile & other, unsigned long profile_size);
        unsigned long get_profile_comparison_size(Profile & other);
        void build_interpolated_profile(unsigned long profile_size);
        static void interpolate_flat(ProfileMetricVectorType & subprofile, ProfileMetricValueType value, unsigned long num_points);
        static void interpolate_linear(ProfileMetricVectorType & subprofile, ProfileMetricValueType x1, ProfileMetricValueType y1, ProfileMetricValueType y2, unsigned long num_points, unsigned long max_points=0);
        static ProfileMetricValueType calc_euclidean_distance(ProfileMetricVectorType & v1, ProfileMetricVectorType & v2);

    private:
        unsigned long                                       fixed_size_;
        Profile::InterpolationMethod                        interpolation_method_;
        ProfileMetricVectorType                             raw_data_;
        std::map<unsigned long, ProfileMetricVectorType>    interpolated_profiles_;
        unsigned long                                       last_profile_comparison_size_;

}; // Profile

} // namespace pstrudel

#endif
