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
        std::vector<double> & get_profile(unsigned long profile_size=0);
        double get_distance(Profile & other, bool normalize_by_profile_size=false);
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
        double calc_distance(Profile & other, unsigned long profile_size);
        unsigned long get_profile_comparison_size(Profile & other);
        void build_interpolated_profile(unsigned long profile_size);
        static void interpolate_flat(std::vector<double> & subprofile, double value, unsigned long num_points);
        static void interpolate_linear(std::vector<double> & subprofile, double x1, double y1, double y2, unsigned long num_points, unsigned long max_points=0);
        static double calc_euclidean_distance(std::vector<double> & v1, std::vector<double> & v2);

    private:
        unsigned long                                   fixed_size_;
        Profile::InterpolationMethod                    interpolation_method_;
        std::vector<double>                             raw_data_;
        std::map<unsigned long, std::vector<double>>    interpolated_profiles_;
        unsigned long                                   last_profile_comparison_size_;

}; // Profile

// class Profile {

//     public:
//         typedef std::vector<double>                 std::vector<double>;
//         enum class InterpolationMethod {
//             STAIRCASE,
//             PIECEWISE_LINEAR,
//         };

//     public:
//         Profile(Profile::InterpolationMethod interpolation_method=Profile::InterpolationMethod::PIECEWISE_LINEAR);

//         template <class iterator>
//         void add_subprofile_source(
//                 const std::string& raw_data_id,
//                 iterator src_begin,
//                 iterator src_end,
//                 unsigned long num_interpolated_points=1000,
//                 bool sort=true) {
//             this->raw_data_[raw_data_id] = std::vector<double>(src_begin, src_end);
//             this->num_interpolated_points_[raw_data_id] = num_interpolated_points;
//             if (sort) {
//                 auto & raw_data = this->raw_data_[raw_data_id];
//                 std::sort(raw_data.begin(), raw_data.end());
//             }
//         }

//         template <class iterator>
//         void insert_subprofile(
//                 const std::string& subprofile_id,
//                 iterator src_begin,
//                 iterator src_end) {
//             this->subprofiles_[subprofile_id] = std::vector<double>(src_begin, src_end);
//         }

//         void build_subprofiles();
//         const std::vector<double>& get_subprofile(const std::string& subprofile_id);
//         double calc_subprofile_distance(const Profile& other, const std::string& subprofile_id) const;
//         double calc_distance(const Profile& other) const;
//         void clear();

//         void poll_max_num_interpolated_points(std::map<const std::string, unsigned long>& num_interpolated_points);
//         void set_num_interpolated_points(std::map<const std::string, unsigned long>& num_interpolated_points);

//     private:
//         void interpolate_flat(std::vector<double>& subprofile, double value, unsigned long num_points);
//         void interpolate_linear(std::vector<double>& subprofile, double x1, double y1, double y2, unsigned long num_points, unsigned long max_points=0);
//         double sum_of_squares(const Profile& other, const std::string& subprofile_id) const;

//     private:
//         Profile::InterpolationMethod                        interpolation_method_;
//         std::map<const std::string, std::vector<double>>     raw_data_;
//         std::map<const std::string, unsigned long>          num_interpolated_points_;
//         std::map<const std::string, std::vector<double>>     subprofiles_;

// }; // Profile

} // namespace pstrudel

#endif
