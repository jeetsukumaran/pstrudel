
#ifndef PSTRUDEL_PROFILE_PROFILE_HPP
#define PSTRUDEL_PROFILE_PROFILE_HPP

#include <algorithm>
#include <map>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <iterator>

namespace pstrudel {

class Profile {

    public:
        typedef std::vector<double>                 SubprofileDataType;
        enum class InterpolationMethod {
            STAIRCASE,
            PIECEWISE_LINEAR,
        };

    public:
        Profile(Profile::InterpolationMethod interpolation_method=Profile::InterpolationMethod::PIECEWISE_LINEAR);

        template <class iterator>
        void add_subprofile_source(
                const std::string& raw_data_id,
                iterator src_begin,
                iterator src_end,
                unsigned long num_interpolated_points=1000,
                bool sort=true) {
            this->raw_data_[raw_data_id] = SubprofileDataType(src_begin, src_end);
            this->num_interpolated_points_[raw_data_id] = num_interpolated_points;
            if (sort) {
                auto & raw_data = this->raw_data_[raw_data_id];
                std::sort(raw_data.begin(), raw_data.end());
            }
        }

        template <class iterator>
        void insert_subprofile(
                const std::string& subprofile_id,
                iterator src_begin,
                iterator src_end) {
            this->subprofiles_[subprofile_id] = SubprofileDataType(src_begin, src_end);
        }

        void build_subprofiles();
        const SubprofileDataType& get_subprofile(const std::string& subprofile_id);
        double calc_subprofile_distance(const Profile& other, const std::string& subprofile_id) const;
        double calc_distance(const Profile& other) const;
        void clear();

        void poll_max_num_interpolated_points(std::map<const std::string, unsigned long>& num_interpolated_points);
        void set_num_interpolated_points(std::map<const std::string, unsigned long>& num_interpolated_points);

    private:
        void interpolate_flat(SubprofileDataType& subprofile, double value, unsigned long num_points);
        void interpolate_linear(SubprofileDataType& subprofile, double x1, double y1, double y2, unsigned long num_points, unsigned long max_points=0);
        double sum_of_squares(const Profile& other, const std::string& subprofile_id) const;

    private:
        Profile::InterpolationMethod                        interpolation_method_;
        std::map<const std::string, SubprofileDataType>     raw_data_;
        std::map<const std::string, unsigned long>          num_interpolated_points_;
        std::map<const std::string, SubprofileDataType>     subprofiles_;

}; // Profile

} // namespace pstrudel

#endif
