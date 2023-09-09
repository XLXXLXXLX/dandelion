#include <random>

#include <catch2/catch_amalgamated.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>

#include "../src/scene/object.h"
#include "../src/utils/math.hpp"
#include "../src/utils/formatter.hpp"

using Eigen::AngleAxisf;
using Eigen::Matrix4f;
using Eigen::Scaling;
using Eigen::Translation3f;
using Eigen::Vector3f;
using Eigen::Vector4f;
using std::default_random_engine;
using std::random_device;
using std::uniform_real_distribution;

constexpr float threshold = 1e-2f;

TEST_CASE("Transformation", "[basic]")
{
    // Generate random positions, scaling factors and rotation angles.
    random_device seed;
    default_random_engine engine(seed());
    uniform_real_distribution<float> position_coord(-100.0f, 100.0f);
    uniform_real_distribution<float> scaling_factor(-100.0f, 100.0f);
    uniform_real_distribution<float> x_angle(-180.0f, 180.0f);
    uniform_real_distribution<float> y_angle(-90.0f, 90.0f);
    uniform_real_distribution<float> z_angle(-180.0f, 180.0f);
    Object test_object("test object");

    for (int i = 0; i < 10; ++i) {
        test_object.center =
            Vector3f(position_coord(engine), position_coord(engine), position_coord(engine));
        test_object.scaling =
            Vector3f(scaling_factor(engine), scaling_factor(engine), scaling_factor(engine));
        const float x_rad    = radians(x_angle(engine));
        const float y_rad    = radians(y_angle(engine));
        const float z_rad    = radians(z_angle(engine));
        test_object.rotation = AngleAxisf(x_rad, Vector3f::UnitX()) *
                               AngleAxisf(y_rad, Vector3f::UnitY()) *
                               AngleAxisf(z_rad, Vector3f::UnitZ());
        INFO(fmt::format("center is {:.2f}", test_object.center));
        INFO(fmt::format("scaling is {:.2f}", test_object.scaling));
        INFO(fmt::format("rotation is {:.2f}",
                         Vector3f(degrees(x_rad), degrees(y_rad), degrees(z_rad))));

        // Use Eigen's transform type as reference.
        Matrix4f reference = (Translation3f(test_object.center) * test_object.rotation *
                              Scaling(test_object.scaling))
                                 .matrix();
        // The answer is given by Object::model.
        Matrix4f answer = test_object.model();
        INFO(fmt::format("the answer is{:>5.1f}", answer));
        INFO(fmt::format("but it should be{:>5.1f}", reference));
        // Difference between the answer matrix and the reference one is measured
        // by L2 norm.
        REQUIRE((reference - answer).norm() < threshold);
    }
}

TEST_CASE("Perspective Projection", "[basic]")
{
    const Vector3f position(3.0f, 4.0f, 5.0f);
    const Vector3f target(0.0f, 0.0f, 0.0f);
    constexpr float near_planes[5]   = {0.1f, 0.01f, 0.93f, 0.82f, 0.4f};
    constexpr float far_planes[5]    = {10.0f, 131.0f, 240.0f, 79.4f, 372.0f};
    constexpr float fov_ys[5]        = {45.0f, 23.0f, 68.2f, 70.5f, 81.8f};
    constexpr float aspect_ratios[5] = {1.33f, 1.67f, 0.98f, 0.72f, 1.24f};

    // References are generated by GLM's perspective projection (glm::perspective).
    Matrix4f references[5];
    references[0] << 1.8152, 0.0000, 0.0000, 0.0000, 0.0000, 2.4142, 0.0000, 0.0000, 0.0000, 0.0000,
        -1.0202, -0.2020, 0.0000, 0.0000, -1.0000, 0.0000;
    references[1] << 2.9432, 0.0000, 0.0000, 0.0000, 0.0000, 4.9152, 0.0000, 0.0000, 0.0000, 0.0000,
        -1.0002, -0.0200, 0.0000, 0.0000, -1.0000, 0.0000;
    references[2] << 1.5071, 0.0000, 0.0000, 0.0000, 0.0000, 1.4770, 0.0000, 0.0000, 0.0000, 0.0000,
        -1.0078, -1.8672, 0.0000, 0.0000, -1.0000, 0.0000;
    references[3] << 1.9652, 0.0000, 0.0000, 0.0000, 0.0000, 1.4150, 0.0000, 0.0000, 0.0000, 0.0000,
        -1.0209, -1.6571, 0.0000, 0.0000, -1.0000, 0.0000;
    references[4] << 0.9310, 0.0000, 0.0000, 0.0000, 0.0000, 1.1544, 0.0000, 0.0000, 0.0000, 0.0000,
        -1.0022, -0.8009, 0.0000, 0.0000, -1.0000, 0.0000;

    for (int i = 0; i < 5; ++i) {
        Camera camera(position, target, near_planes[i], far_planes[i], fov_ys[i], aspect_ratios[i]);
        Matrix4f answer           = camera.projection();
        const Matrix4f& reference = references[i];
        INFO(fmt::format("near: {:.3f}, far: {:.1f}, fov_y: {:.1f} deg, aspect ratio: {:.2f}",
                         near_planes[i], far_planes[i], fov_ys[i], aspect_ratios[i]));
        INFO(fmt::format("the answer is{:>5.1f}", answer));
        INFO(fmt::format("but it should be{:>5.1f}", reference));
        REQUIRE((reference - answer).norm() < threshold);
    }
}
