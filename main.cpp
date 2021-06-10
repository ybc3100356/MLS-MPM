#include <cstddef>
#include "algebra.h"
#include "particle.h"
#include <iostream>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

using namespace cv;

using std::vector;
using std::max;
using std::min;

const size_t width = 512;

const size_t dim = 2;
const size_t numParticles = 1024;
const size_t numGrid = 128;

const Real dx = 1.0 / numGrid;
const Real inv_dx = 1.0 / dx;
const Real dt = 2.0e-4;
const Real p_vol = (dx * 0.5) * (dx * 0.5);
const Real p_rho = 0.5;
const Real p_mass = p_vol * p_rho;
const Real E = 10;

//template<typename T>
//vec2<T> clamp_pos(vec2<T> &&pos) {
//    return vec2<T>(max(min(0.9, pos[0]), 0.1), max(min(0.9, pos[1]), 0.1));
//}

void step(vector<Particle> &particles, vector<vector<vec2<Real>>> &grid_v, vector<vector<Real>> &grid_m) {
    // p2g
    for (auto &particle : particles) {
        auto base = (particle.position * inv_dx - 0.5).cast<int>();
        auto fx = (particle.position * inv_dx) - base.cast<Real>();
        // quadratic B-spline weights
        vector<vec2<Real>> w = {0.5 * ((1.5 - fx) * (1.5 - fx)),
                                0.75 - ((fx - 1) * (fx - 1)),
                                0.5 * (fx - 0.5) * (fx - 0.5)};
        Real stress = -dt * p_vol * (particle.J - 1) * 4 * inv_dx * inv_dx * E;
        auto affine = mat2<Real>(vec2<Real>(stress, 0), vec2<Real>(0, stress)) + p_mass * particle.C;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                auto offset = vec2<int>(i, j);
                auto dpos = (offset.cast<Real>() - fx) * dx;
                auto weight = w[i][0] * w[j][1];
                auto index = base + offset;
                assert (index[0] < numGrid && index[1] < numGrid);
                grid_v[index[0]][index[1]] += weight * (p_mass * particle.velocity + affine.mul(dpos));
                grid_m[index[0]][index[1]] += weight * p_mass;
            }
        }
    }

    for (int i = 0; i < grid_m.size(); i++) {
        for (int j = 0; j < grid_m[0].size(); j++) {
            if (grid_m[i][j] > 0) {
                auto inv_m = 1.0 / grid_m[i][j];
                grid_v[i][j] = inv_m * grid_v[i][j];
                grid_v[i][j][1] -= dt * 9.8;
                auto bound = 3;
                if (i < bound && grid_v[i][j][0] < 0)
                    grid_v[i][j][0] = 0;
                if (i > numGrid - bound && grid_v[i][j][0] > 0)
                    grid_v[i][j][0] = 0;
                if (j < bound && grid_v[i][j][1] < 0)
                    grid_v[i][j][1] = 0;
                if (j > numGrid - bound && grid_v[i][j][1] > 0)
                    grid_v[i][j][1] = 0;
            }
        }
    }

    // g2p
    for (auto &particle : particles) {
        auto base = (particle.position * inv_dx - 0.5).cast<int>();
        auto fx = (particle.position * inv_dx) - base.cast<Real>();
        // quadratic B-spline weights
        vector<vec2<Real>> w = {0.5 * ((1.5 - fx) * (1.5 - fx)),
                                0.75 - ((fx - 1) * (fx - 1)),
                                0.5 * (fx - 0.5) * (fx - 0.5)};

        auto new_v = vec2<Real>(0, 0);
        auto new_C = mat2<Real>();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                auto offset = vec2<int>(i, j);
                auto index = base + offset;
                assert (index[0] < numGrid && index[1] < numGrid);
                auto weight = w[i][0] * w[j][1];
                auto dpos = (offset.cast<Real>() - fx) * dx;
                auto g_v = grid_v[index[0]][index[1]];
                new_v += weight * g_v;
                new_C += 4 * weight * g_v.outer_product(dpos) * inv_dx;
            }
        }

        particle.velocity = new_v;
        particle.position += particle.velocity * dt;// boundary
        particle.J *= 1 + dt * new_C.trace();
        particle.C = new_C;
    }
}

void test() {
    std::cout << "run test!" << std::endl;
    vec2<Real> v1{3, 7}, v2{17, 13};
    mat2<Real> m1(v1, v2), m2(v2, v1);
    v1 = m1.mul(v2);
    std::cout << v1[0] << "\t" << v1[1] << std::endl;
    std::cout << m1[0][0] << "\t" << m1[0][1] << std::endl;
    std::cout << m1[1][0] << "\t" << m1[1][1] << std::endl;
    exit(0);
}

int main() {
//    test();
    vector<Particle> particles(numParticles);
    vector<vector<vec2<Real>>> grid_v(numGrid, vector<vec2<Real>>(numGrid, vec2<Real>(0, 0)));
    vector<vector<Real>> grid_m(numGrid, vector<Real>(numGrid));
    vector<Mat> frames;
    int numFrames = 200;
    while (numFrames--) {
        for (int s = 0; s < 50; s++) {
            for (auto &m : grid_m)
                std::fill(m.begin(), m.end(), 0);
            for (auto &v : grid_v)
                std::fill(v.begin(), v.end(), vec2<Real>(0, 0));
            step(particles, grid_v, grid_m);
        }
        Mat img = Mat::zeros(width, width, CV_8UC3);
        for (auto particle : particles) {
            circle(img,
                   Point2d(particle.position[0] * width, width - particle.position[1] * width),
                   2,
                   Scalar(0, 0, 255),
                   FILLED,
                   LINE_8);
        }
        frames.push_back(img);
    }

    int frame_num = 20;
    namedWindow("test_pic", WINDOW_AUTOSIZE);
    moveWindow("test_pic", 40, 100);
    std::cout << "calc finished, press any key to start " << std::endl;
    waitKey();
    for (const auto &frame:frames) {
        imshow("image", frame);
        auto key = waitKey(int(1000 / frame_num));
        if (key == ('q' & 0xFF))
            break;
    }
    destroyAllWindows();
    std::cout << "finished" << std::endl;
    return 0;
}
