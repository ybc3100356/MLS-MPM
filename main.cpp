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
const size_t numParticles = 512;
const size_t numGrid = 32;
const Real dx = 1.0 / numGrid;
const Real inv_dx = 1.0 / dx;
const Real dt = 2.0e-3;

template<typename T>
vec2<T> clamp_pos(vec2<T> &&pos) {
    return vec2<T>(max(min(0.9, pos[0]), 0.1), max(min(0.9, pos[1]), 0.1));
}

void step(vector<Particle> &particles, vector<vector<vec2<Real>>> &grid_v, vector<vector<Real>> &grid_m) {
    // p2g
    for (auto &particle : particles) {
        auto base = (particle.position * inv_dx - 0.5).cast<int>();
        auto fx = (particle.position * inv_dx) - base.cast<Real>();
        // quadratic B-spline weights
        vector<vec2<Real>> w = {0.5 * ((1.5 - fx) * (1.5 - fx)),
                                0.75 - ((fx - 1) * (fx - 1)),
                                0.5 * (fx - 0.5) * (fx - 0.5)};
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                auto offset = vec2<int>(i, j);
                auto weight = w[i][0] * w[j][1];
                auto index = base + offset;
                if (index[0] < numGrid && index[1] < numGrid) {
                    grid_v[index[0]][index[1]] += weight * particle.velocity;
                    grid_m[index[0]][index[1]] += weight;
                }
            }
        }
    }

    // normalization
    for (int i = 0; i < grid_m.size(); i++) {
        for (int j = 0; j < grid_m[0].size(); j++) {
            if (grid_m[i][j] > 0) {
                auto inv_m = 1.0 / grid_m[i][j];
                grid_v[i][j] = inv_m * grid_v[i][j];
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
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                auto offset = vec2<int>(i, j);
                auto weight = w[i][0] * w[j][1];
                auto index = base + offset;
                if (index[0] < numGrid && index[1] < numGrid) {
                    new_v += weight * grid_v[index[0]][index[1]];
                }
            }
        }

        particle.position = clamp_pos(particle.position + particle.velocity * dt);// boundary
        particle.velocity = new_v;
    }
}

int main() {
    vector<Particle> particles(numParticles);
    vector<vector<vec2<Real>>> grid_v(numGrid, vector<vec2<Real>>(numGrid, vec2<Real>(0, 0)));
    vector<vector<Real>> grid_m(numGrid, vector<Real>(numGrid));
    vector<Mat> frames;
    int numFrames = 1000;
    while (numFrames--) {
        for (int s = 0; s < 10; s++) {
            for (auto &m : grid_m)
                std::fill(m.begin(), m.end(), 0);
            for (auto &v : grid_v)
                std::fill(v.begin(), v.end(), vec2<Real>(0, 0));
            step(particles, grid_v, grid_m);
        }
        Mat img = Mat::zeros(width, width, CV_8UC3);
        for (auto particle : particles) {
            circle(img,
                   Point2d(particle.position[0] * width, particle.position[1] * width),
                   3,
                   Scalar(0, 0, 255),
                   FILLED,
                   LINE_8);
        }
        frames.push_back(img);
    }

    int frame_num = 60;
    namedWindow("test_pic", WINDOW_AUTOSIZE);
    moveWindow("test_pic", 40, 100);
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
