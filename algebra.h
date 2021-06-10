//
// Created by ybc on 2021/6/9.
//

#ifndef MLSMPM_ALGEBRA_H
#define MLSMPM_ALGEBRA_H

template<typename T>
class vec2 {
    T element[2];
public:
    // plus
    friend vec2 operator+(const vec2 &v1, const vec2 &v2) {
        return vec2(v1.element[0] + v2.element[0], v1.element[1] + v2.element[1]);
    }

    friend vec2 operator+(const vec2 &v, T scalar) {
        return vec2(v.element[0] + scalar, v.element[1] + scalar);
    }

    friend vec2 operator+(T scalar, const vec2 &v) {
        return vec2(v.element[0] + scalar, v.element[1] + scalar);
    }

    vec2 &operator+=(const vec2 &rhs) {
        element[0] += rhs.element[0];
        element[1] += rhs.element[1];
        return *this;
    }

    // minus
    friend vec2 operator-(const vec2 &v1, const vec2 &v2) {
        return vec2(v1.element[0] - v2.element[0], v1.element[1] - v2.element[1]);
    }

    friend vec2 operator-(const vec2 &v, T scalar) {
        return vec2(v.element[0] - scalar, v.element[1] - scalar);
    }

    friend vec2 operator-(T scalar, const vec2 &v) {
        return vec2(scalar - v.element[0], scalar - v.element[1]);
    }

    // mul
    friend vec2 operator*(const vec2 &v, T scalar) {
        return vec2(v.element[0] * scalar, v.element[1] * scalar);
    }

    friend vec2 operator*(T scalar, const vec2 &v) {
        return vec2(v.element[0] * scalar, v.element[1] * scalar);
    }

    friend vec2 operator*(const vec2 &v1, const vec2 &v2) {
        return vec2(v1.element[0] * v2.element[0], v1.element[1] * v2.element[1]);
    }

    explicit vec2(T x = 0, T y = 0) : element{x, y} {}

    template<typename U>
    vec2<U> cast() {
        return vec2<U>((U) this->element[0], (U) this->element[1]);
    }

    T &operator[](size_t i) { return element[i]; }

    vec2(const vec2 &v) : element{v.element[0], v.element[1]} {}

};


#endif //MLSMPM_ALGEBRA_H
