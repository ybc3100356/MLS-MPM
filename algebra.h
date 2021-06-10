//
// Created by ybc on 2021/6/9.
//

#ifndef MLSMPM_ALGEBRA_H
#define MLSMPM_ALGEBRA_H

template<typename T>
class mat2;

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
    friend T dot(const vec2 &v1, const vec2 &v2) {
        return v1.element[0] * v2.element[0] + v1.element[1] * v2.element[1];
    }

    friend vec2 operator*(const vec2 &v, T scalar) {
        return vec2(v.element[0] * scalar, v.element[1] * scalar);
    }

    friend vec2 operator*(T scalar, const vec2 &v) {
        return vec2(v.element[0] * scalar, v.element[1] * scalar);
    }

    friend vec2 operator*(const vec2 &v1, const vec2 &v2) {
        return vec2(v1.element[0] * v2.element[0], v1.element[1] * v2.element[1]);
    }

    mat2<T> outer_product(vec2<double> &vec2) {
        return mat2<T>(element[0] * vec2[0], element[0] * vec2[1], element[1] * vec2[0], element[1] * vec2[1]);
    }

    explicit vec2(T x = 0, T y = 0) : element{x, y} {}

    template<typename U>
    vec2<U> cast() {
        return vec2<U>((U) this->element[0], (U) this->element[1]);
    }

    T &operator[](size_t i) { return element[i]; }

    vec2(const vec2 &v) : element{v.element[0],
                                  v.element[1]} {}

};

template<typename T>
class mat2 {
    vec2<T> element[2]; // row vector
public:
    vec2<T> &operator[](size_t i) { return element[i]; }

    // plus
    friend mat2 operator+(const mat2 &m1, const mat2 &m2) {
        return mat2(m1.element[0] + m2.element[0], m1.element[1] + m2.element[1]);
    }

    friend mat2 operator+(const mat2 &m, T scalar) {
        return mat2(m.element[0] + scalar, m.element[1] + scalar);
    }

    friend mat2 operator+(T scalar, const mat2 &m) {
        return mat2(m.element[0] + scalar, m.element[1] + scalar);
    }

    mat2 &operator+=(const mat2 &rhs) {
        element[0] += rhs.element[0];
        element[1] += rhs.element[1];
        return *this;
    }

    // minus
    friend mat2 operator-(const mat2 &m1, const mat2 &m2) {
        return mat2(m1.element[0] - m2.element[0], m1.element[1] - m2.element[1]);
    }

    friend mat2 operator-(const mat2 &m, T scalar) {
        return mat2(m.element[0] - scalar, m.element[1] - scalar);
    }

    friend mat2 operator-(T scalar, const mat2 &m) {
        return mat2(scalar - m.element[0], scalar - m.element[1]);
    }

    // mul

    vec2<T> mul(vec2<T> &v) {
        return vec2<T>(element[0][0] * v[0] + element[0][1] * v[1], element[1][0] * v[0] + element[1][1] * v[1]);
    }

    friend mat2 operator*(const mat2 &m, T scalar) {
        return mat2(m.element[0] * scalar, m.element[1] * scalar);
    }

    friend mat2 operator*(T scalar, const mat2 &m) {
        return mat2(m.element[0] * scalar, m.element[1] * scalar);
    }

    friend mat2 operator*(const mat2 &m1, const mat2 &m2) {
        return mat2(m1.element[0] * m2.element[0], m1.element[1] * m2.element[1]);
    }

    explicit mat2(vec2<T> row1, vec2<T> row2) : element{row1,
                                                        row2} {}

    explicit mat2(T e11 = 0, T e12 = 0, T e21 = 0, T e22 = 0) : element{vec2<T>{e11, e12},
                                                                        vec2<T>{e21, e22}} {}

    T trace() { return element[0][0] + element[1][1]; }
};


#endif //MLSMPM_ALGEBRA_H
