# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, nonecheck=False
from libc.stdlib cimport malloc, free
from libc.string cimport memset
import numpy as np
cimport numpy as cnp

cdef long long g_m_ops = 0
cdef long long g_as_ops = 0

cpdef void reset_ops():
    global g_m_ops, g_as_ops
    g_m_ops = 0
    g_as_ops = 0

cpdef tuple get_ops():
    return (g_m_ops, g_as_ops, g_m_ops + g_as_ops)

cdef inline int next_pow2(int n) noexcept nogil:
    cdef int p = 1
    while p < n:
        p <<= 1
    return p

cdef inline int next_pow3(int n) noexcept nogil:
    cdef int p = 1
    while p < n:
        p *= 3
    return p

cdef void _add_c(const double * A, int sa, const double * B, int sb, double * C, int sc, int r, int c) noexcept nogil:
    global g_as_ops
    cdef int i, j
    for i in range(r):
        for j in range(c):
            C[i * sc + j] = A[i * sa + j] + B[i * sb + j]
    g_as_ops += <long long> r * c

cdef void _sub_c(const double * A, int sa, const double * B, int sb, double * C, int sc, int r, int c) noexcept nogil:
    global g_as_ops
    cdef int i, j
    for i in range(r):
        for j in range(c):
            C[i * sc + j] = A[i * sa + j] - B[i * sb + j]
    g_as_ops += <long long> r * c

cdef void _binet_add_c(const double * A, int sa, const double * B, int sb, double * C, int sc, int n) noexcept nogil:
    global g_m_ops, g_as_ops
    cdef int k
    if n == 1:
        C[0] += A[0] * B[0]
        g_m_ops += 1
        g_as_ops += 1
        return
    if n == 2:
        C[0] += A[0] * B[0] + A[1] * B[sb]
        C[1] += A[0] * B[1] + A[1] * B[sb + 1]
        C[sc] += A[sa] * B[0] + A[sa + 1] * B[sb]
        C[sc + 1] += A[sa] * B[1] + A[sa + 1] * B[sb + 1]
        g_m_ops += 8
        g_as_ops += 8
        return
    k = n >> 1
    _binet_add_c(A, sa, B, sb, C, sc, k)
    _binet_add_c(A + k, sa, B + k * sb, sb, C, sc, k)
    _binet_add_c(A, sa, B + k, sb, C + k, sc, k)
    _binet_add_c(A + k, sa, B + k * sb + k, sb, C + k, sc, k)
    _binet_add_c(A + k * sa, sa, B, sb, C + k * sc, sc, k)
    _binet_add_c(A + k * sa + k, sa, B + k * sb, sb, C + k * sc, sc, k)
    _binet_add_c(A + k * sa, sa, B + k, sb, C + k * sc + k, sc, k)
    _binet_add_c(A + k * sa + k, sa, B + k * sb + k, sb, C + k * sc + k, sc, k)

cdef void _strassen_c(const double * A, int sa, const double * B, int sb, double * C, int sc, int n) noexcept nogil:
    global g_m_ops, g_as_ops
    cdef double a11, a12, a21, a22, b11, b12, b21, b22
    cdef double m1, m2, m3, m4, m5, m6, m7
    if n == 1:
        C[0] = A[0] * B[0]
        g_m_ops += 1
        return
    if n == 2:
        a11 = A[0];
        a12 = A[1]
        a21 = A[sa];
        a22 = A[sa + 1]
        b11 = B[0];
        b12 = B[1]
        b21 = B[sb];
        b22 = B[sb + 1]
        m1 = (a11 + a22) * (b11 + b22)
        m2 = (a21 + a22) * b11
        m3 = a11 * (b12 - b22)
        m4 = a22 * (b21 - b11)
        m5 = (a11 + a12) * b22
        m6 = (a21 - a11) * (b11 + b12)
        m7 = (a12 - a22) * (b21 + b22)
        C[0] = m1 + m4 - m5 + m7
        C[1] = m3 + m5
        C[sc] = m2 + m4
        C[sc + 1] = m1 - m2 + m3 + m6
        g_m_ops += 7
        g_as_ops += 18
        return

    cdef int k = n >> 1
    cdef int k2 = k * k
    cdef double * work = <double *> malloc(9 * k2 * sizeof(double))
    cdef double * T1 = work
    cdef double * T2 = work + k2
    cdef double * M1 = work + 2 * k2
    cdef double * M2 = work + 3 * k2
    cdef double * M3 = work + 4 * k2
    cdef double * M4 = work + 5 * k2
    cdef double * M5 = work + 6 * k2
    cdef double * M6 = work + 7 * k2
    cdef double * M7 = work + 8 * k2

    cdef const double * A11 = A
    cdef const double * A12 = A + k
    cdef const double * A21 = A + k * sa
    cdef const double * A22 = A + k * sa + k

    cdef const double * B11 = B
    cdef const double * B12 = B + k
    cdef const double * B21 = B + k * sb
    cdef const double * B22 = B + k * sb + k

    _add_c(A11, sa, A22, sa, T1, k, k, k)
    _add_c(B11, sb, B22, sb, T2, k, k, k)
    _strassen_c(T1, k, T2, k, M1, k, k)

    _add_c(A21, sa, A22, sa, T1, k, k, k)
    _strassen_c(T1, k, B11, sb, M2, k, k)

    _sub_c(B12, sb, B22, sb, T2, k, k, k)
    _strassen_c(A11, sa, T2, k, M3, k, k)

    _sub_c(B21, sb, B11, sb, T2, k, k, k)
    _strassen_c(A22, sa, T2, k, M4, k, k)

    _add_c(A11, sa, A12, sa, T1, k, k, k)
    _strassen_c(T1, k, B22, sb, M5, k, k)

    _sub_c(A21, sa, A11, sa, T1, k, k, k)
    _add_c(B11, sb, B12, sb, T2, k, k, k)
    _strassen_c(T1, k, T2, k, M6, k, k)

    _sub_c(A12, sa, A22, sa, T1, k, k, k)
    _add_c(B21, sb, B22, sb, T2, k, k, k)
    _strassen_c(T1, k, T2, k, M7, k, k)

    cdef int i, j
    for i in range(k):
        for j in range(k):
            C[i * sc + j] = M1[i * k + j] + M4[i * k + j] - M5[i * k + j] + M7[i * k + j]
            C[i * sc + (j + k)] = M3[i * k + j] + M5[i * k + j]
            C[(i + k) * sc + j] = M2[i * k + j] + M4[i * k + j]
            C[(i + k) * sc + (j + k)] = M1[i * k + j] - M2[i * k + j] + M3[i * k + j] + M6[i * k + j]
    g_as_ops += 8 * <long long> k2

    free(work)

cdef void _laderman_c(const double * X, int sx, const double * Y, int sy, double * C, int sc, int n) noexcept nogil:
    global g_m_ops, g_as_ops
    cdef double a11, a12, a13, a21, a22, a23, a31, a32, a33
    cdef double b11, b12, b13, b21, b22, b23, b31, b32, b33
    cdef double m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12
    cdef double m13, m14, m15, m16, m17, m18, m19, m20, m21, m22, m23
    if n == 1:
        C[0] = X[0] * Y[0]
        g_m_ops += 1
        return
    if n == 3:
        a11 = X[0];
        a12 = X[1];
        a13 = X[2]
        a21 = X[sx];
        a22 = X[sx + 1];
        a23 = X[sx + 2]
        a31 = X[2 * sx];
        a32 = X[2 * sx + 1];
        a33 = X[2 * sx + 2]

        b11 = Y[0];
        b12 = Y[1];
        b13 = Y[2]
        b21 = Y[sy];
        b22 = Y[sy + 1];
        b23 = Y[sy + 2]
        b31 = Y[2 * sy];
        b32 = Y[2 * sy + 1];
        b33 = Y[2 * sy + 2]

        m1 = (a11 + a12 + a13 - a21 - a22 - a32 - a33) * b22
        m2 = (a11 - a21) * (-b12 + b22)
        m3 = a22 * (-b11 + b12 + b21 - b22 - b23 - b31 + b33)
        m4 = (-a11 + a21 + a22) * (b11 - b12 + b22)
        m5 = (a21 + a22) * (-b11 + b12)
        m6 = a11 * b11
        m7 = (-a11 + a31 + a32) * (b11 - b13 + b23)
        m8 = (-a11 + a31) * (b13 - b23)
        m9 = (a31 + a32) * (-b11 + b13)
        m10 = (a11 + a12 + a13 - a22 - a23 - a31 - a32) * b23
        m11 = a32 * (-b11 + b13 + b21 - b22 - b23 - b31 + b32)
        m12 = (-a13 + a32 + a33) * (b22 + b31 - b32)
        m13 = (a13 - a33) * (b22 - b32)
        m14 = a13 * b31
        m15 = (a32 + a33) * (-b31 + b32)
        m16 = (-a13 + a22 + a23) * (b23 + b31 - b33)
        m17 = (a13 - a23) * (b23 - b33)
        m18 = (a22 + a23) * (-b31 + b33)
        m19 = a12 * b21
        m20 = a23 * b32
        m21 = a21 * b13
        m22 = a31 * b12
        m23 = a33 * b33

        C[0] = m6 + m14 + m19
        C[1] = m1 + m4 + m5 + m6 + m12 + m14 + m15
        C[2] = m6 + m7 + m9 + m10 + m14 + m16 + m18
        C[sc] = m2 + m3 + m4 + m6 + m14 + m16 + m17
        C[sc + 1] = m2 + m4 + m5 + m6 + m20
        C[sc + 2] = m14 + m16 + m17 + m18 + m21
        C[2 * sc] = m6 + m7 + m8 + m11 + m12 + m13 + m14
        C[2 * sc + 1] = m12 + m13 + m14 + m15 + m22
        C[2 * sc + 2] = m6 + m7 + m8 + m9 + m23

        g_m_ops += 23
        g_as_ops += 98
        return

    cdef int k = n // 3
    cdef int k2 = k * k
    cdef double * work = <double *> malloc(25 * k2 * sizeof(double))
    cdef double * T1 = work
    cdef double * T2 = work + k2
    cdef double * M = work + 2 * k2

    cdef const double * A11 = X
    cdef const double * A12 = X + k
    cdef const double * A13 = X + 2 * k
    cdef const double * A21 = X + k * sx
    cdef const double * A22 = X + k * sx + k
    cdef const double * A23 = X + k * sx + 2 * k
    cdef const double * A31 = X + 2 * k * sx
    cdef const double * A32 = X + 2 * k * sx + k
    cdef const double * A33 = X + 2 * k * sx + 2 * k

    cdef const double * B11 = Y
    cdef const double * B12 = Y + k
    cdef const double * B13 = Y + 2 * k
    cdef const double * B21 = Y + k * sy
    cdef const double * B22 = Y + k * sy + k
    cdef const double * B23 = Y + k * sy + 2 * k
    cdef const double * B31 = Y + 2 * k * sy
    cdef const double * B32 = Y + 2 * k * sy + k
    cdef const double * B33 = Y + 2 * k * sy + 2 * k

    cdef int i, j, idx
    for i in range(k):
        for j in range(k):
            T1[i * k + j] = A11[i * sx + j] + A12[i * sx + j] + A13[i * sx + j] - A21[i * sx + j] - A22[i * sx + j] - \
                            A32[i * sx + j] - A33[i * sx + j]
    g_as_ops += 6 * <long long> k2
    _laderman_c(T1, k, B22, sy, M + 0 * k2, k, k)

    _sub_c(A11, sx, A21, sx, T1, k, k, k)
    _sub_c(B22, sy, B12, sy, T2, k, k, k)
    _laderman_c(T1, k, T2, k, M + 1 * k2, k, k)

    for i in range(k):
        for j in range(k):
            T2[i * k + j] = B12[i * sy + j] - B11[i * sy + j] + B21[i * sy + j] - B22[i * sy + j] - B23[i * sy + j] - \
                            B31[i * sy + j] + B33[i * sy + j]
    g_as_ops += 6 * <long long> k2
    _laderman_c(A22, sx, T2, k, M + 2 * k2, k, k)

    for i in range(k):
        for j in range(k):
            T1[i * k + j] = A21[i * sx + j] - A11[i * sx + j] + A22[i * sx + j]
            T2[i * k + j] = B11[i * sy + j] - B12[i * sy + j] + B22[i * sy + j]
    g_as_ops += 4 * <long long> k2
    _laderman_c(T1, k, T2, k, M + 3 * k2, k, k)

    _add_c(A21, sx, A22, sx, T1, k, k, k)
    _sub_c(B12, sy, B11, sy, T2, k, k, k)
    _laderman_c(T1, k, T2, k, M + 4 * k2, k, k)

    _laderman_c(A11, sx, B11, sy, M + 5 * k2, k, k)

    for i in range(k):
        for j in range(k):
            T1[i * k + j] = A31[i * sx + j] - A11[i * sx + j] + A32[i * sx + j]
            T2[i * k + j] = B11[i * sy + j] - B13[i * sy + j] + B23[i * sy + j]
    g_as_ops += 4 * <long long> k2
    _laderman_c(T1, k, T2, k, M + 6 * k2, k, k)

    _sub_c(A31, sx, A11, sx, T1, k, k, k)
    _sub_c(B13, sy, B23, sy, T2, k, k, k)
    _laderman_c(T1, k, T2, k, M + 7 * k2, k, k)

    _add_c(A31, sx, A32, sx, T1, k, k, k)
    _sub_c(B13, sy, B11, sy, T2, k, k, k)
    _laderman_c(T1, k, T2, k, M + 8 * k2, k, k)

    for i in range(k):
        for j in range(k):
            T1[i * k + j] = A11[i * sx + j] + A12[i * sx + j] + A13[i * sx + j] - A22[i * sx + j] - A23[i * sx + j] - \
                            A31[i * sx + j] - A32[i * sx + j]
    g_as_ops += 6 * <long long> k2
    _laderman_c(T1, k, B23, sy, M + 9 * k2, k, k)

    for i in range(k):
        for j in range(k):
            T2[i * k + j] = B13[i * sy + j] - B11[i * sy + j] + B21[i * sy + j] - B22[i * sy + j] - B23[i * sy + j] - \
                            B31[i * sy + j] + B32[i * sy + j]
    g_as_ops += 6 * <long long> k2
    _laderman_c(A32, sx, T2, k, M + 10 * k2, k, k)

    for i in range(k):
        for j in range(k):
            T1[i * k + j] = A32[i * sx + j] - A13[i * sx + j] + A33[i * sx + j]
            T2[i * k + j] = B22[i * sy + j] + B31[i * sy + j] - B32[i * sy + j]
    g_as_ops += 4 * <long long> k2
    _laderman_c(T1, k, T2, k, M + 11 * k2, k, k)

    _sub_c(A13, sx, A33, sx, T1, k, k, k)
    _sub_c(B22, sy, B32, sy, T2, k, k, k)
    _laderman_c(T1, k, T2, k, M + 12 * k2, k, k)

    _laderman_c(A13, sx, B31, sy, M + 13 * k2, k, k)

    _add_c(A32, sx, A33, sx, T1, k, k, k)
    _sub_c(B32, sy, B31, sy, T2, k, k, k)
    _laderman_c(T1, k, T2, k, M + 14 * k2, k, k)

    for i in range(k):
        for j in range(k):
            T1[i * k + j] = A22[i * sx + j] - A13[i * sx + j] + A23[i * sx + j]
            T2[i * k + j] = B23[i * sy + j] + B31[i * sy + j] - B33[i * sy + j]
    g_as_ops += 4 * <long long> k2
    _laderman_c(T1, k, T2, k, M + 15 * k2, k, k)

    _sub_c(A13, sx, A23, sx, T1, k, k, k)
    _sub_c(B23, sy, B33, sy, T2, k, k, k)
    _laderman_c(T1, k, T2, k, M + 16 * k2, k, k)

    _add_c(A22, sx, A23, sx, T1, k, k, k)
    _sub_c(B33, sy, B31, sy, T2, k, k, k)
    _laderman_c(T1, k, T2, k, M + 17 * k2, k, k)

    _laderman_c(A12, sx, B21, sy, M + 18 * k2, k, k)
    _laderman_c(A23, sx, B32, sy, M + 19 * k2, k, k)
    _laderman_c(A21, sx, B13, sy, M + 20 * k2, k, k)
    _laderman_c(A31, sx, B12, sy, M + 21 * k2, k, k)
    _laderman_c(A33, sx, B33, sy, M + 22 * k2, k, k)

    for i in range(k):
        for j in range(k):
            idx = i * k + j
            C[i * sc + j] = M[5 * k2 + idx] + M[13 * k2 + idx] + M[18 * k2 + idx]
            C[i * sc + (j + k)] = M[0 * k2 + idx] + M[3 * k2 + idx] + M[4 * k2 + idx] + M[5 * k2 + idx] + M[
                11 * k2 + idx] + M[13 * k2 + idx] + M[14 * k2 + idx]
            C[i * sc + (j + 2 * k)] = M[5 * k2 + idx] + M[6 * k2 + idx] + M[8 * k2 + idx] + M[9 * k2 + idx] + M[
                13 * k2 + idx] + M[15 * k2 + idx] + M[17 * k2 + idx]
            C[(i + k) * sc + j] = M[1 * k2 + idx] + M[2 * k2 + idx] + M[3 * k2 + idx] + M[5 * k2 + idx] + M[
                13 * k2 + idx] + M[15 * k2 + idx] + M[16 * k2 + idx]
            C[(i + k) * sc + (j + k)] = M[1 * k2 + idx] + M[3 * k2 + idx] + M[4 * k2 + idx] + M[5 * k2 + idx] + M[
                19 * k2 + idx]
            C[(i + k) * sc + (j + 2 * k)] = M[13 * k2 + idx] + M[15 * k2 + idx] + M[16 * k2 + idx] + M[17 * k2 + idx] + \
                                            M[20 * k2 + idx]
            C[(i + 2 * k) * sc + j] = M[5 * k2 + idx] + M[6 * k2 + idx] + M[7 * k2 + idx] + M[10 * k2 + idx] + M[
                11 * k2 + idx] + M[12 * k2 + idx] + M[13 * k2 + idx]
            C[(i + 2 * k) * sc + (j + k)] = M[11 * k2 + idx] + M[12 * k2 + idx] + M[13 * k2 + idx] + M[14 * k2 + idx] + \
                                            M[21 * k2 + idx]
            C[(i + 2 * k) * sc + (j + 2 * k)] = M[5 * k2 + idx] + M[6 * k2 + idx] + M[7 * k2 + idx] + M[8 * k2 + idx] + \
                                                M[22 * k2 + idx]
    g_as_ops += 42 * <long long> k2

    free(work)

cdef void _mul_c(const double * A, int sa, int rA, int cA,
                 const double * B, int sb, int rB, int cB,
                 double * C, int sc, int method) noexcept nogil:
    global g_m_ops, g_as_ops
    cdef int max_dim = rA
    if cA > max_dim: max_dim = cA
    if rB > max_dim: max_dim = rB
    if cB > max_dim: max_dim = cB

    cdef int n, i, j
    if method == 0 or method == 1:
        n = next_pow2(max_dim)
    else:
        n = next_pow3(max_dim)

    if n == 1:
        C[0] = A[0] * B[0]
        g_m_ops += 1
        return

    cdef int n2 = n * n
    cdef double * pad = <double *> malloc(3 * n2 * sizeof(double))
    memset(pad, 0, 3 * n2 * sizeof(double))
    cdef double * A_pad = pad
    cdef double * B_pad = pad + n2
    cdef double * C_pad = pad + 2 * n2

    for i in range(rA):
        for j in range(cA):
            A_pad[i * n + j] = A[i * sa + j]
    for i in range(rB):
        for j in range(cB):
            B_pad[i * n + j] = B[i * sb + j]

    if method == 0:
        _binet_add_c(A_pad, n, B_pad, n, C_pad, n, n)
        g_as_ops -= <long long> n2
    elif method == 1:
        _strassen_c(A_pad, n, B_pad, n, C_pad, n, n)
    else:
        _laderman_c(A_pad, n, B_pad, n, C_pad, n, n)

    for i in range(rA):
        for j in range(cB):
            C[i * sc + j] = C_pad[i * n + j]

    free(pad)

cdef void _inv_c(const double * A, int sa, double * B, int sb, int n, int method) noexcept nogil:
    global g_m_ops, g_as_ops
    if n == 1:
        B[0] = 1.0 / A[0]
        g_m_ops += 1
        return

    cdef int k = n >> 1
    cdef int m = n - k
    cdef int i, j

    cdef int sz = 2 * k * k + 4 * k * m + 3 * m * m
    cdef double * work = <double *> malloc(sz * sizeof(double))
    cdef double * A11_inv = work
    cdef double * U = A11_inv + k * k
    cdef double * V = U + k * m
    cdef double * T_mm = V + m * k
    cdef double * S22 = T_mm + m * m
    cdef double * S22_inv = S22 + m * m
    cdef double * US_inv = S22_inv + m * m
    cdef double * T_kk = US_inv + k * m
    cdef double * T_mk = T_kk + k * k

    cdef const double * A11 = A
    cdef const double * A12 = A + k
    cdef const double * A21 = A + k * sa
    cdef const double * A22 = A + k * sa + k

    _inv_c(A11, sa, A11_inv, k, k, method)

    _mul_c(A11_inv, k, k, k, A12, sa, k, m, U, m, method)
    _mul_c(A21, sa, m, k, A11_inv, k, k, k, V, k, method)

    _mul_c(A21, sa, m, k, U, m, k, m, T_mm, m, method)
    _sub_c(A22, sa, T_mm, m, S22, m, m, m)

    _inv_c(S22, m, S22_inv, m, m, method)

    _mul_c(U, m, k, m, S22_inv, m, m, m, US_inv, m, method)
    for i in range(k):
        for j in range(m):
            B[i * sb + (j + k)] = -US_inv[i * m + j]
    g_as_ops += <long long> k * m

    _mul_c(US_inv, m, k, m, V, k, m, k, T_kk, k, method)
    for i in range(k):
        for j in range(k):
            B[i * sb + j] = A11_inv[i * k + j] + T_kk[i * k + j]
    g_as_ops += <long long> k * k

    _mul_c(S22_inv, m, m, m, V, k, m, k, T_mk, k, method)
    for i in range(m):
        for j in range(k):
            B[(i + k) * sb + j] = -T_mk[i * k + j]
    g_as_ops += <long long> m * k

    for i in range(m):
        for j in range(m):
            B[(i + k) * sb + (j + k)] = S22_inv[i * m + j]

    free(work)

cdef void _lu_c(const double * A, int sa, double * L, int sl, double * U, int su, int n, int method) noexcept nogil:
    global g_m_ops, g_as_ops
    cdef int i, j
    for i in range(n):
        for j in range(n):
            L[i * sl + j] = 0.0
            U[i * su + j] = 0.0
    if n == 1:
        L[0] = 1.0
        U[0] = A[0]
        return

    cdef int k = n >> 1
    cdef int m = n - k

    cdef int sz = 2 * k * k + 2 * m * m
    cdef double * work = <double *> malloc(sz * sizeof(double))
    cdef double * U11_inv = work
    cdef double * L11_inv = U11_inv + k * k
    cdef double * T_mm = L11_inv + k * k
    cdef double * S = T_mm + m * m

    cdef const double * A11 = A
    cdef const double * A12 = A + k
    cdef const double * A21 = A + k * sa
    cdef const double * A22 = A + k * sa + k

    cdef double * L11 = L
    cdef double * L21 = L + k * sl
    cdef double * LS = L + k * sl + k

    cdef double * U11 = U
    cdef double * U12 = U + k
    cdef double * US = U + k * su + k

    _lu_c(A11, sa, L11, sl, U11, su, k, method)
    _inv_c(U11, su, U11_inv, k, k, method)
    _mul_c(A21, sa, m, k, U11_inv, k, k, k, L21, sl, method)
    _inv_c(L11, sl, L11_inv, k, k, method)
    _mul_c(L11_inv, k, k, k, A12, sa, k, m, U12, su, method)
    _mul_c(L21, sl, m, k, U12, su, k, m, T_mm, m, method)
    _sub_c(A22, sa, T_mm, m, S, m, m, m)
    _lu_c(S, m, LS, sl, US, su, m, method)

    free(work)

cpdef tuple binet_mul(cnp.ndarray[double, ndim=2] A, cnp.ndarray[double, ndim=2] B):
    reset_ops()
    cdef int rA = A.shape[0], cA = A.shape[1]
    cdef int rB = B.shape[0], cB = B.shape[1]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] B_c = np.ascontiguousarray(B, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] C = np.empty((rA, cB), dtype=np.float64)
    _mul_c(&A_c[0, 0], cA, rA, cA, &B_c[0, 0], cB, rB, cB, &C[0, 0], cB, 0)
    return C, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple strassen_mul(cnp.ndarray[double, ndim=2] A, cnp.ndarray[double, ndim=2] B):
    reset_ops()
    cdef int rA = A.shape[0], cA = A.shape[1]
    cdef int rB = B.shape[0], cB = B.shape[1]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] B_c = np.ascontiguousarray(B, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] C = np.empty((rA, cB), dtype=np.float64)
    _mul_c(&A_c[0, 0], cA, rA, cA, &B_c[0, 0], cB, rB, cB, &C[0, 0], cB, 1)
    return C, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple laderman_mul(cnp.ndarray[double, ndim=2] A, cnp.ndarray[double, ndim=2] B):
    reset_ops()
    cdef int rA = A.shape[0], cA = A.shape[1]
    cdef int rB = B.shape[0], cB = B.shape[1]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] B_c = np.ascontiguousarray(B, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] C = np.empty((rA, cB), dtype=np.float64)
    _mul_c(&A_c[0, 0], cA, rA, cA, &B_c[0, 0], cB, rB, cB, &C[0, 0], cB, 2)
    return C, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple inverse_binet(cnp.ndarray[double, ndim=2] A):
    reset_ops()
    cdef int n = A.shape[0]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] B = np.empty((n, n), dtype=np.float64)
    _inv_c(&A_c[0, 0], n, &B[0, 0], n, n, 0)
    return B, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple inverse_strassen(cnp.ndarray[double, ndim=2] A):
    reset_ops()
    cdef int n = A.shape[0]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] B = np.empty((n, n), dtype=np.float64)
    _inv_c(&A_c[0, 0], n, &B[0, 0], n, n, 1)
    return B, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple inverse_laderman(cnp.ndarray[double, ndim=2] A):
    reset_ops()
    cdef int n = A.shape[0]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] B = np.empty((n, n), dtype=np.float64)
    _inv_c(&A_c[0, 0], n, &B[0, 0], n, n, 2)
    return B, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple lu_binet(cnp.ndarray[double, ndim=2] A):
    reset_ops()
    cdef int n = A.shape[0]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] L = np.empty((n, n), dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] U = np.empty((n, n), dtype=np.float64)
    _lu_c(&A_c[0, 0], n, &L[0, 0], n, &U[0, 0], n, n, 0)
    return L, U, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple lu_strassen(cnp.ndarray[double, ndim=2] A):
    reset_ops()
    cdef int n = A.shape[0]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] L = np.empty((n, n), dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] U = np.empty((n, n), dtype=np.float64)
    _lu_c(&A_c[0, 0], n, &L[0, 0], n, &U[0, 0], n, n, 1)
    return L, U, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple lu_laderman(cnp.ndarray[double, ndim=2] A):
    reset_ops()
    cdef int n = A.shape[0]
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] L = np.empty((n, n), dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] U = np.empty((n, n), dtype=np.float64)
    _lu_c(&A_c[0, 0], n, &L[0, 0], n, &U[0, 0], n, n, 2)
    return L, U, g_m_ops, g_as_ops, g_m_ops + g_as_ops

cpdef tuple det_lu(cnp.ndarray[double, ndim=2] A, int method=0):
    global g_m_ops, g_as_ops
    reset_ops()
    cdef int n = A.shape[0], i
    cdef cnp.ndarray[double, ndim=2] A_c = np.ascontiguousarray(A, dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] L = np.empty((n, n), dtype=np.float64)
    cdef cnp.ndarray[double, ndim=2] U = np.empty((n, n), dtype=np.float64)
    _lu_c(&A_c[0, 0], n, &L[0, 0], n, &U[0, 0], n, n, method)
    cdef double d = 1.0
    for i in range(n):
        d *= U[i, i]
    if n > 1:
        g_m_ops += n - 1
    return d, g_m_ops, g_as_ops, g_m_ops + g_as_ops
