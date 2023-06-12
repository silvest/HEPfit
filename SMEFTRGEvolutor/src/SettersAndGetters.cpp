inline void RGESolver::WC1_set(double * c, int i, int j, double val) {
    c[3 * i + j] = val;
}

inline double RGESolver::WC1(double * c, int i, int j) {
    return (c[3 * i + j]);
}

inline void RGESolver::Yukawa_set(gslpp::matrix<double> *y,
        int i, int j, double val) {
    y->assign(i, j, val);
}

inline double RGESolver::Yukawa(gslpp::matrix<double> *y,
            int i, int j) {
        return y->operator()(i,j);
}





inline double RGESolver::WC2R(double * c, int i, int j) {
    return (j >= i ? c[3 * i + j] : c[3 * j + i]);
}

inline void RGESolver::WC2R_set(double * c, int i, int j, double val) {
    if (j >= i)
        c[3 * i + j] = val;
    else
        c[3 * j + i] = val;
}

inline double RGESolver::WC2I(double * c, int i, int j) {
    return (j > i ? c[3 * i + j] : (j == i ? 0. : -c[3 * j + i]));
}

inline void RGESolver::WC2I_set(double * c, int i, int j, double val) {
    if (j > i)
        c[3 * i + j] = val;
    else if (j < i)
        c[3 * j + i] = - val;
}

inline double RGESolver::WC3(double * c, int i, int j) {
    return (j >= i ? c[3 * i + j] : c[3 * j + i]);
}

inline void RGESolver::WC3_set(double * c, int i, int j, double val) {
    if (j >= i)
        c[3 * i + j] = val;
    else
        c[3 * j + i] = val;
}

inline void RGESolver::WC5_set(double * c, int i, int j, int k, int l, double val) {
    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
}

inline double RGESolver::WC5(double * c, int i, int j, int k, int l) {
    return c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l];
}

inline double RGESolver::WC6R(double * c, int i, int j, int k, int l) {
    return (i < j ? (k < l ? (i <= k ? (j <= l ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : (k == l ? (j <= k ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (i < k ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]))
            : (k == j ? (i == l ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (i > l ? c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]))
            : (k > j ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]))))
            : (i == j ? (k == l ? (i <= k ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : (k > l ? (l >= i ? c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]
            : c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i])
            : (k >= j ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])))
            : (k < l ? (i == l ? (k < j ? c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]
            : (k == j ? c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]
            : c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]))
            : (i > l ? c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]
            : c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]))
            : (k == l ? (j >= k ? c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k])
            : (l >= j ? (i > k ? c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k])
            : c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i])))));
}

inline void RGESolver::WC6R_set(double * c, int i, int j, int k, int l, double
        val) {
    if (i < j) {
        if (k < l) {
            if (i <= k) {
                if (j <= l) {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                } else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
            } else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
        } else if (k == l) {
            if (j <= k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else if (i < k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            }
            else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
        } else { //k>l
            if (k == j) {
                if (i == l) {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                } else if (i > l) c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
                else {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                }
            } else if (k > j) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
        }
    } else if (i == j) {
        if (k == l) {
            if (i <= k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
        } else if (k > l) {
            if (l >= i) c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
            else c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
        } else { //k<l
            if (k >= j) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
        }
    } else { //i>j
        if (k < l) {
            if (i == l) {
                if (k < j) c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                else if (k == j) c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
            } else if (i > l) c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
            else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
        } else if (k == l)
            if (j >= k) c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
            else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
        else
            if (l >= j)
            if (i > k) c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
            else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
        else c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
    }
}

inline double RGESolver::WC6I(double * c, int i, int j, int k, int l) {
    return (i < j ? (k < l ? (i <= k ? (j <= l ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : (k == l ? (j <= k ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (i < k ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]))
            : (k == j ? (i == l ? 0.
            : (i > l ? -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]))
            : (k > j ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]))))
            : (i == j ? (k == l ? 0.
            : (k > l ? (l >= i ? -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]
            : -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i])
            : (k >= j ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])))
            : (k < l ? (i == l ? (k < j ? c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]
            : (k == j ? 0.
            : -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]))
            : (i > l ? c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]
            : -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]))
            : (k == l ? (j >= k ? -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k])
            : (l >= j ? (i > k ? -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k])
            : -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i])))));
}

inline void RGESolver::WC6I_set(double * c, int i, int j, int k, int l,
        double val) {
    if (i < j) {
        if (k < l) {
            if (i <= k)
                if (j <= l) {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                } else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val; //j>l
            else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val; //i>k
        } else if (k == l) {
            if (j <= k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else if (i < k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val; //j>k
        } else { //k>l
            if (k == j) {
                if (i > l)
                    c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
                else if (i < l) {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                }
            } else if (k > j) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val; //k<j
        }
    } else if (i == j) {
        if (k > l)
            if (l >= i)
                c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
            else c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val; //l<i
        else if (k < l) {
            if (k >= j) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val; //k<i
        }
    } else { //i>j
        if (k < l) {
            if (i == l) {
                if (k < j) c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                else if (k > j) c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
            } else if (i > l) c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
            else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
        } else if (k == l)
            if (j >= k) c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
            else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
        else
            if (l >= j)
            if (i > k) c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
            else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
        else c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
    }
}

inline double RGESolver::WC7R(double * c, int i, int j, int k, int l) {
    return (i < j ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (i == j ? (k == l ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (k > l ? c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]
            : c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]))
            : c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]));
}

inline void RGESolver::WC7R_set(double * c, int i, int j, int k, int l, double
        val) {
    if (i < j) c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
    else if (i == j)
        if (k == l) c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
        else if (k > l) c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
        else c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
    else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
}

inline double RGESolver::WC7I(double * c, int i, int j, int k, int l) {
    return (i < j ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (i == j ? (k == l ? 0.
            : (k > l ? -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]
            : c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]))
            : -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]));
}

inline void RGESolver::WC7I_set(double * c, int i, int j, int k, int l, double
        val) {
    if (i < j) c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
    else if (i == j)
        if (k == l) c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = 0.;
        else if (k > l) c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
        else c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
    else c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
}

inline double RGESolver::WC8R(double * c, int i, int j, int k, int l) {
    return (i < j ? (k < l ? (i <= k ? (j <= l ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : (k == l ? (j <= k ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (i < k ? c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]))
            : (k == j ? (i == l ? c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j]
            : (i > l ? c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k]
            : (j > l ? c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j]
            : c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l])))
            : (k > j ? (j > l ? c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k]
            : c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l])
            : (j > l ? (i > l ? c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k]
            : c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j])
            : c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i])))))
            : (i == j ? (k == l ? (i <= k ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : (k > l ? (l >= i ? c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]
            : (l < j ? (k < j ? c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k])
            : c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]))
            : (k >= j ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (k < i ? (l < i ? c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]
            : c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l])
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]))))
            : (k < l ? (i == l ? (k < j ? c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l]
            : (k == j ? (k < i ? c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : c[3 * 3 * 3 * j + 3 * 3 * k + 3 * i + l]))
            : (i > l ? (j < l ? c[3 * 3 * 3 * j + 3 * 3 * k + 3 * l + i]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : (k < i ? c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l]
            : c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k])))
            : (k == l ? (j >= k ? c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : (i > k ? c[3 * 3 * 3 * j + 3 * 3 * k + 3 * l + i]
            : c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]))
            : (l >= j ? (i > k ? c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k])
            : c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i])))));
}

inline void RGESolver::WC8R_set(double * c, int i, int j, int k, int l, double
        val) {
    if (i < j) {
        if (k < l) {
            if (i <= k) {
                if (j <= l) {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                } else {
                    c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                }
            } else {
                c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
            }
        } else if (k == l) {
            if (j <= k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else if (i < k) {
                c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << i+1 << ", " << l+1 << ", " << k+1 << ", " << j+1 << std::endl;
            } else {
                c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
            }
        } else { //k>l
            if (k == j) {
                if (i == l) {
                    c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << i+1 << ", " << l+1 << ", " << k+1 << ", " << j+1 << std::endl;
                } else if (i > l) {
                    c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << i+1 << ", " << j+1 << ", " << k+1 << std::endl;
                } else if (j > l) {
                    c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << i+1 << ", " << l+1 << ", " << k+1 << ", " << j+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                }
            } else if (k > j)
                if (j > l) {
                    c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << i+1 << ", " << j+1 << ", " << k+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                } else { //k<j
                if (j > l) {
                    if (i > l) {
                        c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << i+1 << ", " << j+1 << ", " << k+1 << std::endl;
                    } else {
                        c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << i+1 << ", " << l+1 << ", " << k+1 << ", " << j+1 << std::endl;
                    }
                } else {
                    c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
                }
            }
        }
    } else if (i == j) {
        if (k == l) {
            if (i <= k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else {
                c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
            }
        } else if (k > l) {
            if (l >= i) {
                c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << i+1 << ", " << l+1 << ", " << k+1 << std::endl;
            } else {
                if (l < j) {
                    if (k < j) {
                        c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
                    } else {
                        c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << i+1 << ", " << j+1 << ", " << k+1 << std::endl;
                    }
                } else {
                    c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
                }
            }
        } else { //k<l
            if (k >= j) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else {
                if (k < i) {
                    if (l < i) {
                        c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                    } else {
                        c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << j+1 << ", " << i+1 << ", " << l+1 << std::endl;
                    }
                } else {
                    c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                }
            }
        }
    } else { //i>j
        if (k < l) {
            if (i == l) {
                if (k < j) {
                    c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << j+1 << ", " << i+1 << ", " << l+1 << std::endl;
                } else if (k == j) {
                    if (k < i) {
                        c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << j+1 << ", " << i+1 << ", " << l+1 << std::endl;
                    } else {
                        c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                    }
                } else {
                    c[3 * 3 * 3 * j + 3 * 3 * k + 3 * i + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << k+1 << ", " << i+1 << ", " << l+1 << std::endl;
                }
            } else if (i > l) {
                if (j < l) {
                    c[3 * 3 * 3 * j + 3 * 3 * k + 3 * l + i] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << k+1 << ", " << l+1 << ", " << i+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                }
            } else {
                if (k < i) {
                    c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << j+1 << ", " << i+1 << ", " << l+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << i+1 << ", " << l+1 << ", " << k+1 << std::endl;
                }
            }
        } else if (k == l)
            if (j >= k) {
                c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
            } else {
                if (i > k) {
                    c[3 * 3 * 3 * j + 3 * 3 * k + 3 * l + i] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << k+1 << ", " << l+1 << ", " << i+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << i+1 << ", " << l+1 << ", " << k+1 << std::endl;
                }
            } else
            if (l >= j)
            if (i > k) {
                c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
            } else {
                c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << i+1 << ", " << l+1 << ", " << k+1 << std::endl;
            } else {
            c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = val;
            //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
        }
    }
}

inline double RGESolver::WC8I(double * c, int i, int j, int k, int l) {
    return (i < j ? (k < l ? (i <= k ? (j <= l ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : (k == l ? (j <= k ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (i < k ? c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]))
            : (k == j ? (i == l ? 0.
            : (i > l ? -c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k]
            : (j > l ? c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j]
            : c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l])))
            : (k > j ? (j > l ? -c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k]
            : c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l])
            : (j > l ? (i > l ? -c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k]
            : c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j])
            : -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i])))))
            : (i == j ? (k == l ? 0.
            : (k > l ? (l >= i ? -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]
            : (l < j ? (k < j ? -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : -c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k])
            : -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]))
            : (k >= j ? c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l]
            : (k < i ? (l < i ? c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]
            : c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l])
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]))))
            : (k < l ? (i == l ? (k < j ? c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] //if i>j
            //: (k < l ? (i == l ? (k < j ? c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j]
            : (k == j ? (k < i ? 0.
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : -c[3 * 3 * 3 * j + 3 * 3 * k + 3 * i + l]))
            : (i > l ? (j < l ? -c[3 * 3 * 3 * j + 3 * 3 * k + 3 * l + i]
            : c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j])
            : (k < i ? c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l]
            : -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k])))
            : (k == l ? (j >= k ? -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : (i > k ? -c[3 * 3 * 3 * j + 3 * 3 * k + 3 * l + i]

            : -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k]))
            : (l >= j ? (i > k ? -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i]
            : -c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k])
            : -c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i])))));
}

inline void RGESolver::WC8I_set(double * c, int i, int j, int k, int l, double
        val) {
    if (i < j) {
        if (k < l) {
            if (i <= k) {
                if (j <= l) {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                } else {
                    c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                }
            } else {
                c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
            }
        } else if (k == l) {
            if (j <= k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else if (i < k) {
                c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << i+1 << ", " << l+1 << ", " << k+1 << ", " << j+1 << std::endl;
            } else {
                c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
            }
        } else { //k>l
            if (k == j) {
                if (i == l) {
                    c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j] = 0.;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << i+1 << ", " << l+1 << ", " << k+1 << ", " << j+1 << " and real" << std::endl;
                } else if (i > l) {
                    c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << i+1 << ", " << j+1 << ", " << k+1 << std::endl;
                } else if (j > l) {
                    c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << i+1 << ", " << l+1 << ", " << k+1 << ", " << j+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                }
            } else if (k > j)
                if (j > l) {
                    c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << i+1 << ", " << j+1 << ", " << k+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
                } else { //k<j
                if (j > l) {
                    if (i > l) {
                        c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k] = - val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << i+1 << ", " << j+1 << ", " << k+1 << std::endl;
                    } else {
                        c[3 * 3 * 3 * i + 3 * 3 * l + 3 * k + j] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << i+1 << ", " << l+1 << ", " << k+1 << ", " << j+1 << std::endl;
                    }
                } else {
                    c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
                }
            }
        }
    } else if (i == j) {
        if (k == l) {
            if (i <= k) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = 0.;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent and real" << std::endl;
            } else {
                c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = 0.;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << " and real" << std::endl;
            }
        } else if (k > l) {
            if (l >= i) {
                c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << i+1 << ", " << l+1 << ", " << k+1 << std::endl;
            } else {
                if (l < j) {
                    if (k < j) {
                        c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
                    } else {
                        c[3 * 3 * 3 * l + 3 * 3 * i + 3 * j + k] = - val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << i+1 << ", " << j+1 << ", " << k+1 << std::endl;
                    }
                } else {
                    c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
                }
            }
        } else { //k<l
            if (k >= j) {
                c[3 * 3 * 3 * i + 3 * 3 * j + 3 * k + l] = val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is independent" << std::endl;
            } else {
                if (k < i) {
                    if (l < i) {
                        c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                    } else {
                        c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << j+1 << ", " << i+1 << ", " << l+1 << std::endl;
                    }
                } else {
                    c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                }
            }
        }
    } else { //i>j
        if (k < l) {
            if (i == l) {
                if (k < j) {
                    c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << j+1 << ", " << i+1 << ", " << l+1 << std::endl;
                } else if (k == j) {
                    if (k < i) {
                        c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] = 0.;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << j+1 << ", " << i+1 << ", " << l+1 << " and real" << std::endl;
                    } else {
                        c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                        //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                    }
                } else {
                    c[3 * 3 * 3 * j + 3 * 3 * k + 3 * i + l] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << k+1 << ", " << i+1 << ", " << l+1 << std::endl;
                }
            } else if (i > l) {
                if (j < l) {
                    c[3 * 3 * 3 * j + 3 * 3 * k + 3 * l + i] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << k+1 << ", " << l+1 << ", " << i+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * k + 3 * 3 * l + 3 * i + j] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << l+1 << ", " << i+1 << ", " << j+1 << std::endl;
                }
            } else {
                if (k < i) {
                    c[3 * 3 * 3 * k + 3 * 3 * j + 3 * i + l] = val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to "  << k+1 << ", " << j+1 << ", " << i+1 << ", " << l+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << i+1 << ", " << l+1 << ", " << k+1 << std::endl;
                }
            }
        } else if (k == l)
            if (j >= k) {
                c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
            } else {
                if (i > k) {
                    c[3 * 3 * 3 * j + 3 * 3 * k + 3 * l + i] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << k+1 << ", " << l+1 << ", " << i+1 << std::endl;
                } else {
                    c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
                    //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << i+1 << ", " << l+1 << ", " << k+1 << std::endl;
                }
            } else
            if (l >= j)
            if (i > k) {
                c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
            } else {
                c[3 * 3 * 3 * j + 3 * 3 * i + 3 * l + k] = - val;
                //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << j+1 << ", " << i+1 << ", " << l+1 << ", " << k+1 << std::endl;
            } else {
            c[3 * 3 * 3 * l + 3 * 3 * k + 3 * j + i] = - val;
            //std::cout << i+1 << ", " << j+1 << ", " << k+1 << ", " << l+1 << " is equal to the conjugate of "  << l+1 << ", " << k+1 << ", " << j+1 << ", " << i+1 << std::endl;
        }
    }
}

