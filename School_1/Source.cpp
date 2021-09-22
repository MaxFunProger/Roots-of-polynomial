#include <iostream>
#include <algorithm>
#include <vector>
#include <string.h>
#include <string>
#include <map>
#include <set>
#include <iomanip> 
#include <math.h>


const double kEps = 1e-9;
const double kEpsFunc = 1e-9;
const double kInf = 1e9;
int nn, mm;
double aa, bb, cc, dd, pp;
std::vector<double> input_vector;

double func(double xx, std::vector<double>& polynom) {
    double res = polynom[0];
    for (int i = 1; i < polynom.size(); ++i) {
        res *= xx;
        res += polynom[i];
    }
    return res;
}

bool check_sign(double aa, double bb) {
    if (aa > 0 && bb > 0 || aa < 0 && bb < 0)
        return true;
    return false;
}

double root_finder(double aa, double bb, std::vector<double>& polynom) {
    double right_val = bb;
    if (std::fabs(func(aa, polynom)) < kEps) {
        return aa;
    }
    if (std::fabs(func(bb, polynom)) < kEps) {
        return bb;
    }
    if (check_sign(func(aa, polynom), func(bb, polynom)))
        return kInf;
    double cc = (aa + bb) / 2.0;
    while (std::fabs(aa - bb) > kEps) {
        cc = (aa + bb) / 2;
        if (std::fabs(func(aa, polynom)) < kEps)
            return aa;
        if (std::fabs(func(bb, polynom)) < kEps)
            return bb;
        if (check_sign(func(aa, polynom), func(cc, polynom))) {
            aa = cc;
        } else {
            bb = cc;
        }
    }
    cc = (aa + bb) / 2;
    return cc;
}

std::vector<double> solve(std::vector<double> polynom) {
    if (polynom.size() == 2) {
        return { -polynom[1] / polynom[0] };
    }
    std::vector<double> derivative;
    int sz = polynom.size();
    for (int i = 0; i < sz - 1; ++i) {
        derivative.push_back(polynom[i] * (sz - 1 - i));
    }
    std::vector<double> res = solve(derivative);
    std::vector<double> roots;
    if (res.empty()) {
        roots.push_back(root_finder(aa, bb, polynom));
        if (roots.back() == kInf) {
            roots.pop_back();
        }
    } else {
        if (res[0] != aa) {
            roots.push_back(root_finder(aa, res[0], polynom));
            if (roots.back() == kInf) {
                roots.pop_back();
            }
        }
        for (int i = 0; i < res.size() - 1; ++i) {
            roots.push_back(root_finder(res[i], res[i + 1], polynom));
            if (roots.back() == kInf) {
                roots.pop_back();
            }
        }
        if (res.back() != bb) {
            roots.push_back(root_finder(res.back(), bb, polynom));
            if (roots.back() == kInf) {
                roots.pop_back();
            }
        }
    }
    return roots;
}


int main() {
    std::cin >> nn;
    std::cout.precision(10);
    for (int i = 0; i < nn + 1; ++i) {
        std::cin >> pp;
        input_vector.push_back(pp);
    }
    std::cin >> aa >> bb;
    std::vector<double> res;
    res = solve(input_vector);
    for (double& i : res) {
        if (std::fabs(i - bb) > kEps && i < bb) {
            std::cout << std::fixed << i << " ";
        }
    }
    
    
    return 0;
}
