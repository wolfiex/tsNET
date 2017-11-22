// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2017 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef GRAPH_BLOCKMODEL_WEIGHTS_HH
#define GRAPH_BLOCKMODEL_WEIGHTS_HH

namespace graph_tool
{

// Weighted entropy terms
// ======================

enum weight_type
{
    NONE,
    COUNT,
    REAL_EXPONENTIAL,
    REAL_NORMAL,
    DISCRETE_GEOMETRIC,
    DISCRETE_POISSON,
    DISCRETE_BINOMIAL,
    DELTA_T
};

// exponential
template <class DT>
double positive_w_log_P(DT N, double x, double alpha, double beta,
                        double epsilon)
{
    if (N == 0)
        return 0.;
    if (std::isnan(alpha) && std::isnan(beta))
    {
        if (x < epsilon || N == 1)
            return 0.;
        else
            return lgamma(N) - (N - 1) * log(x);
    }
    return lgamma(N + alpha) - lgamma(alpha) + alpha * log(beta) -
        (alpha + N) * log(beta + x);
}

// normal
template <class DT>
double signed_w_log_P(DT N, double x, double x2, double m0, double k0, double v0,
                      double nu0, double epsilon)
{
    if (N == 0)
        return 0.;
    if (std::isnan(m0) && std::isnan(k0))
    {
        auto smu1 = x * (x / N);

        if (N < 2 || smu1 >= x2 || (x2 - smu1) < std::pow(epsilon, 2))
            return 0.;
        else
            return (lgamma((N - 1) / 2.) + log(N) / 2.
                    - ((int(N) - 3) / 2.) * log(x2 - smu1)
                    - ((N - 1) / 2.) * log(M_PI));
    }
    auto v = x2 - x * (x / N);
    auto k_n = k0 + N;
    auto nu_n = nu0 + N;
    auto v_n = (v0 * nu0 + v + ((N * k0)/(k0 + N)) * pow(m0 - x/N, 2)) / nu_n;
    return lgamma(nu_n / 2.) - lgamma(nu0 / 2.) + (log(k0) - log(k_n)) / 2. +
        (nu0 / 2.) * log(nu0 * v0) - (nu_n / 2.) * log(nu_n * v_n) - (N / 2.) *
        log(M_PI);
}

// discrete: geometric
template <class DT>
double geometric_w_log_P(DT N, double x, double alpha, double beta)
{
    if (N == 0)
        return 0.;
    if (std::isnan(alpha) && std::isnan(beta))
        return -lbinom((N - 1) + x, x);
    return lbeta(N + alpha, x + beta) - lbeta(alpha, beta);
}

// discrete: binomial
template <class DT>
double binomial_w_log_P(DT N, double x, int n, double alpha, double beta)
{
    if (N == 0)
        return 0.;
    if (std::isnan(alpha) && std::isnan(beta))
        return -lbinom(N * n, x);
    return lbeta(x + alpha, N * n - x + beta) - lbeta(alpha, beta);
}

// discrete: Poisson
template <class DT>
double poisson_w_log_P(DT N, double x, double alpha, double beta)
{
    if (N == 0)
        return 0.;
    if (std::isnan(alpha) && std::isnan(beta))
        return lgamma(x+1) - x * log(N);
    return lgamma(x + alpha) - (x + alpha) * log(N + beta) - lgamma(alpha) +
        alpha * log(beta);
}

} //namespace graph_tool

#endif // GRAPH_BLOCKMODEL_WEIGHTS_HH
