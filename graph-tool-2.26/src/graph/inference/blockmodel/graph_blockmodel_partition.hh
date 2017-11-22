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

#ifndef GRAPH_BLOCKMODEL_PARTITION_HH
#define GRAPH_BLOCKMODEL_PARTITION_HH

#include "../support/util.hh"
#include "../support/int_part.hh"

namespace graph_tool
{

// ===============
// Partition stats
// ===============

constexpr size_t null_group = std::numeric_limits<size_t>::max();

typedef vprop_map_t<std::vector<std::tuple<size_t, size_t, size_t>>>::type
    degs_map_t;

struct simple_degs_t {};

template <class Graph, class Vprop, class Eprop, class F>
void degs_op(size_t v, Vprop& vweight, Eprop& eweight, const simple_degs_t&,
             Graph& g, F&& f)
{
    f(in_degreeS()(v, g, eweight), out_degreeS()(v, g, eweight), vweight[v]);
}

template <class Graph, class Vprop, class Eprop, class F>
void degs_op(size_t v, Vprop& vweight, Eprop& eweight,
             typename degs_map_t::unchecked_t& degs, Graph& g, F&& f)
{
    auto& ks = degs[v];
    if (ks.empty())
    {
        degs_op(v, vweight, eweight, simple_degs_t(), g, std::forward<F>(f));
    }
    else
    {
        for (auto& k : ks)
            f(get<0>(k), get<1>(k), get<2>(k));
    }
}

class partition_stats_t
{
public:

    typedef gt_hash_map<pair<size_t,size_t>, int> map_t;

    template <class Graph, class Vprop, class VWprop, class Eprop, class Degs,
              class Mprop, class Vlist>
    partition_stats_t(Graph& g, Vprop& b, Vlist& vlist, size_t E, size_t B,
                      VWprop& vweight, Eprop& eweight, Degs& degs,
                      const Mprop& ignore_degree, std::vector<size_t>& bmap,
                      bool allow_empty)
        : _bmap(bmap), _N(0), _E(E), _total_B(B), _allow_empty(allow_empty)
    {
        for (auto v : vlist)
        {
            if (vweight[v] == 0)
                continue;

            auto r = get_r(b[v]);
            if (v >= _ignore_degree.size())
                _ignore_degree.resize(v + 1, 0);
            _ignore_degree[v] = ignore_degree[v];

            degs_op(v, vweight, eweight, degs, g,
                    [&](auto kin, auto kout, auto n)
                    {
                        if (_ignore_degree[v] == 2)
                            kout = 0;
                        if (_ignore_degree[v] != 1)
                        {
                            _hist[r][make_pair(kin, kout)] += n;
                            _em[r] += kin * n;
                            _ep[r] += kout * n;
                        }
                        _total[r] += n;
                        _N += n;
                    });
        }

        _actual_B = 0;
        for (auto n : _total)
        {
            if (n > 0)
                _actual_B++;
        }
    }

    size_t get_r(size_t r)
    {
        constexpr size_t null =
            std::numeric_limits<size_t>::max();
        if (r >= _bmap.size())
            _bmap.resize(r + 1, null);
        size_t nr = _bmap[r];
        if (nr == null)
            nr = _bmap[r] = _hist.size();
        if (nr >= _hist.size())
        {
            _hist.resize(nr + 1);
            _total.resize(nr + 1);
            _ep.resize(nr + 1);
            _em.resize(nr + 1);
        }
        return nr;
    }

    double get_partition_dl()
    {
        double S = 0;
        if (_allow_empty)
            S += lbinom(_total_B + _N - 1, _N);
        else
            S += lbinom(_N - 1, _actual_B - 1);
        S += lgamma_fast(_N + 1);
        for (auto nr : _total)
            S -= lgamma_fast(nr + 1);
        S += safelog(_N);
        return S;
    }

    double get_deg_dl_ent()
    {
        double S = 0;
        for (size_t r = 0; r < _ep.size(); ++r)
        {
            size_t total = 0;
            for (auto& k_c : _hist[r])
            {
                S -= xlogx(k_c.second);
                total += k_c.second;
            }
            S += xlogx(total);
        }
        return S;
    }

    double get_deg_dl_uniform()
    {
        double S = 0;
        for (size_t r = 0; r < _ep.size(); ++r)
        {
            S += lbinom(_total[r] + _ep[r] - 1, _ep[r]);
            S += lbinom(_total[r] + _em[r] - 1, _em[r]);
        }
        return S;
    }

    double get_deg_dl_dist()
    {
        double S = 0;
        for (size_t r = 0; r < _ep.size(); ++r)
        {
            S += log_q(_ep[r], _total[r]);
            S += log_q(_em[r], _total[r]);

            size_t total = 0;
            for (auto& k_c : _hist[r])
            {
                S -= lgamma_fast(k_c.second + 1);
                total += k_c.second;
            }
            S += lgamma_fast(total + 1);
        }
        return S;
    }

    double get_deg_dl(int kind)
    {
        switch (kind)
        {
        case deg_dl_kind::ENT:
            return get_deg_dl_ent();
        case deg_dl_kind::UNIFORM:
            return get_deg_dl_uniform();
        case deg_dl_kind::DIST:
            return get_deg_dl_dist();
        default:
            return numeric_limits<double>::quiet_NaN();
        }
    }

    template <class VProp>
    double get_delta_partition_dl(size_t v, size_t r, size_t nr, VProp& vweight)
    {
        if (r == nr)
            return 0;

        if (r != null_group)
            r = get_r(r);

        if (nr != null_group)
            nr = get_r(nr);

        int n = vweight[v];
        if (n == 0)
        {
            if (r == null_group)
                n = 1;
            else
                return 0;
        }

        double S_b = 0, S_a = 0;

        if (r != null_group)
        {
            S_b += -lgamma_fast(_total[r] + 1);
            S_a += -lgamma_fast(_total[r] - n + 1);
        }

        if (nr != null_group)
        {
            S_b += -lgamma_fast(_total[nr] + 1);
            S_a += -lgamma_fast(_total[nr] + n + 1);
        }

        int dN = 0;
        if (r == null_group)
            dN += n;
        if (nr == null_group)
            dN -= n;

        S_b += lgamma_fast(_N + 1);
        S_a += lgamma_fast(_N + dN + 1);

        int dB = 0;
        if (r != null_group && _total[r] == n)
            dB--;
        if (nr != null_group && _total[nr] == 0)
            dB++;

        if ((dN != 0 || dB != 0) && !_allow_empty)
        {
            S_b += lbinom_fast(_N - 1, _actual_B - 1);
            S_a += lbinom_fast(_N - 1 + dN, _actual_B + dB - 1);
        }

        if (dN != 0)
        {
            S_b += safelog(_N);
            S_a += safelog(_N + dN);
        }

        return S_a - S_b;
    }

    template <class VProp, class Graph>
    double get_delta_edges_dl(size_t v, size_t r, size_t nr, VProp& vweight,
                              size_t actual_B, Graph&)
    {
        if (r == nr || _allow_empty)
            return 0;

        if (r != null_group)
            r = get_r(r);
        if (nr != null_group)
            nr = get_r(nr);

        double S_b = 0, S_a = 0;

        int n = vweight[v];

        if (n == 0)
        {
            if (r == null_group)
                n = 1;
            else
                return 0;
        }

        int dB = 0;
        if (r != null_group && _total[r] == n)
            dB--;
        if (nr != null_group && _total[nr] == 0)
            dB++;

        if (dB != 0)
        {
            auto get_x = [](size_t B)
                {
                    if (is_directed::apply<Graph>::type::value)
                        return B * B;
                    else
                        return (B * (B + 1)) / 2;
                };

            S_b += lbinom(get_x(actual_B) + _E - 1, _E);
            S_a += lbinom(get_x(actual_B + dB) + _E - 1, _E);
        }

        return S_a - S_b;
    }

    template <class Graph, class VProp, class EProp, class Degs>
    double get_delta_deg_dl(size_t v, size_t r, size_t nr, VProp& vweight,
                            EProp& eweight, Degs& degs, Graph& g, int kind)
    {
        if (r == nr || _ignore_degree[v] == 1 || vweight[v] == 0)
            return 0;
        if (r != null_group)
            r = get_r(r);
        if (nr != null_group)
            nr = get_r(nr);

        auto dop =
            [&](auto&& f)
            {
                if (_ignore_degree[v] == 2)
                {
                    degs_op(v, vweight, eweight, degs, g,
                            [&](auto kin, auto, auto n)
                            {
                                f(kin, 0, n);
                            });
                }
                else
                {
                    degs_op(v, vweight, eweight, degs, g,
                            [&](auto... k) { f(k...); });
                }
            };

        double dS = 0;
        switch (kind)
        {
        case deg_dl_kind::ENT:
            if (r != null_group)
                dS += get_delta_deg_dl_ent_change(r,  dop, -1);
            if (nr != null_group)
                dS += get_delta_deg_dl_ent_change(nr, dop, +1);
            break;
        case deg_dl_kind::UNIFORM:
            if (r != null_group)
                dS += get_delta_deg_dl_uniform_change(v, r,  dop, -1);
            if (nr != null_group)
                dS += get_delta_deg_dl_uniform_change(v, nr, dop, +1);
            break;
        case deg_dl_kind::DIST:
            if (r != null_group)
                dS += get_delta_deg_dl_dist_change(v, r,  dop, -1);
            if (nr != null_group)
                dS += get_delta_deg_dl_dist_change(v, nr, dop, +1);
            break;
        default:
            dS = numeric_limits<double>::quiet_NaN();
        }
        return dS;
    }

    template <class DegOP>
    double get_delta_deg_dl_ent_change(size_t r, DegOP&& dop, int diff)
    {
        int nr = _total[r];
        auto get_Sk = [&](size_t s, pair<size_t, size_t>& deg, int delta)
            {
                int nd = 0;
                auto iter = _hist[s].find(deg);
                if (iter != _hist[s].end())
                    nd = iter->second;
                assert(nd + delta >= 0);
                return -xlogx(nd + delta);
            };

        double S_b = 0, S_a = 0;
        int dn = 0;

        dop([&](size_t kin, size_t kout, int nk)
            {
                dn += diff * nk;
                auto deg = make_pair(kin, kout);
                S_b += get_Sk(r, deg,         0);
                S_a += get_Sk(r, deg, diff * nk);
            });

        S_b += xlogx(nr);
        S_a += xlogx(nr + dn);

        return S_a - S_b;
    }

    template <class DegOP>
    double get_delta_deg_dl_uniform_change(size_t v, size_t r, DegOP&& dop,
                                           int diff)
    {
        auto get_Se = [&](int dn, int dkin, int dkout)
            {
                double S = 0;
                S += lbinom(_total[r] + dn + _ep[r] - 1 + dkout, _ep[r] + dkout);
                S += lbinom(_total[r] + dn + _em[r] - 1 + dkin,  _em[r] + dkin);
                return S;
            };

        double S_b = 0, S_a = 0;
        int tkin = 0, tkout = 0, n = 0;
        dop([&](auto kin, auto kout, int nk)
            {
                tkin += kin * nk;
                if (_ignore_degree[v] != 2)
                    tkout += kout * nk;
                n += nk;
            });

        S_b += get_Se(       0,           0,            0);
        S_a += get_Se(diff * n, diff * tkin, diff * tkout);
        return S_a - S_b;
    }

    template <class DegOP>
    double get_delta_deg_dl_dist_change(size_t v, size_t r, DegOP&& dop,
                                        int diff)
    {
        auto get_Se = [&](int delta, int kin, int kout)
            {
                double S = 0;
                assert(_total[r] + delta >= 0);
                assert(_em[r] + kin >= 0);
                assert(_ep[r] + kout >= 0);
                S += log_q(_em[r] + kin, _total[r] + delta);
                S += log_q(_ep[r] + kout, _total[r] + delta);
                return S;
            };

        auto get_Sr = [&](int delta)
            {
                assert(_total[r] + delta + 1 >= 0);
                return lgamma_fast(_total[r] + delta + 1);
            };

        auto get_Sk = [&](pair<size_t, size_t>& deg, int delta)
            {
                int nd = 0;
                auto iter = _hist[r].find(deg);
                if (iter != _hist[r].end())
                    nd = iter->second;
                assert(nd + delta >= 0);
                return -lgamma_fast(nd + delta + 1);
            };

        double S_b = 0, S_a = 0;
        int tkin = 0, tkout = 0, n = 0;
        dop([&](size_t kin, size_t kout, int nk)
            {
                tkin += kin * nk;
                if (_ignore_degree[v] != 2)
                    tkout += kout * nk;
                n += nk;

                auto deg = make_pair(kin, kout);
                S_b += get_Sk(deg,         0);
                S_a += get_Sk(deg, diff * nk);
            });

        S_b += get_Se(       0,           0,            0);
        S_a += get_Se(diff * n, diff * tkin, diff * tkout);

        S_b += get_Sr(       0);
        S_a += get_Sr(diff * n);

        return S_a - S_b;
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void change_vertex(size_t v, size_t r, bool deg_corr, Graph& g,
                       VWeight& vweight, EWeight& eweight, Degs& degs,
                       int diff)
    {
        degs_op(v, vweight, eweight, degs, g,
                [&](auto kin, auto kout, auto n)
                {
                    this->change_k(v, r, deg_corr, n, kin, kout, diff);
                });
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void remove_vertex(size_t v, size_t r, bool deg_corr, Graph& g,
                       VWeight& vweight, EWeight& eweight, Degs& degs)
    {
        if (r == null_group || vweight[v] == 0)
            return;
        r = get_r(r);
        change_vertex(v, r, deg_corr, g, vweight, eweight, degs, -1);
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void add_vertex(size_t v, size_t nr, bool deg_corr, Graph& g,
                    VWeight& vweight, EWeight& eweight, Degs& degs)
    {
        if (nr == null_group || vweight[v] == 0)
            return;
        nr = get_r(nr);
        change_vertex(v, nr, deg_corr, g, vweight, eweight, degs, 1);
    }

    void change_k(size_t v, size_t r, bool deg_corr, int vweight,
                  int kin, int kout, int diff)
    {
        if (_total[r] == 0 && diff * vweight > 0)
            _actual_B++;

        if (_total[r] == vweight && diff * vweight < 0)
            _actual_B--;

        _total[r] += diff * vweight;
        _N += diff * vweight;

        assert(_total[r] >= 0);

        if (deg_corr && _ignore_degree[v] != 1)
        {
            if (_ignore_degree[v] == 2)
                kout = 0;
            auto deg = make_pair(kin, kout);
            _hist[r][deg] += diff * vweight;
            _em[r] += diff * deg.first * vweight;
            _ep[r] += diff * deg.second * vweight;
        }
    }

    size_t get_N()
    {
        return _N;
    }

    size_t get_actual_B()
    {
        return _actual_B;
    }

    void add_block()
    {
        _total_B++;
    }


private:
    vector<size_t>& _bmap;
    size_t _N;
    size_t _E;
    size_t _actual_B;
    size_t _total_B;
    bool _allow_empty;
    vector<map_t> _hist;
    vector<int> _total;
    vector<int> _ep;
    vector<int> _em;
    vector<uint8_t> _ignore_degree;
};

} //namespace graph_tool
#endif // GRAPH_BLOCKMODEL_PARTITION_HH
