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

#ifndef URN_SAMPLER_HH
#define URN_SAMPLER_HH

#include "random.hh"

namespace graph_tool
{
using namespace std;
using namespace boost;

// Sampling with and without replacement from an urn

template <class Value, bool replacement>
class UrnSampler
{
public:
    UrnSampler()  {}

    template <class Int>
    UrnSampler(const vector<Value>& items,
               const vector<Int>& counts)
    {
        for (size_t i = 0; i < items.size(); ++i)
        {
            const auto& x = items[i];
            size_t n = counts[i];
            for (size_t j = 0; j < n; ++j)
                _urn.push_back(x);
        }
    }

    template <class RNG>
    typename std::conditional<replacement,
                              const Value&,
                              Value>::type
    sample(RNG& rng)
    {
        uniform_int_distribution<size_t> sample(0, _urn.size() - 1);
        size_t i = sample(rng);
        if (replacement)
        {
            return _urn[i];
        }
        else
        {
            std::swap(_urn[i], _urn.back());
            Value temp = _urn.back();
            _urn.pop_back();
            return temp;
        }
    }

    bool empty()
    {
        return _urn.empty();
    }

private:
    vector<Value> _urn;
};

} // namespace graph_tool

#endif // URN_SAMPLER_HH
