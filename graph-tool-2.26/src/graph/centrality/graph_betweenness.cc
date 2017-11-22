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

#include "graph_filtering.hh"

#include <boost/python.hpp>
#include <boost/graph/betweenness_centrality.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

template <class Graph, class EdgeBetweenness, class VertexBetweenness>
void normalize_betweenness(const Graph& g,
                           std::vector<size_t>& pivots,
                           EdgeBetweenness edge_betweenness,
                           VertexBetweenness vertex_betweenness,
                           size_t n)
{
    size_t p = pivots.size();
    double pfactor = (p > 1 && n > 2) ? ((p - 1) * (n - 2)) : .0;
    double vfactor = (p > 0 && n > 2) ? (p * (n - 2)) : .0;
    double efactor = (p > 0 && n > 1) ? (p * (n - 1)) : .0;
    if (std::is_convertible<typename graph_traits<Graph>::directed_category,
                            undirected_tag>::value)
    {
        pfactor /= 2;
        vfactor /= 2;
        efactor /= 2;
    }

    pfactor = (pfactor > 0) ? 1./pfactor : 0;
    vfactor = (vfactor > 0) ? 1./vfactor : 0;
    efactor = (efactor > 0) ? 1./efactor : 0;

    typedef typename vprop_map_t<bool>::type::unchecked_t vprop_t;
    vprop_t is_pivot(get(vertex_index, g), num_vertices(g));

    parallel_loop(pivots, [&](size_t, auto v){ is_pivot[v] = true;});

    parallel_vertex_loop
        (g,
         [&](auto v)
         {
             if (is_pivot[v])
                 put(vertex_betweenness, v, pfactor * get(vertex_betweenness, v));
             else
                 put(vertex_betweenness, v, vfactor * get(vertex_betweenness, v));
         });

    parallel_edge_loop
        (g,
         [&](const auto& e)
         {
             put(edge_betweenness, e, efactor * get(edge_betweenness, e));
         });
}

struct get_betweenness
{
    typedef void result_type;
    template <class Graph, class EdgeBetweenness, class VertexBetweenness>
    void operator()(Graph& g,
                    std::vector<size_t>& pivots,
                    GraphInterface::vertex_index_map_t index_map,
                    EdgeBetweenness edge_betweenness,
                    VertexBetweenness vertex_betweenness,
                    bool normalize, size_t n) const
    {
        vector<vector<typename graph_traits<Graph>::edge_descriptor> >
            incoming_map(num_vertices(g));
        vector<size_t> distance_map(num_vertices(g));
        vector<typename property_traits<VertexBetweenness>::value_type>
            dependency_map(num_vertices(g));
        vector<size_t> path_count_map(num_vertices(g));
        brandes_betweenness_centrality
            (g, pivots, vertex_betweenness, edge_betweenness,
             make_iterator_property_map(incoming_map.begin(), index_map),
             make_iterator_property_map(distance_map.begin(), index_map),
             make_iterator_property_map(dependency_map.begin(), index_map),
             make_iterator_property_map(path_count_map.begin(), index_map),
             index_map);
        if (normalize)
            normalize_betweenness(g, pivots, edge_betweenness, vertex_betweenness, n);
    }
};

struct get_weighted_betweenness
{
    typedef void result_type;
    template <class Graph, class EdgeBetweenness, class VertexBetweenness,
              class VertexIndexMap>
    void operator()(Graph& g, std::vector<size_t>& pivots,
                    VertexIndexMap vertex_index,
                    EdgeBetweenness edge_betweenness,
                    VertexBetweenness vertex_betweenness,
                    boost::any weight_map, bool normalize,
                    size_t n, size_t max_eindex) const
    {
        vector<vector<typename graph_traits<Graph>::edge_descriptor> >
            incoming_map(num_vertices(g));
        vector<typename property_traits<EdgeBetweenness>::value_type>
            distance_map(num_vertices(g));
        vector<typename property_traits<VertexBetweenness>::value_type>
            dependency_map(num_vertices(g));
        vector<size_t> path_count_map(num_vertices(g));

        typename EdgeBetweenness::checked_t weight =
            any_cast<typename EdgeBetweenness::checked_t>(weight_map);

        brandes_betweenness_centrality
            (g, pivots, vertex_betweenness, edge_betweenness,
             make_iterator_property_map(incoming_map.begin(), vertex_index),
             make_iterator_property_map(distance_map.begin(), vertex_index),
             make_iterator_property_map(dependency_map.begin(), vertex_index),
             make_iterator_property_map(path_count_map.begin(), vertex_index),
             vertex_index, weight.get_unchecked(max_eindex+1));
        if (normalize)
            normalize_betweenness(g, pivots, edge_betweenness, vertex_betweenness, n);
    }
};

void betweenness(GraphInterface& g, std::vector<size_t>& pivots,
                 boost::any weight,
                 boost::any edge_betweenness,
                 boost::any vertex_betweenness,
                 bool normalize)
{
    if (!belongs<edge_floating_properties>()(edge_betweenness))
        throw ValueException("edge property must be of floating point value"
                             " type");

    if (!belongs<vertex_floating_properties>()(vertex_betweenness))
        throw ValueException("vertex property must be of floating point value"
                             " type");

    if (!weight.empty())
    {
        run_action<>()
            (g, std::bind<>(get_weighted_betweenness(),
                            std::placeholders::_1,
                            std::ref(pivots),
                            g.get_vertex_index(),
                            std::placeholders::_2,
                            std::placeholders::_3, weight, normalize,
                            g.get_num_vertices(), g.get_edge_index_range()),
             edge_floating_properties(),
             vertex_floating_properties())
            (edge_betweenness, vertex_betweenness);
    }
    else
    {
        run_action<>()
            (g, std::bind<void>(get_betweenness(), std::placeholders::_1,
                                std::ref(pivots),
                                g.get_vertex_index(), std::placeholders::_2,
                                std::placeholders::_3, normalize,
                                g.get_num_vertices()),
             edge_floating_properties(),
             vertex_floating_properties())
            (edge_betweenness, vertex_betweenness);
    }
}

struct get_central_point_dominance
{
    template <class Graph, class VertexBetweenness>
    void operator()(Graph& g, VertexBetweenness vertex_betweenness, double& c)
        const
    {
        c = double(central_point_dominance(g, vertex_betweenness));
    }
};

double central_point(GraphInterface& g,
                     boost::any vertex_betweenness)
{
    double c = 0.0;
    run_action<graph_tool::detail::never_reversed>()
        (g, std::bind<>(get_central_point_dominance(), std::placeholders::_1,
                        std::placeholders::_2, std::ref(c)),
         vertex_scalar_properties()) (vertex_betweenness);
    return c;
}

void export_betweenness()
{
    using namespace boost::python;
    def("get_betweenness", &betweenness);
    def("get_central_point_dominance", &central_point);
}
