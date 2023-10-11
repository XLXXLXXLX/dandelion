#include "halfedge.h"

#include <cstddef>
#include <set>
#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>

using Eigen::Matrix3f;
using Eigen::Matrix4f;
using Eigen::Vector3f;
using Eigen::Vector4f;
using std::optional;
using std::set;
using std::size_t;
using std::string;
using std::unordered_map;
using std::vector;

HalfedgeMesh::EdgeRecord::EdgeRecord(unordered_map<Vertex*, Matrix4f>& vertex_quadrics, Edge* e)
    : edge(e)
{
    (void)vertex_quadrics;
    optimal_pos = Vector3f(0.0f, 0.0f, 0.0f);
    cost        = 0.0f;
}

bool operator<(const HalfedgeMesh::EdgeRecord& a, const HalfedgeMesh::EdgeRecord& b)
{
    return a.cost < b.cost;
}

optional<Edge*> HalfedgeMesh::flip_edge(Edge* e)
{
    // 要用到的半边
    Halfedge* h12 = e->halfedge;
    Halfedge* h21 = h12->inv;
    Halfedge* h23 = h12->next;
    Halfedge* h31 = h23->next;
    Halfedge* h14 = h21->next;
    Halfedge* h42 = h14->next;
    auto h43      = h12;
    auto h34      = h21;
    // 要用到的顶点
    // v1 and v2 are vertices along the edge
    Vertex* v1 = h12->from;
    Vertex* v2 = h21->from;
    // v3 and v4 are vertices opposite the edge
    Vertex* v3 = h31->from;
    Vertex* v4 = h42->from;
    // 要用到的面片
    Face* f1 = h12->face;
    Face* f2 = h21->face;

    f1 = new_face(f1->is_boundary || f2->is_boundary);
    f2 = new_face(f1->is_boundary || f2->is_boundary);

    f1->halfedge = h12;
    f2->halfedge = h21;
    e->halfedge  = h12;
    // 重新连接各基本元素
    v1->halfedge = h14;
    v2->halfedge = h23;
    h43->set_neighbors(h31, h14, h21, v4, e, f1);
    h34->set_neighbors(h42, h23, h12, v3, e, f2);
    h31->set_neighbors(h14, h43, h31->inv, v3, h31->edge, f1);
    h14->set_neighbors(h43, h31, h14->inv, v1, h14->edge, f1);
    h42->set_neighbors(h23, h34, h42->inv, v4, h42->edge, f2);
    h23->set_neighbors(h34, h42, h23->inv, v2, h23->edge, f2);

    return e;
    // }
}

optional<Vertex*> HalfedgeMesh::split_edge(Edge* e)
{
    // 要用到的半边
    auto h12 = e->halfedge;
    auto h21 = h12->inv;
    auto h23 = h12->next;
    auto h31 = h23->next;
    auto h14 = h21->next;
    auto h42 = h14->next;
    auto v1  = h12->from;
    auto v2  = h21->from;
    auto v3  = h31->from;
    auto v4  = h42->from;
    auto f1  = h12->face;
    auto f2  = h21->face;

    auto vn = new_vertex();
    vn->pos = e->center();

    auto h1n  = h12;
    auto h2n  = h21;
    auto h3n  = new_halfedge();
    auto h4n  = new_halfedge();
    auto hn1  = new_halfedge();
    auto hn2  = new_halfedge();
    auto hn3  = new_halfedge();
    auto hn4  = new_halfedge();
    auto e1n  = e;
    auto e2n  = new_edge();
    auto e3n  = new_edge();
    auto e4n  = new_edge();
    auto fn31 = new_face(f1->is_boundary);
    auto fn23 = new_face(f1->is_boundary);
    auto fn14 = new_face(f2->is_boundary);
    auto fn42 = new_face(f2->is_boundary);
    // 初始化edge
    e1n->is_new   = false;
    e2n->is_new   = false;
    e3n->is_new   = true;
    e4n->is_new   = true;
    e1n->halfedge = hn1;
    e2n->halfedge = hn2;
    e3n->halfedge = hn3;
    e4n->halfedge = hn4;
    // 初始化face
    fn31->halfedge = h31;
    fn14->halfedge = h14;
    fn23->halfedge = h23;
    fn42->halfedge = h42;
    // 连接halfedge
    // h->set_netghbors(next, prev, inv, vertex, edge, face);
    hn1->set_neighbors(h14, h4n, h1n, vn, e1n, fn14);
    hn2->set_neighbors(h23, h3n, h2n, vn, e2n, fn23);
    hn3->set_neighbors(h31, h1n, h3n, vn, e3n, fn31);
    hn4->set_neighbors(h42, h2n, h4n, vn, e4n, fn42);
    h1n->set_neighbors(hn3, h31, hn1, v1, e1n, fn31);
    h2n->set_neighbors(hn4, h42, hn2, v2, e2n, fn42);
    h3n->set_neighbors(hn2, h23, hn3, v3, e3n, fn23);
    h4n->set_neighbors(hn1, h14, hn4, v4, e4n, fn14);
    h14->set_neighbors(h4n, hn1, h14->inv, v1, h14->edge, fn14);
    h23->set_neighbors(h3n, hn2, h23->inv, v2, h23->edge, fn23);
    h31->set_neighbors(h1n, hn3, h31->inv, v3, h31->edge, fn31);
    h42->set_neighbors(h2n, hn4, h42->inv, v4, h42->edge, fn42);

    vn->halfedge = hn1;
    vn->is_new   = true;
    return vn;
}

optional<Vertex*> HalfedgeMesh::collapse_edge(Edge* e)
{
    // 要用到的半边
    auto h12 = e->halfedge;
    auto h21 = h12->inv;
    auto h23 = h12->next, h32 = h23->inv;
    auto h31 = h12->prev, h13 = h31->inv;
    auto h14 = h21->next, h41 = h14->inv;
    auto h42 = h21->prev, h24 = h42->inv;
    erase(e);
    erase(h12);
    erase(h21);
    erase(h23);
    erase(h31);
    erase(h14);
    erase(h42);
    auto h1 = h13, h2 = h24;
    auto hinv1 = h32, hinv2 = h41;
    auto hfrom2 = hinv1->next;
    auto hfrom1 = hinv2->next;
    // 要用到的边
    auto e1n = h1->edge;
    auto e2n = h2->edge;

    e1n->halfedge = h1;
    e2n->halfedge = h2;
    // 要用到的顶点
    auto v1 = h13->from;
    auto v2 = h24->from;
    auto v3 = h32->from;
    auto v4 = h41->from;
    erase(v1);
    erase(v2);
    v3->halfedge = hinv1;
    v4->halfedge = hinv2;
    auto vn      = new_vertex();
    vn->pos      = (v3->pos + v4->pos) / 2;
    vn->halfedge = h1;
    h1->set_neighbors(h1->next, h1->prev, hinv1, vn, e1n, h1->face);
    h2->set_neighbors(h2->next, h2->prev, hinv2, vn, e2n, h2->face);
    hinv1->set_neighbors(hinv1->next, hinv1->prev, h1, v3, e1n, hinv1->face);
    hinv2->set_neighbors(hinv2->next, hinv2->prev, h2, v4, e2n, hinv2->face);
    while (hfrom1 != h1) {
        hfrom1->from = vn;
        hfrom1       = hfrom1->inv->next;
    }
    while (hfrom2 != h2) {
        hfrom2->from = vn;
        hfrom2       = hfrom2->inv->next;
    }
    return vn;
}

void HalfedgeMesh::loop_subdivide()
{
    optional<HalfedgeMeshFailure> check_result = validate();
    if (check_result.has_value()) {
        return;
    }
    logger->info("subdivide object {} (ID: {}) with Loop Subdivision strategy", object.name,
                 object.id);
    logger->info("original mesh: {} vertices, {} faces in total", vertices.size, faces.size);

    // Each vertex and edge of the original mesh can be associated with a vertex
    // in the new (subdivided) mesh.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions using the connectivity of the original
    // (coarse) mesh. Navigating this mesh will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.
    // We will then assign vertex positions in the new mesh based on the values
    // we computed for the original mesh.

    // Compute new positions for all the vertices   in the input mesh using
    // the Loop subdivision rule and store them in Vertex::new_pos.
    //    At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.
    for (auto v = vertices.head; v != nullptr; v = v->next_node) {
        auto n = v->degree();
        double u;
        if (n == 3) {
            u = 3.0f / 16.0f;
        } else {
            u = 3.0f / (8.0f * n);
        }
        v->new_pos = (1 - u * n) * v->pos + u * (v->neighborhood_center() * n);
    }
    // Next, compute the subdivided vertex positions associated with edges, and
    // store them in Edge::new_pos.
    for (auto e = edges.head; e != nullptr; e = e->next_node) {
        vector<Edge*> original_edges;
        original_edges.push_back(e);
        e->new_pos = 3 / 8 * (e->halfedge->from->pos + e->halfedge->inv->from->pos) +
                     1 / 8 * (e->halfedge->prev->from->pos + e->halfedge->inv->prev->from->pos);
    }

    // Next, we're going to split every edge in the mesh, in any order.
    // We're also going to distinguish subdivided edges that came from splitting
    // an edge in the original mesh from new edges by setting the boolean Edge::is_new.
    // Note that in this loop, we only want to iterate over edges of the original mesh.
    // Otherwise, we'll end up splitting edges that we just split (and the
    // loop will never end!)
    // I use a vector to store iterators of original because there are three kinds of
    // edges: original edges, edges split from original edges and newly added edges.
    // The newly added edges are marked with Edge::is_new property, so there is not
    // any other property to mark the edges I just split.

    for (auto f = faces.head; f != nullptr; f = f->next_node) {

        auto e1 = f->halfedge->edge;
        auto e2 = f->halfedge->next->edge;
        auto e3 = f->halfedge->prev->edge;

        Vertex* vn = nullptr;
        if (!e1->is_new) {
            vn = split_edge(e1).value_or(nullptr);
        }
        if (!e2->is_new) {
            split_edge(e2);
        }
        if (!e3->is_new) {
            split_edge(e3);
        }
        if (vn != nullptr) {
            flip_edge(vn->halfedge->edge);
        }
    }

    // Now flip any new edge that connects an old and new vertex.

    // Finally, copy new vertex positions into the Vertex::pos.

    // Once we have successfully subdivided the mesh, set global_inconsistent
    // to true to trigger synchronization with GL::Mesh.

    global_inconsistent = true;
    logger->info("subdivided mesh: {} vertices, {} faces in total", vertices.size, faces.size);
    logger->info("Loop Subdivision done");
    logger->info("");
    validate();
}

void HalfedgeMesh::simplify()
{
    optional<HalfedgeMeshFailure> check_result = validate();
    if (check_result.has_value()) {
        return;
    }
    logger->info("simplify object {} (ID: {})", object.name, object.id);
    logger->info("original mesh: {} vertices, {} faces", vertices.size, faces.size);
    unordered_map<Vertex*, Matrix4f> vertex_quadrics;
    unordered_map<Face*, Matrix4f> face_quadrics;
    unordered_map<Edge*, EdgeRecord> edge_records;
    set<EdgeRecord> edge_queue;

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics

    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics

    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.

    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    logger->info("simplified mesh: {} vertices, {} faces", vertices.size, faces.size);
    logger->info("simplification done\n");
    global_inconsistent = true;
    validate();
}

void HalfedgeMesh::isotropic_remesh()
{
    optional<HalfedgeMeshFailure> check_result = validate();
    if (check_result.has_value()) {
        return;
    }
    logger->info("remesh the object {} (ID: {}) with strategy Isotropic Remeshing", object.name,
                 object.id);
    logger->info("original mesh: {} vertices, {} faces", vertices.size, faces.size);
    // Compute the mean edge length.

    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions
    static const size_t iteration_limit = 5;
    set<Edge*> selected_edges;
    for (size_t i = 0; i != iteration_limit; ++i) {
        // Split long edges.

        // Collapse short edges.

        // Flip edges.

        // Vertex averaging.
    }
    logger->info("remeshed mesh: {} vertices, {} faces\n", vertices.size, faces.size);
    global_inconsistent = true;
    validate();
}
