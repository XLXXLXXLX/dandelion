#include "halfedge.h"

#include <cstddef>
#include <iterator>
#include <set>
#include <map>
#include <sys/types.h>
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
    if (e->on_boundary())
        return optional<Edge*>{};
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

    auto fn1 = new_face(f1->is_boundary || f2->is_boundary);
    auto fn2 = new_face(f1->is_boundary || f2->is_boundary);
    erase(f1), erase(f2);
    fn1->halfedge = h12;
    fn2->halfedge = h21;
    e->halfedge   = h12;
    // 重新连接各基本元素
    v1->halfedge = h14;
    v2->halfedge = h23;
    h43->set_neighbors(h31, h14, h34, v4, e, fn1);
    h34->set_neighbors(h42, h23, h43, v3, e, fn2);
    h31->set_neighbors(h14, h43, h31->inv, v3, h31->edge, fn1);
    h14->set_neighbors(h43, h31, h14->inv, v1, h14->edge, fn1);
    h42->set_neighbors(h23, h34, h42->inv, v4, h42->edge, fn2);
    h23->set_neighbors(h34, h42, h23->inv, v2, h23->edge, fn2);

    return e;
    // }
}

optional<Vertex*> HalfedgeMesh::split_edge(Edge* e)
{
    auto h12      = e->halfedge;
    auto h21      = h12->inv;
    auto v1       = h12->from;
    auto v2       = h21->from;
    auto vn       = new_vertex();
    vn->pos       = e->center();
    vn->is_new    = true;
    auto h1n      = h12;
    auto h2n      = h21;
    auto hn1      = new_halfedge();
    auto hn2      = new_halfedge();
    auto e1n      = e;
    auto e2n      = new_edge();
    e1n->is_new   = false;
    e2n->is_new   = false;
    e1n->halfedge = hn1;
    e2n->halfedge = hn2;
    if (h21->is_boundary()) {
        logger->info("boundary......");
        // 如果h21在边上，那么h21->face即为虚假面
        hn1->set_neighbors(h21->next, h2n, h1n, vn, e1n, h21->face);
        h2n->set_neighbors(hn1, h21->prev, hn2, v2, e2n, h21->face);
        h21->face->halfedge = hn1;
    } else {
        // 如果h21不在边上，则正常处理
        auto h14      = h21->next;
        auto h42      = h14->next;
        auto v4       = h42->from;
        auto f2       = h21->face;
        auto h4n      = new_halfedge();
        auto hn4      = new_halfedge();
        auto fn14     = new_face(f2->is_boundary);
        auto fn42     = new_face(f2->is_boundary);
        auto e4n      = new_edge();
        e4n->is_new   = true;
        e4n->halfedge = hn4;
        erase(f2);
        fn14->halfedge = h14;
        fn42->halfedge = h42;
        hn1->set_neighbors(h14, h4n, h1n, vn, e1n, fn14);
        hn4->set_neighbors(h42, h2n, h4n, vn, e4n, fn42);
        h2n->set_neighbors(hn4, h42, hn2, v2, e2n, fn42);
        h4n->set_neighbors(hn1, h14, hn4, v4, e4n, fn14);
        h14->set_neighbors(h4n, hn1, h14->inv, v1, h14->edge, fn14);
        h42->set_neighbors(h2n, hn4, h42->inv, v4, h42->edge, fn42);
    }
    if (h12->is_boundary()) {
        logger->info("boundary......");
        h1n->set_neighbors(hn2, h12->prev, hn1, v1, e1n, h12->face);
        hn2->set_neighbors(h12->next, h1n, h2n, vn, e2n, h12->face);
        h12->face->halfedge = h1n;
    } else {
        auto h23      = h12->next;
        auto h31      = h23->next;
        auto v3       = h31->from;
        auto f1       = h12->face;
        auto h3n      = new_halfedge();
        auto hn3      = new_halfedge();
        auto fn31     = new_face(f1->is_boundary);
        auto fn23     = new_face(f1->is_boundary);
        auto e3n      = new_edge();
        e3n->is_new   = true;
        e3n->halfedge = hn3;
        erase(f1);
        fn31->halfedge = h31;
        fn23->halfedge = h23;
        hn2->set_neighbors(h23, h3n, h2n, vn, e2n, fn23);
        hn3->set_neighbors(h31, h1n, h3n, vn, e3n, fn31);
        h1n->set_neighbors(hn3, h31, hn1, v1, e1n, fn31);
        h3n->set_neighbors(hn2, h23, hn3, v3, e3n, fn23);
        h23->set_neighbors(h3n, hn2, h23->inv, v2, h23->edge, fn23);
        h31->set_neighbors(h1n, hn3, h31->inv, v3, h31->edge, fn31);
    }

    // 初始化face
    // 连接halfedge
    vn->halfedge = hn1;
    return vn;
}

optional<Vertex*> HalfedgeMesh::collapse_edge(Edge* e)
{
    logger->info("collapse begins");
    // values to be destroyed
    auto e12  = e;
    auto h12  = e->halfedge;
    auto h21  = h12->inv;
    auto h23  = h12->next;
    auto h31  = h12->prev;
    auto h14  = h21->next;
    auto h42  = h21->prev;
    auto v1   = h12->from;
    auto v2   = h21->from;
    auto f123 = h12->face;
    auto f124 = h21->face;
    logger->info("values to be adjusted");
    // values to be adjusted
    auto h1    = h31->inv;
    auto h2    = h42->inv;
    auto hinv1 = h23->inv;
    auto hinv2 = h14->inv;
    auto e1    = h1->edge;
    auto e2    = h2->edge;
    // values to be destroyed
    auto ed1 = hinv1->edge;
    auto ed2 = hinv2->edge;
    // values to be created
    logger->info("new_vertex");
    auto vn    = new_vertex();
    vn->is_new = true;
    // values to be used
    auto v3      = hinv1->from;
    auto v4      = hinv2->from;
    auto hfromv1 = v1->halfedge;
    auto hfromv2 = v2->halfedge;

    // adjustments
    logger->info("adjustments");
    auto v1degree = v1->degree();
    auto v2degree = v2->degree();
    logger->info("{} {}", v1degree, v2degree);
    for (size_t i = 0; i < v1degree; i++) {
        logger->info("hfromv1");
        hfromv1->from = vn;
        hfromv1       = hfromv1->inv->next;
    }
    logger->info("done");
    for (size_t i = 0; i < v2degree; i++) {
        logger->info("hfromv2");
        hfromv2->from = vn;
        hfromv2       = hfromv2->inv->next;
    }
    logger->info("done");
    if (h12->is_boundary()) {
        logger->info("boundary......");
        h12->face->halfedge = h12->next;
        vn->pos             = e->center();
        vn->halfedge        = h2;
        h2->set_neighbors(h2->next, h2->prev, hinv2, vn, e2, h2->face);
        hinv2->set_neighbors(hinv2->next, hinv2->prev, h2, hinv2->from, e2, hinv2->face);
        erase(h21), erase(h14), erase(h42);
        erase(e12), erase(ed2);
        erase(f124);
    } else if (h21->is_boundary()) {
        logger->info("boundary......");
        h21->face->halfedge = h21->next;
        vn->pos             = e->center();
        vn->halfedge        = h1;
        h1->set_neighbors(h1->next, h1->prev, hinv1, vn, e1, h1->face);
        hinv1->set_neighbors(hinv1->next, hinv1->prev, h1, hinv1->from, e1, hinv1->face);
        erase(h12), erase(h23), erase(h31);
        erase(e12), erase(ed1);
        erase(f123);
    } else {
        h1->set_neighbors(h1->next, h1->prev, hinv1, vn, e1, h1->face);
        h2->set_neighbors(h2->next, h2->prev, hinv2, vn, e2, h2->face);
        hinv1->set_neighbors(hinv1->next, hinv1->prev, h1, hinv1->from, e1, hinv1->face);
        hinv2->set_neighbors(hinv2->next, hinv2->prev, h2, hinv2->from, e2, hinv2->face);
        e1->halfedge = h1;
        e2->halfedge = h2;
        vn->pos      = (v3->pos + v4->pos) / 2;
        vn->halfedge = h1;
        // destroy values
        logger->info("destroying...");
        erase(h12), erase(h21), erase(h23), erase(h31), erase(h14), erase(h42);
        erase(ed1), erase(ed2), erase(e12);
        erase(f123), erase(f124);
    }
    erase(v1), erase(v2);
    logger->info("return");
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
            u = 3.0f / 16;
        } else {
            u = 3.0f / (8 * n);
        }

        v->new_pos = (1 - u * n) * v->pos + u * (v->neighborhood_center() * n);
    }
    // Next, compute the subdivided vertex positions associated with edges, and
    // store them in Edge::new_pos.
    vector<Edge*> original_edges;
    for (auto e = edges.head; e != nullptr; e = e->next_node) {
        original_edges.push_back(e);
        e->new_pos = 3.0f / 8 * (e->halfedge->from->pos + e->halfedge->inv->from->pos) +
                     1.0f / 8 * (e->halfedge->prev->from->pos + e->halfedge->inv->prev->from->pos);
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
    for (auto e : original_edges) {
        auto new_pos = e->new_pos;
        auto v       = split_edge(e).value_or(nullptr);
        if (v) {
            v->new_pos = new_pos;
        }
    }
    // Now flip any new edge that connects an old and new vertex.
    for (auto e = edges.head; e != nullptr; e = e->next_node) {
        if (e->is_new && (e->halfedge->from->is_new != e->halfedge->inv->from->is_new)) {
            flip_edge(e);
        }
    }

    // Finally, copy new vertex positions into the Vertex::pos.
    for (auto v = vertices.head; v != nullptr; v = v->next_node) v->pos = v->new_pos;
    for (auto e = edges.head; e != nullptr; e = e->next_node) e->is_new = false;
    for (auto v = vertices.head; v != nullptr; v = v->next_node) v->is_new = false;

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
    // 计算平均边长
    // 重复以下步骤5～6次
    // - 将比目标长度长得多的边分开（写循环时要注意
    // - 将比目标长度短得多的边坍缩掉。这里循环需要非常小心，因为许多边可能已经被摧毁掉了
    // - 翻转每个会增加结点度数的边
    // - 最后对节点位置作优化
    static const size_t iteration_limit = 5;
    set<Edge*> selected_edges;
    double average_edge_length = 0;
    for (auto e = edges.head; e != nullptr; e = e->next_node) {
        selected_edges.insert(e);
        average_edge_length += e->length();
        logger->info("adding edges...:{}", (e)->length());
    }
    average_edge_length /= edges.size;
    auto up_lim   = average_edge_length * 4.0f / 3;
    auto down_lim = average_edge_length * 4.0f / 5;
    for (size_t i = 0; i != iteration_limit; ++i) {
        vector<Edge*> save_delete_edges;
        vector<Edge*> save_add_edges;
        for (auto pe = selected_edges.begin(); pe != selected_edges.end(); ++pe) {
            // 分开长边
            logger->info("split?:{}", (*pe)->length());
            if (((*pe)->length() > up_lim)) {
                logger->info("spliting...");
                auto vn = split_edge(*pe).value_or(nullptr);
                auto e1 = vn->halfedge->edge;
                auto e2 = vn->halfedge->prev->edge;
                auto e3 = vn->halfedge->prev->inv->next->edge;
                auto e4 = vn->halfedge->inv->next->edge;
                save_delete_edges.push_back(*pe);
                save_add_edges.push_back(e1);
                save_add_edges.push_back(e2);
                save_add_edges.push_back(e3);
                save_add_edges.push_back(e4);
            }
        }
        for (auto pe : save_add_edges) {
            selected_edges.erase(pe);
        }
        for (auto pe : save_delete_edges) {
            selected_edges.insert(pe);
        }
        for (auto pe = selected_edges.rbegin(); pe != selected_edges.rend(); ++pe) {
            // 摧毁短边
            logger->info("colloapse?:{}", (*pe)->length());
            if (((*pe)->length() < down_lim)) {
                auto hinv1 = (*pe)->halfedge->next;
                auto hinv2 = (*pe)->halfedge->inv->next;
                auto e14   = hinv1->edge;
                auto e42   = (*pe)->halfedge->prev->edge;
                auto e23   = hinv2->edge;
                auto e31   = (*pe)->halfedge->inv->prev->edge;
                save_delete_edges.push_back(*pe);
                save_delete_edges.push_back(e14);
                save_delete_edges.push_back(e42);
                save_delete_edges.push_back(e23);
                save_delete_edges.push_back(e31);
                logger->info("colloapsing...");
                collapse_edge(*pe);
                auto e1 = hinv1->edge;
                auto e2 = hinv2->edge;
                save_add_edges.push_back(e1);
                save_add_edges.push_back(e2);
            }
        }
        for (auto pe : save_add_edges) {
            selected_edges.erase(pe);
        }
        for (auto pe : save_delete_edges) {
            selected_edges.erase(pe);
        }
        save_add_edges.clear();
        save_delete_edges.clear();

        // 翻转边
        for (auto e = edges.head; e != nullptr; e = e->next_node) {
            auto v1 = e->halfedge->from;
            auto v2 = e->halfedge->inv->from;
            auto v3 = e->halfedge->next->from;
            auto v4 = e->halfedge->inv->next->from;
            if (std::abs(v1->degree() - 6.0f) + std::abs(v2->degree() - 6.0f) +
                    std::abs(v3->degree() - 6.0f) + std::abs(v4->degree() - 6.0f) >
                std::abs(v1->degree() - 7.0f) + std::abs(v2->degree() - 7.0f) +
                    std::abs(v3->degree() - 5.0f) + std::abs(v4->degree() - 5.0f)) {
                flip_edge(e);
            }
        }
        for (auto v = vertices.head; v != nullptr; v = v->next_node) {
            // 将节点平均化
        }
    }
    logger->info("remeshed mesh: {} vertices, {} faces\n", vertices.size, faces.size);
    global_inconsistent = true;
    validate();
}
