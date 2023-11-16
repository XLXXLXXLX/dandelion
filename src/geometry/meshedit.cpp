#include "halfedge.h"
#include <cmath>
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
    optimal_pos  = Vector3f::Zero();
    cost         = 0.0f;
    Vector3f x   = e->center();
    Halfedge* h1 = e->halfedge;
    Halfedge* h2 = h1->inv;
    Matrix4f Kf  = Matrix4f::Zero();

    Kf = vertex_quadrics[h1->from] + vertex_quadrics[h2->from];

    Matrix3f A = Kf.topLeftCorner(3, 3);
    Vector3f b = -Kf.topRightCorner(3, 1);
    if (std::abs(A.determinant()) < 1e-9) {
        optimal_pos = A.partialPivLu().solve(b);
    } else {
        optimal_pos = x;
    }
}

bool operator<(const HalfedgeMesh::EdgeRecord& a, const HalfedgeMesh::EdgeRecord& b)
{
    return a.cost < b.cost;
}

optional<Edge*> HalfedgeMesh::flip_edge(Edge* e)
{
    if (e->on_boundary()) {
        return optional<Edge*>{};
    }
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
    // 插入一种特殊情况：如果删去的边是一个三角锥的底边，
    // 则会产生一个三角形片，它由两个重合的三角形面片构成
    // 遇到这种情况，放弃坍缩该边
    auto v1d = v1->degree();
    auto v2d = v2->degree();
    if (v1d == 3 || v2d == 3) {
        return optional<Edge*>{};
    }

    auto fn1 = new_face(f1->is_boundary || f2->is_boundary);
    auto fn2 = new_face(f1->is_boundary || f2->is_boundary);
    erase(f1);
    erase(f2);
    fn1->halfedge = h12;
    fn2->halfedge = h21;
    e->halfedge   = h12;
    if (!(e->halfedge)) {
        logger->error("halfedge is null");
    }
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
        logger->debug("flip boundary");
        // 如果h21在边上，那么h21->face即为虚假面
        hn1->set_neighbors(h21->next, h2n, h1n, vn, e1n, h21->face);
        h2n->set_neighbors(hn1, h21->prev, hn2, v2, e2n, h21->face);
        h21->face->halfedge = hn1;
    } else if (h12->is_boundary()) {
        // 不认为有可能出现两边的面都是虚假面的情况（那样不就成一条线段了啊喂！！
        logger->debug("flip boundary");
        h1n->set_neighbors(hn2, h12->prev, hn1, v1, e1n, h12->face);
        hn2->set_neighbors(h12->next, h1n, h2n, vn, e2n, h12->face);
        h12->face->halfedge = h1n;
    } else {
        // 如果不在边上，则正常处理
        // 先不考虑有重叠面片的情况
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
        if ((!hn4)) {
            logger->error("halfedge is null");
        }
        erase(f2);
        fn14->halfedge = h14;
        fn42->halfedge = h42;
        hn1->set_neighbors(h14, h4n, h1n, vn, e1n, fn14);
        hn4->set_neighbors(h42, h2n, h4n, vn, e4n, fn42);
        h2n->set_neighbors(hn4, h42, hn2, v2, e2n, fn42);
        h4n->set_neighbors(hn1, h14, hn4, v4, e4n, fn14);
        h14->set_neighbors(h4n, hn1, h14->inv, v1, h14->edge, fn14);
        h42->set_neighbors(h2n, hn4, h42->inv, v4, h42->edge, fn42);

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
        if ((!hn3)) {
            logger->error("halfedge is null");
        }
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
    // 变量初始化：
    //  要摧毁的对象：
    //      坍缩边本边
    auto e12 = e;
    //      坍缩三角形面的半边
    auto h12 = e->halfedge;
    auto h21 = h12->inv;
    auto h23 = h12->next;
    auto h31 = h12->prev;
    auto h14 = h21->next;
    auto h42 = h21->prev;
    //      坍缩边的两个端点
    auto v1 = h12->from;
    auto v2 = h21->from;
    //      坍缩三角形面
    auto f123 = h12->face;
    auto f124 = h21->face;
    //      删去右上、左下两个边
    auto ed1 = h23->edge;
    auto ed2 = h14->edge;
    //-------------------------
    //  需要调整的对象
    //      保留菱形外侧四个半边与新节点相连
    auto h1    = h31->inv;
    auto h2    = h42->inv;
    auto hinv1 = h23->inv;
    auto hinv2 = h14->inv;
    //      保留左上、右下两个边与新节点相连
    auto e1 = h1->edge;
    auto e2 = h2->edge;
    //      上下两个节点的连接关系需要调整
    auto v3 = hinv1->from;
    auto v4 = hinv2->from;
    if (v1->id == v2->id) {
        logger->critical("v1->id=v2->id");
        return optional<Vertex*>{};
    }
    if (v3->id == v4->id) {
        logger->critical("v3->id=v4->id");
        return optional<Vertex*>{};
    }
    //-------------------------
    //  其他对象
    auto hfromv1 = v1->halfedge;
    auto hfromv2 = v2->halfedge;

    // 插入一种特殊情况：如果删去的边是一个三角锥的底边，
    // 则会产生一个三角形片，它由两个重合的三角形面片构成
    // 遇到这种情况，放弃坍缩该边
    auto v3d = v3->degree();
    auto v4d = v4->degree();
    if (v3d == 3 || v4d == 3) {
        return optional<Vertex*>{};
    }

    //  新创建的对象
    auto vn    = new_vertex();
    vn->is_new = true;
    //=========================
    // 调整
    auto v1d = v1->degree();
    auto v2d = v2->degree();
    // 将所有原先指向v1的半边指向vn
    for (size_t i = 0; i < v1d; i++) {
        hfromv1->from = vn;
        hfromv1       = hfromv1->inv->next;
    }
    // 将所有原先指向v2的半边指向vn
    for (size_t i = 0; i < v2d; i++) {
        hfromv2->from = vn;
        hfromv2       = hfromv2->inv->next;
    }
    // 如果上三角面为虚假面
    if (h12->is_boundary()) {
        logger->info("collapse boundary");
        h12->face->halfedge = h12->next;
        vn->pos             = e->center();
        vn->halfedge        = h2;
        v4->halfedge        = v4->halfedge->inv->next;
        h2->set_neighbors(h2->next, h2->prev, hinv2, vn, e2, h2->face);
        hinv2->set_neighbors(hinv2->next, hinv2->prev, h2, hinv2->from, e2, hinv2->face);
        erase(h21);
        erase(h14);
        erase(h42);
        erase(e12);
        erase(ed2);
        erase(f124);
    } // 如果下三角面为虚假面
    else if (h21->is_boundary()) {
        logger->info("collapse boundary");
        h21->face->halfedge = h21->next;
        vn->pos             = e->center();
        vn->halfedge        = h1;
        v3->halfedge        = v3->halfedge->inv->next;
        h1->set_neighbors(h1->next, h1->prev, hinv1, vn, e1, h1->face);
        hinv1->set_neighbors(hinv1->next, hinv1->prev, h1, hinv1->from, e1, hinv1->face);
        erase(h12);
        erase(h23);
        erase(h31);
        erase(e12);
        erase(ed1);
        erase(f123);
    } else { // 如果没在边上，正常处理
        // 调整四个剩余半边的连接关系
        h1->inv     = hinv1;
        h1->edge    = e1;
        hinv1->inv  = h1;
        hinv1->edge = e1;
        h2->inv     = hinv2;
        h2->edge    = e2;
        hinv2->inv  = h2;
        hinv2->edge = e2;
        // 调整两个边的连接关系
        e1->halfedge = h1;
        e2->halfedge = h2;
        // 调整三个节点的连接关系
        vn->halfedge = h1;
        v3->halfedge = hinv1;
        v4->halfedge = hinv2;
        // 调整vn的位置
        vn->pos = (v1->pos + v2->pos) / 2;
        //==================
        // 摧毁的对象
        //  六个内部的半边
        erase(h12);
        erase(h21);
        erase(h23);
        erase(h31);
        erase(h14);
        erase(h42);
        //  中心边和右上、左下两个边
        erase(ed1);
        erase(ed2);
        erase(e12);
        //  两个面
        erase(f123);
        erase(f124);
    }
    // 摧毁两个端点
    erase(v1), erase(v2);
    v3->degree(), v4->degree();

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

    // 下面遍历v1相邻的所有面
    // do {
    //     Face* f1 = h1->face;
    //     // 指向下一个面的半边
    //     h1 = h1->inv->next;
    //     // 面的法向量
    //     Vector3f N = f1->normal();
    //     // x到面的投影
    //     Vector3f p = x - (x.transpose() * N) * N;
    //     // Vector4f u(x(0), x(1), x(2), 1.0f);
    //     Vector4f v(N(0), N(1), N(2), -N.transpose() * p);
    //     Matrix4f Kfi = v * v.transpose();
    //     Kf += Kfi;
    // } while (h1->edge != e);
    // // 下面遍历v2相邻的所有面
    // do {
    //     Face* f2 = h2->face;
    //     // 指向下一个面的半边
    //     h2 = h2->inv->next;
    //     // 面的法向量
    //     Vector3f N = f2->normal();
    //     // x到面的投影
    //     Vector3f p = x - (x.transpose() * N) * N;
    //     // Vector4f u(x(0), x(1), x(2), 1.0f);
    //     Vector4f v(N(0), N(1), N(2), -N.transpose() * p);
    //     Matrix4f Kfj = v * v.transpose();
    //     Kf += Kfj;
    // } while (h2->edge != e);

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics
    for (auto iter = faces.head; iter != nullptr; iter = iter->next_node) {

        // Compute the plane equation for the face
        Vector3f plane_normal = iter->normal();
        // Vector3f plane_point           = iter->point();
        // float plane_d                  = -plane_normal.dot(plane_point);
        Matrix4f face_quadric          = Matrix4f::Zero();
        face_quadric.block<1, 3>(3, 0) = plane_normal.transpose();
        // face_quadric(3, 3)             = plane_d;
        face_quadrics[iter] = face_quadric;
    }

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
    double average_edge_length = 0;
    for (auto e = edges.head; e != nullptr; e = e->next_node) {
        selected_edges.insert(e);
        average_edge_length += e->length();
    }
    average_edge_length /= edges.size;
    auto up_lim   = average_edge_length * 4.0f / 3;
    auto down_lim = average_edge_length * 4.0f / 5;
    for (size_t i = 0; i != iteration_limit; ++i) {
        logger->info("###Iteration {}###", i + 1);
        logger->info("...splits begin...");
        // int count = 0;
        for (auto pe = selected_edges.begin(); pe != selected_edges.end(); ++pe) {
            // 分开长边
            if (erased_edges.find((*pe)->id) != erased_edges.end()) {
                selected_edges.erase(pe++);
            } else if (((*pe)->length() > up_lim)) {
                auto vn = split_edge(*pe).value();
                // 将新增的四个边也加入被选中的边中（因为它们很可能过短，需要坍缩
                auto e1 = vn->halfedge->edge;
                auto e2 = vn->halfedge->prev->edge;
                auto e3 = vn->halfedge->prev->inv->next->edge;
                auto e4 = vn->halfedge->inv->next->edge;
                selected_edges.erase(pe++);
                selected_edges.insert(e1);
                selected_edges.insert(e2);
                selected_edges.insert(e3);
                selected_edges.insert(e4);
            }
        }
        logger->info("...splits end...");
        logger->info("...collapses begin...");
        for (auto pe = selected_edges.begin(); pe != selected_edges.end(); ++pe) {
            // 摧毁短边
            if (erased_edges.find((*pe)->id) != erased_edges.end()) {
                logger->trace("e erased, not collapse");
                selected_edges.erase(pe++);
            } else if ((*pe)->length() < down_lim) {
                bool collapsable = true;
                auto h12         = (*pe)->halfedge;
                auto hfromv1     = h12->inv->next;
                auto hfromv2     = h12->inv;
                auto ecenter     = (*pe)->center();
                do {
                    hfromv1 = hfromv1->inv->next;
                    auto v_ = hfromv1->inv->from;
                    if ((v_->pos - ecenter).norm() > up_lim) {
                        collapsable = false;
                    }
                } while (hfromv1 != h12->prev->inv);
                do {
                    hfromv2 = hfromv2->inv->next;
                    auto v_ = hfromv2->inv->from;
                    if ((v_->pos - ecenter).norm() > up_lim) {
                        collapsable = false;
                    }
                } while (hfromv1 != h12->prev->inv);
                if (collapsable) {
                    logger->trace("[collapse begins]");
                    if (!collapse_edge(*pe).has_value()) {
                        logger->trace("[collapse failed]");
                    }
                    logger->trace("[collapse ends]");
                } else {
                    logger->trace("[collapse not needed]");
                }
                selected_edges.erase(pe++);
            }
        }
        logger->info("...collapses end...");
        // 翻转边
        logger->info("...flips begin...");
        for (auto e = edges.head; e != nullptr; e = e->next_node) {
            auto v1  = e->halfedge->from;
            auto v2  = e->halfedge->inv->from;
            auto v3  = e->halfedge->prev->from;
            auto v4  = e->halfedge->inv->prev->from;
            auto v1d = v1->degree();
            auto v2d = v2->degree();
            auto v3d = v3->degree();
            auto v4d = v4->degree();
            logger->trace("v1d:{} v2d:{} v3d:{} v4d:{}", v1d, v2d, v3d, v4d);
            auto d0 = std::abs((int)(v1d - 6)) + std::abs((int)(v2d - 6)) +
                      std::abs((int)(v3d - 6)) + std::abs((int)(v4d - 6));
            auto d1 = std::abs((int)(v1d - 7)) + std::abs((int)(v2d - 7)) +
                      std::abs((int)(v3d - 5)) + std::abs((int)(v4d - 5));
            logger->trace("d0:{}  d1:{}", d0, d1);
            if (d0 > d1) {
                logger->trace("[flip begins]");
                if (!flip_edge(e)) {
                    logger->trace("[not flip]");
                }
                logger->trace("[flip ends]");
            }
        }
        logger->info("...flips end...");
        // 将节点平均化
        logger->info("...smooth begins...");
        for (auto v = vertices.head; v != nullptr; v = v->next_node) {
            using Eigen::Vector3f;
            Vector3f p  = v->pos;
            Vector3f c  = v->neighborhood_center();
            Vector3f N  = v->normal();
            Vector3f V  = c - p;
            Vector3f V_ = V - (N.dot(V)) * N;
            float w     = 1.0f / 5;
            v->new_pos  = p + w * V_;
        }
        for (auto v = vertices.head; v != nullptr; v = v->next_node) {
            v->pos = v->new_pos;
        }
        logger->info("...smooth ends...");
    }
    logger->info("remeshed mesh: {} vertices, {} faces\n", vertices.size, faces.size);
    global_inconsistent = true;
    validate();
}
