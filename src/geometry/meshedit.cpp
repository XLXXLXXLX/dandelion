#include "halfedge.h"
#include "spdlog/logger.h"
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <set>
#include <map>
#include <sys/types.h>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>

#include <spdlog/formatter.hpp>
#include "../utils/logger.h"

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
    auto tmplogger = get_logger("EdgeRecord");
    (void)vertex_quadrics;
    {
        Halfedge* h1 = e->halfedge;
        Halfedge* h2 = e->halfedge->inv;
        Vertex* v1   = h1->from;
        Vertex* v2   = h2->from;
        Vector3f x   = e->center();
        for (Halfedge* h = h1->inv->next; h != h1; h = h->inv->next) {
            if (h->from != v1) {
                break;
            }
            Face* f     = h->face;
            Vector3f N  = f->normal();
            Vector3f p  = (x - f->center()).norm() * N / N.norm();
            Vector4f v  = {N.x(), N.y(), N.z(), -N.dot(p)};
            Matrix4f Kf = v * v.transpose();
            if (vertex_quadrics.find(v1) == vertex_quadrics.end())
                vertex_quadrics[v1] = Matrix4f::Zero();
            vertex_quadrics[v1] += Kf;
        }
        for (Halfedge* h = h2->inv->next; h != h2; h = h->inv->next) {
            if (h->from != v2) {
                break;
            }
            Face* f     = h->face;
            Vector3f N  = f->normal();
            Vector3f p  = (x - f->center()).norm() * N / N.norm();
            Vector4f v  = {N.x(), N.y(), N.z(), -N.dot(p)};
            Matrix4f Kf = v * v.transpose();
            if (vertex_quadrics.find(v2) == vertex_quadrics.end())
                vertex_quadrics[v2] = Matrix4f::Zero();
            vertex_quadrics[v2] += Kf;
        }
    }

    Halfedge* h1 = e->halfedge;
    Halfedge* h2 = h1->inv;
    Vertex* v1   = h1->from;
    Vertex* v2   = h2->from;
    Matrix4f Kf  = vertex_quadrics[v1] + vertex_quadrics[v2];
    Vector4f u   = Vector4f::Zero();

    Vector4f b(0.0f, 0.0f, 0.0f, 1.0f);
    Matrix4f Kf_copy = Kf;
    Kf_copy(3, 0)    = 0.0f;
    Kf_copy(3, 1)    = 0.0f;
    Kf_copy(3, 2)    = 0.0f;
    Kf_copy(3, 3)    = 1.0f;
    Eigen::FullPivLU<Matrix4f> lu(Kf_copy);
    if (lu.isInjective()) {
        u           = (Kf_copy.inverse() * b);
        optimal_pos = {u.x(), u.y(), u.z()};
    } else {
        optimal_pos = e->center();
        u           = {optimal_pos.x(), optimal_pos.y(), optimal_pos.z(), 1.0f};
    }

    if (optimal_pos == Vector3f::Zero()) {
        optimal_pos = e->center();
        u           = {optimal_pos.x(), optimal_pos.y(), optimal_pos.z(), 1.0f};
        cost        = u.transpose() * Kf * u;
    } else {
        cost = u.transpose() * Kf * u;
    }
    tmplogger->trace("[EdgeRecord] u: {}", u);
    tmplogger->trace("[EdgeRecord] optimal_pos: {}", optimal_pos);
    tmplogger->trace("[EdgeRecord] cost: {}", cost);
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
    Halfedge* h43 = h12;
    Halfedge* h34 = h21;
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
    size_t v1d = v1->degree();
    size_t v2d = v2->degree();
    if (v1d == 3 || v2d == 3) {
        return optional<Edge*>{};
    }

    Face* fn1 = new_face(f1->is_boundary || f2->is_boundary);
    Face* fn2 = new_face(f1->is_boundary || f2->is_boundary);
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
    Halfedge* h12 = e->halfedge;
    Halfedge* h21 = h12->inv;
    Vertex* v1    = h12->from;
    Vertex* v2    = h21->from;
    Vertex* vn    = new_vertex();
    vn->pos       = e->center();
    vn->is_new    = true;
    Halfedge* h1n = h12;
    Halfedge* h2n = h21;
    Halfedge* hn1 = new_halfedge();
    Halfedge* hn2 = new_halfedge();
    Edge* e1n     = e;
    Edge* e2n     = new_edge();
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
        Halfedge* h14 = h21->next;
        Halfedge* h42 = h14->next;
        Vertex* v4    = h42->from;
        Face* f2      = h21->face;
        Halfedge* h4n = new_halfedge();
        Halfedge* hn4 = new_halfedge();
        Face* fn14    = new_face(f2->is_boundary);
        Face* fn42    = new_face(f2->is_boundary);
        Edge* e4n     = new_edge();
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

        Halfedge* h23 = h12->next;
        Halfedge* h31 = h23->next;
        Vertex* v3    = h31->from;
        Face* f1      = h12->face;
        Halfedge* h3n = new_halfedge();
        Halfedge* hn3 = new_halfedge();
        Face* fn31    = new_face(f1->is_boundary);
        Face* fn23    = new_face(f1->is_boundary);
        Edge* e3n     = new_edge();
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
    logger->trace("[collapse] edge to collapse:{}", e->id);
    // 变量初始化：
    //  要摧毁的对象：
    //      坍缩边本边
    Edge* e12 = e;
    //      坍缩三角形面的半边
    Halfedge* h12 = e->halfedge;
    Halfedge* h21 = h12->inv;
    Halfedge* h23 = h12->next;
    Halfedge* h31 = h12->prev;
    Halfedge* h14 = h21->next;
    Halfedge* h42 = h21->prev;
    //      坍缩边的两个端点
    Vertex* v1 = h12->from;
    Vertex* v2 = h21->from;
    //      坍缩三角形面
    Face* f123 = h12->face;
    Face* f124 = h21->face;
    //      删去右上、左下两个边
    Edge* ed1 = h23->edge;
    Edge* ed2 = h14->edge;
    //-------------------------
    //  需要调整的对象
    //      保留菱形外侧四个半边与新节点相连
    Halfedge* h1    = h31->inv;
    Halfedge* h2    = h42->inv;
    Halfedge* hinv1 = h23->inv;
    Halfedge* hinv2 = h14->inv;
    //      保留左上、右下两个边与新节点相连
    Edge* e1 = h1->edge;
    Edge* e2 = h2->edge;
    //      上下两个节点的连接关系需要调整
    Vertex* v3 = hinv1->from;
    Vertex* v4 = hinv2->from;
    logger->trace("[collapse] edge to collapse:{}", e->id);
    logger->trace("[collapse] v1->id={}", v1->id);
    logger->trace("[collapse] v2->id={}", v2->id);
    logger->trace("[collapse] v3->id={}", v3->id);
    logger->trace("[collapse] v4->id={}", v4->id);
    if (v1->id == v2->id) {
        logger->critical("[collapse] v1->id=v2->id");
        return optional<Vertex*>{};
    }
    if (v3->id == v4->id) {
        logger->critical("[collapse] v3->id=v4->id");
        return optional<Vertex*>{};
    }
    //-------------------------

    // 插入一种特殊情况：如果删去的边是一个三角锥的底边，
    // 则会产生一个三角形片，它由两个重合的三角形面片构成
    // 遇到这种情况，放弃坍缩该边
    size_t v3d = v3->degree();
    size_t v4d = v4->degree();
    if (v3d == 3 || v4d == 3) {
        return optional<Vertex*>{};
    }

    //  新创建的对象
    Vertex* vn = new_vertex();
    vn->is_new = true;
    //=========================
    // 调整
    // 将所有原先指向v1的半边指向vn
    h12->from = vn;
    for (Halfedge* hfromv1 = h12->inv->next; hfromv1 != h12; hfromv1 = hfromv1->inv->next) {
        if (hfromv1->from != v1) {
            logger->critical("[collapse] hfromv1->from!=v1");
            break;
        }
        hfromv1->from = vn;
    }
    // 将所有原先指向v2的半边指向vn
    h21->from = vn;
    for (Halfedge* hfromv2 = h21->inv->next; hfromv2 != h21; hfromv2 = hfromv2->inv->next) {
        if (hfromv2->from != v2) {
            logger->critical("[collapse] hfromv2->from!=v2");
            break;
        }
        hfromv2->from = vn;
    }
    // 如果上三角面为虚假面
    if (h12->is_boundary()) {
        logger->info("[collapse] collapse boundary");
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
        logger->info("[collapse] collapse boundary");
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
        if (vn->id != h1->from->id) {
            logger->critical("[collapse] fucked up");
        }
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
    for (Vertex* v = vertices.head; v != nullptr; v = v->next_node) {
        size_t n = v->degree();
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
    for (Edge* e = edges.head; e != nullptr; e = e->next_node) {
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
    for (Edge* e : original_edges) {
        Vector3f new_pos = e->new_pos;
        Vertex* v        = split_edge(e).value_or(nullptr);
        if (v) {
            v->new_pos = new_pos;
        }
    }
    // Now flip any new edge that connects an old and new vertex.
    for (Edge* e = edges.head; e != nullptr; e = e->next_node) {
        if (e->is_new && (e->halfedge->from->is_new != e->halfedge->inv->from->is_new)) {
            flip_edge(e);
        }
    }

    // Finally, copy new vertex positions into the Vertex::pos.
    for (Vertex* v = vertices.head; v != nullptr; v = v->next_node) v->pos = v->new_pos;
    for (Edge* e = edges.head; e != nullptr; e = e->next_node) e->is_new = false;
    for (Vertex* v = vertices.head; v != nullptr; v = v->next_node) v->is_new = false;

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
    logger->info("[simplify] initialization start");
    for (Edge* e = edges.head; e != nullptr; e = e->next_node) {
        edge_records[e] = EdgeRecord(vertex_quadrics, e);
    }
    for (auto pair : edge_records) {
        edge_queue.insert(pair.second);
    }
    logger->info("[simplify] initialization done");
    logger->info("[simplify] collapses begin");
    unsigned int faces_to_erase      = faces.size / 4;
    unsigned int erased_faces_amount = 0;
    while (true) {
        // -> Remove edges from the queue that have already been processed
        if (erased_faces_amount >= faces_to_erase) {
            logger->trace("[simplify] erased_faces_amount >= faces_to_erase");
            break;
        }
        if (edge_queue.size() <= 0) {
            logger->trace("[simplify] edge_queue.size() <= 0");
            break;
        }
        EdgeRecord er = *edge_queue.begin();
        if (edge_queue.find(er) == edge_queue.end()) {
            logger->critical("[simplify] edge_queue.find(er) == edge_queue.end()");
            continue;
        }
        Edge* e = er.edge;
        edge_queue.erase(er);
        Vector3f new_pos = edge_records[e].optimal_pos;
        Halfedge* h1 = e->halfedge;
        for (Halfedge* h = h1->inv->next; h != h1; h = h->inv->next) {
            if (h->from != h1->from) {
                logger->critical("[simplify] hfromv1->from!=v1");
                break;
            }
            Edge* e = h->edge;
            logger->trace("[simplify] Rerecord edges beside the collapse edge, id={}", e->id);
            EdgeRecord er1 = edge_records[e];
            edge_queue.erase(er1);
        }
        Halfedge* h2 = e->halfedge->inv;
        for (Halfedge* h = h2->inv->next; h != h2; h = h->inv->next) {
            if (h->from != h2->from) {
                logger->critical("[simplify] hfromv2->from!=v2");
                break;
            }
            Edge* e = h->edge;
            logger->trace("[simplify] Rerecord edges beside the collapse edge, id={}", e->id);
            EdgeRecord er1 = edge_records[e];
            edge_queue.erase(er1);
        }

        logger->trace("[simplify] collapse begins");
        optional<Vertex*> res = collapse_edge(e);
        if (res.has_value()) {
            Vertex* vn = res.value();
            vn->pos      = new_pos;
            Halfedge* h1 = vn->halfedge;
            EdgeRecord er = edge_records[e] = {vertex_quadrics, e};
            edge_queue.insert(er);
            for (Halfedge* h = h1->inv->next; h != h1; h = h->inv->next) {
                EdgeRecord er = edge_records[e] = {vertex_quadrics, e};
                edge_queue.insert(er);
            }
            erased_faces_amount += 2;
        } else {
            logger->info("[simplify] collapse failed");
            for (Halfedge* h = h1->inv->next; h != h1; h = h->inv->next) {
                Edge* e        = h->edge;
                EdgeRecord er1 = edge_records[e];
                edge_queue.insert(er1);
            }
            for (Halfedge* h = h2->inv->next; h != h2; h = h->inv->next) {
                Edge* e        = h->edge;
                EdgeRecord er1 = edge_records[e];
                edge_queue.insert(er1);
            }
        }
        logger->trace("[simplify] collapse ends");
    }
    logger->info("[simplify] collapses end");
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
    double average_edge_length = 0;
    for (Edge* e = edges.head; e != nullptr; e = e->next_node) {
        selected_edges.insert(e);
        average_edge_length += e->length();
    }
    average_edge_length /= edges.size;
    double up_lim   = average_edge_length * 4.0f / 3;
    double down_lim = average_edge_length * 4.0f / 5;
    for (size_t i = 0; i != iteration_limit; ++i) {
        logger->info("[isotropic_remesh]###Iteration {}###", i + 1);
        logger->info("[isotropic_remesh]...splits begin...");
        // int count = 0;
        for (auto pe = selected_edges.begin(); pe != selected_edges.end(); ++pe) {
            // 分开长边
            if (erased_edges.find((*pe)->id) != erased_edges.end()) {
                selected_edges.erase(pe++);
            } else if (((*pe)->length() > up_lim)) {
                Vertex* vn = split_edge(*pe).value();
                // 将新增的四个边也加入被选中的边中（因为它们很可能过短，需要坍缩
                Edge* e1 = vn->halfedge->edge;
                Edge* e2 = vn->halfedge->prev->edge;
                Edge* e3 = vn->halfedge->prev->inv->next->edge;
                Edge* e4 = vn->halfedge->inv->next->edge;
                selected_edges.erase(pe++);
                selected_edges.insert(e1);
                selected_edges.insert(e2);
                selected_edges.insert(e3);
                selected_edges.insert(e4);
            }
        }
        logger->info("[isotropic_remesh]...splits end...");
        logger->info("[isotropic_remesh]...collapses begin...");
        for (auto pe = selected_edges.begin(); pe != selected_edges.end(); ++pe) {
            // 摧毁短边
            if (erased_edges.find((*pe)->id) != erased_edges.end()) {
                logger->trace("[isotropic_remesh]e erased, not collapse");
                selected_edges.erase(pe++);
            } else if ((*pe)->length() < down_lim) {
                bool collapsable  = true;
                Halfedge* h12     = (*pe)->halfedge;
                Halfedge* hfromv1 = h12;
                Halfedge* hfromv2 = h12->inv;
                Vector3f ecenter  = (*pe)->center();
                do {
                    hfromv1    = hfromv1->inv->next;
                    Vertex* v_ = hfromv1->inv->from;
                    if ((v_->pos - ecenter).norm() > up_lim) {
                        collapsable = false;
                    }
                } while (hfromv1 != h12->prev->inv);
                do {
                    hfromv2    = hfromv2->inv->next;
                    Vertex* v_ = hfromv2->inv->from;
                    if ((v_->pos - ecenter).norm() > up_lim) {
                        collapsable = false;
                    }
                } while (hfromv2 != h12->inv->prev->inv);
                if (collapsable) {
                    logger->trace("[isotropic_remesh][collapse begins]");
                    if (!collapse_edge(*pe).has_value()) {
                        logger->trace("[isotropic_remesh][collapse failed]");
                    }
                    logger->trace("[isotropic_remesh][collapse ends]");
                } else {
                    logger->trace("[isotropic_remesh][collapse not needed]");
                }
                selected_edges.erase(pe++);
            }
        }
        logger->info("[isotropic_remesh]...collapses end...");
        // 翻转边
        logger->info("[isotropic_remesh]...flips begin...");
        for (Edge* e = edges.head; e != nullptr; e = e->next_node) {
            Vertex* v1 = e->halfedge->from;
            Vertex* v2 = e->halfedge->inv->from;
            Vertex* v3 = e->halfedge->prev->from;
            Vertex* v4 = e->halfedge->inv->prev->from;
            size_t v1d = v1->degree();
            size_t v2d = v2->degree();
            size_t v3d = v3->degree();
            size_t v4d = v4->degree();
            logger->trace("[isotropic_remesh]v1d:{} v2d:{} v3d:{} v4d:{}", v1d, v2d, v3d, v4d);
            int d0 = std::abs((int)(v1d - 6)) + std::abs((int)(v2d - 6)) +
                     std::abs((int)(v3d - 6)) + std::abs((int)(v4d - 6));
            int d1 = std::abs((int)(v1d - 7)) + std::abs((int)(v2d - 7)) +
                     std::abs((int)(v3d - 5)) + std::abs((int)(v4d - 5));
            logger->trace("[isotropic_remesh]d0:{}  d1:{}", d0, d1);
            if (d0 > d1) {
                logger->trace("[isotropic_remesh][flip begins]");
                if (!flip_edge(e)) {
                    logger->trace("[isotropic_remesh][not flip]");
                }
                logger->trace("[isotropic_remesh][flip ends]");
            }
        }
        logger->info("[isotropic_remesh]...flips end...");
        // 将节点平均化
        logger->info("[isotropic_remesh]...smooth begins...");
        for (Vertex* v = vertices.head; v != nullptr; v = v->next_node) {
            using Eigen::Vector3f;
            Vector3f p  = v->pos;
            Vector3f c  = v->neighborhood_center();
            Vector3f N  = v->normal();
            Vector3f V  = c - p;
            Vector3f V_ = V - (N.dot(V)) * N;
            float w     = 1.0f / 5;
            v->new_pos  = p + w * V_;
        }
        for (Vertex* v = vertices.head; v != nullptr; v = v->next_node) {
            v->pos = v->new_pos;
        }
        logger->info("[isotropic_remesh]...smooth ends...");
    }
    logger->info("remeshed mesh: {} vertices, {} faces\n", vertices.size, faces.size);
    global_inconsistent = true;
    validate();
}
