
Matrix4f Object::model()
{
    const Quaternionf& r             = rotation;
    auto [x_angle, y_angle, z_angle] = quaternion_to_ZYX_euler<double>(r.w(), r.x(), r.y(), r.z());

    Matrix4f rotation_x_matrix, rotation_y_matrix, rotation_z_matrix;

    rotation_x_matrix << 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, cos(-radians<double>(x_angle)),
        sin(-radians<double>(x_angle)), 0.0f, 0.0f, -sin(-radians<double>(x_angle)),
        cos(-radians<double>(x_angle)), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f;
    rotation_y_matrix << cos(-radians<double>(y_angle)), 0.0f, -sin(-radians<double>(y_angle)),
        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, sin(-radians<double>(y_angle)), 0.0f,
        cos(-radians<double>(y_angle)), 0.0f, 0.0f, 0.0f, 0.0f, 1.0f;
    rotation_z_matrix << cos(-radians<double>(z_angle)), sin(-radians<double>(z_angle)), 0.0f, 0.0f,
        -sin(-radians<double>(z_angle)), cos(-radians<double>(z_angle)), 0.0f, 0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f;
    Matrix4f rotation_matrix = rotation_x_matrix * rotation_y_matrix * rotation_z_matrix;

    Matrix4f scaling_matrix = Matrix4f::Identity();
    scaling_matrix(0, 0)    = scaling[0];
    scaling_matrix(1, 1)    = scaling[1];
    scaling_matrix(2, 2)    = scaling[2];
    scaling_matrix(3, 3)    = 1.0f;

    Matrix4f center_matrix;
    center_matrix << 1.0f, 0.0f, 0.0f, center[0], 0.0f, 1.0f, 0.0f, center[1], 0.0f, 0.0f, 1.0f,
        center[2], 0.0f, 0.0f, 0.0f, 1.0f;

    return center_matrix * rotation_matrix * scaling_matrix;
}

Matrix4f Camera::projection()
{
    const float fov_y = radians(fov_y_degrees); // Y方向视角大小（弧度制）
    // 使用平行投影时，用户并不能从画面上直观地感受到相机的位置，
    // 因而会产生近处物体裁剪过多的错觉。为了产程更好的观察效果，
    // 这里没有使用相机本身的 near 而是取 near = -far 来让相机能看到“背后”的物体。
    Matrix4f projection = Matrix4f::Zero();
    projection(1, 1)    = 1 / std::tan(fov_y / 2.0f);
    projection(0, 0)    = projection(1, 1) / aspect_ratio;
    projection(2, 2)    = -(far + near) / (far - near);
    projection(3, 2)    = -1.0f;
    projection(2, 3)    = -(2.0f * far * near) / (far - near);

    return projection;
}
// changes {
Eigen::Quaternionf MyAngleAxisf(float radian_angle, Vector3f axis)
{
    double axisLength = axis.norm();
    axis /= axisLength;
    Eigen::Quaternionf q(
        1.0f * std::cos(radian_angle / 2.0f), axis(0) * std::sin(radian_angle / 2.0f),
        axis(1) * std::sin(radian_angle / 2.0f), axis(2) * std::sin(radian_angle / 2.0f));
    return q;
}
//}

void Toolbar::layout_mode(Scene& scene)
{
    if (ImGui::BeginTabItem("Layout")) {
        if (mode != WorkingMode::LAYOUT) {
            on_selection_canceled();
            mode = WorkingMode::LAYOUT;
        }
        scene_hierarchies(scene);

        Object* selected_object = scene.selected_object;
        if (selected_object != nullptr) {
            material_editor(selected_object->mesh.material);
            ImGui::SeparatorText("Transform");
            ImGui::Text("Translation");
            ImGui::PushID("Translation##");
            Vector3f& center = selected_object->center;
            xyz_drag(&center.x(), &center.y(), &center.z(), POSITION_UNIT);
            ImGui::PopID();

            ImGui::Text("Scaling");
            ImGui::PushID("Scaling##");
            Vector3f& scaling = selected_object->scaling;
            xyz_drag(&scaling.x(), &scaling.y(), &scaling.z(), SCALING_UNIT);
            ImGui::PopID();

            const Quaternionf& r             = selected_object->rotation;
            auto [x_angle, y_angle, z_angle] = quaternion_to_ZYX_euler(r.w(), r.x(), r.y(), r.z());
            ImGui::Text("Rotation (ZYX Euler)");
            ImGui::PushID("Rotation##");
            ImGui::PushItemWidth(0.3f * ImGui::CalcItemWidth());
            ImGui::DragFloat("pitch", &x_angle, ANGLE_UNIT, -180.0f, 180.0f, "%.1f deg",
                             ImGuiSliderFlags_AlwaysClamp);
            ImGui::SameLine();
            ImGui::DragFloat("yaw", &y_angle, ANGLE_UNIT, -90.0f, 90.0f, "%.1f deg",
                             ImGuiSliderFlags_AlwaysClamp);
            ImGui::SameLine();
            ImGui::DragFloat("roll", &z_angle, ANGLE_UNIT, -180.0f, 180.0f, "%.1f deg",
                             ImGuiSliderFlags_AlwaysClamp);
            ImGui::PopItemWidth();
            ImGui::PopID();
            // changes{
            // selected_object->rotation = AngleAxisf(radians(x_angle), Vector3f::UnitX()) *
            //                             AngleAxisf(radians(y_angle), Vector3f::UnitY()) *
            //                             AngleAxisf(radians(z_angle), Vector3f::UnitZ());
            selected_object->rotation = MyAngleAxisf(radians(x_angle), Vector3f::UnitX()) *
                                        MyAngleAxisf(radians(y_angle), Vector3f::UnitY()) *
                                        MyAngleAxisf(radians(z_angle), Vector3f::UnitZ());
            //}
        }
        ImGui::EndTabItem();
    }
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

    auto fn1 = new_face(f1->is_boundary || f2->is_boundary);
    auto fn2 = new_face(f1->is_boundary || f2->is_boundary);
    erase(f1);
    erase(f2);
    fn1->halfedge = h12;
    fn2->halfedge = h21;
    e->halfedge   = h12;
    // 重新连接各基本元素
    v1->halfedge = h14;
    v2->halfedge = h23;
    h43->set_neighbors(h31, h14, h21, v4, e, fn1);
    h34->set_neighbors(h42, h23, h12, v3, e, fn2);
    h31->set_neighbors(h14, h43, h31->inv, v3, h31->edge, fn1);
    h14->set_neighbors(h43, h31, h14->inv, v1, h14->edge, fn1);
    h42->set_neighbors(h23, h34, h42->inv, v4, h42->edge, fn2);
    h23->set_neighbors(h34, h42, h23->inv, v2, h23->edge, fn2);

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

    auto vn    = new_vertex();
    vn->pos    = e->center();
    vn->is_new = true;

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
    erase(f1);
    erase(f2);
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
