
Matrix4f Object::model()
{
    const Quaternionf& r             = rotation;
    auto [x_angle, y_angle, z_angle] = quaternion_to_ZYX_euler<double>(r.w(), r.x(), r.y(), r.z());

    Matrix4f rotation_x_matrix, rotation_y_matrix, rotation_z_matrix;

    rotation_x_matrix <<  
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, cos(-radians<double>(x_angle)), sin(-radians<double>(x_angle)), 0.0f,
        0.0f, -sin(-radians<double>(x_angle)), cos(-radians<double>(x_angle)), 0.0f, 
        0.0f, 0.0f, 0.0f, 1.0f;
    rotation_y_matrix << 
        cos(-radians<double>(y_angle)), 0.0f, -sin(-radians<double>(y_angle)), 0.0f, 
        0.0f, 1.0f, 0.0f, 0.0f, 
        sin(-radians<double>(y_angle)), 0.0f, cos(-radians<double>(y_angle)), 0.0f, 
                         0.0f, 0.0f, 0.0f, 1.0f;
    rotation_z_matrix << cos(-radians<double>(z_angle)), sin(-radians<double>(z_angle)), 0.0f, 0.0f, 
                         -sin(-radians<double>(z_angle)), cos(-radians<double>(z_angle)), 0.0f, 0.0f,
                         0.0f, 0.0f, 1.0f, 0.0f,
                         0.0f, 0.0f, 0.0f, 1.0f;
    Matrix4f rotation_matrix = rotation_x_matrix * rotation_y_matrix * rotation_z_matrix;

    Matrix4f scaling_matrix = Matrix4f::Identity();
    scaling_matrix(0, 0)    = scaling[0];
    scaling_matrix(1, 1)    = scaling[1];
    scaling_matrix(2, 2)    = scaling[2];
    scaling_matrix(3, 3)    = 1.0f;

    Matrix4f center_matrix;
    center_matrix << 1.0f, 0.0f, 0.0f, center[0], 
                     0.0f, 1.0f, 0.0f, center[1],
                     0.0f, 0.0f, 1.0f, center[2], 
                     0.0f, 0.0f, 0.0f, 1.0f;

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

//changes {
Vector4f MyAngleAxisf(float radian_angle, Vector3f axis)
{
    Vector4f q;

    // 将旋转角度转换为弧度

    // 计算旋转轴的长度
    double axisLength = axis.norm();

    // 将旋转轴标准化为单位向量
    axis /= axisLength;

    // 计算四元数的实部和虚部
    q(0)                = std::cos(radian_angle);
    double sinHalfAngle = std::sin(radian_angle);
    q(1)                = axis(0) * sinHalfAngle;
    q(2)                = axis(1) * sinHalfAngle;
    q(3)                = axis(2) * sinHalfAngle;
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
            selected_object->rotation = MyAngleAxisf(radians(x_angle), Vector3f::UnitX()) *
                                        MyAngleAxisf(radians(y_angle), Vector3f::UnitY()) *
                                        MyAngleAxisf(radians(z_angle), Vector3f::UnitZ());
            //}
        }
        ImGui::EndTabItem();
    }
}