#ifndef CGOGN_MODULE_SPH_MODULE_H_
#define CGOGN_MODULE_SPH_MODULE_H_

#include <cgogn/core/ui_modules/mesh_provider.h>

#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/rendering/shape_drawer.h>

#include <cgogn/simulation/algos/SPH.h>

#include <cgogn/geometry/algos/picking.h>



#include <iomanip>


namespace cgogn {

namespace ui {

template <typename MESH>

class SimSPH : public ViewModule
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
    using Vertex = typename mesh_traits<MESH>::Vertex;
    using Volume = typename mesh_traits<MESH>::Volume;
    using Edge = typename mesh_traits<MESH>::Edge;
    using Vec3 = geometry::Vec3;
    using Mat3d = geometry::Mat3d;
    using Scalar = geometry::Scalar;

    struct Parameters
    {
        Parameters()
            : vertex_forces_(nullptr), vertex_position_(nullptr), have_selected_vertex_(false), init_vertex_position_(nullptr),
              fixed_vertex(nullptr), show_frame_manipulator_(nullptr), vertex_masse_(nullptr)
        {
            param_move_vertex_ = rendering::ShaderPointSprite::generate_param();
            param_move_vertex_->color_ = rendering::GLColor(1, 1, 0, 0.65);
            param_move_vertex_->set_vbos({&move_vertex_vbo_});

            param_edge_ = rendering::ShaderBoldLine::generate_param();
            param_edge_->color_ = rendering::GLColor(1, 0, 1, 0.65);
            param_edge_->width_ = 2.0f;
            param_edge_->set_vbos({&edges_vbo_});
        }

        CGOGN_NOT_COPYABLE_NOR_MOVABLE(Parameters);

    public:
        void update_move_vertex_vbo()
        {
            if(have_selected_vertex_)
            {
                std::vector<Vec3> vertices_position;
                vertices_position.push_back(move_vertex_);
                vertices_position.push_back(value<Vec3>(*mesh_, vertex_position_.get(), selected_vertex_));

                rendering::update_vbo(vertices_position, &move_vertex_vbo_);
                rendering::update_vbo(vertices_position, &edges_vbo_);
            }
        }

        MESH* mesh_;
        std::shared_ptr<Attribute<Vec3>> vertex_position_;
        std::shared_ptr<Attribute<Vec3>> vertex_forces_;
        std::shared_ptr<Attribute<double>> vertex_masse_;
        Vertex selected_vertex_;
        bool have_selected_vertex_;
        Vec3 move_vertex_;
        rendering::VBO move_vertex_vbo_;
        rendering::VBO edges_vbo_;

        rendering::FrameManipulator frame_manipulator_;


        std::shared_ptr<Attribute<Vec3>> init_vertex_position_;
        float32 vertex_base_size_;
        float32 vertex_scale_factor_;


        std::shared_ptr<Attribute<bool>> fixed_vertex;
        bool show_frame_manipulator_;

        std::unique_ptr<rendering::ShaderPointSprite::Param> param_move_vertex_;
        std::unique_ptr<rendering::ShaderBoldLine::Param> param_edge_;
    };

public:
    SimSPH(const App& app)
        : ViewModule(app, "SPH Simulation"), selected_view_(app.current_view()), running_(false),
        apply_gravity(true), shape_(nullptr), draw_cube(true), size_cube(250, 5, 250), pos_cube(120, -200, 100), Zaxis_cube(0, 1, 0),
        c_t(0.000f), show_particles_(true)
    {
        param_particle_ = rendering::ShaderPointSprite::generate_param();
        param_particle_->color_ = rendering::GLColor(1, 1, 0, 0.65);
        param_particle_->set_vbos({&particle_vbo_});
    }


    ~SimSPH() {}

private:
    void init_mesh(MESH* m)
    {
        Parameters& p = parameters_[m];
        p.mesh_ = m;
        mesh_connections_[m].push_back(
            boost::synapse::connect<typename MeshProvider<MESH>::template attribute_changed_t<Vec3>>(
                m, [this, m](Attribute<Vec3>* attribute) {
                    Parameters& p = parameters_[m];
                    if (p.vertex_position_.get() == attribute)
                    {
                        p.vertex_base_size_ = geometry::mean_edge_length(*m, p.vertex_position_.get()) / 6.0;
                    }
                }));
        mesh_connections_[m].push_back(
            boost::synapse::connect<typename MeshProvider<MESH>::connectivity_changed>(m, [this, m]() {
            Parameters& p = parameters_[m];
            if (p.vertex_position_ && p.init_vertex_position_ && p.vertex_forces_ && p.vertex_masse_)
            {
                    simu_solver.update_topo(*m, {});
            }

            }));
    }

public:
    void find_box(const MESH& mesh, Vec3& pt_min, Vec3& pt_max, const std::shared_ptr<Attribute<Vec3>>& vertex_position_)
    {
        pt_min = Vec3(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        pt_max = Vec3(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());

        foreach_cell(mesh, [&](Vertex v) -> bool {

            Vec3 pos = value<Vec3>(mesh, vertex_position_, v);
            for (int i = 0; i < 3; ++i)
            {
                if (pos[i] < pt_min[i])
                    pt_min[i] = pos[i];
                if (pos[i] > pt_max[i])
                    pt_max[i] = pos[i];
            }
            return true;
        });
    }

    Vec3 scaling_factor(const Vec3& pt_min, const Vec3& pt_max, const Vec3& new_box)
    {
        Vec3 box_size = pt_max - pt_min;
        double dim_max = std::max({box_size[0], box_size[1], box_size[2]});
        double new_box_min_dim = std::min({new_box[0], new_box[1], new_box[2]});
        double scale_factor = new_box_min_dim / dim_max;
        return Vec3(scale_factor, scale_factor, scale_factor);
    }

    void resize_mesh(const MESH& mesh, const Vec3& scaling_factor, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
    {
        foreach_cell(mesh, [&](Vertex v) -> bool {
            Vec3& pos = value<Vec3>(mesh, vertex_position, v);
            for (int i = 0; i < 3; ++i)
            {
                pos[i] *= scaling_factor[i];
            }
            return true;
        });
    }

    void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
    {
          Parameters& p = parameters_[&m];
//        Vec3 pt_min;
//        Vec3 pt_max;

//        find_box(m, pt_min, pt_max, vertex_position);

//        Vec3 new_box_(300, 300, 300);
//        Vec3 scaling_factor_ = scaling_factor(pt_min, pt_max, new_box_);

//        resize_mesh(m, scaling_factor_, vertex_position);

        //simu_solver.init_solver(*selected_mesh_, vertex_position.get());
        p.vertex_position_ = vertex_position;
        //p.vertex_forces_ = simu_solver.particule_vertex_;

        if (p.vertex_position_)
        {
            p.vertex_base_size_ = geometry::mean_edge_length(m, p.vertex_position_.get()) / 6.0;
        }
    }


    void set_init_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& init_vertex_position)
    {
        Parameters& p = parameters_[&m];

        p.init_vertex_position_ = init_vertex_position;
    }

    void set_vertex_force(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_forces)
    {
        Parameters& p = parameters_[&m];

        p.vertex_forces_ = vertex_forces;
    }

    void set_vertex_masse(const MESH& m, const std::shared_ptr<Attribute<double>>& vertex_masse)
    {
        Parameters& p = parameters_[&m];

        p.vertex_masse_ = vertex_masse;
    }


//    void startSimulation()
//    {
//        mesh_provider_->foreach_mesh([this](MESH& m, const std::string& name) {
//            selected_mesh_ = &m;
//        });

//        Parameters& p = parameters_[selected_mesh_];

//        if (!p.vertex_position_)
//        {
//            foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
//                set_vertex_position(*selected_mesh_, attribute);
//                return;
//            });
//        }

//        if (p.vertex_position_)
//        {
//            Parameters& p = parameters_[selected_mesh_];

//            if (!running_)
//            {
//                start();
//            }

//            if (need_update_) {
//                p.update_move_vertex_vbo();
//                mesh_provider_->emit_attribute_changed(*selected_mesh_, p.vertex_position_.get());
//                need_update_ = false;
//            }
//        }
//    }


protected:

    void init() override
    {
        mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
            app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
        mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
        connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
            mesh_provider_, this, &SimSPH<MESH>::init_mesh));
        shape_ = rendering::ShapeDrawer::instance();
        shape_->color(rendering::ShapeDrawer::CUBE) = rendering::GLColor(1, 0.5, 0, 1);
    }

    void mouse_press_event(View *view, int32 button, int32 x, int32 y) override
    {
        Parameters& p = parameters_[selected_mesh_];

        if(selected_mesh_ && view->shift_pressed())
        {
            if(p.vertex_position_)
            {
              rendering::GLVec3d near = view->unproject(x, y, 0.0);
              rendering::GLVec3d far = view->unproject(x, y, 1.0);
              Vec3 A{near.x(), near.y(), near.z()};
              Vec3 B{far.x(), far.y(), far.z()};
              std::vector<Vertex> picked;
              cgogn::geometry::picking(*selected_mesh_, p.vertex_position_.get(), A, B, picked);
              if(!picked.empty())
              {
                  p.selected_vertex_ = picked[0];
                  p.have_selected_vertex_ = true;
                  p.move_vertex_ = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), picked[0]);
                  p.update_move_vertex_vbo();
                  view->request_update();
              }
            }
        }
    }

    void mouse_release_event(View *view, int32, int32, int32) override
    {
        if(selected_mesh_)
        {
            Parameters& p = parameters_[selected_mesh_];
            p.frame_manipulator_.release();
            view->request_update();
        }
    }



#define TIME_STEP 0.015f
    float incline_angle = 2*M_PI / 6; // 30 degrees in radians
    float incline_angle_p2 = -M_PI / 6;

    void start()
    {
        running_ = true;

        launch_thread([this]() {
            while(this->running_)
            {
                c_t+=TIME_STEP;
                std::cout<<"Simulation time = "<<c_t<<std::endl;

                Parameters& p = parameters_[selected_mesh_];


        for (int i = 0; i < 1; i++)
        {
            if(apply_gravity)
            {
                simu_solver.set_particule_forces(
                    [](simulation::Particule_SPH& p) -> Vec3 { return p.masse_ * Vec3(0, -9.81, 0); });
            }

            simu_solver.solve_constraint(*selected_mesh_, p.vertex_position_.get(), nullptr, TIME_STEP);

        }

        if(draw_cube){

            Vec3 position;
            Vec3 axis_x;
            Vec3 axis_y;
            Vec3 axis_z;
            p.frame_manipulator_.get_position(position);
            p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Xt, axis_x);
            p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Yt, axis_y);
            p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Zt, axis_z);


                    parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {

                        Vec3& pos_particule =
                                value<simulation::Particule_SPH*>(*selected_mesh_, simu_solver.particule_vertex_.get(), v)->current_position_;

        //************Plan 0degree************//
                        Vec3 pos_vertex_cube = pos_particule - pos_cube.cast<double>();
                        double dist0 = pos_vertex_cube.dot(axis_y)-20;
                        double d0 = position.dot(axis_y);
        //************END************//

        //******************Collision Plan 0*************************//
                        if(dist0<d0)
                        {
                            value<simulation::Particule_SPH*>(*selected_mesh_, simu_solver.particule_vertex_.get(), v)->current_position_
                                    += (d0  - dist0) * axis_y;
                            Vec3& speed_ = value<simulation::Particule_SPH*>(*selected_mesh_, simu_solver.particule_vertex_.get(), v)->speed_;
                            speed_ -= axis_y.dot(speed_) * axis_y;
                        }
        //*******************END****************************//

                         return true;
                 });

}
        //********//


                   need_update_ = true;
            }

         });
    }


    void stop()
    {
        running_ = false;
    }

    void draw(View* view) override
    {

        const rendering::GLMat4& proj_matrix = view->projection_matrix();
        const rendering::GLMat4& view_matrix = view->modelview_matrix();
        if (selected_mesh_)
        {
            auto& m = selected_mesh_;
            auto& p = parameters_[selected_mesh_];

            MeshData<MESH>& md = mesh_provider_->mesh_data(*m);

            if (p.have_selected_vertex_ && p.param_move_vertex_->attributes_initialized())
            {
                p.param_move_vertex_->point_size_ = p.vertex_base_size_ * p.vertex_scale_factor_;
                p.param_move_vertex_->bind(proj_matrix, view_matrix);
                glDrawArrays(GL_POINTS, 0, 2);
                p.param_move_vertex_->release();
            }

            if (show_particles_)
            {
                param_particle_->point_size_ = 0.1;
                param_particle_->bind(proj_matrix, view_matrix);
                glDrawArrays(GL_POINTS, 0, simu_solver.particules_.size());
                param_particle_->release();
            }

            if (p.have_selected_vertex_ && p.param_edge_->attributes_initialized())
            {
                p.param_edge_->bind(proj_matrix, view_matrix);
                glDrawArrays(GL_LINES, 0, 2);
                p.param_edge_->release();
            }

            if (p.show_frame_manipulator_)
            {
                Scalar size = (md.bb_max_ - md.bb_min_).norm() / 10;
                p.frame_manipulator_.set_size(size);
                p.frame_manipulator_.draw(true, true, proj_matrix, view_matrix);
            }
        }

        if (draw_cube)
        {
            Eigen::Affine3f transfo0 = Eigen::Translation3f(pos_cube)*Eigen::Scaling(size_cube);
            shape_->draw(rendering::ShapeDrawer::CUBE, proj_matrix, view_matrix * transfo0.matrix());
        }
    }

    void left_panel() override
    {

        std::stringstream ss;
        ss << std::setw(6) << std::fixed << std::setprecision(2) << App::fps();
        std::string str_fps = ss.str() + " fps";
        ImGui::Text(str_fps.c_str());

        if (ImGui::BeginCombo("View", selected_view_->name().c_str()))
        {
            for (View* v : linked_views_)
            {
                bool is_selected = v == selected_view_;
                if (ImGui::Selectable(v->name().c_str(), is_selected))
                    selected_view_ = v;
                if (is_selected)
                    ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        if (ImGui::ListBoxHeader("Mesh"))
        {
            mesh_provider_->foreach_mesh([this](MESH& m, const std::string& name) {
                if (ImGui::Selectable(name.c_str(), &m == selected_mesh_))
                {
                    selected_mesh_ = &m;

                    Parameters& p = parameters_[selected_mesh_];
                    p.fixed_vertex = get_attribute<bool, Vertex>(m, "fixed_vertex");
                    if (p.fixed_vertex == nullptr)
                        p.fixed_vertex = add_attribute<bool, Vertex>(m, "fixed_vertex");

                }
            });
            ImGui::ListBoxFooter();
        }

        if (selected_mesh_)
        {
            double X_button_width = ImGui::CalcTextSize("X").x + ImGui::GetStyle().FramePadding.x * 2;

            Parameters& p = parameters_[selected_mesh_];

            need_update_ |= ImGui::Checkbox("Show plane", &p.show_frame_manipulator_);

            if (ImGui::BeginCombo("Position", p.vertex_position_ ? p.vertex_position_->name().c_str() : "-- select --"))
            {
                foreach_attribute<Vec3, Vertex>(*selected_mesh_,
                                                [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
                                                    bool is_selected = attribute == p.vertex_position_;
                                                    if (ImGui::Selectable(attribute->name().c_str(), is_selected))
                                                        set_vertex_position(*selected_mesh_, attribute);
                                                    if (is_selected)
                                                        ImGui::SetItemDefaultFocus();
                                                });
                ImGui::EndCombo();
            }
            if (p.vertex_position_)
            {
                ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
                if (ImGui::Button("X##position"))
                    set_vertex_position(*selected_mesh_, nullptr);
            }

//            if (ImGui::BeginCombo("Forces", p.vertex_forces_ ? p.vertex_forces_->name().c_str() : "-- select --"))
//            {
//                foreach_attribute<Vec3, Vertex>(*selected_mesh_,
//                                                [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
//                                                    bool is_selected = attribute == p.vertex_forces_;
//                                                    if (ImGui::Selectable(attribute->name().c_str(), is_selected))
//                                                        set_vertex_force(*selected_mesh_, attribute);
//                                                    if (is_selected)
//                                                        ImGui::SetItemDefaultFocus();
//                                                });
//                ImGui::EndCombo();
//            }
//            if (p.vertex_forces_)
//            {
//                ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
//                if (ImGui::Button("X##force"))
//                    set_vertex_force(*selected_mesh_, nullptr);
//            }

//            if (ImGui::BeginCombo("vertex masse", p.vertex_masse_ ? p.vertex_masse_->name().c_str() : "-- select --"))
//            {
//                foreach_attribute<double, Vertex>(*selected_mesh_,
//                                                  [&](const std::shared_ptr<Attribute<double>>& attribute) {
//                                                      bool is_selected = attribute == p.vertex_masse_;
//                                                      if (ImGui::Selectable(attribute->name().c_str(), is_selected))
//                                                      {
//                                                          set_vertex_masse(*selected_mesh_, attribute);
//                                                      }
//                                                      if (is_selected)
//                                                          ImGui::SetItemDefaultFocus();
//                                                  });
//                ImGui::EndCombo();
//            }
//            if (p.vertex_masse_)
//            {
//                ImGui::SameLine(ImGui::GetWindowContentRegionMax().x - X_button_width);
//                if (ImGui::Button("X##masse"))
//                    set_vertex_masse(*selected_mesh_, nullptr);
//            }
//            if (p.vertex_masse_)
//            {
//                if (ImGui::Button("init masse"))
//                {
//                    foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
//                        value<double>(*selected_mesh_, p.vertex_masse_.get(), v) = 1.0f;
//                        return true;
//                    });
//                }
//            }


            if (p.vertex_position_)
            {
                if (ImGui::Button("init solver"))
                {
                    simu_solver.init_solver(*selected_mesh_, p.vertex_position_.get());
                    std::vector<Vec3> vertices_position;

                    for (auto& p : simu_solver.particules_)
                    {
                        vertices_position.push_back(p.current_position_);
                    }
                    rendering::update_vbo(vertices_position, &particle_vbo_);
                }
            }
            if (p.vertex_position_)
            {
                ImGui::Separator();

                MeshData<MESH>& md = mesh_provider_->mesh_data(*selected_mesh_);
                Parameters& p = parameters_[selected_mesh_];

                if (!running_)
                {
                    if (ImGui::Button("play"))
                    {
                        start();
                    }

                }
                else
                {
                    if (ImGui::Button("stop"))
                    {
                        stop();
                    }
                }
                if (need_update_)
                {
                    p.update_move_vertex_vbo();
                    mesh_provider_->emit_attribute_changed(*selected_mesh_, p.vertex_position_.get());
                    need_update_ = false;
                }
            }
        }
    }

public:
    MESH* selected_mesh_;

    std::unordered_map<const MESH*, Parameters> parameters_;
    cgogn::simulation::SPH_Particule_constraint_solver<MESH> simu_solver;

    bool show_particles_;
    View* selected_view_;
    bool running_;
    bool apply_gravity;
    bool draw_cube;
    Eigen::Vector3f size_cube;
    Eigen::Vector3f pos_cube;
    Vec3 Zaxis_cube;
    rendering::ShapeDrawer* shape_;
    float32 c_t;
    rendering::VBO particle_vbo_;
    std::unique_ptr<rendering::ShaderPointSprite::Param> param_particle_;
    std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
    std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
    MeshProvider<MESH>* mesh_provider_;
    bool need_update_;

};

} //ui

} //cgogn

#endif // SPH_H
