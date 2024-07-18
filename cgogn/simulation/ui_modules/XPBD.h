/*******************************************************************************
 * CGoGN                                                                        *
 * Copyright (C) 2019, IGG Group, ICube, University of Strasbourg, France       *
 *                                                                              *
 * This library is free software; you can redistribute it and/or modify it      *
 * under the terms of the GNU Lesser General Public License as published by the *
 * Free Software Foundation; either version 2.1 of the License, or (at your     *
 * option) any later version.                                                   *
 *                                                                              *
 * This library is distributed in the hope that it will be useful, but WITHOUT  *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
 * for more details.                                                            *
 *                                                                              *
 * You should have received a copy of the GNU Lesser General Public License     *
 * along with this library; if not, write to the Free Software Foundation,      *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
 *                                                                              *
 * Web site: http://cgogn.unistra.fr/                                           *
 * Contact information: cgogn@unistra.fr                                        *
 *                                                                              *
 *******************************************************************************/

#ifndef CGOGN_MODULE_XPBD_MODULE_H_
#define CGOGN_MODULE_XPBD_MODULE_H_

#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/app.h>
#include <cgogn/ui/view.h>

#include <cgogn/rendering/shape_drawer.h>
#include <cgogn/rendering/vbo_update.h>
#include <cgogn/rendering/frame_manipulator.h>
#include <cgogn/rendering/ui_modules/volume_render.h>

#include <cgogn/core/types/mesh_traits.h>

#include <cgogn/geometry/algos/selection.h> // parallel_foreach_cell
#include <cgogn/geometry/algos/picking.h>
#include <cgogn/geometry/algos/length.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/simulation/algos/XPBD/XPBD.h>

#include <boost/synapse/connect.hpp>
#include <iomanip> //setprecision

#include <cgogn/modeling/algos/subdivision.h>
#include <stdlib.h>

#include <fstream>

namespace cgogn
{

namespace ui
{

template <typename MESH>

class SimBall : public ViewModule
{
    template <typename T>
    using Attribute = typename mesh_traits<MESH>::template Attribute<T>;
    using Vertex = typename mesh_traits<MESH>::Vertex;
    using Edge = typename mesh_traits<MESH>::Edge;
    using Vec3 = geometry::Vec3;
    using Mat3d = geometry::Mat3d;

    using Volume = typename mesh_traits<MESH>::Volume;
    std::vector<Vertex> vertices_inferior;
    std::vector<Vertex> vertices_superior;

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
    SimBall(const App& app) : ViewModule(app, "Cube Simulation"), selected_view_(app.current_view()), running_(false),
    apply_gravity(false), apply_force(true), shape_(nullptr), draw_cube(true),
    size_cube_ref(350, 1, 350), pos_cube_ref(0,0,0), Zaxis_cube(0, 1, 0), c_t(0.000f)
    {

    }

    ~SimBall()
    {

    }

private:
 //********************************************************************************************************************************//
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
                    // sm_solver_.update_topo(*m, {});
                }
            }));
    }
 //********************************************************************************************************************************//

public:

//************************resize mesh***************************************************//
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
//************************END resize mesh***************************************************//

bool t=false;
    void set_vertex_position(const MESH& m, const std::shared_ptr<Attribute<Vec3>>& vertex_position)
    {
        Parameters& p = parameters_[&m];
        Vec3 pt_min;
        Vec3 pt_max;

        find_box(m, pt_min, pt_max, vertex_position);

        Vec3 new_box_(300, 300, 300);
        Vec3 scaling_factor_ = scaling_factor(pt_min, pt_max, new_box_);

        resize_mesh(m, scaling_factor_, vertex_position);
        t=true;


        simu_solver.init_solver(*selected_mesh_, vertex_position);
        p.vertex_position_ = vertex_position;
        p.vertex_forces_ = simu_solver.f_ext_;

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


    void call_set_vertex_position()
    {
        mesh_provider_->foreach_mesh([this](MESH& m, const std::string& name) {
            selected_mesh_ = &m;
        });

        //need_update_ = true;

        Parameters& p = parameters_[selected_mesh_];

        if (!p.vertex_position_)
        {
            foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
                set_vertex_position(*selected_mesh_, attribute);
                return;
            });
        }

        if (p.vertex_position_)
        {
            Parameters& p = parameters_[selected_mesh_];



            if (need_update_) {
                p.update_move_vertex_vbo();
                mesh_provider_->emit_attribute_changed(*selected_mesh_, p.vertex_position_.get());
                need_update_ = false;
            }
        }
    }
    void startSimulation() {
        mesh_provider_->foreach_mesh([this](MESH& m, const std::string& name) {
            selected_mesh_ = &m;
        });

        //need_update_ = true;

        Parameters& p = parameters_[selected_mesh_];

        if (!p.vertex_position_)
        {
            foreach_attribute<Vec3, Vertex>(*selected_mesh_, [&](const std::shared_ptr<Attribute<Vec3>>& attribute) {
                set_vertex_position(*selected_mesh_, attribute);
                return;
            });
        }

        if (p.vertex_position_)
        {
            Parameters& p = parameters_[selected_mesh_];

            if (!running_)
            {
                start();
            }

            if (need_update_) {
                p.update_move_vertex_vbo();
                mesh_provider_->emit_attribute_changed(*selected_mesh_, p.vertex_position_.get());
                need_update_ = false;
            }
        }
    }

protected:
    //
    //********************************************************************************************************************************//
    void init() override
    {
        mesh_provider_ = static_cast<ui::MeshProvider<MESH>*>(
            app_.module("MeshProvider (" + std::string{mesh_traits<MESH>::name} + ")"));
        mesh_provider_->foreach_mesh([this](MESH& m, const std::string&) { init_mesh(&m); });
        connections_.push_back(boost::synapse::connect<typename MeshProvider<MESH>::mesh_added>(
            mesh_provider_, this, &SimBall<MESH>::init_mesh));
        shape_ = rendering::ShapeDrawer::instance();
        shape_->color(rendering::ShapeDrawer::CUBE) = rendering::GLColor(0.24f, 0.0f, 0.5f, 1.0f);
    }
    //--------------------------------------------------------------------------------------------------------------------------------//
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

float32 dt = 0.000f;

bool k = true;
float incline_angle = M_PI / 6; // 30 degrees in radians
float incline_angle_p2 = -1*M_PI / 6;
Vec3 speed_;



#define TIME_STEP 0.00002f

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
                        parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
                            value<Vec3>(*selected_mesh_, p.vertex_forces_, v)+=
                                    value<double>(*selected_mesh_, simu_solver.masse_, v) *Vec3(0, -9.81, 0);

                            return true;
                        });
                    }


//                    parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
//                        Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);
//                        std::cout<<"pos = "<<pos<<std::endl;
//                        return true;
//                    });


                        if(apply_force)
                        {

                            if(p.fixed_vertex == nullptr){
                                p.fixed_vertex = get_or_add_attribute<bool,Vertex>(*selected_mesh_,"fixed_vertex_XPBD");
                                simu_solver.fixed_vertex = p.fixed_vertex;
                            }


                        for (auto v : vertices_inferior) {
                            Parameters& p = parameters_[selected_mesh_];
                            Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);
                            value<bool>(*selected_mesh_, p.fixed_vertex.get(), v) = true;
                            pos.y() -= 30;
//                            value<Vec3>(*selected_mesh_, p.vertex_forces_, v)+=
//                                    value<double>(*selected_mesh_, simu_solver.masse_, v) *Vec3(0, -600, 0);

                        }

                        for (auto v : vertices_superior) {
                            Parameters& p = parameters_[selected_mesh_];
                            Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);
                            value<bool>(*selected_mesh_, p.fixed_vertex.get(), v) = true;
                            pos.y() += 30;
//                            value<Vec3>(*selected_mesh_, p.vertex_forces_, v)+=
//                                    value<double>(*selected_mesh_, simu_solver.masse_, v) *Vec3(0, 600, 0);

                        }



                        if(c_t>=0.00008)
                        {
                            apply_force = false;
                        }

                        }
                    simu_solver.solver(*selected_mesh_, TIME_STEP);

                    parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
                        value<Vec3>(*selected_mesh_, p.vertex_forces_, v) = Vec3(0, 0, 0);
                        return true;
                    });

                }

                    if(draw_cube)
                    {
                        Vec3 position;
                        Vec3 axis_x = Vec3(pos_cube_ref.x(),0,0);
                        Vec3 axis_y = Vec3(0, pos_cube_ref.y(),0);
                        Vec3 axis_z = Vec3(0,0,pos_cube_ref.z());
                        //p.frame_manipulator_.get_position(position);

                        position = Vec3(pos_cube_ref.x(), pos_cube_ref.y(), pos_cube_ref.z());


                        //p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Xt, axis_x);
                        //p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Yt, axis_y);
                        //p.frame_manipulator_.get_axis(cgogn::rendering::FrameManipulator::Zt, axis_z);

                        int total_vertex_count = 0;
                        int contact_with_planes_count = 0;

                        parallel_foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
                            Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);
 //************Plan 0degree************//
                            Vec3 pos_vertex_cube = pos - pos_cube.cast<double>();
                            double dist0 = pos_vertex_cube.dot(axis_y)-20;
                            double d0 = position.dot(axis_y);
 //************END************//

 //************Plan 30degree (inclinaison yz, x fixe)************//
                            Vec3 normal(0, std::cos(incline_angle), std::sin(incline_angle));
                            double dist = pos.dot(normal) + 280 / normal.norm(); // 280 = 300 - 20
                            double d = position.dot(normal) / normal.norm();
//************END************//

//************Plan 150degree (inclinaison yz, x fixe)************//
                             Vec3 normal_p2(0, std::cos(incline_angle_p2), std::sin(incline_angle_p2));
                             double dist_p2 = pos.dot(normal_p2) + 280 / normal_p2.norm(); // 280 = 300 - 20
                             double d_p2 = position.dot(normal_p2) / normal_p2.norm();
//************END************//

                            speed_ = value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v);
                            //std::cout<<"Speed = "<<speed_<<std::endl;


                            if((fabs(speed_.dot(axis_x))>0.4 || fabs(speed_.dot(axis_y))>0.4 || fabs(speed_.dot(axis_z))>0.4) || (dt<1000))
                            {
                                dt+=TIME_STEP;
                            }

//******************Collision Plan 0*************************//
//                            if(dist0<d0)
//                            {
//                                value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v) += (d0-dist0) * axis_y;
//                                value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v) -=
//                                    axis_y.dot(value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v)) * axis_y;
//                            }
//*******************END****************************//


//*******************Collision Deux Plans************************//


//                            total_vertex_count++;
//                            bool collision_p1 = false;
//                            bool collision_p2 = false;

//                            double friction_coefficient = 0.5;
//                            double minimal_speed = 0.01;

//                            Vec3& speed_ = value<Vec3>(*selected_mesh_, simu_solver.speed_.get(), v);

//                            if (dist < d)
//                            {
//                                collision_p1 = true;
//                                value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v) += (d - dist) * normal;

//                                Vec3 speed_normal1 = normal.dot(speed_) * normal;
//                                Vec3 speed_tangent1 = speed_ - speed_normal1;

//                                Vec3 friction_force_normal1 = friction_coefficient * speed_normal1;
//                                Vec3 friction_force_tangent1 = friction_coefficient * speed_tangent1;

//                                speed_ -= (friction_force_normal1 + friction_force_tangent1);
//                            }


//                            if (dist_p2 < d_p2)
//                            {
//                                collision_p2 = true;
//                                value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v) += (d_p2 - dist_p2) * normal_p2;

//                                Vec3 speed_normal2 = normal_p2.dot(speed_) * normal_p2;
//                                Vec3 speed_tangent2 = speed_ - speed_normal2;

//                                Vec3 friction_force_normal2 = friction_coefficient * speed_normal2;
//                                Vec3 friction_force_tangent2 = friction_coefficient * speed_tangent2;

//                                speed_ -= (friction_force_normal2 + friction_force_tangent2);
//                            }

//                            if (collision_p1 || collision_p2)
//                            {
//                                  contact_with_planes_count++;
//                            }

//                            if (collision_p1 && collision_p2 && speed_.norm() < minimal_speed) {
//                                Vec3 speed_normal1 = normal.dot(speed_) * normal;
//                                Vec3 speed_normal2 = normal_p2.dot(speed_) * normal_p2;

//                                Vec3 friction_force_normal1 = friction_coefficient * speed_normal1;
//                                Vec3 friction_force_normal2 = friction_coefficient * speed_normal2;

//                                speed_ -= (friction_force_normal1 + friction_force_normal2);

//                                if (speed_.norm() < minimal_speed) {
//                                    speed_ = Vec3(0, 0, 0);
//                                }
//                            }

                            //*******************END************************//

                            return true;
                         });//END parallel_foreach_cell

//                        double percentage_in_contact = (static_cast<double>(contact_with_planes_count) / total_vertex_count) * 100.0;
//                        std::cout << "Pourcentage des vertex en contact avec les deux plans: " << percentage_in_contact << "%" << std::endl;


//**************find energy********************//
                        double PSI_NEO =0;
                        double h = NUM_SUBSTEP / TIME_STEP;
                        foreach_cell(*selected_mesh_, [&](Volume v) -> bool {
                        PSI_NEO += simu_solver.energy_calculation(*selected_mesh_, v, h);
                        return true;
                          });
                        //std::cout << "Energie = " << PSI_NEO << std::endl;
                        std::ofstream csv_file(csv_data, std::ios::app);
                        if (csv_file.is_open())
                        {
                            csv_file << c_t << "," << PSI_NEO << std::endl;
                        }
                         csv_file.close();
//**************END********************//


//**************find centre de masse position********************//
                         double total_mass = 0.0;
                         Vec3 weighted_sum = Vec3::Zero();

                         foreach_cell(*selected_mesh_, [&](Vertex w) -> bool {
                         double vertex_mass = value<double>(*selected_mesh_, simu_solver.masse_, w);
                         Vec3 vertex_position = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), w);

                         total_mass += vertex_mass;
                         weighted_sum += vertex_mass * vertex_position;

                         return true;
                         });

                         Vec3 cm = weighted_sum / total_mass;

                         std::ofstream csv_file_(csv_data_position, std::ios::app);
                         if (csv_file_.is_open())
                         {
                             csv_file_ << c_t << "," << cm.dot(axis_y)-20 << std::endl;
                         }
                          csv_file_.close();
//**************END********************//

//**************stop condition********************//
                           if(static_cast<int>(c_t)==400)
                           {
                                std::cout<<"Temps de convergence = "<<dt<<std::endl;

                                running_ = false;
                                exit(EXIT_SUCCESS);

                                dt = 0;
                                c_t = 0;
                            }
                     }
                need_update_ = true;
            }
        });
    }

    void stop()
    {
        running_ = false;
    }


//************Draw Plan************//
    void draw(ui::View *view) override
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

            if (p.have_selected_vertex_ && p.param_edge_->attributes_initialized())
            {
                p.param_edge_->bind(proj_matrix, view_matrix);
                glDrawArrays(GL_LINES, 0, 2);
                p.param_edge_->release();
            }

            if (p.show_frame_manipulator_)
            {
                double size = (md.bb_max_ - md.bb_min_).norm() / 10;
                p.frame_manipulator_.set_size(size);
                p.frame_manipulator_.draw(true, true, proj_matrix, view_matrix);
            }


//                auto& m = selected_mesh_;
//                MeshData<MESH>& md = mesh_provider_->mesh_data(*m);
                Eigen::Vector3f bb_min = Eigen::Vector3f(md.bb_min_.x()+md.bb_max_.x()/2,
                                                         md.bb_min_.y()-10,
                                                         md.bb_min_.z()+md.bb_max_.z()/2);
                Eigen::Vector3f bb_max = Eigen::Vector3f(md.bb_max_.x()-md.bb_max_.x()/2,
                                                         md.bb_max_.y()+10,
                                                         md.bb_max_.z()-md.bb_max_.z()/2);
                bbmin = bb_min;
                bbmax = bb_max;

//                std::cout<<"bb min = "<<bbmin<<std::endl;
//                std::cout<<"bb max = "<<bbmax<<std::endl;
//                std::cout<<"min ="<<md.bb_min_.x()<<", "<<md.bb_min_.y()<<", "<<md.bb_min_.z()<<std::endl;
//                std::cout<<"max ="<<md.bb_max_.x()<<", "<<md.bb_max_.y()<<", "<<md.bb_max_.z()<<std::endl;

           }
            Eigen::Vector3f pos_cube_ (0, 30, 0) ;
            Eigen::Affine3f transfo0 = Eigen::Translation3f(bbmin+pos_cube_)*Eigen::Scaling(size_cube);
            shape_->draw(rendering::ShapeDrawer::CUBE, proj_matrix, view_matrix * transfo0.matrix());

            Eigen::Affine3f transfo1 = Eigen::Translation3f(bbmax-pos_cube_)*Eigen::Scaling(size_cube);
            shape_->draw(rendering::ShapeDrawer::CUBE, proj_matrix, view_matrix * transfo1.matrix());

            foreach_cell(*selected_mesh_, [&](Vertex v) -> bool {
                Parameters& p = parameters_[selected_mesh_];
                Vec3& pos = value<Vec3>(*selected_mesh_, p.vertex_position_.get(), v);

                if (pos.y() <= (bbmin.y()+10)+20 && pos.y() >= bbmin.y()+10)
                {
                    vertices_inferior.push_back(v);
                }

                if (pos.y() <= bbmax.y()-10 && pos.y() >= (bbmax.y()-10)-20)
                {
                    vertices_superior.push_back(v);
                }

                return true;
            });


//************Draw Plan ref************//
//                Eigen::Affine3f transfo0 = Eigen::Translation3f(pos_cube_ref)*Eigen::Scaling(size_cube_ref);
//                shape_->draw(rendering::ShapeDrawer::CUBE, proj_matrix, view_matrix * transfo0.matrix());
//************END************//

//************Draw Plan 0degree************//
//               Eigen::Affine3f transfo0 = Eigen::Translation3f(pos_cube)*Eigen::Scaling(size_cube);
//               shape_->draw(rendering::ShapeDrawer::CUBE, proj_matrix, view_matrix * transfo0.matrix());
//************END************//

//************Draw Plan 30degree************//
//            Eigen::Matrix3f rotation_matrix;
//            rotation_matrix = Eigen::AngleAxisf(incline_angle, Eigen::Vector3f::UnitX());
//            Eigen::Affine3f incline_transformation = Eigen::Affine3f::Identity();
//            incline_transformation.rotate(rotation_matrix);
//            Eigen::Affine3f transfo = Eigen::Translation3f(pos_cube)*Eigen::Scaling(size_cube);
//            Eigen::Affine3f final_transformation = incline_transformation * transfo;
//            shape_->draw(rendering::ShapeDrawer::CUBE, proj_matrix, view_matrix * final_transformation.matrix());
//************END************//

//************Draw Plan 150degree************//
//            Eigen::Matrix3f rotation_matrix_p2;
//            rotation_matrix_p2 = Eigen::AngleAxisf(incline_angle_p2, Eigen::Vector3f::UnitX());
//            Eigen::Affine3f incline_transformation_p2 = Eigen::Affine3f::Identity();
//            incline_transformation_p2.rotate(rotation_matrix_p2);
//            Eigen::Affine3f transfo_p2 = Eigen::Translation3f(pos_cube)*Eigen::Scaling(size_cube);
//            Eigen::Affine3f final_transformation_p2 = incline_transformation_p2 * transfo_p2;
//            shape_->draw(rendering::ShapeDrawer::CUBE, proj_matrix, view_matrix * final_transformation_p2.matrix());
//************END************//

    }
//************END************//


    //
    //********************************************************************************************************************************//


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
                    simu_solver.fixed_vertex = p.fixed_vertex;
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
  //--------------------------------------------------------------------------------------------------------------------------------//

public:
    MESH* selected_mesh_;

    std::vector<Volume> vec_volume;

    std::unordered_map<const MESH*, Parameters> parameters_;
    cgogn::simulation::XPBD<MESH> simu_solver;
    View* selected_view_;
    bool running_;
    bool apply_gravity;
    bool apply_force;
    bool draw_cube;
//  Eigen::Vector3f pos_cube;
//  Eigen::Vector3f size_cube;
    Eigen::Vector3f size_cube_ref;
    Eigen::Vector3f pos_cube_ref;
    Vec3 Zaxis_cube;

    rendering::ShapeDrawer* shape_;

    float32 c_t;
    //
    //********************************************************************************************************************************//

    std::vector<std::shared_ptr<boost::synapse::connection>> connections_;
    std::unordered_map<const MESH*, std::vector<std::shared_ptr<boost::synapse::connection>>> mesh_connections_;
    MeshProvider<MESH>* mesh_provider_;
    bool need_update_;

    //--------------------------------------------------------------------------------------------------------------------------------//

};

} // namespace ui

} // namespace cgogn

#endif // CGOGN_MODULE_XPBD_MODULE_H_
