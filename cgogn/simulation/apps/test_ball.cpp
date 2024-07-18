/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C), IGG Group, ICube, University of Strasbourg, France            *
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
#include <cgogn/core/types/maps/cmap/cmap3.h>

#include <cgogn/core/types/mesh_traits.h>
#include <cgogn/geometry/types/vector_traits.h>

#include <cgogn/ui/app.h>
#include <cgogn/ui/module.h>
#include <cgogn/ui/view.h>

#include <cgogn/core/types/maps/cmap/cmap3.h>

#include <cgogn/rendering/shape_drawer.h>

#include <GLFW/glfw3.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/volume.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>


//#include <cgogn/simulation/ui_modules/XPBD.h>

using Mesh = cgogn::CMap3;

using namespace ::cgogn;
using namespace ::cgogn::rendering;


class ShapeRender : public ui::ViewModule
{
public:
    ShapeRender(const ui::App & app)
        : ViewModule(app, "Test Shape"), app_(app), shape_(nullptr)
    {}

    ~ShapeRender()	{}

public:
    void init() override
    {
        app_.current_view()->set_scene_radius(5);
        app_.current_view()->set_scene_center(GLVec3(0,-2,0));
        shape_ = ShapeDrawer::instance();
        shape_->update_subdivision(24);
        shape_->color(ShapeDrawer::SPHERE) = GLColor(1, 0, 0, 1);
        shape_->color(ShapeDrawer::CUBE) = GLColor(1, 0, 0, 1);
    }

    void draw(ui::View * view) override
    {
        const GLMat4& proj_matrix = view->projection_matrix();
        const GLMat4& view_matrix = view->modelview_matrix();

        Eigen::Affine3f transfo = Eigen::Translation3f(GLVec3(4, 0, 0)) * Eigen::Scaling(2.0f);
        transfo = Eigen::Translation3f(GLVec3(0, 4, 0)) * Eigen::Scaling(2.0f);
        shape_->draw(ShapeDrawer::SPHERE, proj_matrix, view_matrix * transfo.matrix());


//      float incline_angle = M_PI / 4; // 30 degrees in radians

//      Eigen::Matrix3f rotation_matrix;
//      rotation_matrix = Eigen::AngleAxisf(incline_angle, Eigen::Vector3f::UnitX());
//      Eigen::Affine3f incline_transformation = Eigen::Affine3f::Identity();
//      incline_transformation.rotate(rotation_matrix);

        transfo = Eigen::Translation3f(GLVec3(0, -8, 0)) * Eigen::Scaling(GLVec3(5, 0.2f, 5));
        //Eigen::Affine3f final_transformation = incline_transformation * transfo;

        shape_->draw(ShapeDrawer::CUBE, proj_matrix, view_matrix * transfo.matrix());


        std::cout << "DRAW" << std::endl;
    }

        private:
            const ui::App& app_;
            ShapeDrawer* shape_;
        };


int main(int argc, char** argv)
{
    cgogn::thread_start();

    cgogn::ui::App app;
    app.set_window_title("Simple Shape Test");
    app.set_window_size(1000, 800);


    ShapeRender sr(app);
    //cgogn::ui::XPBD_Module<Mesh> xpbd(app);
    app.current_view()->link_module(&sr);

    //app.current_view()->link_module(&xpbd);
    app.init_modules();




    std::srand(std::time(nullptr));
    return app.launch();
}

