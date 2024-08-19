/*******************************************************************************
 * CGoGN: Combinatorial and Geometric modeling with Generic N-dimensional Maps  *
 * Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
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
#include <cgogn/ui/view.h>
#include <GLFW/glfw3.h>

#include <cgogn/core/functions/attributes.h>
#include <cgogn/core/functions/traversals/edge.h>
#include <cgogn/core/functions/traversals/volume.h>
#include <cgogn/core/ui_modules/mesh_provider.h>
#include <cgogn/modeling/algos/subdivision.h>
#include <cgogn/rendering/ui_modules/surface_render.h>
#include <cgogn/rendering/ui_modules/volume_render.h>
#include <cgogn/simulation/ui_modules/XPBD.h>
#include <cgogn/simulation/algos/XPBD/XPBD.h>

#include <fstream>
#include <cstring>
#include <iostream>

#include <cgogn/io/reflect/reflect.hpp>
#include <cgogn/io/reflect/reflect.std.map.hpp>
#include <cgogn/io/reflect/reflect.std.unordered_map.hpp>
#include <cgogn/io/reflect/reflect.std.vector.hpp>

#include <cgogn/io/reflect/codecs/json/decoder.hpp>
#include <cgogn/io/reflect/codecs/json/encoder.hpp>

#define DEFAULT_MESH_PATH CGOGN_STR(CGOGN_DATA_PATH)

using Mesh = cgogn::CMap3;
using MAP = cgogn::CMap3;

template <typename T>
using Attribute = typename cgogn::mesh_traits<Mesh>::Attribute<T>;
using Vertex = typename cgogn::mesh_traits<Mesh>::Vertex;
using Edge = typename cgogn::mesh_traits<Mesh>::Edge;
using Face = typename cgogn::mesh_traits<Mesh>::Face;
using Volume = typename cgogn::mesh_traits<Mesh>::Volume;

using Vec3 = cgogn::geometry::Vec3;

struct parameters_a
{
    reflect_fields(
            ((int), density),
            ((float), poisson_ratio),
            ((int), young_modulus),
            ((float), time_step)
            );
};

struct parameters_b
{
    reflect_fields(
            ((int), position_x1),
            ((int), position_y1),
            ((int), position_z1),
            ((int), dimension_x1),
            ((int), dimension_y1),
            ((int), dimension_z1),
            ((int), angle1),
            ((float), friction_coef1)
            );
};

struct parameters_c
{
    reflect_fields(
            ((int), position_x2),
            ((int), position_y2),
            ((int), position_z2),
            ((int), dimension_x2),
            ((int), dimension_y2),
            ((int), dimension_z2),
            ((int), angle2),
            ((float), friction_coef2)
            );
};

struct json_data
{
    reflect_fields(
            ((parameters_a), param_sim),
            ((parameters_b), plan1),
            ((parameters_c), plan2)
            );
};

int main(int argc, char** argv)
{
    std::string filename;
    std::string filename_csv_energy;
    std::string filename_csv_position;

    if (argc < 5)
    {
        filename = std::string("/home/yehya/Desktop/JSON/cube86.mesh");
    }
    filename = std::string(argv[1]);
    const char* json_file = argv[2];
    filename_csv_energy = std::string(argv[3]);
    filename_csv_position = std::string(argv[4]);

    cgogn::thread_start();

    cgogn::ui::App app;
    app.set_window_title("Shape Matching");
    app.set_window_size(1000, 800);

    cgogn::ui::MeshProvider<Mesh> mp(app);
    cgogn::ui::VolumeRender<Mesh> mrsr(app);
    cgogn::ui::SimXPBD<Mesh> xpbd(app);

    cgogn::ui::View* v1 = app.current_view();
    v1->link_module(&mp);
    v1->link_module(&mrsr);
    v1->link_module(&xpbd);

    app.init_modules();

    Mesh* m = mp.load_volume_from_file(filename);
    if (!m)
    {
        std::cout << "File could not be loaded" << std::endl;
        return 1;
    }

    std::shared_ptr<Attribute<Vec3>> position = cgogn::get_attribute<Vec3, Vertex>(*m, "position");

    cgogn::foreach_cell(*m, [&](Vertex v) -> bool {
        cgogn::value<Vec3>(*m, position, v) *= 100;
        return true;
    });

    cgogn::index_cells<Mesh::Volume>(*m);
    cgogn::index_cells<Mesh::Edge>(*m);
    cgogn::index_cells<Mesh::Face>(*m);

    mp.emit_attribute_changed(*m, position.get());

    mrsr.set_vertex_position(*v1, *m, position);

    //xpbd.set_vertex_position(*m, position);

    //*****************JSON**********************//
    std::cout<<"test"<<std::endl;

    reflect::codecs::json::preferences prefs;
    prefs.indent = "    ";
    prefs.newline = "\n";
    prefs.newline_at_eof = true;

    reflect::stream_writer writer(std::cout);
    reflect::codecs::json::encoder encode(writer,prefs);

    parameters_a param_sim{1800, 0.46f, 100000000, 0.03f};
    parameters_b plan1{0, -300, 0, 300, 20, 300, 30, 0.5f};
    parameters_c plan2{0, -300, 0, 300, 20, 300, -30, 0.5f};

    json_data data{param_sim, plan1, plan2};

    encode(data);

    std::ifstream file(json_file);

    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << json_file << std::endl;
        return 1;
    }

    std::string buffer;
    std::string line;
    while (std::getline(file, line)) {
        buffer += line;
    }
    file.close();

    reflect::string_reader reader(buffer);
    reflect::codecs::json::decoder decode(reader);

    decode(data);

    int density_ = data.param_sim.density;
    float poisson_ratio_ = data.param_sim.poisson_ratio;
    int young_modulus_ = data.param_sim.young_modulus;
    float time_step_ = data.param_sim.time_step;


    int position_x_1 = data.plan1.position_x1;
    int position_y_1 = data.plan1.position_y1;
    int position_z_1 = data.plan1.position_z1;
    int dimension_x_1 = data.plan1.dimension_x1;
    int dimension_y_1 = data.plan1.dimension_y1;
    int dimension_z_1 = data.plan1.dimension_z1;
    int angle_1 = data.plan1.angle1;
    float friction_coef_1 = data.plan1.friction_coef1;

    int position_x_2 = data.plan2.position_x2;
    int position_y_2 = data.plan2.position_y2;
    int position_z_2 = data.plan2.position_z2;
    int dimension_x_2 = data.plan2.dimension_x2;
    int dimension_y_2 = data.plan2.dimension_y2;
    int dimension_z_2 = data.plan2.dimension_z2;
    int angle_2 = data.plan2.angle2;
    float friction_coef_2 = data.plan2.friction_coef2;

    setDensity(density_);
    setPoissonRatio(poisson_ratio_);
    setYoungModule(young_modulus_);
    setTimeStep(time_step_);

//    setPosition(position_x_1, position_y_1, position_z_1);
//    setDimension(dimension_x_1, dimension_y_1, dimension_z_1);
//    setAngle(angle_1);
//    setFrictionCoef(friction_coef_1);

    setPlan1(position_x_1, position_y_1, position_z_1, dimension_x_1, dimension_y_1, dimension_z_1, angle_1, friction_coef_1);
    setPlan2(position_x_2, position_y_2, position_z_2, dimension_x_2, dimension_y_2, dimension_z_2, angle_2, friction_coef_2);

    encode(data);

    Setcsv_energy(filename_csv_energy);
    Setcsv_position(filename_csv_position);
    //*******************END JSON******************************//

    //    mp.emit_attribute_changed(*m1, position1.get());

    std::srand(std::time(nullptr));

    //xpbd.set_vertex_position(*m, position);
    //xpbd.startSimulation();
    xpbd.call_set_vertex_position();

    return app.launch();
}
