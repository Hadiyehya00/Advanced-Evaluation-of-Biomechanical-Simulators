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

#include <cgogn/core/functions/attributes.h>

#include <GLFW/glfw3.h>
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

struct parameters_ {
    reflect_fields(
        ((int), density),
        ((float), poisson_ratio),
        ((int), young_modulus),
        ((float), time_step),
        ((int), position_x),
        ((int), position_y),
        ((int), position_z),
        ((int), dimension_x),
        ((int), dimension_y),
        ((int), dimension_z),
        ((int), angle),
        ((float), friction_coef)
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
    cgogn::ui::SimBall<Mesh> xpbd(app);

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

        parameters_ p {1800, 0.46f, 100000000, 0.02f, 0, -300, 0, 400, 20, 400, 0, 1.0f};

        encode(p);

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

            decode(p);

            int density_ = p.density;
            float poisson_ratio_ = p.poisson_ratio;
            int young_modulus_ = p.young_modulus;
            float time_step_ = p.time_step;
            int position_x_ = p.position_x;
            int position_y_ = p.position_y;
            int position_z_ = p.position_z;
            int dimension_x_ = p.dimension_x;
            int dimension_y_ = p.dimension_y;
            int dimension_z_ = p.dimension_z;
            int angle_ = p.angle;
            float friction_coef_ = p.friction_coef;


            setDensity(density_);
            setPoissonRatio(poisson_ratio_);
            setYoungModule(young_modulus_);
            setTimeStep(time_step_);

            setPosition(position_x_, position_y_, position_z_);
            setDimension(dimension_x_, dimension_y_, dimension_z_);
            setAngle(angle_);
            setFrictionCoef(friction_coef_);

            encode(p);

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
