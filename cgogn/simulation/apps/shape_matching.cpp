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
        ((float), cp),
        ((int), ym)

    );
};

struct pos_cube_parameters{
    reflect_fields(
        ((int), x),
        ((int), y),
        ((int), z),
        ((int), width),
        ((int), height),
        ((int), depth)
    );
};

//struct size_cube_parameters {
//    reflect_fields(

//    );
//};

int main(int argc, char** argv)
{

//*******************************************************************************************************//

    std::string filename;
    std::string filename_csv_energy;
    std::string filename_csv_position;

    //std::string filename1;

    if (argc < 6)
    {
        filename = std::string("/home/yehya/Desktop/JSON/cube86.mesh");
        //filename1 = std::string("/home/yehya/Downloads/cube86.mesh");
    }
        filename = std::string(argv[1]);
        const char* json_file = argv[2];
        filename_csv_energy = std::string(argv[3]);
        filename_csv_position = std::string(argv[4]);
        const char* json_file_plan = argv[5];


    cgogn::thread_start();

    cgogn::ui::App app;
    app.set_window_title("Shape Matching");
    app.set_window_size(1000, 800);

    cgogn::ui::MeshProvider<Mesh> mp(app);
    cgogn::ui::VolumeRender<Mesh> mrsr(app);
    cgogn::ui::SimBall<Mesh> xpbd(app);
    //cgogn::ui::SimBall::Parameters<Mesh> p;
    //cgogn::simulation::XPBD<MAP> xp(app);

    cgogn::ui::View* v1 = app.current_view();
    v1->link_module(&mp);
    v1->link_module(&mrsr);
    v1->link_module(&xpbd);
    //v1->link_module(&xp);

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


//    Vec3 shiftVec = cgogn::rendering::GLVec3d(1,200,1);

//if(xpbd.c_t ==5000){
//    cgogn::foreach_cell(*m, [&](Vertex v) -> bool {
//        cgogn::value<Vec3>(*m, position, v) = cgogn::value<Vec3>(*m, position, v)+shiftVec ;
//        return true;
//    });
//}

//xpbd.startSimulation();

//if (xpbd.need_update_)
//{

//    p.update_move_vertex_vbo();

//    xpbd.mesh_provider_->emit_attribute_changed(*m, position.get());


//    xpbd.need_update_ = false;
//}


    //************************************************
//    Mesh* m1 = mp.load_volume_from_file(filename1);
//    if (!m1)
//    {
//        std::cout << "File could not be loaded" << std::endl;
//        return 1;
//    }

//    std::shared_ptr<Attribute<Vec3>> position1 = cgogn::get_attribute<Vec3, Vertex>(*m1, "position");

//    cgogn::foreach_cell(*m1, [&](Vertex v) -> bool {
//        cgogn::value<Vec3>(*m1, position1, v) *= 100;
//        return true;
//    });

//    cgogn::index_cells<Mesh::Volume>(*m1);
//    cgogn::index_cells<Mesh::Edge>(*m1);
//    cgogn::index_cells<Mesh::Face>(*m1);

//    mrsr.set_vertex_position(*v1, *m1, position1);
//    Vec3 shiftVec = cgogn::rendering::GLVec3d(0,200,0);
//    cgogn::foreach_cell(*m, [&](Vertex v) -> bool {
//        cgogn::value<Vec3>(*m, position1, v) = cgogn::value<Vec3>(*m, position1, v) +shiftVec;
//        return true;
//    });



//            xpbd.selected_mesh_ = m;
//            xpbd.start();








//***********************************************************************************************************************//



    std::cout<<"test"<<std::endl;

        reflect::codecs::json::preferences prefs;
        prefs.indent = "    ";
        prefs.newline = "\n";
        prefs.newline_at_eof = true;

        reflect::stream_writer writer(std::cout);
        reflect::codecs::json::encoder encode(writer,prefs);

        parameters_ p {500, 0.46f};

        pos_cube_parameters pos_params {0, -200, 100, 250, 20, 250} ;

//        size_cube_parameters size_params {250, 20, 250};


        encode(p);
        encode(pos_params);
//        encode(size_params);


//        if (argc < 3) {
//                std::cerr << "Usage: " << argv[0] << " <json_file>" << std::endl;
//                return 1;
//            }

//            const char* json_file = argv[2];
            std::ifstream file(json_file);

            if (!file.is_open()) {
                std::cerr << "Failed to open file: " << json_file << std::endl;
                return 1;
            }

        std::string line;

        std::getline(file, line);
            //std::istringstream lineStream(line);
            reflect::string_reader reader(line);
            reflect::codecs::json::decoder decode(reader);
            decode(p);
            int density_ = p.density;
            float poisson_ratio_ = p.cp;
            int young_modulus_ = p.ym;
            setDensity(density_);
            setPoissonRatio(poisson_ratio_);
            setYoungModule(young_modulus_);

            encode(p);


            Setcsv_energy(filename_csv_energy);
            Setcsv_position(filename_csv_position);


            std::ifstream file_plan(json_file_plan);

            if (!file_plan.is_open()) {
                std::cerr << "Failed to open file: " << json_file_plan << std::endl;
                return 1;
            }

        std::string line_plan;

        std::getline(file_plan, line_plan);

            reflect::string_reader reader_plan(line_plan);
            reflect::codecs::json::decoder decode_plan(reader_plan);
            decode_plan(pos_params);
//            decode_plan(size_params);

            int x = pos_params.x;
            int y = pos_params.y;
            int z = pos_params.z;
            int width = pos_params.width;
            int height = pos_params.height;
            int depth = pos_params.depth;


//            setCP(poisson_ratio_);

            setPlanPosition(x, y, z, width, height, depth);
//            SetPlanSize(width, height, depth);
            encode(pos_params);
//            encode(size_params);


//***********************************************************************************************************************//

//    mp.emit_attribute_changed(*m1, position1.get());

    std::srand(std::time(nullptr));

    //xpbd.set_vertex_position(*m, position);
    //xpbd.startSimulation();
    xpbd.call_set_vertex_position();

    return app.launch();
}
