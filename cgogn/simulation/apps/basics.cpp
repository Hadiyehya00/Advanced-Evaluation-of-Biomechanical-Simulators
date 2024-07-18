////usr/bin/env $(dirname $0)/cxx -c++17 -I $(dirname $0)/../.. -r $0; exit $?


#include <fstream>
#include <cstring>

#include <iostream>

#include <cgogn/io/reflect/reflect.hpp>
#include <cgogn/io/reflect/reflect.std.map.hpp>
#include <cgogn/io/reflect/reflect.std.unordered_map.hpp>
#include <cgogn/io/reflect/reflect.std.vector.hpp>

#include <cgogn/io/reflect/codecs/json/decoder.hpp>
#include <cgogn/io/reflect/codecs/json/encoder.hpp>


struct parameters_ {
    reflect_fields(
        ((int),density),
        ((float),cp)
    );
};
int i = 0;
int main(int argc, char** argv) {

std::cout<<"test"<<std::endl;
std::cout<<i<<std::endl;
i++;
std::cout<<i<<std::endl;
    reflect::codecs::json::preferences prefs;
    prefs.indent = "    ";
    prefs.newline = "\n";
    prefs.newline_at_eof = true;

    reflect::stream_writer writer(std::cout);
    reflect::codecs::json::encoder encode(writer,prefs);

    parameters_ p {500, 0.46f};
    encode(p);


    if (argc != 2) {
            std::cerr << "Usage: " << argv[0] << " <json_file>" << std::endl;
            return 1;
        }

        const char* json_file = argv[1];
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
        encode(p);

return 0;
}
