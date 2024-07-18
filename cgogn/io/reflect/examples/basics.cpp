////usr/bin/env $(dirname $0)/cxx -c++17 -I $(dirname $0)/../.. -r $0; exit $?

//------------------------------------------------------------------------------
// Declare a struct or class with some reflected fields:

#include <cstring>
#include <cgogn/io/reflect/reflect.hpp>

#include <fstream>


//------------------------------------------------------------------------------
// Additional headers provide reflection for standard library types,
// std::vector, std::map, and std::unordered_map:

#include <cgogn/io/reflect/reflect.std.map.hpp>
#include <cgogn/io/reflect/reflect.std.unordered_map.hpp>
#include <cgogn/io/reflect/reflect.std.vector.hpp>

struct example_struct_2 {
    // use the reflect_fields macro to declare reflected fields
    reflect_fields(
        ((float),masse),
        ((float),cp)

    )

};

//------------------------------------------------------------------------------
// Serialize to/from JSON:

#include <cgogn/io/reflect/codecs/json/decoder.hpp>
#include <cgogn/io/reflect/codecs/json/encoder.hpp>

int main(int,char**) {
    reflect::codecs::json::preferences prefs;
    prefs.indent = "    ";
    prefs.newline = "\n";
    prefs.newline_at_eof = true;

    reflect::stream_writer writer(std::cout);
    reflect::codecs::json::encoder encode(writer,prefs);

    example_struct_2 s{55.5f, 0.2f};
    encode(s);
    // {
    //     "i":1,
    //     "f":2.3,
    //     "s":"hello"dhadi yehya
    // }


    std::ifstream file("/home/yehya/Desktop/file.json");
    if (!file.is_open()){
        std::cerr<<"erreur";
        return 1;
    }


std::string line;
std::getline(file, line);
    reflect::string_reader reader(R"({ "massse":2.0, "cf":3.4})");
    reflect::codecs::json::decoder decode(reader);
    decode(s);
    encode(s);
    // {
    //     "i":2,
    //     "f":3.5,
    //     "s":"world"
    // }

//    example_struct_2 s2 {
//        /*numbers*/{{1,2.3f},{4,5.6f}},
//        /*people*/{
//            {"me",{"Me","123 My St."}},
//            {"you",{"You","456 Your St."}},
//        },
//        /*vector*/{
//            {1,2,3},
//            {4,5,6},
//        },
//        /*bitmap*/{
//            false,
//            true,
//        },
//    };
//    encode(s2);
    // {
    //     "numbers":{
    //         "1":2.3,
    //         "4":5.6
    //     },
    //     "people":{
    //         "you":{
    //             "name":"You",
    //             "address":"456 Your St."
    //         },
    //         "me":{
    //             "name":"Me",
    //             "address":"123 My St."
    //         }
    //     },
    //     "vector":[
    //         {
    //             "x":1,
    //             "y":2,
    //             "z":3
    //         },
    //         {
    //             "x":4,
    //             "y":5,
    //             "z":6
    //         }
    //     ],
    //     "bitmap":[
    //         false,
    //         true
    //     ]
    // }

    return 0;
}
