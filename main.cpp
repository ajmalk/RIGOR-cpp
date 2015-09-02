#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/contrib/contrib.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iterator>
#include <climits>

#include "Parameters.h"
#include "PropGenerator.h"

namespace fs = boost::filesystem;
namespace po = boost::program_options;

using namespace cv;
using namespace std;



int main(int argc, char* argv[])
{
    po::options_description desc("Allowed options and parameters");
    desc.add_options()
    ("help,h",                                                         "produce help message")
    ("imfilepath,i",      po::value<fs::path>()->required(),           "file path to the input image")
    ("spxnum,s",          po::value<int>()->default_value(1000),       "desired number of superpixels")
    ("spxcompact,c",      po::value<double>()->default_value(20),      "compactness factor. use a value ranging from 10 to 40 - higher value would give more regular-shaped superpixels")
    ("seeds,g", po::value<int>()->default_value(0), "number of seeds to generate")
    ("radius,r", po::value<int>()->default_value(1), "seed radius")
    ;
    
    try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        
        if (vm.count("help")) {
            std::cerr << desc << std::endl;
            return 0;
        }
        
        po::notify(vm);
        
        fs::path im_filepath = vm["imfilepath"].as<fs::path>();
        
        Parameters params = Parameters::getDebug();
        PropGenerator generator = PropGenerator(params);
        
        generator.generate(im_filepath);
        
        
    } catch (po::error& e) {
        std::cerr << "Program Options error: " << e.what() << std::endl;
        std::cerr << desc << std::endl;
        return -1;
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "Unknown error!" << std::endl;
        return -1;
    }
    
    return 0;
}