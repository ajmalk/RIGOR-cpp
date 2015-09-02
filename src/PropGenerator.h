//
//  PropGenerator.h
//  Rigor
//
//  Created by Ajmal Kunnummal on 2/20/15.
//  Copyright (c) 2015 ajmal. All rights reserved.
//

#ifndef __Rigor__PropGenerator__
#define __Rigor__PropGenerator__

#include <boost/filesystem.hpp>

#include "Parameters.h"

namespace fs = boost::filesystem;

class PropGenerator {
private:
   Parameters params;
public:
   PropGenerator();
   PropGenerator(Parameters &params);
   void generate(fs::path filename);
};

#endif /* defined(__Rigor__PropGenerator__) */
