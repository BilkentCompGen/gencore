#ifndef RBAM_H
#define RBAM_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "args.h"
#include "lps.h"
#include "chtslib.h"
#include "logging.h"

void read_bam( struct targs& arguments, const struct pargs& program_arguments );

#endif