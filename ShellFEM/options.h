#pragma once

#include <string>

#include "utils.h"
#include "s3element.h"
#include "fem.h"

void parse_arguments(int argc, char** argv, SimulationProperties &simProps, JobProperties &jobProps);