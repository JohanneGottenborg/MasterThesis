#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('/faststorage/project/ostrich_thermal/people/johanneg/workflows/JohanneG_popgen_w_GLs/workflow_source/'))
from popgen_w_GLs_workflow_source import *

gwf = running_GL_population_genomics_analyses() 