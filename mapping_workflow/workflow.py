#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('/faststorage/project/ostrich_thermal/people/johanneg/workflows/mapping_and_filtering/workflow_source/'))
from mapping_filtering_workflow_source import *

gwf = filtering_raw_reads_mapping_filtering_bams()