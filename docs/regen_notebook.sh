#!/bin/sh
# disable AVX-512 features by default, but let an empty value override it
NPY_DISABLE_CPU_FEATURES=${NPY_DISABLE_CPU_FEATURES-AVX512F AVX512CD AVX512_SKX}
export NPY_DISABLE_CPU_FEATURES
jupyter nbconvert --to=notebook --execute --ExecutePreprocessor.record_timing=False --inplace "$@"
