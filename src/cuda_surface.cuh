#ifndef CUDA_SURFACE_H
#define CUDA_SURFACE_H

#include <vector>
//#include <cuda_runtime.h>
#include <vector_types.h>

namespace CUDA_Surface {
    struct Points3D {
        int size;
        float *x;
        float *y;
        float *z;

        Points3D() {}
        // Points3D(const int _size, float4 *_x) : size(_size), x(_x) {}
        Points3D(const int _size, float *_x, float *_y, float *_z) : size(_size), x(_x), y(_y), z(_z) {}

        void allocate(const int _size);
        void free();
        void load_from(const Points3D &cpu) const;
        // void save_to(const Points3D &cpu) const;
    };

    struct Points1D {
        int size;
        float *s;

        Points1D() {}
        Points1D(const int _size, float *_s) : size(_size), s(_s) {}

        void allocate(const int _size);
        void free();
        // void load_from(const Points1D &cpu) const;
        void save_to(const Points1D &gpu) const;
    };

    float compute_surface(Points1D &gpu_surface, Points1D &gpu_max_time, 
                          const Points3D &gpu_centroids, const Points3D &gpu_pos); 
}

#endif