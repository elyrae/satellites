#include "cuda_surface.cuh"

#include <cstdio>
#include <cuda.h>

void CUDA_Surface::Points3D::allocate()
{
    cudaMalloc((void**) &x, size*sizeof(float));
    cudaMalloc((void**) &y, size*sizeof(float));
    cudaMalloc((void**) &z, size*sizeof(float));
}

void CUDA_Surface::Points3D::allocate(const int _size)
{
    size = _size;
    cudaMalloc((void**) &x, size*sizeof(float));
    cudaMalloc((void**) &y, size*sizeof(float));
    cudaMalloc((void**) &z, size*sizeof(float));
}

void CUDA_Surface::Points3D::free()
{
    cudaFree(x);
    cudaFree(y);
    cudaFree(z);
}

void CUDA_Surface::Points3D::load_from(const CUDA_Surface::Points3D &cpu) const
{
    cudaMemcpy(x, cpu.x, size*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(y, cpu.y, size*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(z, cpu.z, size*sizeof(float), cudaMemcpyHostToDevice);
}

void CUDA_Surface::Points3D::save_to(const CUDA_Surface::Points3D &cpu) const
{
    cudaMemcpy(cpu.x, x, size*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cpu.y, y, size*sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(cpu.z, z, size*sizeof(float), cudaMemcpyDeviceToHost);
}

void CUDA_Surface::Points1D::allocate()
{
    cudaMalloc((void**) &s, size*sizeof(float));
    cudaMemset(s, 0, size);
}

void CUDA_Surface::Points1D::allocate(const int _size)
{
    size = _size;
    cudaMalloc((void**) &s, size*sizeof(float));
    cudaMemset(s, 0, size);
}

void CUDA_Surface::Points1D::free()
{
    cudaFree(s);
}

void CUDA_Surface::Points1D::load_from(const CUDA_Surface::Points1D &cpu) const
{
    cudaMemcpy(s, cpu.s, size*sizeof(float), cudaMemcpyHostToDevice);
}

void CUDA_Surface::Points1D::save_to(const CUDA_Surface::Points1D &cpu) const
{
    cudaMemcpy(cpu.s, s, size*sizeof(float), cudaMemcpyDeviceToHost);
}


__global__ void compute_surface_layer(float *surface, const float *cx, const float *cy, const float *cz, 
                                                      const float *x, const float *y, const float *z, const float hor)
{
    const int SURFACE_SIZE       = 5120;
    const int CONFIGURATION_SIZE = 20;
    const int CONFIGURATIONS     = 100;
    const int TIMESTEPS          = 240;

    const int global_conf  = blockDim.x*blockIdx.x + threadIdx.x; // configuration 
    const int global_point = blockDim.y*blockIdx.y + threadIdx.y; // point
    const int global_time  = blockDim.z*blockIdx.z + threadIdx.z; // time

    float surface_result = 0.0, sat_x = 0.0, sat_y = 0.0, sat_z = 0.0; 
    for (int iorb = 0; iorb < CONFIGURATION_SIZE; iorb++)
    {
        sat_x = x[(global_time*CONFIGURATIONS + global_conf)*CONFIGURATION_SIZE + iorb];
        sat_y = y[(global_time*CONFIGURATIONS + global_conf)*CONFIGURATION_SIZE + iorb];
        sat_z = z[(global_time*CONFIGURATIONS + global_conf)*CONFIGURATION_SIZE + iorb];
        // sat_h = h[(global_time*CONFIGURATIONS + global_conf)*CONFIGURATION_SIZE + iorb]    
        surface_result |= (sat_x*cx[global_point] + sat_y*cy[global_point] + sat_z*cz[global_point] > sat_h);
    } 
    surface[(global_time*SURFACE_SIZE + global_point)*CONFIGURATIONS + global_conf] = surface_result;
}