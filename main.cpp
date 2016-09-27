#include <iostream>
#include <sstream>
#include <fstream>
#include <CL/cl.hpp>
#include <ctime>

const int screenWidth=640;
const int screenHeight=480;
 
typedef struct {
    cl_float3 P,D;
} Ray;
Ray cons_Ray(cl_float3 p, cl_float3 d){
    Ray ray; ray.P=p; ray.D=d; return ray;
}
Ray init_Ray(){
    return cons_Ray((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 1.0f});
}

typedef struct{
    cl_float t;
    cl_float3 P,N;
} Hit;
Hit cons_Hit(cl_float t, cl_float3 p, cl_float3 n){
    Hit hit; hit.t=t; hit.P=p; hit.N=n; return hit;
}
Hit init_Hit(){
    return cons_Hit((cl_float){-1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.0f, 0.0f, 0.0f});
}

typedef struct {
    cl_float3 r1,r2,r3,N;
} Triangle;
Triangle cons_Triangle(cl_float3 r1, cl_float3 r2, cl_float3 r3, cl_float3 n){
    Triangle tri; tri.r1=r1; tri.r2=r2; tri.r3=r3; tri.N=n; return tri;
}
Triangle init_Triangle(){
    return cons_Triangle((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, -1.0f});
}

class Scene{
private:
    std::vector<Triangle> tris;
    int tris_size;
    int rays_size=screenWidth*screenHeight;
    
    cl::Context context;
    cl::Program program;
    cl::Buffer buffer_tris;
    cl::Buffer buffer_rays;
    cl::Buffer buffer_hits;
    cl::CommandQueue queue;
public:
    void init_Scene(){
        //get all platforms (drivers)
        std::vector<cl::Platform> all_platforms;
        cl::Platform::get(&all_platforms);
        if(all_platforms.size()==0){
            std::cout<<" No platforms found. Check OpenCL installation!\n";
            exit(1);
        }
        cl::Platform default_platform=all_platforms[0];
        std::cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";
        
        //get default device of the default platform
        std::vector<cl::Device> all_devices;
        default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
        if(all_devices.size()==0){
            std::cout<<" No devices found. Check OpenCL installation!\n";
            exit(1);
        }
        cl::Device default_device=all_devices[0];
        std::cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";
        
        //read the source file
        std::ifstream inFile("prog.cl");
        std::stringstream strStream;
        strStream << inFile.rdbuf();
        std::string str = strStream.str();
        std::string kernel_code=str;
        
        cl::Program::Sources sources;
        sources.push_back({kernel_code.c_str(),kernel_code.length()});
        
        //build the source file
        context=cl::Context({default_device});
        program=cl::Program(context,sources);
        if(program.build({default_device})!=CL_SUCCESS){
            std::cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";
            exit(1);
        }
        
        //create queue to which we will push commands for the device.
        queue=cl::CommandQueue(context,default_device);
    }
    void add_Triangle(Triangle tri){
        tris.push_back(tri);
    }
    void upload_Triangles(){
        tris_size=tris.size();
        buffer_tris=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(Triangle)*tris_size);
        queue.enqueueWriteBuffer(buffer_tris,CL_TRUE,0,sizeof(Triangle)*tris_size,&tris[0]);
    }
    void trace_rays(){
        cl::Buffer buffer_rays(context,CL_MEM_READ_ONLY,sizeof(Ray)*rays_size);
        cl::Buffer buffer_hits(context,CL_MEM_WRITE_ONLY ,sizeof(Hit)*rays_size);
        
        Ray* rays=new Ray[rays_size];
        Hit* hits=new Hit[rays_size];
        
        for(int i=0;i<rays_size;++i)
        rays[i]=cons_Ray((cl_float3){-0.99f+i/(float)rays_size, -0.99f+i/(float)rays_size, -10.0f}, (cl_float3){0.0f, 0.0f, 1.0f});
        
        //write data to the device
        queue.enqueueWriteBuffer(buffer_rays,CL_TRUE,0,sizeof(Ray)*rays_size,rays);
        
        clock_t begin=clock();
        //run the kernel
        cl::Kernel kernel_trace_ray=cl::Kernel(program,"trace_ray");
        kernel_trace_ray.setArg(0,buffer_tris);
        kernel_trace_ray.setArg(1,tris_size);
        kernel_trace_ray.setArg(2,buffer_rays);
        kernel_trace_ray.setArg(3,buffer_hits);
        
        queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(rays_size),cl::NullRange);

        queue.enqueueReadBuffer(buffer_hits,CL_TRUE,0,sizeof(Hit)*rays_size,hits);
        
        clock_t end=clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        printf("%f\n", elapsed_secs);
        //for(int i=0;i<rays_size;++i){
        for(int i=0;i<50;++i){
            Hit hit=hits[i];
            printf("hits[%03d]=\tt=%06.2f \tP=[%06.2f %06.2f %06.2f] \tN=[%06.2f %06.2f %06.2f]\n", i, hit.t, hit.P.s[0], hit.P.s[1], hit.P.s[2], hit.N.s[0], hit.N.s[1], hit.N.s[2]);
        }
    }
    void finish(){
        queue.finish();
    }
};

int main(){
    Scene scene;
    scene.init_Scene();
    
    for(int i=0;i<50;++i){
        scene.add_Triangle(cons_Triangle((cl_float3){-1.0f, -1.0f, 0.0f+i}, (cl_float3){0.0f, 1.0f, 0.0f+i}, (cl_float3){1.0f, -1.0f, 0.0f+i}, (cl_float3){0.0f, 0.0f, -1.0f}));
    }
    scene.upload_Triangles();
    
    scene.trace_rays();
    scene.finish();
    
    return 0;
}