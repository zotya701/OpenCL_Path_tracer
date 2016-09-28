#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <CL/cl.hpp>

const int screenWidth=1;
const int screenHeight=10;

typedef struct{
    cl_float3 kd,ks,emission,F0;    //diffuse, specular, emission, Fresnel
    float n, shininess, glossiness;
    int type;   //0-diffuse, 1-x, 2-x, 3-Emitter
} Material;

typedef struct {
    cl_float3 P,D;  //origo and direction
} Ray;
Ray cons_Ray(cl_float3 p, cl_float3 d){
    Ray ray; ray.P=p; ray.D=d; return ray;
}

typedef struct{
    cl_float t;         //time
    cl_float3 P,N;      //hitposition and normal vector in hitposition
    Material mat;  //material of the triangle
} Hit;

typedef struct {
    cl_float3 r1,r2,r3,N;   //vertices of the triangle and it's normal vector
    Material mat;
} Triangle;
Triangle cons_Triangle(cl_float3 r1, cl_float3 r2, cl_float3 r3, cl_float3 n){
    Triangle tri; tri.r1=r1; tri.r2=r2; tri.r3=r3; tri.N=n; return tri;
}

class Color{
public:
    float r,g,b;
    Color(){
        r=g=b=0.0f;
    }
    Color(float r, float g, float b){
        this->r=r; this->g=g; this->b;
    }
};

Color color_image[screenWidth*screenHeight];

class Scene{
private:
    std::vector<Triangle> tris;
    int tris_size;
    int rays_size=screenWidth*screenHeight;
    Ray* rays;
    cl_float3* cl_float3_image;
    
    cl::Context context;
    cl::Program program;
    cl::Buffer buffer_tris;
    cl::Buffer buffer_rays;
    cl::Buffer buffer_colors;
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
        
        buffer_rays=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(Ray)*rays_size);
        buffer_colors=cl::Buffer(context,CL_MEM_WRITE_ONLY ,sizeof(cl_float3)*rays_size);
        
        rays=new Ray[rays_size];
        cl_float3_image=new cl_float3[rays_size];
    }
    void add_Triangle(Triangle tri){
        tris.push_back(tri);
    }
    void upload_Triangles(){
        tris_size=tris.size();
        buffer_tris=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(Triangle)*tris_size);
        queue.enqueueWriteBuffer(buffer_tris,CL_TRUE,0,sizeof(Triangle)*tris_size,&tris[0]);
    }
    void generate_rays(){
        for(int i=0;i<rays_size;++i)
            rays[i]=cons_Ray((cl_float3){1.0f+i/(float)rays_size, 1.0f+i/(float)rays_size, -10.0f}, (cl_float3){0.0f, 0.0f, 1.0f});
        
        //write data to the device
        queue.enqueueWriteBuffer(buffer_rays,CL_TRUE,0,sizeof(Ray)*rays_size,rays);
    }
    void trace_rays(){
        clock_t begin=clock();
        //run the kernel
        cl::Kernel kernel_trace_ray=cl::Kernel(program,"trace_ray");
        kernel_trace_ray.setArg(0,buffer_tris);
        kernel_trace_ray.setArg(1,tris_size);
        kernel_trace_ray.setArg(2,buffer_rays);
        kernel_trace_ray.setArg(3,buffer_colors);
        
        queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(rays_size),cl::NullRange);

        queue.enqueueReadBuffer(buffer_colors,CL_TRUE,0,sizeof(cl_float3)*rays_size,cl_float3_image);
        
        clock_t end=clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        printf("%f\n", elapsed_secs);
        for(int i=0;i<10;++i){
            cl_float3 c=cl_float3_image[i];
            printf("r=%6.2f g=%6.2f b=%6.2f\n", c.s[0], c.s[1], c.s[2]);
        }
    }
    void finish(){
        queue.finish();
    }
};

int main(){
    Scene scene;
    scene.init_Scene();
    
    for(int i=0;i<25;++i){
        scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 0.0f, 1000.0f+i}, (cl_float3){0.0f, 1000.0f, 1000.0f+i}, (cl_float3){1000.0f, 1000.0f, 1000.0f+i}, (cl_float3){0.0f, 0.0f, -1.0f}));
        scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 1000.0f+i}, (cl_float3){1000.0f, 0.0f, 1000.0f+i}, (cl_float3){0.0f, 0.0f, 1000.0f+i}, (cl_float3){0.0f, 0.0f, -1.0f}));
    }
    scene.upload_Triangles();
    scene.generate_rays();
    scene.trace_rays();
    
//    for(int i=0;i<500;++i){
//        clock_t begin=clock();
//        scene.trace_rays();
//        clock_t end=clock();
//        float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
//        float fps=1.0f/elapsed_secs;
//        printf("\r                                                                                                                                       \r");
//        printf("%f fps",fps);
//        fflush(stdout);
//    }
    
    scene.finish();
    
    return 0;
}