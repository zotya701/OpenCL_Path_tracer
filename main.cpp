#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include "tiny_obj_loader.h"
#include <windows.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <random>
#include <vector>
#include <queue>
//#include <GL/gl.h>
//#include <GL/glu.h>
#include <GL/glut.h>
#include <CL/cl.hpp>

#define clear_line() printf("\r                                                                                                                             \r");

int WINID=0;
const int screen_width=192*8;
const int screen_height=108*8;
//const int screen_width=600;
//const int screen_height=400;
//const int screen_width=20;
//const int screen_height=20;
const int max_iterations=50;
int iterations=1;
int current_sample=0;
int old_sample=0;
float global_fov=75.0f;
float global_yaw=-13.800002-50;
float global_pitch=5.599997+10;
//float global_fov=60.0f;
//float global_yaw=0;
//float global_pitch=0;
float global_forward=0;
float global_rightward=0;
float global_upward=0;
cl_float3 global_shift=(cl_float3){265.055481, 162.305969, 360.414001};
//cl_float3 global_shift=(cl_float3){0,0,0};
enum ControllKeys {W, A, S, D, Q, Y, E, C, keys_num};
bool keys_down[keys_num];
bool real_time=true;

std::minstd_rand0 minstd_rand0;

cl_float3 rotate_z(cl_float3 v, float alpha){
    alpha=alpha/180.0f*3.141593f;
    float r[3];
    r[0]=v.s[0]*cos(alpha)-v.s[1]*sin(alpha);
    r[1]=v.s[0]*sin(alpha)+v.s[1]*cos(alpha);
    r[2]=v.s[2];
    return (cl_float3){r[0], r[1], r[2]};
}
cl_float3 rotate_y(cl_float3 v, float beta){
    beta=beta/180.0f*3.141593f;
    float r[3];
    r[0]=v.s[0]*cos(beta)+v.s[2]*sin(beta);
    r[1]=v.s[1];
    r[2]=-v.s[0]*sin(beta)+v.s[2]*cos(beta);
    return (cl_float3){r[0], r[1], r[2]};
}
cl_float3 rotate_x(cl_float3 v, float gamma){
    gamma=gamma/180.0f*3.141593f;
    float r[3];
    r[0]=v.s[0];
    r[1]=v.s[1]*cos(gamma)-v.s[2]*sin(gamma);
    r[2]=v.s[1]*sin(gamma)+v.s[2]*cos(gamma);
    return (cl_float3){r[0], r[1], r[2]};
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    std::vector<std::string> elems;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

cl_float3 str_to_cl_float3(std::string str){
    std::vector<std::string> arr=split(str, ' ');
    return (cl_float3){atof(arr[0].c_str()), atof(arr[1].c_str()), atof(arr[2].c_str())};
}

cl_int str_to_cl_int(std::string str){
    return atoi(str.c_str());
}

class Material{
private:
    cl_float3 kd,ks,emission,F0;    //diffuse, specular, emission, Fresnel
    cl_float n, shininess;
    cl_int type;   //0-diffuse, 1-x, 2-x, 3-Emitter
public:
    Material(){
        this->type=-1;
    }
    Material(cl_float3 kd, cl_float3 ks, cl_float3 emission, cl_float3 N, cl_float3 K, cl_float shininess, cl_int type){
        this->kd=kd; this->ks=ks; this->emission=emission; this->shininess=shininess; this->type=type;
        this->n=(cl_float){(N.s[0]+N.s[1]+N.s[2])/3.0f};
        float F0[3];
        for(int i=0;i<3;++i){
            float a=(N.s[i]-1)*(N.s[i]-1);
            float b=(N.s[i]+1)*(N.s[i]+1);
            F0[i]=(K.s[i]*K.s[i]+a)/(K.s[i]*K.s[i]+b);
        }
        this->F0=(cl_float3){F0[0], F0[1], F0[2]};
    }
};

class Ray{
private:
    cl_float3 P,D;  //origo and direction
};

class BBox{
private:
    cl_float3 bl,tr;
public:
    BBox(){
        bl=(cl_float3){0,0,0};
        tr=(cl_float3){0,0,0};
    }
    BBox(cl_float3 bl, cl_float3 tr){
        this->bl=bl;
        this->tr=tr;
    }
    void expand(BBox box){
        for(int i=0;i<3;++i){
            if (box.bl.s[i] < bl.s[i]) bl.s[i] = box.bl.s[i];
            if (box.tr.s[i] > tr.s[i]) tr.s[i] = box.tr.s[i];
        }
    }
};

class Triangle{
private:
    cl_float3 r1,r2,r3,N;   //vertices of the triangle and it's normal vector
    cl_ushort mati;
public:
    Triangle(cl_float3 r1, cl_float3 r2, cl_float3 r3, cl_ushort mati){
        this->r1=r1; this->r2=r2; this->r3=r3; this->mati=mati;
        float v1[3],v2[3],n[3];

        //calculate (r2-r1) and (r3-r1)
        for(int i=0;i<3;++i){
            v1[i]=r2.s[i]-r1.s[i];
            v2[i]=r3.s[i]-r1.s[i];
        }

        //cross product of (r2-r1) and (r3-r1)
        n[0]=v1[1]*v2[2] - v1[2]*v2[1];
        n[1]=v1[2]*v2[0] - v1[0]*v2[2];
        n[2]=v1[0]*v2[1] - v1[1]*v2[0];

        //normalize the normal vector
        float length=sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        for(int i=0;i<3;++i){
            n[i]=n[i]/length;
        }

        this->N=(cl_float3){n[0], n[1], n[2]};
    }
    BBox bbox(){
        cl_float3 bl, tr;
        for(int i=0;i<3;++i){
            bl.s[i]=std::min(std::min(r1.s[i], r2.s[i]),r3.s[i]);
            tr.s[i]=std::max(std::max(r1.s[i], r2.s[i]),r3.s[i]);
        }
        return BBox(bl, tr);
    }
    cl_float3 midpoint(){
        cl_float3 mp;
        for(int i=0;i<3;++i){
            mp.s[i]=(r1.s[i] + r2.s[i] + r3.s[i])/3.0f;
        }
        return mp;
    }
};

class Node{
private:
    cl_int2 trii;
    BBox bbox;
public:
    Node(cl_int2 trii, BBox bbox){
        this->trii=trii;
        this->bbox=bbox;
    }
};

class NodeOnHost{
private:
    NodeOnHost* left;
    NodeOnHost* right;
    BBox box;
    bool leaf;
    std::vector<Triangle> triangles;
public:
    NodeOnHost(){
        left=NULL;
        right=NULL;
        box=BBox();
        leaf=false;
        triangles=std::vector<Triangle>();
    };
    NodeOnHost* build(std::vector<Triangle> &tris, int depth){
        NodeOnHost *node=new NodeOnHost();
        if(tris.size()<=6){
            node->left=node->right=NULL;
            node->leaf=true;
            node->triangles=tris;
            node->box=tris[0].bbox();
            for(long i=1;i<tris.size();i++){
                node->box.expand(tris[i].bbox());
            }
            return node;
        }
        
        node->box=tris[0].bbox();
        cl_float3 midpoint=tris[0].midpoint();
        for(long i=1;i<tris.size();i++){
            node->box.expand(tris[i].bbox());
            cl_float3 mp=tris[i].midpoint();
            for(int j=0;j<3;++j){
                midpoint.s[j]=midpoint.s[j] + mp.s[j];
            }
        }
        for(int i=0;i<3;++i){
            midpoint.s[i]=midpoint.s[i]/tris.size();
        }
        
        int axis=depth%3;
        std::vector<Triangle> left_tris;
        std::vector<Triangle> right_tris;
        for(long i=0;i<tris.size();i++){
            cl_float3 mp=tris[i].midpoint();
            if(midpoint.s[axis] >= mp.s[axis])
                right_tris.push_back(tris[i]);
            else
                left_tris.push_back(tris[i]);
        }
        while(left_tris.size()==0 || right_tris.size()==0){
            left_tris.clear();
            right_tris.clear();
            axis=(axis+1)%3;
            for(long i=0;i<tris.size();i++){
                cl_float3 mp=tris[i].midpoint();
                if(midpoint.s[axis] >= mp.s[axis])
                    right_tris.push_back(tris[i]);
                else
                    left_tris.push_back(tris[i]);
            }
        }

        node->left=build(left_tris, depth+1);
        node->right=build(right_tris, depth+1);
        return node;
    }
    void convert(NodeOnHost* root, std::vector<Node>& kdarr, std::vector<Triangle>& neworder, int tri_shift, int kd_shift){
        if (root==NULL) return;
        
        struct tuple{
            NodeOnHost* node;
            int ptr;
            tuple(NodeOnHost* node, int i){
                this->node=node; this->ptr=i;
            }
        };
        
        std::queue<tuple> q;
        q.push(tuple(root, 1));
        int from=tri_shift;
        int to=tri_shift;
        while(q.empty()==false){
            tuple pair=q.front();
            if(pair.node->leaf){
                to=to+pair.node->triangles.size();
                while(kdarr.size()<=pair.ptr+kd_shift)
                        kdarr.push_back(Node((cl_int2){-1, -1}, pair.node->box));
                kdarr[pair.ptr+kd_shift]=Node((cl_int2){from, to}, pair.node->box);
                for(int i=0;i<pair.node->triangles.size();++i){
                    neworder.push_back(pair.node->triangles[i]);
                }
                from=from+pair.node->triangles.size();
            }else{
                while(kdarr.size()<=pair.ptr+kd_shift)
                        kdarr.push_back(Node((cl_int2){-1, -1}, pair.node->box));
                kdarr[pair.ptr+kd_shift]=Node((cl_int2){-1, -1}, pair.node->box);
            }
            q.pop();

            if(pair.node->left!= NULL){
                q.push(tuple(pair.node->left, pair.ptr*2));
            }
            if(pair.node->right != NULL){
                q.push(tuple(pair.node->right, pair.ptr*2+1));
            }
        }
    }
};

class Camera{
private:
    cl_float3 eye, lookat, up, right;
    cl_float XM, YM;
public:
    Camera(){
        float fov=global_fov;
        float yaw=global_yaw;
        float pitch=global_pitch;

        this->XM=(cl_float){(float)screen_width};
        this->YM=(cl_float){(float)screen_height};
        
        float up_length=this->YM/2.0f;
        float right_length=this->XM/2.0f;
        float ahead_length=right_length/tan(fov/2.0f/180.0f*3.141593f);

        cl_float3 up=(cl_float3){0.0f, 1.0f, 0.0f};
        cl_float3 right=(cl_float3){1.0f, 0.0f, 0.0f};
        cl_float3 ahead=(cl_float3){0.0f, 0.0f, 1.0f};

        up=rotate_x(up, pitch);
        up=rotate_y(up, yaw);
        right=rotate_x(right, pitch);
        right=rotate_y(right, yaw);
        ahead=rotate_x(ahead, pitch);
        ahead=rotate_y(ahead, yaw);

        for(int i=0;i<3;++i){
            global_shift.s[i]=global_shift.s[i] + ahead.s[i]*global_forward + right.s[i]*global_rightward + up.s[i]*global_upward;
        }
        for (int i=0;i<3;++i){
            up.s[i]=up.s[i]*up_length;
            right.s[i]=right.s[i]*right_length;
            ahead.s[i]=ahead.s[i]*ahead_length;
        }

        this->eye=(cl_float3){500.0f+global_shift.s[0], 500.0f+global_shift.s[1], -1299.037842f+global_shift.s[2]};
        this->up=(cl_float3){up.s[0], up.s[1], up.s[2]};
        this->right=(cl_float3){right.s[0], right.s[1], right.s[2]};
        this->lookat=(cl_float3){eye.s[0]+ahead.s[0], eye.s[1]+ahead.s[1], eye.s[2]+ahead.s[2]};
    }
};

class Color{
public:
    float r,g,b;
    Color(){
        r=g=b=0.0f;
    }
    Color(float r, float g, float b){
        this->r=r; this->g=g; this->b=b;
    }
};

//Color color_image[screen_width*screen_height];

class Scene{
private:
    Camera camera;
    std::vector<Triangle> tris;
    int tris_size;
    int tri_shift=0;
    std::vector<cl_int> kd_tree_shift;
    std::vector<Material> mats;
    std::vector<Node> kd_tree;
    int kd_tree_shift_size;
    int rays_size=screen_width*screen_height;
    
    cl::Context context;
    cl::Program program;
    cl::Buffer buffer_tris;
    cl::Buffer buffer_mats;
    cl::Buffer buffer_kd_tree;
    cl::Buffer buffer_kd_tree_shift;
    cl::Buffer buffer_rays;
    cl::Buffer buffer_rnds;
    cl::Buffer buffer_colors;
    
    cl::ImageGL imageFromGL;
    cl::ImageGL cl_screen;
    cl::CommandQueue queue;
public:
    void list_info(){
        int i, j;
        char* value;
        size_t valueSize;
        cl_uint platformCount;
        cl_platform_id* platforms;
        cl_uint deviceCount;
        cl_device_id* devices;
        cl_uint maxComputeUnits;

        // get all platforms
        clGetPlatformIDs(0, NULL, &platformCount);
        platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
        clGetPlatformIDs(platformCount, platforms, NULL);

        for (i = 0; i < platformCount; i++) {
            // get all devices
            clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
            devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
            clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);

            // for each device print critical attributes
            for (j = 0; j < deviceCount; j++) {
                // print device name
                clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
                value = (char*) malloc(valueSize);
                clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
                printf("%d. Device: %s\n", j+1, value);
                free(value);

                // print hardware device version
                clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
                value = (char*) malloc(valueSize);
                clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
                printf(" %d.%d Hardware version: %s\n", j+1, 1, value);
                free(value);

                // print software driver version
                clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
                value = (char*) malloc(valueSize);
                clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
                printf(" %d.%d Software version: %s\n", j+1, 2, value);
                free(value);

                // print c version supported by compiler for device
                clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
                value = (char*) malloc(valueSize);
                clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
                printf(" %d.%d OpenCL C version: %s\n", j+1, 3, value);
                free(value);

                // print parallel compute units
                clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                        sizeof(maxComputeUnits), &maxComputeUnits, NULL);
                printf(" %d.%d Parallel compute units: %d\n", j+1, 4, maxComputeUnits);
                
                // extensions
                clGetDeviceInfo(devices[j], CL_DEVICE_EXTENSIONS, 0, NULL, &valueSize);
                value = (char*) malloc(valueSize);
                clGetDeviceInfo(devices[j], CL_DEVICE_EXTENSIONS, valueSize, value, NULL);
                printf(" %d.%d OpenCL C extensions: %s\n", j+1, 3, value);
                free(value);
            }
            free(devices);
        }
        free(platforms);
    }
    void init_Scene(){
        list_info();
        
        //get all platforms (drivers)
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
//        if(platforms.size()==0){
//            std::cout<<" No platforms found. Check OpenCL installation!\n";
//            exit(1);
//        }
        cl::Platform platform=platforms[0];
//        std::cout << "Using platform: "<<platform.getInfo<CL_PLATFORM_NAME>()<<"\n";
        
        //get default device of the default platform
        std::vector<cl::Device> devices;
        platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
//        if(devices.size()==0){
//            std::cout<<" No devices found. Check OpenCL installation!\n";
//            exit(1);
//        }
        cl::Device device=devices[0];
//        std::cout<< "Using device: "<<device.getInfo<CL_DEVICE_NAME>()<<"\n";
        
        //read the source file
        std::ifstream stream("prog.cl");
        std::stringstream str_stream;
        str_stream << stream.rdbuf();
        std::string kernel_code=str_stream.str();
        
        cl::Program::Sources sources;
        sources.push_back({kernel_code.c_str(),kernel_code.length()});
        
        
        cl_context_properties properties[]={ 
            CL_CONTEXT_PLATFORM, (cl_context_properties)(platform)(), 
            CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(), 
            CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(), 0
        };
        
        //build the source file
        context=cl::Context({device}, properties);
        program=cl::Program(context,sources);
        if(program.build({device})!=CL_SUCCESS){
            std::string nfo=program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device);
            printf("%s\n\r",nfo.c_str());
            //std::cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(device)<<"\n";
            exit(1);
        }
        
        //create queue to which we will push commands for the device.
        queue=cl::CommandQueue(context,device);
        
        buffer_rays=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(Ray)*rays_size);
        buffer_rnds=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(cl_int)*rays_size);
        
        // create the texture which we will use in kernel
        glEnable(GL_TEXTURE_2D);
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screen_width, screen_height, 0, GL_RGBA, GL_FLOAT, NULL);
        //cl_screen=clCreateFromGLTexture(context(), CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, texture, NULL);
        cl_screen=cl::ImageGL(context, CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, texture, NULL);
        buffer_colors=cl::Buffer(context,CL_MEM_READ_WRITE ,sizeof(cl_float3)*rays_size);
        
        cl_int* RNDS=new cl_int[rays_size];
        for(int i=0;i<rays_size;++i){
            RNDS[i]=minstd_rand0();
        }
        queue.enqueueWriteBuffer(buffer_rnds,CL_TRUE,0,sizeof(cl_int)*rays_size,&RNDS[0]);
        delete[] RNDS;
    }
    void add_Triangle(Triangle tri){
        tris.push_back(tri);
    }
    int add_Material(Material mat){
        mats.push_back(mat);
        return (mats.size()-1);
    }
    void end_Obj(){
        if(kd_tree_shift.empty())
            kd_tree_shift.push_back((cl_int){0});
        else
            kd_tree_shift.push_back((cl_int){kd_tree.size()-1});
        NodeOnHost* root=new NodeOnHost();
        std::vector<Triangle> inorder;
        for(int i=tri_shift;i<tris.size();++i)
            inorder.push_back(tris[i]);
        root=root->build(inorder,0);
        inorder.clear();
        root->convert(root, kd_tree, inorder, tri_shift, kd_tree_shift[kd_tree_shift.size()-1]);
        for(int i=tri_shift;i<tris.size();++i)
            tris[i]=inorder[i-tri_shift];
        tri_shift=tris.size();
    }
    void add_Obj(std::string file, cl_float3 pos, cl_float3 scale, float pitch, float yaw){
        std::string matpath=file.substr(0,file.find_last_of("/")+1);
        tinyobj::attrib_t attrib;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;
        std::string err;
        bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, file.c_str(), matpath.c_str());
        if(!err.empty()) printf("Error msg: %s\n\r",err.c_str());
        if(!ret) exit(1);
        
        int mat_offset=this->mats.size();
        
        for(int i=0;i<materials.size();++i){
            cl_float3 diffuse_color=(cl_float3){materials[i].diffuse[0], materials[i].diffuse[1], materials[i].diffuse[2]};
            cl_float3 specular_color=(cl_float3){materials[i].specular[0], materials[i].specular[1], materials[i].specular[2]};
            cl_float3 emission=(cl_float3){materials[i].emission[0], materials[i].emission[1], materials[i].emission[2]};
            cl_float3 refractive_index=str_to_cl_float3(materials[i].unknown_parameter.at("Kn"));
            cl_float3 extinction_coefficient=str_to_cl_float3(materials[i].unknown_parameter.at("Kk"));
            cl_float shininess=materials[i].shininess;
            cl_int type=str_to_cl_int(materials[i].unknown_parameter.at("Tp"));
            this->add_Material(Material(diffuse_color, specular_color, emission, refractive_index, extinction_coefficient, shininess, type));
            printf("name = %s\n\r",materials[i].name.c_str());
            printf("  diffuse_color          = %4.2f %4.2f %4.2f\n\r", diffuse_color.s[0], diffuse_color.s[1], diffuse_color.s[2]);
            printf("  specular_color         = %4.2f %4.2f %4.2f\n\r", specular_color.s[0], specular_color.s[1], specular_color.s[2]);
            printf("  emission               = %4.2f %4.2f %4.2f\n\r", emission.s[0], emission.s[1], emission.s[2]);
            printf("  refractive_index       = %4.2f %4.2f %4.2f\n\r", refractive_index.s[0], refractive_index.s[1], refractive_index.s[2]);
            printf("  extinction_coefficient = %4.2f %4.2f %4.2f\n\r", extinction_coefficient.s[0], extinction_coefficient.s[1], extinction_coefficient.s[2]);
            printf("  shininess              = %4.2f\n\r", shininess);
            printf("  type                   = %d\n\r", type);
        }
//                                                           diffuse_color                  specular_color                     emission                      refractive_index              extinction_coefficient        shininess       type
        //int mat=this->.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){60.0f*2, 50.0f*2, 40.0f*2}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){  0}, (cl_int){3}));
        
        
        //loop over shapes
        for (size_t s = 0; s < shapes.size(); s++) {
            // Loop over faces(polygon)
            size_t index_offset = 0;
            for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++){
                int fv = shapes[s].mesh.num_face_vertices[f];

                cl_float3 vertices[fv];
                // Loop over vertices in the face.
                for (size_t v = 0; v < fv; v++) {
                    // access to vertex
                    tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
                    float vx = -attrib.vertices[3*idx.vertex_index+0];
                    float vy = attrib.vertices[3*idx.vertex_index+1];
                    float vz = attrib.vertices[3*idx.vertex_index+2];
                    vertices[v]=(cl_float3){vx, vy, vz};
                    vertices[v]=rotate_x(vertices[v], pitch);
                    vertices[v]=rotate_y(vertices[v], yaw);
                    for(int i=0;i<3;++i){
                        vertices[v].s[i]=vertices[v].s[i]*scale.s[i]+pos.s[i];
                    }
                }
                index_offset += fv;

                // per-face material
                shapes[s].mesh.material_ids[f];
                
                this->add_Triangle(Triangle(vertices[0], vertices[1], vertices[2], mat_offset+shapes[s].mesh.material_ids[f]));
            }
            this->end_Obj();
        }
    }
    void upload_Triangles(){
        buffer_kd_tree=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(Node)*kd_tree.size());
        queue.enqueueWriteBuffer(buffer_kd_tree,CL_TRUE,0,sizeof(Node)*kd_tree.size(),&kd_tree[0]);
        
        kd_tree_shift_size=kd_tree_shift.size();
        buffer_kd_tree_shift=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(cl_int)*kd_tree_shift_size);
        queue.enqueueWriteBuffer(buffer_kd_tree_shift,CL_TRUE,0,sizeof(cl_int)*kd_tree_shift_size,&kd_tree_shift[0]);
        
        tris_size=tris.size();
        buffer_tris=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(Triangle)*tris_size);
        queue.enqueueWriteBuffer(buffer_tris,CL_TRUE,0,sizeof(Triangle)*tris_size,&tris[0]);
        printf("%d %d\n\r",kd_tree.size(), tris.size());
    }
    void upload_Materials(){
        buffer_mats=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(Material)*mats.size());
        queue.enqueueWriteBuffer(buffer_mats,CL_TRUE,0,sizeof(Material)*mats.size(),&mats[0]);
    }
    void generate_rays(){
        camera=Camera();
        cl::Kernel kernel_gen_ray=cl::Kernel(program,"gen_ray");
        kernel_gen_ray.setArg(0,buffer_rays);
        kernel_gen_ray.setArg(1,camera);
        kernel_gen_ray.setArg(2,buffer_rnds);
        
        queue.enqueueNDRangeKernel(kernel_gen_ray,cl::NullRange,cl::NDRange(rays_size),cl::NullRange);
        //queue.finish();
    }
    void trace_rays(){
        //run the kernel
        cl::Kernel kernel_trace_ray=cl::Kernel(program,"trace_ray");
        kernel_trace_ray.setArg(0,cl_screen);
        kernel_trace_ray.setArg(1,buffer_tris);
        kernel_trace_ray.setArg(2,tris_size);
        kernel_trace_ray.setArg(3,buffer_mats);
        kernel_trace_ray.setArg(4,buffer_kd_tree);
        kernel_trace_ray.setArg(5,buffer_kd_tree_shift);
        kernel_trace_ray.setArg(6,kd_tree_shift_size);
        kernel_trace_ray.setArg(7,buffer_rays);
        kernel_trace_ray.setArg(8,buffer_rnds);
        kernel_trace_ray.setArg(9,iterations);
        kernel_trace_ray.setArg(10,current_sample);
        kernel_trace_ray.setArg(11,camera);
        kernel_trace_ray.setArg(12,buffer_colors);
        
        //queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(screen_width, screen_height),cl::NDRange(16, 16));
        //queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(screen_width, screen_height),cl::NullRange);
        
//        cl::Kernel kernel_filt_im=cl::Kernel(program,"filt_im");
//        kernel_filt_im.setArg(0,cl_screen);
//        kernel_filt_im.setArg(1,buffer_colors);
//        queue.enqueueNDRangeKernel(kernel_filt_im,cl::NullRange,cl::NDRange(screen_width, screen_height),cl::NullRange);
        
        if(real_time){
            glFinish();
            std::vector<cl::Memory> objs{cl_screen};
            queue.enqueueAcquireGLObjects(&objs, NULL, NULL);
            queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(screen_width, screen_height),cl::NullRange);
            queue.finish();
            queue.enqueueReleaseGLObjects(&objs, NULL, NULL);
        }else{
            queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(screen_width, screen_height),cl::NullRange);
            if(current_sample%3==0)
                queue.finish();
        }
    }
    void render(){
        this->generate_rays();
        this->trace_rays();
        current_sample++;
    }
//    void uniform_test(){
//        cl_float2* RNDS=new cl_float2[rays_size];
//        queue.enqueueReadBuffer(buffer_rnds,CL_TRUE,0,sizeof(cl_float2)*rays_size,RNDS);
//        int parts[10];
//        for(int i=0;i<10;++i)
//            parts[i]=0;
//        for(int i=0;i<rays_size;++i){
//            for(int j=0;j<2;++j){
//                if(RNDS[i].s[j]>0.9)
//                    parts[9]++;
//                else if(RNDS[i].s[j]>0.8)
//                    parts[8]++;
//                else if(RNDS[i].s[j]>0.7)
//                    parts[7]++;
//                else if(RNDS[i].s[j]>0.6)
//                    parts[6]++;
//                else if(RNDS[i].s[j]>0.5)
//                    parts[5]++;
//                else if(RNDS[i].s[j]>0.4)
//                    parts[4]++;
//                else if(RNDS[i].s[j]>0.3)
//                    parts[3]++;
//                else if(RNDS[i].s[j]>0.2)
//                    parts[2]++;
//                else if(RNDS[i].s[j]>0.1)
//                    parts[1]++;
//                else
//                    parts[0]++;
//            }
//        }
//        delete[] RNDS;
//        float sum=0;
//        for(int i=0;i<10;++i){
//            float percent=(float)parts[i]/(rays_size*2)*100;
//            sum+=percent;
//            printf("Quantity in range: 0.%d-%d.%d = %07.4f%\n\r", i, (i+1)/10, (i+1)%10, percent);
//        }
//        printf("%.4f\n\r",sum);
//    }
//    void download_image(){
//        //queue.enqueueReadBuffer(buffer_colors,CL_TRUE,0,sizeof(cl_float3)*rays_size,cl_float3_image);
//        
//        for(int i=0;i<rays_size;++i){
//            cl_float3 c=cl_float3_image[i];
////            cl_float3 c3=(cl_float3){fmax(0.0f, c.s[0]-0.004f), fmax(0.0f, c.s[1]-0.004f), fmax(0.0f, c.s[2]-0.004f)};
////            cl_float3 c2;
////            for(int i=0;i<3;++i){
////                c2.s[i]=(c3.s[i]*(c3.s[i]*6.2f+0.5f))/(c3.s[i]*(c3.s[i]*6.2f+1.7f)+0.06f);
////                c.s[i]=pow(c2.s[i], 2.2f);
////            }
//            color_image[i]=Color(c.s[0], c.s[1], c.s[2]);
//        }
//        glutPostRedisplay();
//    }
};

Scene scene;
void onInitialization( ) { 
    srand(time(0));
    glViewport(0, 0, screen_width, screen_height);

    scene.init_Scene();

    unsigned short LAMP, SUN, WHITE_DIFFUSE, RED_DIFFUSE, GREEN_DIFFUSE, PUPRLE_SPECULAR, BLACK_SPECULAR, CHROMIUM, GLASS, GOLD;
    //                                                           diffuse_color                  specular_color                     emission                      refractive_index              extinction_coefficient        shininess       type
    LAMP=0;             scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){60.0f*2, 50.0f*2, 40.0f*2}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){  0}, (cl_int){3}));
    SUN=1;              scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){60.0f*5, 50.0f*5, 40.0f*5}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){  0}, (cl_int){3}));
    WHITE_DIFFUSE=2;    scene.add_Material(Material((cl_float3){0.3f, 0.3f, 0.3f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){ 50}, (cl_int){0}));
    RED_DIFFUSE=3;      scene.add_Material(Material((cl_float3){0.3f, 0.1f, 0.1f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){ 50}, (cl_int){0}));
    GREEN_DIFFUSE=4;    scene.add_Material(Material((cl_float3){0.1f, 0.3f, 0.1f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){ 50}, (cl_int){0}));
    PUPRLE_SPECULAR=5;  scene.add_Material(Material((cl_float3){0.3f, 0.0f, 0.0f}, (cl_float3){0.3f, 0.3f, 0.3f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){200}, (cl_int){0}));
    BLACK_SPECULAR=6;   scene.add_Material(Material((cl_float3){0.05f, 0.05f, 0.05f}, (cl_float3){0.3f, 0.3f, 0.3f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){200}, (cl_int){0}));
    CHROMIUM=7;         scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){3.10f, 3.05f, 2.05f}, (cl_float3){3.3f, 3.3f, 2.9f}, (cl_float){  0}, (cl_int){1}));
    GOLD=8;             scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.17f, 0.35f, 1.50f}, (cl_float3){3.1f, 2.7f, 1.9f}, (cl_float){  0}, (cl_int){1}));
    GLASS=9;            scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){1.50f, 1.50f, 1.50f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){  0}, (cl_int){2}));
    
    //lámpa
//    scene.add_Triangle(Triangle((cl_float3){300.0f, 999.9f, 700.0f}, (cl_float3){300.0f, 999.9f, 300.0f}, (cl_float3){700.0f, 999.9f, 700.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){700.0f, 999.9f, 700.0f}, (cl_float3){300.0f, 999.9f, 300.0f}, (cl_float3){700.0f, 999.9f, 300.0f}, LAMP));
    
//    float asd=0;
//    for(int i=0;i<6;++i){
//        scene.add_Triangle(Triangle((cl_float3){-100.0f+asd, 700.9f, 600.0f}, (cl_float3){-100.0f+asd, 700.9f, 500.0f}, (cl_float3){0.0f+asd, 700.9f, 600.0f}, LAMP));
//        scene.add_Triangle(Triangle((cl_float3){0.0f+asd, 700.9f, 600.0f}, (cl_float3){-100.0f+asd, 700.9f, 500.0f}, (cl_float3){0.0f+asd, 700.9f, 500.0f}, LAMP));
//        asd=asd+220.00f;
//    }
    
//    scene.add_Triangle(Triangle((cl_float3){-50.0f, 999.9f, 750.0f}, (cl_float3){-50.0f, 999.9f, 950.0f}, (cl_float3){150.0f, 999.9f, 750.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){150.0f, 999.9f, 750.0f}, (cl_float3){-50.0f, 999.9f, 950.0f}, (cl_float3){150.0f, 999.9f, 950.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-50.0f+900, 999.9f, 750.0f}, (cl_float3){-50.0f+900, 999.9f, 950.0f}, (cl_float3){150.0f+900, 999.9f, 750.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){150.0f+900, 999.9f, 750.0f}, (cl_float3){-50.0f+900, 999.9f, 950.0f}, (cl_float3){150.0f+900, 999.9f, 950.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-50.0f, 999.9f, 750.0f-700}, (cl_float3){-50.0f, 999.9f, 950.0f-700}, (cl_float3){150.0f, 999.9f, 750.0f-700}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){150.0f, 999.9f, 750.0f-700}, (cl_float3){-50.0f, 999.9f, 950.0f-700}, (cl_float3){150.0f, 999.9f, 950.0f-700}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-50.0f+900, 999.9f, 750.0f-700}, (cl_float3){-50.0f+900, 999.9f, 950.0f-700}, (cl_float3){150.0f+900, 999.9f, 750.0f-700}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){150.0f+900, 999.9f, 750.0f-700}, (cl_float3){-50.0f+900, 999.9f, 950.0f-700}, (cl_float3){150.0f+900, 999.9f, 950.0f-700}, LAMP));
    
//    scene.add_Triangle(Triangle((cl_float3){-99.999f, 750.0f, 750.0f}, (cl_float3){-99.999f, 950.0f, 750.0f}, (cl_float3){-99.999f, 750.0f, 950.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-99.999f, 750.0f, 950.0f}, (cl_float3){-99.999f, 950.0f, 750.0f}, (cl_float3){-99.999f, 950.0f, 950.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-99.999f+1199.9f, 750.0f, 750.0f}, (cl_float3){-99.999f+1199.9f, 950.0f, 750.0f}, (cl_float3){-99.999f+1199.9f, 750.0f, 950.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-99.999f+1199.9f, 750.0f, 950.0f}, (cl_float3){-99.999f+1199.9f, 950.0f, 750.0f}, (cl_float3){-99.999f+1199.9f, 950.0f, 950.0f}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-99.999f, 750.0f, 750.0f-1000}, (cl_float3){-99.999f, 950.0f, 750.0f-1000}, (cl_float3){-99.999f, 750.0f, 950.0f-1000}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-99.999f, 750.0f, 950.0f-1000}, (cl_float3){-99.999f, 950.0f, 750.0f-1000}, (cl_float3){-99.999f, 950.0f, 950.0f-1000}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-99.999f+1199.9f, 750.0f, 750.0f-1000}, (cl_float3){-99.999f+1199.9f, 950.0f, 750.0f-1000}, (cl_float3){-99.999f+1199.9f, 750.0f, 950.0f-1000}, LAMP));
//    scene.add_Triangle(Triangle((cl_float3){-99.999f+1199.9f, 750.0f, 950.0f-1000}, (cl_float3){-99.999f+1199.9f, 950.0f, 750.0f-1000}, (cl_float3){-99.999f+1199.9f, 950.0f, 950.0f-1000}, LAMP));
    
    //elől
//    scene.add_Triangle(Triangle((cl_float3){-100.0f, 0.0f, 1000.0f}, (cl_float3){-100.0f, 1000.0f, 1000.0f}, (cl_float3){1100.0f, 1000.0f, 1000.0f}, WHITE_DIFFUSE));
//    scene.add_Triangle(Triangle((cl_float3){1100.0f, 1000.0f, 1000.0f}, (cl_float3){1100.0f, 0.0f, 1000.0f}, (cl_float3){-100.0f, 0.0f, 1000.0f}, WHITE_DIFFUSE));
//    
//    //balra
//    scene.add_Triangle(Triangle((cl_float3){-100.0f, 0.0f, 1000.0f}, (cl_float3){-100.0f, 0.0f, -1000.0f}, (cl_float3){-100.0f, 1000.0f, 1000.0f}, RED_DIFFUSE));
//    scene.add_Triangle(Triangle((cl_float3){-100.0f, 1000.0f, 1000.0f}, (cl_float3){-100.0f, 0.0f, -1000.0f}, (cl_float3){-100.0f, 1000.0f, -1000.0f}, RED_DIFFUSE));
//    
//    //jobbra
//    scene.add_Triangle(Triangle((cl_float3){1100.0f, 1000.0f, 1000.0f}, (cl_float3){1100.0f, 0.0f, -1000.0f}, (cl_float3){1100.0f, 0.0f, 1000.0f}, GREEN_DIFFUSE));
//    scene.add_Triangle(Triangle((cl_float3){1100.0f, 1000.0f, -1000.0f}, (cl_float3){1100.0f, 0.0f, -1000.0f}, (cl_float3){1100.0f, 1000.0f, 1000.0f}, GREEN_DIFFUSE));
//    
//    //felül
//    scene.add_Triangle(Triangle((cl_float3){-100.0f, 1000.0f, 1000.0f}, (cl_float3){-100.0f, 1000.0f, -1000.0f}, (cl_float3){1100.0f, 1000.0f, 1000.0f}, WHITE_DIFFUSE));
//    scene.add_Triangle(Triangle((cl_float3){1100.0f, 1000.0f, 1000.0f}, (cl_float3){-100.0f, 1000.0f, -1000.0f}, (cl_float3){1100.0f, 1000.0f, -1000.0f}, WHITE_DIFFUSE));
//    scene.end_Obj();
    //hátul
//    scene.add_Triangle(Triangle((cl_float3){-100.0f-1000, 0.0f, -1000.0f-1000}, (cl_float3){-100.0f-1000, 1000.0f, -1000.0f-1000}, (cl_float3){1100.0f+1000, 1000.0f, -1000.0f-1000}, WHITE_DIFFUSE));
//    scene.add_Triangle(Triangle((cl_float3){1100.0f+1000, 1000.0f, -1000.0f-1000}, (cl_float3){1100.0f+1000, 0.0f, -1000.0f-1000}, (cl_float3){-100.0f-1000, 0.0f, -1000.0f-1000}, WHITE_DIFFUSE));
    
    //alul
    scene.add_Triangle(Triangle((cl_float3){-10000.0f, 0.0f, -10000.0f}, (cl_float3){-10000.0f, 0.0f, 10000.0f}, (cl_float3){10000.0f, 0.0f, 10000.0f}, WHITE_DIFFUSE));
    scene.add_Triangle(Triangle((cl_float3){10000.0f, 0.0f, 10000.0f}, (cl_float3){10000.0f, 0.0f, -10000.0f}, (cl_float3){-10000.0f, 0.0f, -10000.0f}, WHITE_DIFFUSE));
    scene.end_Obj();
    
    cl_float3 sun_pos[6];
    sun_pos[0]=(cl_float3){-20000.0f, 300000.0f, 20000.0f}; sun_pos[1]=(cl_float3){-20000.0f, 300000.0f, -20000.0f}; sun_pos[2]=(cl_float3){20000.0f, 300000.0f, 20000.0f};
    sun_pos[3]=(cl_float3){20000.0f, 300000.0f, 20000.0f}; sun_pos[4]=(cl_float3){-20000.0f, 300000.0f, -20000.0f}; sun_pos[5]=(cl_float3){20000.0f, 300000.0f, -20000.0f};
    for(int i=0;i<6;++i){
        sun_pos[i]=rotate_x(sun_pos[i], 45);
        sun_pos[i]=rotate_y(sun_pos[i], 150);
    }
//    scene.add_Triangle(Triangle(sun_pos[0], sun_pos[1], sun_pos[2], SUN));
//    scene.add_Triangle(Triangle(sun_pos[3], sun_pos[4], sun_pos[5], SUN));
//    scene.end_Obj();

    //arany cucc
//    cl_float3 move=(cl_float3){750.0f, 0.001f, 300.0f};
    //cl_float3 move=(cl_float3){950.0f, 0.001f, 850.0f-500};
//    cl_float3 move=(cl_float3){400.0f, 100.001f, -2400.0f};
    //cl_float3 scale=(cl_float3){1.5f, 0.28284f*0.28284f/40*500, 1.5f};
//    cl_float3 scale=(cl_float3){6.0f, 0.28284f*0.28284f/40*800, 6.0f};
    //felül
    
    cl_float3 move=(cl_float3){50.0f+350-150, 0.000f, 850-350-1300};
    cl_float3 scale=(cl_float3){3.0f, 0.28284f*0.28284f/40*400, 3.0f};
//    cl_float3 move=(cl_float3){-50.0f, 0.001f, 950};
//    cl_float3 scale=(cl_float3){0.0125f, 0.28284f*0.28284f/40*1.25, 0.0125f};
    for(int szar3=0;szar3<1*0;++szar3){
        for(int szar2=0;szar2<1*1;++szar2){
            for(int szar=0;szar<1*1;++szar){
//                move=(cl_float3){-50.0f+szar*100, 100.001f+szar3*20, 950-szar2*100};
//                move=(cl_float3){400+szar*100/20.0, 100.001f+szar3*5, 600-szar2*100/20.0};
                scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GOLD));
                scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GOLD));
                //alul
                scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GOLD));
                scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GOLD));
                //jobb alul
                scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, GOLD));
                scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, GOLD));
                //jobb felül
                scene.add_Triangle(Triangle((cl_float3){100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GOLD));
                scene.add_Triangle(Triangle((cl_float3){100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, GOLD));
                //bal alul
                scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, GOLD));
                scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, GOLD));
                //bal felül
                scene.add_Triangle(Triangle((cl_float3){-100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GOLD));
                scene.add_Triangle(Triangle((cl_float3){-100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, GOLD));
            }
            //scene.start_New_Obj();
        }
        //scene.start_New_Obj();
    }
//    clock_t begin = clock();
//    scene.end_Obj();
//    clock_t end = clock();
//    double elapsed_secs = double(end - begin) / (CLOCKS_PER_SEC/1000);
//    printf("time=%fms\n\r",elapsed_secs);
    
//    move=(cl_float3){400.0f, 200.0f, 600.0f};
//    scale=(cl_float3){1.5f/200.0f*400.0f, 0.28284f*0.28284f/40*400, 1.5f/200.0f*400.0f};
//    move=(cl_float3){400.0f-200, 200.0f, 600.0f-100};
//    scale=(cl_float3){1.5f/200.0f*200.0f, 0.28284f*0.28284f/40*200, 1.5f/200.0f*200.0f};
//    scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GLASS));
//    //alul
//    scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GLASS));
//    //jobb alul
//    scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, GLASS));
//    //jobb felül
//    scene.add_Triangle(Triangle((cl_float3){100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, (cl_float3){100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, GLASS));
//    //bal alul
//    scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], -100.0f*scale.s[2]+move.s[2]}, GLASS));
//    //bal felül
//    scene.add_Triangle(Triangle((cl_float3){-100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){-100.0f*scale.s[0]+move.s[0], 0.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, (cl_float3){0.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 100.0f*scale.s[2]+move.s[2]}, (cl_float3){-100.0f*scale.s[0]+move.s[0], 500.0f*scale.s[1]+move.s[1], 0.0f*scale.s[2]+move.s[2]}, GLASS));
//    scene.start_New_Obj();
//    
//    scene.add_Triangle(Triangle((cl_float3){800.0f, 200.0f, 300.0f}, (cl_float3){800.0f, 200.0f, 700.0f}, (cl_float3){1100.0f, 600.0f, 700.0f}, GOLD));
//    scene.add_Triangle(Triangle((cl_float3){1100.0f, 600.0f, 300.0f}, (cl_float3){1100.0f, 600.0f, 700.0f}, (cl_float3){800.0f, 200.0f, 300.0f}, GOLD));
//    scene.start_New_Obj();
//    scene.add_Triangle(Triangle((cl_float3){0.0f, 800.9f, 800.0f}, (cl_float3){0.0f, 800.9f, 200.0f}, (cl_float3){600.0f, 800.9f, 800.0f}, WHITE_DIFFUSE));
//    scene.add_Triangle(Triangle((cl_float3){600.0f, 800.9f, 800.0f}, (cl_float3){0.0f, 800.9f, 200.0f}, (cl_float3){600.0f, 800.9f, 200.0f}, WHITE_DIFFUSE));
//    scene.start_New_Obj();

    //üveg hasáb
    float movx=200+300;
    float movy=0;
    float movz=-150-300-1200;
//    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 500.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){400.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, GLASS));
//    
//    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){0.0f+movx, 0.0f+movy, 600.0f+movz}, GLASS));
//    
//    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
//    
//    scene.add_Triangle(Triangle((cl_float3){400.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){400.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
//    
//    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
//    
//    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 600.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 600.0f+movz}, GLASS));
//    scene.end_Obj();
    
    
//    movx=50;
//    movy=750;
//    movz=850;
    movx=400;
    movy=400;
    movz=600;
//    GLASS=BLACK_DIFFUSE;
//    scene.add_Triangle(Triangle((cl_float3){-100.0f+movx, 0.0f+movy, -100.0f+movz}, (cl_float3){-100.0f+movx, 0.0f+movy, 100.0f+movz}, (cl_float3){0.0f+movx, 150.0f+movy, 0.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){100.0f+movx, 0.0f+movy, -100.0f+movz}, (cl_float3){100.0f+movx, 0.0f+movy, 100.0f+movz}, (cl_float3){0.0f+movx, 150.0f+movy, 0.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){-100.0f+movx, 0.0f+movy, -100.0f+movz}, (cl_float3){100.0f+movx, 0.0f+movy, -100.0f+movz}, (cl_float3){0.0f+movx, 150.0f+movy, 0.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){-100.0f+movx, 0.0f+movy, 100.0f+movz}, (cl_float3){100.0f+movx, 0.0f+movy, 100.0f+movz}, (cl_float3){0.0f+movx, 150.0f+movy, 0.0f+movz}, GLASS));
//
//    scene.add_Triangle(Triangle((cl_float3){-100.0f+movx, 0.0f+movy, -100.0f+movz}, (cl_float3){-100.0f+movx, 0.0f+movy, 100.0f+movz}, (cl_float3){0.0f+movx, -150.0f+movy, 0.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){100.0f+movx, 0.0f+movy, -100.0f+movz}, (cl_float3){100.0f+movx, 0.0f+movy, 100.0f+movz}, (cl_float3){0.0f+movx, -150.0f+movy, 0.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){-100.0f+movx, 0.0f+movy, -100.0f+movz}, (cl_float3){100.0f+movx, 0.0f+movy, -100.0f+movz}, (cl_float3){0.0f+movx, -150.0f+movy, 0.0f+movz}, GLASS));
//    scene.add_Triangle(Triangle((cl_float3){-100.0f+movx, 0.0f+movy, 100.0f+movz}, (cl_float3){100.0f+movx, 0.0f+movy, 100.0f+movz}, (cl_float3){0.0f+movx, -150.0f+movy, 0.0f+movz}, GLASS));
//    scene.start_New_Obj();
    
//    movx=0;
//    movy=300;
//    movz=500;
//    for(int i=0;i<6;++i){
//        scene.add_Triangle(Triangle((cl_float3){-100.0f+movx, 0.0f+movy, 0.0f+movz}, (cl_float3){100.0f+movx, 0.0f+movy, 0.0f+movz}, (cl_float3){100.0f+movx, 200.0f+movy, 200.0f+movz}, 7));
//        scene.add_Triangle(Triangle((cl_float3){100.0f+movx, 200.0f+movy, 200.0f+movz}, (cl_float3){-100.0f+movx, 0.0f+movy, 0.0f+movz}, (cl_float3){-100.0f+movx, 200.0f+movy, 200.0f+movz}, 7));
//        movx=movx+200;
//    }
    movx=700;
    movy=300;
    movz=200;
    cl_float3 movs=(cl_float3){700, 500, 100};
    cl_float3 v[4];
    v[0]=(cl_float3){-100,-100,0};
    v[1]=(cl_float3){100,-100,0};
    v[2]=(cl_float3){100,100,0};
    v[3]=(cl_float3){-100,100,0};
    for(int i=0;i<4;++i){
        v[i]=rotate_x(v[i], 30);
        v[i]=rotate_y(v[i], 120);
        for(int j=0;j<3;++j)
            v[i].s[j]=v[i].s[j]+movs.s[j];
    }
//    scene.add_Triangle(Triangle(v[0], v[1], v[2], 7));
//    scene.add_Triangle(Triangle(v[0], v[3], v[2], 7));
//    scene.end_Obj();
    
    for(int i=0;i<12;++i){
        cl_float3 asd=(cl_float3){120,90,-900};
        asd=rotate_y(asd, i*(360/12));
//        scene.add_Obj("../models/sphere.obj", asd, (cl_float3){200,200,200}, 0, 0, GOLD);
    }
    
    
    
//    scene.add_Obj("../models/lsphere.obj", (cl_float3){-1000,1000,-50}, (cl_float3){200,200,200}, 0, 0, LAMP);
//    scene.add_Obj("../models/lsphere.obj", (cl_float3){-1000+200,1000-250,-50}, (cl_float3){200,200,200}, 0, 0, LAMP);
//    scene.add_Obj("../models/lsphere.obj", (cl_float3){-400,1000-400,-400}, (cl_float3){200,200,200}, 0, 0, LAMP);
//    scene.add_Obj("../models/lsphere.obj", (cl_float3){-978,965+80,-1120}, (cl_float3){80,80,80}, 0, 0, SUN);
//    scene.add_Obj("../models/Studiolight.obj", (cl_float3){-2800,0+80,-200}, (cl_float3){70,70,70}, 0, 180, BLACK_SPECULAR);
//    scene.add_Obj("../models/templeobj.obj", (cl_float3){0,0,0}, (cl_float3){20,20,20}, 0, 180, RED_DIFFUSE);
//    scene.add_Obj("../models/dragon.obj", (cl_float3){-1000+600,0,0}, (cl_float3){60,60,60}, 0, -20+20+180, GLASS);
//    scene.add_Obj("../models/dragon.obj", (cl_float3){-1000,0,0}, (cl_float3){60,60,60}, 0, -20+180, PUPRLE_SPECULAR);
//    scene.add_Obj("../models/sphere.obj", (cl_float3){-1000,230,-700}, (cl_float3){500,500,500}, 0, 0, GOLD);
//    scene.add_Obj("../models/Wineglass.obj", (cl_float3){-100,0,-600}, (cl_float3){35,35,35}, 0, 0, GLASS);
//    scene.add_Obj("../models/spherical_rind_shape_wide_strips.obj", (cl_float3){0,130,-50}, (cl_float3){10,10,10}, 0, 0, PUPRLE_SPECULAR);
//    scene.add_Obj("../models/sphere.obj", (cl_float3){-400,90,-600}, (cl_float3){200,200,200}, 0, 0, GLASS);
//    scene.add_Obj("../models/sphere.obj", (cl_float3){-400,90,-900}, (cl_float3){200,200,200}, 0, 0, GLASS);
//    scene.add_Obj("../models/sphere.obj", (cl_float3){-700,250,700}, (cl_float3){500,500,500}, 0, 0, GOLD);
//    scene.add_Obj("../models/FordGT402.obj", (cl_float3){-900,0,-900}, (cl_float3){10,10,10}, -90, 180, RED_DIFFUSE);
//    scene.add_Obj("../models/Wineglass.obj", (cl_float3){-100,0,-900}, (cl_float3){35,35,35}, 0, 0, GLASS);
//    scene.add_Obj("../models/champagne_bottle.obj", (cl_float3){-600,0,-900}, (cl_float3){40,40,40}, 0, 0, GLASS);
//    scene.add_Obj("../models/CocaCola_Glass.obj", (cl_float3){200,0,-900}, (cl_float3){0.5,0.5,0.5}, 0, 180, GLASS);
//    scene.add_Obj("../models/Engagement Ring.obj", (cl_float3){300,0,-500}, (cl_float3){3,3,3}, 0, 0, GOLD);
//    scene.add_Obj("../models/egg.obj", (cl_float3){300,0,-500}, (cl_float3){6,6,6}, 0, 0, PUPRLE_SPECULAR);
    
    scene.add_Obj("../models/lsphere.obj", (cl_float3){0,1000,-50}, (cl_float3){200,200,200}, 0, 0);
//    scene.add_Obj("../models/lamp.obj", (cl_float3){-830,0,-700}, (cl_float3){10,10,10}, 0, 0);
//    scene.add_Obj("../models/lamp2.obj", (cl_float3){-200,330,-480}, (cl_float3){4,4,4}, 0, 0);
    scene.add_Obj("../models/chair.obj", (cl_float3){50,0,-150}, (cl_float3){190,190,190}, 0, 0);
    scene.add_Obj("../models/egg.obj", (cl_float3){-350,330,-400}, (cl_float3){0.5,0.5,0.5}, 0, 0);
    scene.add_Obj("../models/dragon.obj", (cl_float3){-670,330,-410}, (cl_float3){10,10,10}, 0, -20+70);
    scene.add_Obj("../models/Wineglass.obj", (cl_float3){-300,330,-270}, (cl_float3){1,1,1}, 0, 0);
    scene.add_Obj("../models/sphere.obj", (cl_float3){-490,377,-400}, (cl_float3){100,100,100}, 0, 0);
    scene.add_Obj("../models/glass-table.obj", (cl_float3){-200,0,-200}, (cl_float3){500,500,500}, 0, 0);
//    scene.add_Obj("../models/building2.obj", (cl_float3){-300,0,-500}, (cl_float3){50,50,50}, 0, 0);
//    scene.add_Obj("../models/FordGT402.obj", (cl_float3){-450,330,-270}, (cl_float3){1.2,1.2,1.2}, -90, 0);
//    scene.add_Obj("../models/sphere.obj", (cl_float3){0,200,0}, (cl_float3){500,500,500}, 0, 0);
    
    scene.upload_Triangles();
    scene.upload_Materials();
}

void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    float h=((screen_height+100.0f)/2.0f-100.0f)/((screen_height+100.0f)/2.0f);
    h=1;
    glBindTexture(GL_TEXTURE_2D, 1);
    glBegin(GL_QUADS);
    glTexCoord2i(0, 1);
    glVertex3f(-1.0f, h, 0.0f);
    glTexCoord2i(1, 1);
    glVertex3f( 1.0f, h, 0.0f);
    glTexCoord2i(1, 0);
    glVertex3f( 1.0f,-1.0f, 0.0f);
    glTexCoord2i(0, 0);
    glVertex3f(-1.0f,-1.0f, 0.0f);
    glEnd();
    glFinish();
    
    glutSwapBuffers();
}

bool space=true;
void onKeyboard(unsigned char key, int x, int y) {
    if(key=='-'){
        if(iterations>1){
            iterations--;
            current_sample=0;
        }
    }
    if(key=='+'){
        if(iterations<max_iterations){
            iterations++;
            current_sample=0;
        }
    }
    if(key==27){//esc key
        glutDestroyWindow(WINID);
        exit(0);
    }
    if(key==' '){
        if(space){
            glutFullScreen();
        }else{
            glutReshapeWindow(screen_width,screen_height);
        }
        space=!space;
    }
    if(key=='r'){
        real_time=!real_time;
    }
    switch(key) {
        case 'w': case 'W':
            keys_down[W] = true;
        break;
        case 's': case 'S':
            keys_down[S] = true;
        break;
        case 'a': case 'A':
            keys_down[A] = true;
        break;
        case 'd': case 'D':
            keys_down[D] = true;
        break;
        case 'q': case 'Q':
            keys_down[Q] = true;
        break;
        case 'y': case 'Y':
            keys_down[Y] = true;
        break;
        case 'e': case 'E':
            keys_down[E] = true;
        break;
        case 'c': case 'C':
            keys_down[C] = true;
        break;
    }
}

void onKeyboardUp(unsigned char key, int x, int y) {
    switch(key) {
        case 'w': case 'W':
            keys_down[W] = false;
            current_sample=0;
        break;
        case 's': case 'S':
            keys_down[S] = false;
            current_sample=0;
        break;
        case 'a': case 'A':
            keys_down[A] = false;
            current_sample=0;
        break;
        case 'd': case 'D':
            keys_down[D] = false;
            current_sample=0;
        break;
        case 'q': case 'Q':
            keys_down[Q] = false;
            current_sample=0;
        break;
        case 'y': case 'Y':
            keys_down[Y] = false;
            current_sample=0;
        break;
        case 'e': case 'E':
            keys_down[E] = false;
            current_sample=0;
        break;
        case 'c': case 'C':
            keys_down[C] = false;
            current_sample=0;
        break;
    }
}

int last_x, last_y;
bool mouse_down=false;
void onMouse(int button, int state, int x, int y) {
    last_x = x;
    last_y = y;
    if ((button == GLUT_LEFT_BUTTON ) && (state == GLUT_DOWN)){
        mouse_down=true;
        current_sample=0;
    }

    if ((button == GLUT_LEFT_BUTTON ) && (state == GLUT_UP)){
        mouse_down=false;
        current_sample=0;
    }
}
 
void onMouseMotion(int x, int y) {
    int dx=x-last_x;
    int dy=y-last_y;
    float speed=0.2f;
    if(global_fov<10)
        speed=0.05f;
    if(global_fov<2)
        speed=0.01f;
    global_yaw=global_yaw+dx*speed;
    global_pitch=global_pitch+dy*speed;
    last_x = x;
    last_y = y;
}

float old=0.0f;
float newTime=0.0f;
float dt=0.0f;
float start=0.0f;
bool reset_timer=true;
clock_t begin,end;
void onIdle( ) {
    if(reset_timer){
        begin=clock();
        old_sample=current_sample;
        reset_timer=false;
    }
    
    int before=iterations;
    if(keys_down[W] || keys_down[A] || keys_down[S] || keys_down[D] || keys_down[Q] || keys_down[Y] || keys_down[E] || keys_down[C] || mouse_down){
        //iterations=1;
        current_sample=0;
        start=glutGet(GLUT_ELAPSED_TIME)/1000.0f;
    }
    
    old = newTime;
    newTime = glutGet(GLUT_ELAPSED_TIME)/1000.0f;
    dt=newTime-old;
    
    float speed=1000.0f;
    if(keys_down[W])
        global_forward=speed*dt;
    else if(keys_down[S])
        global_forward=-speed*dt;
    else
        global_forward=0;
    
    if(keys_down[A])
        global_rightward=-speed*dt;
    else if(keys_down[D])
        global_rightward=speed*dt;
    else
        global_rightward=0;
    
    if(keys_down[Q])
        global_upward=speed*dt;
    else if(keys_down[Y])
        global_upward=-speed*dt;
    else
        global_upward=0;
    
    if(keys_down[E])
        if(global_fov>10)
            global_fov=global_fov-20*dt;
        else if(global_fov>0.1f)
            global_fov=global_fov-2*dt;
        else
            global_fov=0.1f;
    else if(keys_down[C])
        if(global_fov<10)
            global_fov=global_fov+2*dt;
        else if(global_fov<90)
            global_fov=global_fov+20*dt;
        else
            global_fov=90.0f;
    
    scene.render();
    if(real_time){
        glutPostRedisplay();
    }
    end=clock();
    double elapsed_secs=double(end-begin)/CLOCKS_PER_SEC;
    
    if(elapsed_secs>1.0){
        reset_timer=true;
        clear_line();
        printf("Samples=%010d  Samples/sec=%08.3f Render time=%08.3fms  real_time=%d  Iterations=%02d  Elapsed seconds=%f", current_sample, (current_sample-old_sample)/elapsed_secs, elapsed_secs/(current_sample-old_sample)*1000.0f, real_time, iterations, glutGet(GLUT_ELAPSED_TIME)/1000.0f-start);
    }
    
    fflush(stdout);
    iterations=before;
}

int main(int argc, char **argv) {
    glutInit(&argc, argv); 				// GLUT inicializalasa
    glutInitWindowSize(screen_width, screen_height);	// Alkalmazas ablak kezdeti merete 600x600 pixel 
    glutInitWindowPosition(100, 100);			// Az elozo alkalmazas ablakhoz kepest hol tunik fel
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);	// 8 bites R,G,B,A + dupla buffer + melyseg buffer

    WINID=glutCreateWindow("Path tracer");                    // Alkalmazas ablak megszuletik es megjelenik a kepernyon

    glMatrixMode(GL_MODELVIEW);				// A MODELVIEW transzformaciot egysegmatrixra inicializaljuk
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);			// A PROJECTION transzformaciot egysegmatrixra inicializaljuk
    glLoadIdentity();

    onInitialization();					// Az altalad irt inicializalast lefuttatjuk

    glutDisplayFunc(onDisplay);				// Esemenykezelok regisztralasa
    glutMouseFunc(onMouse); 
    glutIdleFunc(onIdle);
    glutKeyboardFunc(onKeyboard);
    glutKeyboardUpFunc(onKeyboardUp);
    glutMotionFunc(onMouseMotion);
    
    glutMainLoop();					// Esemenykezelo hurok
    
    return 0;
}