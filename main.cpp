#include <windows.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <random>
#include <vector>
#include <queue>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <CL/cl.hpp>

#define clear_line() printf("\r                                                                                                                             \r");

int WINID=0;
const int screen_width=192*3;
const int screen_height=108*3;
//const int screen_width=600;
//const int screen_height=600;
//const int screen_width=20;
//const int screen_height=20;
const int max_iterations=15;
int iterations=1;
int current_sample=0;
int old_sample=0;
float global_yaw=-13.800002;
float global_pitch=5.599997;
float global_forward=0;
float global_rightward=0;
cl_float3 global_shift=(cl_float3){265.055481, 162.305969, 360.414001};
enum ControllKeys {W, A, S, D, keys_num};
bool keys_down[keys_num];
bool real_time=true;

//std::minstd_rand0 generator;
//std::uniform_int_distribution<int> distribution(1,2147483646);
//#define RND distribution(generator);

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
    cl_int axis;
    cl_float split;
    BBox bbox;
public:
    Node(cl_int2 trii, cl_int axis, cl_float split, BBox bbox){
        this->trii=trii;
        this->split=split;
        this->bbox=bbox;
    }
};

class KDNode{
private:
    BBox box;
    KDNode* left;
    KDNode* right;
    std::vector<Triangle> triangles;
    bool leaf;
    int axis;
    float split;
public:
    KDNode(){
        box=BBox();
        left=NULL;
        right=NULL;
        triangles=std::vector<Triangle>();
        leaf=false;
        split=0;
    };
    KDNode* build(std::vector<Triangle> &tris, int depth){
        KDNode *node=new KDNode();
        node->leaf=false;
        node->triangles=std::vector<Triangle>();
        node->left=NULL;
        node->right=NULL;
        node->box=BBox();

        if(tris.size()==0) return node;
        
        if(tris.size()<=3){
            node->triangles=tris;
            node->leaf=true;
            node->box=tris[0].bbox();
            for(long i=1;i<tris.size();i++){
                node->box.expand(tris[i].bbox());
            }
            node->left=NULL;
            node->right=NULL;
            //node->left->triangles=std::vector<Triangle>();
            //node->right->triangles=std::vector<Triangle>();
            return node;
        }

        node->box=tris[0].bbox();
        cl_float3 midpt=tris[0].midpoint();

        for(long i=1;i<tris.size();i++){
            node->box.expand(tris[i].bbox());
            cl_float3 mp=tris[i].midpoint();
            for(int i=0;i<3;++i){
                midpt.s[i]=midpt.s[i] + mp.s[i];
            }
        }
        
        for(int i=0;i<3;++i){
            midpt.s[i]=midpt.s[i]/tris.size();
        }

        std::vector<Triangle> left_tris;
        std::vector<Triangle> right_tris;
        int axis=depth%3;

        node->split=midpt.s[axis];
        node->axis=axis;
        for(long i=0;i<tris.size();i++){
            cl_float3 mp=tris[i].midpoint();
            //midpt.s[axis] >= mp.s[axis] ? right_tris.push_back(tris[i]) : left_tris.push_back(tris[i]);
            if(midpt.s[axis] >= mp.s[axis])
                right_tris.push_back(tris[i]);
            else
                left_tris.push_back(tris[i]);
        }

        if(tris.size()==left_tris.size() || tris.size()==right_tris.size()){
            node->triangles=tris;
            node->leaf=true;
            node->box=tris[0].bbox();

            for(long i=1;i<tris.size();i++){
                node->box.expand(tris[i].bbox());
            }

            node->left=NULL;
            node->right=NULL;
            //node->left->triangles=std::vector<Triangle>();
            //node->right->triangles=std::vector<Triangle>();

            return node;
        }
        
        if(!node->leaf){
            node->left=build(left_tris, depth+1);
            node->right=build(right_tris, depth+1);
        }
        return node;
    }
    void convert(KDNode* root, std::vector<Node>& kdarr, std::vector<Triangle>& neworder, int tri_shift, int kd_shift){
        if (root==NULL){
            return;
        }
        
        struct tuple{
            KDNode* node;
            int i;
            tuple(KDNode* node, int i){
                this->node=node; this->i=i;
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
                while(kdarr.size()<=pair.i+kd_shift)
                        kdarr.push_back(Node((cl_int2){-1, -1}, pair.node->axis, pair.node->split, pair.node->box));
                kdarr[pair.i+kd_shift]=Node((cl_int2){from, to}, pair.node->axis, pair.node->split, pair.node->box);
                for(int i=0;i<pair.node->triangles.size();++i){
                    neworder.push_back(pair.node->triangles[i]);
                }
                from=from+pair.node->triangles.size();
            }else{
                while(kdarr.size()<=pair.i+kd_shift)
                        kdarr.push_back(Node((cl_int2){-1, -1}, pair.node->axis, pair.node->split, pair.node->box));
                kdarr[pair.i+kd_shift]=Node((cl_int2){-1, -1}, pair.node->axis, pair.node->split, pair.node->box);
            }
            q.pop();

            if(pair.node->left!= NULL){
                q.push(tuple(pair.node->left, pair.i*2));
            }
            if(pair.node->right != NULL){
                q.push(tuple(pair.node->right, pair.i*2+1));
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
        float fov=75;
        float yaw=0.0f+global_yaw;
        float pitch=0.0f+global_pitch;
        float roll=0.0f;

        float up_length=1.0f;
        float right_length=1.0f*((float)screen_width/screen_height);
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
            global_shift.s[i]=global_shift.s[i] + ahead.s[i]*global_forward + right.s[i]*global_rightward;
        }

        eye=(cl_float3){500.0f+global_shift.s[0], 500.0f+global_shift.s[1], -1299.037842f+global_shift.s[2]};

        this->up=(cl_float3){up.s[0]*up_length, up.s[1]*up_length, up.s[2]*up_length};
        this->right=(cl_float3){right.s[0]*right_length, right.s[1]*right_length, right.s[2]*right_length};
        lookat=(cl_float3){eye.s[0]+ahead.s[0]*ahead_length, eye.s[1]+ahead.s[1]*ahead_length, eye.s[2]+ahead.s[2]*ahead_length};

        XM=(cl_float){(float)screen_width};
        YM=(cl_float){(float)screen_height};
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
    int tri_shift;
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
    cl::Image2D cl_screen;
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
            }
            free(devices);
        }
        free(platforms);
    }
    void init_Scene(){
        list_info();
        
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
        
        
        cl_context_properties properties[] = 
        { 
          CL_CONTEXT_PLATFORM, (cl_context_properties)(default_platform)(), 
          CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(), 
          CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(), 0
        };
        
        //build the source file
        context=cl::Context({default_device}, properties);
        program=cl::Program(context,sources);
        if(program.build({default_device})!=CL_SUCCESS){
            std::cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";
            exit(1);
        }
        
        //create queue to which we will push commands for the device.
        queue=cl::CommandQueue(context,default_device);
        
        buffer_rays=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(Ray)*rays_size);
        buffer_rnds=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(cl_int)*rays_size);
        
        // create the texture which we will use in kernel
        glEnable(GL_TEXTURE_2D);
        GLuint texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, screen_width, screen_height, 0, GL_RGBA, GL_FLOAT, NULL);
        cl_screen=clCreateFromGLTexture(context(), CL_MEM_WRITE_ONLY, GL_TEXTURE_2D, 0, texture, NULL);
        buffer_colors=cl::Buffer(context,CL_MEM_READ_WRITE ,sizeof(cl_float3)*rays_size);
        
        cl_int* RNDS=new cl_int[rays_size];
        for(int i=0;i<rays_size;++i){
            RNDS[i]=i+1;
        }
        queue.enqueueWriteBuffer(buffer_rnds,CL_TRUE,0,sizeof(cl_int)*rays_size,&RNDS[0]);
        delete[] RNDS;
        
        camera=Camera();
        tri_shift=0;
    }
    void add_Triangle(Triangle tri){
        tris.push_back(tri);
    }
    void add_Material(Material mat){
        mats.push_back(mat);
    }
    void start_New_Obj(){
        if(kd_tree_shift.empty())
            kd_tree_shift.push_back((cl_int){0});
        else
            kd_tree_shift.push_back((cl_int){kd_tree.size()-1});
        KDNode* root=new KDNode();
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
        queue.finish();
    }
    void trace_rays(){
        glFinish();
        clEnqueueAcquireGLObjects(queue(), 1, &cl_screen(), 0, NULL, NULL);
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
        queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(screen_width, screen_height),cl::NullRange);
        
//        cl::Kernel kernel_filt_im=cl::Kernel(program,"filt_im");
//        kernel_filt_im.setArg(0,cl_screen);
//        kernel_filt_im.setArg(1,buffer_colors);
//        queue.enqueueNDRangeKernel(kernel_filt_im,cl::NullRange,cl::NDRange(screen_width, screen_height),cl::NullRange);
        
        queue.finish();
        clEnqueueReleaseGLObjects(queue(), 1, &cl_screen(), 0, NULL, NULL);
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

    unsigned short LAMP, SUN, WHITE_DIFFUSE, RED_DIFFUSE, GREEN_DIFFUSE, BLACK_DIFFUSE, CHROMIUM, GLASS, GOLD;
    //                                                           diffuse_color                  specular_color                     emission                      refractive_index              extinction_coefficient        shininess       type
    LAMP=0;             scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){60.0f*2, 50.0f*2, 40.0f*2}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){  0}, (cl_int){3}));
    SUN=1;              scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){60.0f*10, 50.0f*10, 40.0f*10}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){  0}, (cl_int){3}));
    WHITE_DIFFUSE=2;    scene.add_Material(Material((cl_float3){0.3f, 0.3f, 0.3f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){ 50}, (cl_int){0}));
    RED_DIFFUSE=3;      scene.add_Material(Material((cl_float3){0.3f, 0.1f, 0.1f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){ 50}, (cl_int){0}));
    GREEN_DIFFUSE=4;    scene.add_Material(Material((cl_float3){0.1f, 0.3f, 0.1f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){ 50}, (cl_int){0}));
    BLACK_DIFFUSE=5;    scene.add_Material(Material((cl_float3){0.0f, 0.1f, 0.1f}, (cl_float3){0.3f, 0.3f, 0.3f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.00f, 0.00f, 0.00f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){100}, (cl_int){0}));
    CHROMIUM=6;         scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){3.10f, 3.05f, 2.05f}, (cl_float3){3.3f, 3.3f, 2.9f}, (cl_float){  0}, (cl_int){1}));
    GOLD=7;             scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){0.17f, 0.35f, 1.50f}, (cl_float3){3.1f, 2.7f, 1.9f}, (cl_float){  0}, (cl_int){1}));
    GLASS=8;            scene.add_Material(Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){ 0.0f,  0.0f,  0.0f}, (cl_float3){1.50f, 1.50f, 1.50f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){  0}, (cl_int){2}));
    
    //lámpa
    scene.add_Triangle(Triangle((cl_float3){300.0f, 999.9f, 700.0f}, (cl_float3){300.0f, 999.9f, 300.0f}, (cl_float3){700.0f, 999.9f, 700.0f}, LAMP));
    scene.add_Triangle(Triangle((cl_float3){700.0f, 999.9f, 700.0f}, (cl_float3){300.0f, 999.9f, 300.0f}, (cl_float3){700.0f, 999.9f, 300.0f}, LAMP));
    
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
    scene.add_Triangle(Triangle((cl_float3){-100.0f, 0.0f, 1000.0f}, (cl_float3){-100.0f, 1000.0f, 1000.0f}, (cl_float3){1100.0f, 1000.0f, 1000.0f}, WHITE_DIFFUSE));
    scene.add_Triangle(Triangle((cl_float3){1100.0f, 1000.0f, 1000.0f}, (cl_float3){1100.0f, 0.0f, 1000.0f}, (cl_float3){-100.0f, 0.0f, 1000.0f}, WHITE_DIFFUSE));
    
    //balra
    scene.add_Triangle(Triangle((cl_float3){-100.0f, 0.0f, 1000.0f}, (cl_float3){-100.0f, 0.0f, -1000.0f}, (cl_float3){-100.0f, 1000.0f, 1000.0f}, RED_DIFFUSE));
    scene.add_Triangle(Triangle((cl_float3){-100.0f, 1000.0f, 1000.0f}, (cl_float3){-100.0f, 0.0f, -1000.0f}, (cl_float3){-100.0f, 1000.0f, -1000.0f}, RED_DIFFUSE));
    
    //jobbra
    scene.add_Triangle(Triangle((cl_float3){1100.0f, 1000.0f, 1000.0f}, (cl_float3){1100.0f, 0.0f, -1000.0f}, (cl_float3){1100.0f, 0.0f, 1000.0f}, GREEN_DIFFUSE));
    scene.add_Triangle(Triangle((cl_float3){1100.0f, 1000.0f, -1000.0f}, (cl_float3){1100.0f, 0.0f, -1000.0f}, (cl_float3){1100.0f, 1000.0f, 1000.0f}, GREEN_DIFFUSE));
    
    //felül
    scene.add_Triangle(Triangle((cl_float3){-100.0f, 1000.0f, 1000.0f}, (cl_float3){-100.0f, 1000.0f, -1000.0f}, (cl_float3){1100.0f, 1000.0f, 1000.0f}, WHITE_DIFFUSE));
    scene.add_Triangle(Triangle((cl_float3){1100.0f, 1000.0f, 1000.0f}, (cl_float3){-100.0f, 1000.0f, -1000.0f}, (cl_float3){1100.0f, 1000.0f, -1000.0f}, WHITE_DIFFUSE));
    scene.start_New_Obj();
    //hátul
//    scene.add_Triangle(Triangle((cl_float3){-100.0f-1000, 0.0f, -1000.0f-1000}, (cl_float3){-100.0f-1000, 1000.0f, -1000.0f-1000}, (cl_float3){1100.0f+1000, 1000.0f, -1000.0f-1000}, WHITE_DIFFUSE));
//    scene.add_Triangle(Triangle((cl_float3){1100.0f+1000, 1000.0f, -1000.0f-1000}, (cl_float3){1100.0f+1000, 0.0f, -1000.0f-1000}, (cl_float3){-100.0f-1000, 0.0f, -1000.0f-1000}, WHITE_DIFFUSE));
    
    //alul
    scene.add_Triangle(Triangle((cl_float3){-10000.0f, 0.0f, -10000.0f}, (cl_float3){-10000.0f, 0.0f, 10000.0f}, (cl_float3){10000.0f, 0.0f, 10000.0f}, WHITE_DIFFUSE));
    scene.add_Triangle(Triangle((cl_float3){10000.0f, 0.0f, 10000.0f}, (cl_float3){10000.0f, 0.0f, -10000.0f}, (cl_float3){-10000.0f, 0.0f, -10000.0f}, WHITE_DIFFUSE));
    scene.start_New_Obj();
    
    cl_float3 sun_pos[6];
    sun_pos[0]=(cl_float3){-20000.0f, 300000.0f, 20000.0f}; sun_pos[1]=(cl_float3){-20000.0f, 300000.0f, -20000.0f}; sun_pos[2]=(cl_float3){20000.0f, 300000.0f, 20000.0f};
    sun_pos[3]=(cl_float3){20000.0f, 300000.0f, 20000.0f}; sun_pos[4]=(cl_float3){-20000.0f, 300000.0f, -20000.0f}; sun_pos[5]=(cl_float3){20000.0f, 300000.0f, -20000.0f};
    for(int i=0;i<6;++i){
        sun_pos[i]=rotate_x(sun_pos[i], 45);
        sun_pos[i]=rotate_y(sun_pos[i], 150);
    }
    scene.add_Triangle(Triangle(sun_pos[0], sun_pos[1], sun_pos[2], SUN));
    scene.add_Triangle(Triangle(sun_pos[3], sun_pos[4], sun_pos[5], SUN));
    scene.start_New_Obj();

    //arany cucc
//    cl_float3 move=(cl_float3){750.0f, 0.001f, 300.0f};
    cl_float3 move=(cl_float3){950.0f, 0.001f, 850.0f-500};
//    cl_float3 move=(cl_float3){400.0f, 100.001f, -2400.0f};
    cl_float3 scale=(cl_float3){1.5f, 0.28284f*0.28284f/40*500, 1.5f};
//    cl_float3 scale=(cl_float3){6.0f, 0.28284f*0.28284f/40*800, 6.0f};
    //felül
    
//    cl_float3 move=(cl_float3){50.0f+350-200, 0.000f, 850-350-1300};
//    cl_float3 scale=(cl_float3){1.5f, 0.28284f*0.28284f/40*300, 1.5f};
//    cl_float3 move=(cl_float3){-50.0f, 0.001f, 950};
//    cl_float3 scale=(cl_float3){0.0125f, 0.28284f*0.28284f/40*1.25, 0.0125f};
    for(int szar3=0;szar3<0;++szar3){
        for(int szar2=0;szar2<1;++szar2){
            for(int szar=0;szar<1;++szar){
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
    scene.start_New_Obj();
    
    //move=(cl_float3){400.0f, 200.0f, 600.0f};
    //scale=(cl_float3){1.5f/200.0f*400.0f, 0.28284f*0.28284f/40*400, 1.5f/200.0f*400.0f};
    move=(cl_float3){400.0f-200, 200.0f, 600.0f-100};
    scale=(cl_float3){1.5f/200.0f*200.0f, 0.28284f*0.28284f/40*200, 1.5f/200.0f*200.0f};
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
    float movx=200;
    float movy=0;
    float movz=-150-300;
    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 500.0f+movz}, GLASS));
    scene.add_Triangle(Triangle((cl_float3){400.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, GLASS));
    
    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
    scene.add_Triangle(Triangle((cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){0.0f+movx, 0.0f+movy, 600.0f+movz}, GLASS));
    
    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
    
    scene.add_Triangle(Triangle((cl_float3){400.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
    scene.add_Triangle(Triangle((cl_float3){400.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
    
    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 600.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 600.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 600.0f+movy, 600.0f+movz}, GLASS));
    
    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 600.0f+movz}, GLASS));
    scene.add_Triangle(Triangle((cl_float3){0.0f+movx, 0.0f+movy, 500.0f+movz}, (cl_float3){0.0f+movx, 0.0f+movy, 600.0f+movz}, (cl_float3){400.0f+movx, 0.0f+movy, 600.0f+movz}, GLASS));
    scene.start_New_Obj();
    
    
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
//        scene.add_Triangle(Triangle((cl_float3){-100.0f+movx, 0.0f+movy, 0.0f+movz}, (cl_float3){100.0f+movx, 0.0f+movy, 0.0f+movz}, (cl_float3){100.0f+movx, 200.0f+movy, 200.0f+movz}, 8+i));
//        scene.add_Triangle(Triangle((cl_float3){100.0f+movx, 200.0f+movy, 200.0f+movz}, (cl_float3){-100.0f+movx, 0.0f+movy, 0.0f+movz}, (cl_float3){-100.0f+movx, 200.0f+movy, 200.0f+movz}, 8+i));
//        movx=movx+200;
//    }
    
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
    if(keys_down[W] || keys_down[A] || keys_down[S] || keys_down[D] || mouse_down){
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
    
    scene.generate_rays();
    scene.trace_rays();
    if(real_time){
        glutPostRedisplay();
    }
    end=clock();
    double elapsed_secs=double(end-begin)/CLOCKS_PER_SEC;
    
    current_sample++;
    
    if(elapsed_secs>1.0){
        reset_timer=true;
        clear_line();
        printf("Samples=%010d  Samples/sec=%08.3f Render time=%08.3fms  real_time=%d  Iterations=%02d  Million ray/sec=%08.3f  Elapsed seconds=%f", current_sample, (current_sample-old_sample)/elapsed_secs, elapsed_secs/(current_sample-old_sample)*1000.0f, real_time, iterations, (current_sample-old_sample)/elapsed_secs*iterations*screen_width*screen_height/1000000.0f, glutGet(GLUT_ELAPSED_TIME)/1000.0f-start);
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