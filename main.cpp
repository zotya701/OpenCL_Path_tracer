#include <windows.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <random>
#include <thread>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <CL/cl.hpp>

const int screen_width=300;
const int screen_height=300;
const int max_iterations=5;
int iterations=1;
int current_iteration=0;
float global_yaw=0;
float global_pitch=0;
float global_forward=0;;
float global_rightward=0;
cl_float3 global_shift=(cl_float3){0.0f, 0.0f, 0.0f};
enum ControllKeys {W, A, S, D, keys_num};
bool keys_down[keys_num];

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0,1.0);

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

typedef struct{
    cl_float3 kd,ks,emission,F0;    //diffuse, specular, emission, Fresnel
    cl_float n, shininess, glossiness;
    cl_int type;   //0-diffuse, 1-x, 2-x, 3-Emitter
} Material;
Material cons_Material(cl_float3 kd, cl_float3 ks, cl_float3 emission, cl_float3 N, cl_float3 K, cl_float shininess, cl_int type){
    Material mat; mat.kd=kd; mat.ks=ks; mat.emission=emission; mat.shininess=shininess; mat.type=type;
    mat.n=(cl_float){(N.s[0]+N.s[1]+N.s[2])/3.0f};
    float F0[3];
    for(int i=0;i<3;++i){
        float a=(N.s[i]-1)*(N.s[i]-1);
        float b=(N.s[i]+1)*(N.s[i]+1);
        F0[i]=(K.s[i]*K.s[i]+a)/(K.s[i]*K.s[i]+b);
    }
    mat.F0=(cl_float3){F0[0], F0[1], F0[2]};
    return mat;
}

typedef struct{
    cl_float3 P,D;  //origo and direction
} Ray;
Ray cons_Ray(cl_float3 p, cl_float3 d){
    Ray ray; ray.P=p; ray.D=d; return ray;
}

typedef struct{
    cl_float t;         //time
    cl_float3 P,N;      //hitposition and normal vector in hitposition
    Material mat;       //material of the triangle what was hit by the ray
} Hit;

typedef struct{
    cl_float3 r1,r2,r3,N;   //vertices of the triangle and it's normal vector
    Material mat;
} Triangle;
Triangle cons_Triangle(cl_float3 r1, cl_float3 r2, cl_float3 r3, Material mat){
    Triangle tri; tri.r1=r1; tri.r2=r2; tri.r3=r3; tri.mat=mat;
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
    
    tri.N=(cl_float3){n[0], n[1], n[2]};
    return tri;
}

typedef struct{
    cl_float3 eye, lookat, up, right;
    cl_float XM, YM;
} Camera;
Camera cons_Camera(){
    Camera cam;
    float fov=60;
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

    cam.eye=(cl_float3){500.0f+global_shift.s[0], 500.0f+global_shift.s[1], -1299.037842f+global_shift.s[2]};
    
    cam.up=(cl_float3){up.s[0]*up_length, up.s[1]*up_length, up.s[2]*up_length};
    cam.right=(cl_float3){right.s[0]*right_length, right.s[1]*right_length, right.s[2]*right_length};
    cam.lookat=(cl_float3){cam.eye.s[0]+ahead.s[0]*ahead_length, cam.eye.s[1]+ahead.s[1]*ahead_length, cam.eye.s[2]+ahead.s[2]*ahead_length};

    cam.XM=(cl_float){(float)screen_width};
    cam.YM=(cl_float){(float)screen_height};
    return cam;
}

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

Color color_image[screen_width*screen_height];

class Scene{
private:
    Camera camera;
    std::vector<Triangle> tris;
    int tris_size;
    int rays_size=screen_width*screen_height;
    cl_float3* cl_float3_image;
    cl_float2* RNDS;
    
    cl::Context context;
    cl::Program program;
    cl::Buffer buffer_tris;
    cl::Buffer buffer_rays;
    cl::Buffer buffer_rnds;
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
        
        buffer_rays=cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(Ray)*rays_size*max_iterations);
        buffer_rnds=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(cl_float2)*rays_size*(max_iterations-1));
        buffer_colors=cl::Buffer(context,CL_MEM_WRITE_ONLY ,sizeof(cl_float3)*rays_size);
        
        cl_float3_image=new cl_float3[rays_size];
        RNDS=new cl_float2[rays_size*(max_iterations-1)];
        camera=cons_Camera();
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
        camera=cons_Camera();
        cl::Kernel kernel_gen_ray=cl::Kernel(program,"gen_ray");
        kernel_gen_ray.setArg(0,buffer_rays);
        kernel_gen_ray.setArg(1,camera);
        
        queue.enqueueNDRangeKernel(kernel_gen_ray,cl::NullRange,cl::NDRange(rays_size),cl::NullRange);
    }
    void trace_rays(){
        int size=rays_size*(iterations-1);
        for(int i=0;i<size;++i){
            RNDS[i].s[0]=distribution(generator);
            RNDS[i].s[1]=distribution(generator);
        }
        queue.enqueueWriteBuffer(buffer_rnds,CL_TRUE,0,sizeof(cl_float2)*rays_size*(iterations-1),&RNDS[0]);

        //run the kernel
        cl::Kernel kernel_trace_ray=cl::Kernel(program,"trace_ray");
        kernel_trace_ray.setArg(0,buffer_tris);
        kernel_trace_ray.setArg(1,tris_size);
        kernel_trace_ray.setArg(2,buffer_rays);
        kernel_trace_ray.setArg(3,buffer_rnds);
        kernel_trace_ray.setArg(4,iterations);
        kernel_trace_ray.setArg(5,current_iteration);
        kernel_trace_ray.setArg(6,camera);
        kernel_trace_ray.setArg(7,buffer_colors);
        
        queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(rays_size),cl::NullRange);
        
        queue.enqueueReadBuffer(buffer_colors,CL_TRUE,0,sizeof(cl_float3)*rays_size,cl_float3_image);
        
        for(int i=0;i<rays_size;++i){
            cl_float3 c=cl_float3_image[i];
            //printf("cl_float3_image[%03d] r=%06.2f g=%06.2f b=%06.2f\n", i, c.s[0], c.s[1], c.s[2]);
            color_image[i]=Color(c.s[0], c.s[1], c.s[2]);
        }
    }
    void finish(){
        queue.finish();
    }
};

void task1(std::string msg)
{
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    while(1){
        //std::cout << "task1 says: " << msg << std::endl;
        //printf("%d\n",sizeof(cl_float2));
        //fflush(stdout);
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
}

Scene scene;
void onInitialization( ) { 
    srand(time(0));
    glViewport(0, 0, screen_width, screen_height);

    for(int Y = 0; Y < screen_height; Y++)
        for(int X = 0; X < screen_width; X++){
            color_image[Y*screen_width + X] = Color((float)Y/screen_height, (float)X/screen_width, 0);
        }
    
    
    scene.init_Scene();

    Material mat;
    
    //lámpa
    mat=cons_Material((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){50.0f, 50.0f, 50.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){0}, (cl_int){3});
    scene.add_Triangle(cons_Triangle((cl_float3){400.0f, 999.0f, 600.0f}, (cl_float3){400.0f, 999.0f, 400.0f}, (cl_float3){600.0f, 999.0f, 600.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){600.0f, 999.0f, 600.0f}, (cl_float3){400.0f, 999.0f, 400.0f}, (cl_float3){600.0f, 999.0f, 400.0f}, mat));
    
    //elől
    mat=cons_Material((cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 0.0f, 1000.0f}, (cl_float3){0.0f, 1000.0f, 1000.0f}, (cl_float3){1000.0f, 1000.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 1000.0f}, (cl_float3){1000.0f, 0.0f, 1000.0f}, (cl_float3){0.0f, 0.0f, 1000.0f}, mat));
    
    //balra
    mat=cons_Material((cl_float3){1.0f, 0.0f, 0.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 0.0f, 1000.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 1000.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 1000.0f, 1000.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 1000.0f, 0.0f}, mat));
    
    //jobbra
    mat=cons_Material((cl_float3){0.0f, 1.0f, 0.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 1000.0f}, (cl_float3){1000.0f, 0.0f, 0.0f}, (cl_float3){1000.0f, 0.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 0.0f}, (cl_float3){1000.0f, 0.0f, 0.0f}, (cl_float3){1000.0f, 1000.0f, 1000.0f}, mat));
    
    //alul
    mat=cons_Material((cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 1000.0f}, (cl_float3){1000.0f, 0.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 0.0f, 1000.0f}, (cl_float3){1000.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, mat));
    
    //felül
    mat=cons_Material((cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 1000.0f, 1000.0f}, (cl_float3){0.0f, 1000.0f, 0.0f}, (cl_float3){1000.0f, 1000.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 1000.0f}, (cl_float3){0.0f, 1000.0f, 0.0f}, (cl_float3){1000.0f, 1000.0f, 0.0f}, mat));
    
    mat=cons_Material((cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){250.0f, 0.0f, 250.0f}, (cl_float3){500.0f, 500.0f, 500.0f}, (cl_float3){750.0f, 0.0f, 250.0f}, mat));
    //mat=cons_Material((cl_float3){0.0f, 1.0f, 0.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){250.0f, 0.0f, 250.0f}, (cl_float3){250.0f, 0.0f, 750.0f}, (cl_float3){500.0f, 500.0f, 500.0f}, mat));
    //mat=cons_Material((cl_float3){0.0f, 0.0f, 1.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){750.0f, 0.0f, 250.0f}, (cl_float3){500.0f, 500.0f, 500.0f}, (cl_float3){750.0f, 0.0f, 750.0f}, mat));
    //mat=cons_Material((cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){10}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){250.0f, 0.0f, 750.0f}, (cl_float3){500.0f, 500.0f, 500.0f}, (cl_float3){750.0f, 0.0f, 750.0f}, mat));
    
    std::thread t1(task1, "Hello");
    t1.detach();
    printf("%d\n",sizeof(Ray));
}

void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDrawPixels(screen_width, screen_height, GL_RGB, GL_FLOAT, color_image);
    
    glutSwapBuffers();
}

bool space=false;
void onKeyboard(unsigned char key, int x, int y) {
    if(key=='-'){
        if(iterations>1){
            iterations--;
            current_iteration=0;
        }
    }
    if(key=='+'){
        if(iterations<max_iterations){
            iterations++;
            current_iteration=0;
        }
    }
    if(key==' '){
        space=!space;
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
            current_iteration=0;
        break;
        case 's': case 'S':
            keys_down[S] = false;
            current_iteration=0;
        break;
        case 'a': case 'A':
            keys_down[A] = false;
            current_iteration=0;
        break;
        case 'd': case 'D':
            keys_down[D] = false;
            current_iteration=0;
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
        current_iteration=0;
    }

    if ((button == GLUT_LEFT_BUTTON ) && (state == GLUT_UP)){
        mouse_down=false;
        current_iteration=0;
    }
}
 
void onMouseMotion(int x, int y) {
    //camera.updateDir(x-last_x, last_y-y);
    int dx=x-last_x;
    int dy=y-last_y;
    float speed=0.2f;
    global_yaw=global_yaw+dx*speed;
    global_pitch=global_pitch+dy*speed;
    last_x = x;
    last_y = y;
}

float old=0;
float newTime=0;
float dt=0;
void onIdle( ) {
    int before=iterations;
    if(keys_down[W] || keys_down[A] || keys_down[S] || keys_down[D] || mouse_down){
        if(!space)
            iterations=1;
        current_iteration=0;
    }
    
    old = newTime;
    newTime = glutGet(GLUT_ELAPSED_TIME)/1000.0f;
    dt=newTime-old;
    
//    const long nrolls=100000000;  // number of experiments
//    const int nstars=95;     // maximum number of stars to distribute
//    const int nintervals=10; // number of intervals
//
//    std::default_random_engine generator;
//    std::uniform_real_distribution<double> distribution(0.0,1.0);
//    double number = distribution(generator);

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
    
    scene.upload_Triangles();
    scene.generate_rays();
    scene.trace_rays();
    
    fflush(stdout);
    glutPostRedisplay();
    iterations=before;
    current_iteration++;
}

int main(int argc, char **argv) {
    glutInit(&argc, argv); 				// GLUT inicializalasa
    glutInitWindowSize(screen_width, screen_height);	// Alkalmazas ablak kezdeti merete 600x600 pixel 
    glutInitWindowPosition(100, 100);			// Az elozo alkalmazas ablakhoz kepest hol tunik fel
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);	// 8 bites R,G,B,A + dupla buffer + melyseg buffer

    glutCreateWindow("Path tracer");                    // Alkalmazas ablak megszuletik es megjelenik a kepernyon

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