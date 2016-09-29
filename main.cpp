#include <windows.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <cmath>
#include <CL/cl.hpp>

const int screenWidth=600;
const int screenHeight=400;
float sinvar=0;

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
    float h_rotate=0.0f;
    float v_rotate=0.0f;
    
//    float up_length=500.0f;
//    float right_length=500.0f*((float)screenWidth/screenHeight);
//    float ahead_length=right_length/tan(fov/2.0f/180.0f*3.141593f);
//    cam.lookat=(cl_float3){500.0f, 500.0f, 0.0f};
//    cam.eye=(cl_float3){cam.lookat.s[0], cam.lookat.s[1], cam.lookat.s[2]-ahead_length};
    
    
    cam.lookat=(cl_float3){500.0f, 500.0f, 0.0f};
    cam.eye=(cl_float3){cam.lookat.s[0], cam.lookat.s[1], cam.lookat.s[2]-500.0f/tan(fov/2.0f/180.0f*3.141593f)};
    cam.up=(cl_float3){0.0f, 500, 0.0f};
    cam.right=(cl_float3){500*((float)screenWidth/screenHeight), 0.0f, 0.0f};
    cam.XM=(cl_float){(float)screenWidth};
    cam.YM=(cl_float){(float)screenHeight};
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

Color color_image[screenWidth*screenHeight];

class Scene{
private:
    Camera camera;
    std::vector<Triangle> tris;
    int tris_size;
    int rays_size=screenWidth*screenHeight;
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
        
        cl_float3_image=new cl_float3[rays_size];
        camera=cons_Camera();
    }
    void add_Triangle(Triangle tri){
        tris.push_back(tri);
    }
    void upload_Triangles(){
        tris_size=tris.size();
        buffer_tris=cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(Triangle)*tris_size);
        queue.enqueueWriteBuffer(buffer_tris,CL_TRUE,0,sizeof(Triangle)*tris_size,&tris[0]);
        
//        for(int i=0;i<tris_size;++i){
//            printf("[%06.2f %06.2f %06.2f]\n", tris[i].N.s[0], tris[i].N.s[1], tris[i].N.s[2]);
//        }
    }
    void generate_rays(){
        camera=cons_Camera();
        cl::Kernel kernel_gen_ray=cl::Kernel(program,"gen_ray");
        kernel_gen_ray.setArg(0,buffer_rays);
        kernel_gen_ray.setArg(1,camera);
        
        queue.enqueueNDRangeKernel(kernel_gen_ray,cl::NullRange,cl::NDRange(rays_size),cl::NullRange);
    }
    void trace_rays(){
        //clock_t begin=clock();
        //run the kernel
        cl::Kernel kernel_trace_ray=cl::Kernel(program,"trace_ray");
        kernel_trace_ray.setArg(0,buffer_tris);
        kernel_trace_ray.setArg(1,tris_size);
        kernel_trace_ray.setArg(2,buffer_rays);
        kernel_trace_ray.setArg(3,buffer_colors);
        
        queue.enqueueNDRangeKernel(kernel_trace_ray,cl::NullRange,cl::NDRange(rays_size),cl::NullRange);

        queue.enqueueReadBuffer(buffer_colors,CL_TRUE,0,sizeof(cl_float3)*rays_size,cl_float3_image);
        
        clock_t end=clock();
        //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        //printf("%f\n", elapsed_secs);
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

Scene scene;

void onInitialization( ) { 
    srand(time(0));
    glViewport(0, 0, screenWidth, screenHeight);

    for(int Y = 0; Y < screenHeight; Y++)
        for(int X = 0; X < screenWidth; X++){
            color_image[Y*screenWidth + X] = Color((float)Y/screenHeight, (float)X/screenWidth, 0);
        }
    
    
    scene.init_Scene();
    
//    for(int i=0;i<25;++i){
//        Material mat=cons_Material((cl_float3){0.5f, 0.3f, 1.0f+i}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
//        scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 0.0f, 1000.0f+i}, (cl_float3){0.0f, 1000.0f, 1000.0f+i}, (cl_float3){1000.0f, 1000.0f, 1000.0f+i}, (cl_float3){0.0f, 0.0f, -1.0f}, mat));
//        scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 1000.0f+i}, (cl_float3){1000.0f, 0.0f, 1000.0f+i}, (cl_float3){0.0f, 0.0f, 1000.0f+i}, (cl_float3){0.0f, 0.0f, -1.0f}, mat));
//    }
    Material mat;
    
    //elől
    mat=cons_Material((cl_float3){1.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 0.0f, 1000.0f}, (cl_float3){0.0f, 1000.0f, 1000.0f}, (cl_float3){1000.0f, 1000.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 1000.0f}, (cl_float3){1000.0f, 0.0f, 1000.0f}, (cl_float3){0.0f, 0.0f, 1000.0f}, mat));
    
    //balra
    mat=cons_Material((cl_float3){0.0f, 1.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 0.0f, 1000.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 1000.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 1000.0f, 1000.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 1000.0f, 0.0f}, mat));
    
    //jobbra
    mat=cons_Material((cl_float3){0.0f, 0.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 1000.0f}, (cl_float3){1000.0f, 0.0f, 0.0f}, (cl_float3){1000.0f, 0.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 0.0f}, (cl_float3){1000.0f, 0.0f, 0.0f}, (cl_float3){1000.0f, 1000.0f, 1000.0f}, mat));
    
    //alul
    mat=cons_Material((cl_float3){1.0f, 1.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 1000.0f}, (cl_float3){1000.0f, 0.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 0.0f, 1000.0f}, (cl_float3){1000.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, mat));
    
    //felül
    mat=cons_Material((cl_float3){0.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){0.0f, 1000.0f, 1000.0f}, (cl_float3){0.0f, 1000.0f, 0.0f}, (cl_float3){1000.0f, 1000.0f, 1000.0f}, mat));
    scene.add_Triangle(cons_Triangle((cl_float3){1000.0f, 1000.0f, 1000.0f}, (cl_float3){0.0f, 1000.0f, 0.0f}, (cl_float3){1000.0f, 1000.0f, 0.0f}, mat));
    
    mat=cons_Material((cl_float3){1.0f, 1.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){250.0f, 0.0f, 250.0f}, (cl_float3){500.0f, 500.0f, 500.0f}, (cl_float3){750.0f, 0.0f, 250.0f}, mat));
    mat=cons_Material((cl_float3){0.0f, 1.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){250.0f, 0.0f, 250.0f}, (cl_float3){250.0f, 0.0f, 750.0f}, (cl_float3){500.0f, 500.0f, 500.0f}, mat));
    mat=cons_Material((cl_float3){0.0f, 0.0f, 1.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float3){1.5f, 1.5f, 1.5f}, (cl_float3){0.0f, 0.0f, 0.0f}, (cl_float){5}, (cl_int){0});
    scene.add_Triangle(cons_Triangle((cl_float3){750.0f, 0.0f, 250.0f}, (cl_float3){500.0f, 500.0f, 500.0f}, (cl_float3){750.0f, 0.0f, 750.0f}, mat));
    
    fflush(stdout);
}

void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDrawPixels(screenWidth, screenHeight, GL_RGB, GL_FLOAT, color_image);
    
    glutSwapBuffers();
}

void onKeyboard(unsigned char key, int x, int y) {
    if(key==' '){
        printf("Render!\n");fflush(stdout);
        scene.upload_Triangles();
        scene.generate_rays();
        scene.trace_rays();
        glutPostRedisplay();
        printf("Render finished!\n");fflush(stdout);
    }
}

void onKeyboardUp(unsigned char key, int x, int y) {
    
}

void onMouse(int button, int state, int x, int y) {
    
}

void onMouseMotion(int x, int y){
    
}

void onIdle( ) {
    //clock_t begin=clock();
    float time=glutGet(GLUT_ELAPSED_TIME)/1000.0f;
    float sinus=sin(time);
    sinvar=sinus*500.0f;
    scene.upload_Triangles();
    scene.generate_rays();
    scene.trace_rays();
    //clock_t end=clock();
    //float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
    //float fps=1.0f/elapsed_secs;
    //printf("\r                                                                                                                                       \r");
    //printf("%f fps",fps);
    fflush(stdout);
    glutPostRedisplay();
}

int main(int argc, char **argv) {
    glutInit(&argc, argv); 				// GLUT inicializalasa
    glutInitWindowSize(screenWidth, screenHeight);	// Alkalmazas ablak kezdeti merete 600x600 pixel 
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