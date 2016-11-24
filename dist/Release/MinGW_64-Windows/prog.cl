typedef struct{
    float3 bl,tr; //bottom left, top right
} BBox;

typedef struct{
    int2 trii;
    float split;
    BBox bbox;
} Node;

typedef struct{
    float3 kd,ks,emission,F0;
    float n, shininess;
    int type;
} Material;

typedef struct{
    float3 P,D; //origin,direction
} Ray;

bool BBox_intersection(global const BBox* box,
                        Ray* ray,
                        float* tmin,
                        float* tmax){
    float tx1=(box->bl.x-ray->P.x)/(ray->D.x);
    float tx2=(box->tr.x-ray->P.x)/(ray->D.x);
    *tmin=fmin(tx1,tx2);
    *tmax=fmax(tx1,tx2);

    float ty1=(box->bl.y-ray->P.y)/(ray->D.y);
    float ty2=(box->tr.y-ray->P.y)/(ray->D.y);
    *tmin=fmax(*tmin, fmin(ty1,ty2));
    *tmax=fmin(*tmax, fmax(ty1,ty2));

    float tz1=(box->bl.z-ray->P.z)/(ray->D.z);
    float tz2=(box->tr.z-ray->P.z)/(ray->D.z);
    *tmin=fmax(*tmin, fmin(tz1, tz2));
    *tmax=fmin(*tmax, fmax(tz1, tz2));

    return *tmax>=*tmin;
}

typedef struct{
    float t;
    float3 P,N;
    ushort mati;
    Material mat;
} Hit;

typedef struct{
    float3 r1,r2,r3,N;
    ushort mati;
} Triangle;

typedef struct{
    float3 eye, lookat, up, right;
    float XM, YM; //screen width, screen height
} Camera;

void rand(global float2* seed){
    int const a = 48271;
    int const m = 2147483647;

    int s=(int)(((*seed).x*2.0f-1.0f)*m);
    s = ((long)(s * a))%m;
    (*seed).x=(s/2147483647.0f+1.0f) / 2.0f;

    s=(int)(((*seed).y*2.0f-1.0f)*m);
    s = ((long)(s * a))%m;
    (*seed).y=(s/2147483647.0f+1.0f) / 2.0f;
}

void orthonormal_base(const float3* V1, float3* V2, float3* V3){
    const float E=0.001f;
    float3 v1,v2,v3;
    v1=(*V1);v2=(*V2);v3=(*V3);
    if(fabs(v1.x)<=E && fabs(v1.z)<=E){
        float length=sqrt(v1.y*v1.y + v1.z*v1.z);
        length=1/length;
        v2.x=0;
        v2.y=-v1.z*length;
        v2.z=v1.y*length;
    }
    else
    {
        float length=sqrt(v1.x*v1.x + v1.z*v1.z);
        length=1/length;
        v2.x=-v1.z*length;
        v2.y=0;
        v2.z=v1.x*length;
    }
    v3=cross(v1,v2);
    (*V2)=v2;
    (*V3)=v3;
}

float3 Fresnel(Hit* hit, Ray* old_ray){
    float cosa=fabs(dot(hit->N, old_ray->D));
    return hit->mat.F0 + ((float3)(1.0f, 1.0f, 1.0f)-hit->mat.F0)*pow(1-cosa, 5);
}

Ray cons_Ray(float3 p, float3 d){
    Ray ray; ray.P=p; ray.D=d; return ray;
}

Hit cons_Hit(float t, float3 p, float3 n, ushort mati){
    Hit hit; hit.t=t; hit.P=p; hit.N=n; hit.mati=mati; return hit;
}

Hit init_Hit(){
    return cons_Hit(-1.0f, (float3)(0.0f, 0.0f, 0.0f), (float3)(0.0f, 0.0f, 0.0f), 0);
}

Ray camera_get_ray(int id, const Camera* cam, float2 rnds){
    int X=cam->XM;
    int Y=cam->YM;
    float x=id%X+rnds.x;
    float y=id/X+rnds.y;
    float3 right=cam->right*(2.0f*x/X-1);
    float3 up=cam->up*(2*y/Y-1);
    float3 p=cam->lookat + right + up;
    float3 d=normalize(p-cam->eye);
    
    return cons_Ray(cam->eye, d);
}

float3 camera_view_dir(Hit hit, const Camera cam){
    return normalize((cam.eye-hit.P));
}

Hit triangle_intersect(global const Triangle* tri, const Ray* ray){
    Hit hit=init_Hit();
    float3 P=ray->P;
    float3 V=ray->D;
    float3 N=tri->N;
    float t=dot((tri->r1-P),N)/dot(V,N);
    if(t<0){
        return hit;
    }
    float3 p=P+V*t;
    if( dot( cross((tri->r2-tri->r1),(p-tri->r1)) , N) >= 0){
        if( dot( cross((tri->r3-tri->r2),(p-tri->r2)) , N)>=0){
            if( dot( cross((tri->r1-tri->r3),(p-tri->r3)) , N)>=0){
                return cons_Hit(t,p,N,tri->mati);
            }
        }
    }
    return hit;
}

//Hit first_intersect(global const Triangle* tris, const int tris_size, const Ray ray){
Hit first_intersect(global const Triangle* tris, const int from, const int to, const Ray ray){
    Hit best_hit=init_Hit();
    for(int i=from; i<to; ++i){
        Hit hit=triangle_intersect(&tris[i], &ray);
        if(hit.t>0 && (best_hit.t<0 || hit.t<best_hit.t)){
            best_hit=hit;
        }
    }
    return best_hit;
}

void new_ray_diffuse(Hit* hit, float2 rnds, global Ray* ray){
    float3 X,Y,Z;
    Y=hit->N;
    orthonormal_base(&Y,&Z,&X);
    float rnd1,rnd2,r,theta,x,y,z;
    rnd1=rnds.x;
    rnd2=rnds.y;
    r=sqrt(rnd1);
    theta=2*M_PI*rnd2;
    x=r*cos(theta);
    y=r*sin(theta);
    z=sqrt(1-rnd1);
    float3 new_d=normalize(X*x+Y*z+Z*y);
    (*ray)=cons_Ray(hit->P+Y*0.001f, new_d);
}

Ray new_ray_specular(Hit hit, Ray old_ray){
    float3 new_d=normalize(old_ray.D - hit.N*dot(hit.N, old_ray.D)*2.0f);
    return cons_Ray(hit.P+hit.N*0.001f, new_d);
}

Ray new_ray_refractive(Hit hit, Ray old_ray, bool in, float rnd){
    if(in){
        hit.mat.n=1.0f/hit.mat.n;
    }
    hit.mat.n=1.0f/hit.mat.n;
    float cosa=-dot(old_ray.D, hit.N);
    float disc=1.0f - (1.0f - cosa*cosa)*hit.mat.n*hit.mat.n;
    float3 F=hit.mat.F0 + ((float3)(1.0f, 1.0f, 1.0f) - hit.mat.F0)*pow(1-cosa, 5);
    float prob=(F.x+F.y+F.z)/3.0f;
    if(disc>0 && rnd>prob){
        return cons_Ray(hit.P - hit.N*0.001f, normalize(old_ray.D*hit.mat.n + hit.N*(cosa*hit.mat.n - sqrt(disc))));
    }else{
        return cons_Ray(hit.P + hit.N*0.001f, normalize(old_ray.D + hit.N*cosa*2.0f));
    }
}

float4 filmic_tone(float3 c){
    c=max(0.0f, c-0.004f);
    c=(c*(c*6.2f+0.5f))/(c*(c*6.2f+1.7f)+0.06f);
    c=pow(c, 2.2f);
    //c=c/(c+1);
    return (float4)(c, 1.0f);
}

void stack_push(int* stack, int* ptr, int val){
    if((*ptr)<300){
        stack[*ptr]=val;
        (*ptr)=(*ptr)+1;
    }
}
int stack_pop(int* stack, int* ptr){
    if((*ptr)>0){
        *ptr=*ptr-1;;
        return stack[*ptr];
    }
    return stack[0];
}

Hit kd_intersect(global const Triangle* tris, global const Node* kd_tree, const int kd_tree_size, const Ray ray ){
    Hit hit=init_Hit();
    Hit best_hit=init_Hit();
    int ptr=1;
    float tmin=999999;;
    float dist=0;
    float tmax=-999999;
    int stack[300];
    int stack_ptr=0;
    bool fail=false;
    bool empty=false;
    while(!empty && !fail){
        if(ptr<kd_tree_size){
            if(BBox_intersection(&kd_tree[ptr].bbox, &ray, &dist, &tmax)){
                if(tmax>=0){
                    if(dist>tmin){  // check if bbox closer than the best hit, becouse if it's not, it's impossible for tha bbox to have a triangle closer what we already have
                        if(stack_ptr==0){
                            empty=true;
                        }else{
                            ptr=stack_pop(stack, &stack_ptr);
                        }
                    }else if(kd_tree[ptr].trii.x<0){  // check if not leaf
                        stack_push(stack, &stack_ptr, 2*ptr+1);
                        ptr=2*ptr;
                    }else{  // if leaf, it will contain triangles
                        hit=first_intersect(tris, kd_tree[ptr].trii.x, kd_tree[ptr].trii.y, ray);   // x and y is from and to indexes, and tris is an orderes array of triangles
                        if(hit.t>0 && (best_hit.t<0 || hit.t<best_hit.t)){ // keep track of the best hit
                            tmin=hit.t;
                            best_hit=hit;
                        }
                        if(stack_ptr==0){
                            empty=true;
                        }else{  // check other nodes
                            ptr=stack_pop(stack, &stack_ptr);
                        }
                    }
                }else{
                    if(stack_ptr==0){
                        empty=true;
                    }else{
                        ptr=stack_pop(stack, &stack_ptr);
                    }
                }
            }else if(stack_ptr==0){
                empty=true;
            }else{
                ptr=stack_pop(stack, &stack_ptr);
            }
        }else{
            fail=true;
        }
    }
    return best_hit;
}

void kernel trace_ray(write_only image2d_t tex,
                        global const Triangle* tris,
                        const int tris_size,
                        global const Material* materials,
                        global const Node* kd_tree,
                        const int kd_tree_size,
                        global Ray* rays,
                        global float2* rnds,
                        const int iterations,
                        const int current_sample,
                        const Camera cam,
                        global float3* colors){

    int id=get_global_id(1)*get_global_size(0) + get_global_id(0);
    float3 factor_L=(float3)(1.0f, 1.0f, 1.0f);
    float3 factor_B=(float3)(1.0f, 1.0f, 1.0f);
    float3 factor_S=(float3)(1.0f, 1.0f, 1.0f);
    float3 factor_R=(float3)(1.0f, 1.0f, 1.0f);
    float3 color=(float3)(0.05f, 0.05f, 0.05f);
    
    if(current_sample==0){
        colors[id]=color;
    }

    bool in=false;
    for(int current=0; current<iterations; ++current){
        Hit hit=first_intersect(tris, 0, tris_size, rays[id]);
        //Hit hit=kd_intersect(tris, kd_tree, kd_tree_size, rays[id]);

        if(hit.t>0){
            hit.mat=materials[hit.mati];
            if(iterations==1){
                color=hit.mat.kd+hit.mat.emission;
            }
            if(dot(rays[id].D,hit.N)>0){                                                                                                            // hence the angle between D and N will always be less than 90 degree
                hit.N=-hit.N;
            }
            if(hit.mat.type==0){                                                                                    // diffuse
                rand(&rnds[id]);
                new_ray_diffuse(&hit, rnds[id], &rays[id]);

                float cos_theta=dot(rays[id].D, hit.N);
                float intensity_diffuse=fmax(0.0f, cos_theta);
                factor_L=factor_L*(hit.mat.kd*intensity_diffuse);

                float3 halfway=normalize(camera_view_dir(hit, cam) + rays[id].D);
                float cos_delta=dot(hit.N, halfway);
                float intensity_specular=fmax(0.0f, cos_delta);
                factor_B=factor_B*(hit.mat.ks*pow(intensity_specular, hit.mat.shininess));
                color=color+(float3)(0.05f, 0.05f, 0.05f)*hit.mat.kd;
            }
            if(hit.mat.type==1){                                                                                    // specular
                Ray old_ray=rays[id];
                rays[id]=new_ray_specular(hit, old_ray);
                factor_S=factor_S*Fresnel(&hit, &old_ray);
            }
            if(hit.mat.type==2){                                                                                    // refractive
                Ray old_ray=rays[id];
                rand(&rnds[id]);
                rays[id]=new_ray_refractive(hit, old_ray, in, rnds[id].x);
                factor_R=factor_R*(1-Fresnel(&hit, &old_ray));
                in=!in;
            }
            if(hit.mat.type==3){                                                                                    // emitter
                rand(&rnds[id]);
                new_ray_diffuse(&hit, rnds[id], &rays[id]);
                color=color + hit.mat.emission*(factor_L + factor_B)*factor_S*factor_R;
            }
        }else{
            break;
        }
    }

    colors[id]=(colors[id]*current_sample + color)/(current_sample+1);
    write_imagef(tex, (int2)(get_global_id(0), get_global_id(1)), filmic_tone(colors[id]));
}

void kernel gen_ray(global Ray* rays,
                    const Camera camera,
                    global float2* rnds){
    int id=get_global_id(0);
    rand(&rnds[id]);
    rays[id]=camera_get_ray(id, &camera, rnds[id]);
}

void kernel filt_im(write_only image2d_t tex, global float3* colors){
    int x=get_global_id(0);
    int y=get_global_id(1);
    int width=get_global_size(0);
    int height=get_global_size(1);
    float3 arr[9];
    if(x>0 && y>0 && x<width && y<height){
        for(int i=0;i<3;++i){
            for(int j=0;j<3;++j){
                int id=(y-1+i)*width + (x-1+j);
                arr[i*3+j]=colors[id];
            }
        }
        float grey[9];
        for(int i=0;i<9;++i){
            float v=(arr[i].x+arr[i].y+arr[i].z)/3.0f;
            grey[i]=v;
        }
        for(int j=9;j>1;j--){
            int maxi=0;
            for(int i=1;i<j;++i){
                if(grey[i]>grey[maxi]){
                    maxi=i;
                }
            }
            float temp=grey[j-1];
            grey[j-1]=grey[maxi];
            grey[maxi]=temp;

            float3 temp2=arr[j-1];
            arr[j-1]=arr[maxi];
            arr[maxi]=temp2;
        }
        float3 med=(float3)(arr[4].x, arr[4].y, arr[4].z);
        write_imagef(tex, (int2)(x, y), filmic_tone(med));
    }
}
